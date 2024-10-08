#!/usr/bin/env python
'''
Author: Ifsoul
Date: 2020-08-29 16:01:13
LastEditTime: 2022-03-14 21:56:18
LastEditors: Ifsoul
Description: Functions for file IO
'''
import numpy as np
import scipy.sparse as SP
import os
import zlib
import msgpack
# import multiprocessing as mp
# import ujson
from chardet import detect

from .data import Elements, Basis_List
from .parallel import Job, run_functions_thread
# from .parallel import run_functions
from .core import atom_no, atom_type


def save_npz(File, Data, SPMod=False):
    if SPMod:
        SP.save_npz(File, SP.csc_matrix(np.ascontiguousarray(Data)))
    else:
        if isinstance(Data, dict):
            np.savez_compressed(File, **{k: v for k, v in Data.items()})
        else:
            np.savez_compressed(File, np.ascontiguousarray(Data))


def load_npz(File, SPMod=False):
    if SPMod:
        try:
            data = SP.load_npz(File).tocsc()
        except Exception:
            data = SP.csc_matrix(np.load(File)["arr_0"])
    else:
        try:
            fdata = np.load(File)
            if len(fdata.files) > 1:
                data = {k: fdata[k] for k in fdata.files}
            else:
                data = fdata[fdata.files[0]]
        except Exception:
            data = SP.load_npz(File).todense()
    return data


# def loadnpz(file, spmod=False):
#     if spmod:
#         try:
#             data = SP.load_npz(file).tocsc()
#         except Exception:
#             data = SP.csc_matrix(np.load(file)["arr_0"])
#     else:
#         try:
#             data = np.load(file)["arr_0"]
#         except Exception:
#             data = SP.load_npz(file).todense()
#     return data

# def savenpz(file, data, spmod=False):
#     if spmod:
#         SP.save_npz(file, SP.csc_matrix(data))
#     else:
#         np.savez_compressed(file, data)


def load_csv(path):
    return np.loadtxt(path, delimiter=',')


def load_trj(Input_Name):
    assert os.path.isfile(Input_Name)
    with open(Input_Name, 'r') as Input:
        Input.readline()
        Input.readline()
        # Time=int(Input.readline())
        Input.readline()
        AtomNum = int(Input.readline())
        Input.readline()
        BoxX = [float(x) for x in Input.readline().split()]
        BoxY = [float(x) for x in Input.readline().split()]
        BoxZ = [float(x) for x in Input.readline().split()]
        BoxL = np.vstack((BoxX, BoxY, BoxZ))
        Input.readline()
        Data = np.zeros((AtomNum, 5))
        for i in range(AtomNum):
            Data[i, :] = [float(x) for x in Input.readline().split()]
    Data[:, 2:] = Data[:, 2:] * np.dot(np.ones((AtomNum, 1)), np.resize(BoxL[:, 1] - BoxL[:, 0], (1, 3))) + BoxL[:, 0]
    Data = Data[np.argsort(Data[:, 0]), :]
    return Data, BoxL


def save_trj(Output_Name, Output_Data, Box):
    OutLines = [
        "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n",
        "%d\n" % (len(Output_Data)),
        "ITEM: BOX BOUNDS pp pp pp\n",
        "%.16f\t%.16f\n" % (Box[0, 0], Box[0, 1]),
        "%.16f\t%.16f\n" % (Box[1, 0], Box[1, 1]),
        "%.16f\t%.16f\n" % (Box[2, 0], Box[2, 1]),
        "ITEM: ATOMS id type xs ys zs\n",
    ]
    Output_Data[:, 2:] = (Output_Data[:, 2:] - Box[:, 0]) / np.resize(Box[:, 1] - Box[:, 0], (1, 3))
    for Atom in Output_Data:
        OutLines.append("%d %d %.8f %.8f %.8f\n" % tuple(Atom[:5]))
    OutLines.append("\n")
    with open(Output_Name, 'w') as Output:
        Output.writelines(OutLines)


def load_xyz(Input_Name):
    assert os.path.isfile(Input_Name)
    with open(Input_Name, 'r') as Input:
        AtomNum = int(Input.readline())
        Input.readline()
        Data = np.zeros((AtomNum, 6))
        Data[:, 0] = np.arange(1, AtomNum + 1, 1)
        for i in range(AtomNum):
            StrList = Input.readline().split()
            ele = StrList[0][0].upper() + StrList[0][1:].lower() if len(StrList[0]) > 1 else StrList[0].upper()
            assert ele in Elements, "ERROR: Unkown atom type %s\n" % (ele)
            Data[i, 1] = Elements.index(ele) + 1
            Data[i, 3:] = [float(x) for x in StrList[1:]]
    return Data


def save_xyz(Output_Name, Output_Data, Col_Type=1, Col_Coor=3):
    AtomNum = len(Output_Data)
    TypeList = np.array(Elements)[Output_Data[:, Col_Type].astype(int) - 1]
    Data = Output_Data[:, Col_Coor:Col_Coor + 3]
    OutLines = ["%d\n\tGenerated by Ifsoul\n" % AtomNum]
    for i in range(AtomNum):
        OutLines.append('%s\t%16.10f\t%16.10f\t%16.10f\n' % (TypeList[i].upper(), Data[i, 0], Data[i, 1], Data[i, 2]))
    with open(Output_Name, 'w') as outfile:
        outfile.writelines(OutLines)


# import h5py
# def save_h5(File, Data):
#     with h5py.File(File, 'w') as f:
#         if isinstance(Data, dict):
#             for k, v in Data.items():
#                 if isinstance(v, (tuple, list, set, np.ndarray)):
#                     if isinstance(v[0], str):
#                         s = np.array(v, dtype=h5py.string_dtype(encoding='utf-8'))
#                         f.create_dataset(k, data=s, compression='lzf')
#                         continue
#                 f.create_dataset(k, data=v, compression='lzf')
#         else:
#             f.create_dataset('dset1', data=Data, compression='lzf')

# def load_h5(File):
#     with h5py.File(File, 'r') as f:
#         if len(f.keys()) > 1:
#             data = {k: f[k][()] for k in f.keys()}
#         else:
#             data = f[list(f.keys())[0]][()]
#     return data


def dtype_to_typecode(DataType):
    if DataType == int:
        return 'i'
    elif DataType == float:
        return 'd'
    else:
        raise KeyError(DataType)


typecode_to_dtype = {'i': int, 'd': float}

if msgpack.version[0] >= 1:
    msg_load_kwarg = {'strict_map_key': False}
else:
    msg_load_kwarg = {'encoding': 'utf-8'}

# def init_worker(X, X_shape):
#     # Using a dictionary is not strictly necessary. You can also use global variables.
#     global var_dict
#     var_dict = {}
#     var_dict['X'] = X
#     var_dict['X_shape'] = X_shape

# def compress_worker(FragIdx, Compress=None):
#     X_np = np.frombuffer(var_dict['X']).reshape(var_dict['X_shape'])
#     return zlib.compress(X_np[FragIdx[0]:FragIdx[1], :])

# def decompress_worker(Data, TypeCode, FragIdx):
#     X_np = np.frombuffer(var_dict['X'].get_obj()).reshape(var_dict['X_shape'])
#     data = zlib.decompress(Data)
#     X_np[FragIdx[0]:FragIdx[1], :] = np.frombuffer(data, dtype=typecode_to_dtype[TypeCode]).reshape((FragIdx[1] - FragIdx[0], -1))

# def save_fsz_process(Output_Name, Data, PoolSize=Job.CPUAvail):
#     Lx = Data.shape[0] / PoolSize
#     xlist = [(int(i * Lx), int((i + 1) * Lx)) for i in range(PoolSize)]
#     X = mp.RawArray(dtype_to_typecode(Data.dtype), Data.shape[0] * Data.shape[1])
#     X_np = np.frombuffer(X).reshape(Data.shape)
#     np.copyto(X_np, Data)
#     funcs = [(compress_worker, (idx,)) for idx in xlist]
#     data_compressed = run_functions(funcs, PoolSize, PoolKwargs={'initializer': init_worker, 'initargs': (X, Data.shape)})
#     data_pack = {'data': (Data.shape, dtype_to_typecode(Data.dtype), xlist, data_compressed)}
#     with open(Output_Name, 'wb') as f:
#         msgpack.dump(data_pack, f, use_bin_type=True)

# def load_fsz_process(Input_Name, PoolSize=Job.CPUAvail):
#     with open(Input_Name, 'rb') as f:
#         data_pack = msgpack.load(f, use_list=False, **msg_load_kwarg)
#     DataShape, TypeCode, FragIdxs, data_compressed = data_pack[list(data_pack.keys())[0]]
#     NFrags = len(FragIdxs)
#     data = mp.Array(TypeCode, DataShape[0] * DataShape[1])
#     funcs = [(decompress_worker, (data_compressed[i], TypeCode, FragIdxs[i])) for i in range(NFrags)]
#     run_functions(funcs, PoolSize, PoolKwargs={'initializer': init_worker, 'initargs': (data, DataShape)})
#     return np.frombuffer(data.get_obj()).reshape(DataShape)


def save_fsz_thread(Output_Name, Data, PoolSize=Job.CPUAvail):
    """Save .fsz file"""
    if isinstance(Data, dict):
        DataDict = Data
    else:
        DataDict = {'data': Data}
    data_pack = {}
    for k, v in DataDict.items():
        if isinstance(v, np.ndarray) and v.dtype in (int, float):
            DataShape = v.shape
            TypeCode = dtype_to_typecode(v.dtype)
            if PoolSize > 1 and v.size > 10000 * PoolSize:
                LFrags = v.size / PoolSize
                FragIdxs = tuple((int(i * LFrags), int((i + 1) * LFrags)) for i in range(PoolSize))
                funcs = [(zlib.compress, (v.ravel()[idx[0]:idx[1]],)) for idx in FragIdxs]
                data_compressed = run_functions_thread(funcs, PoolSize)
            else:
                FragIdxs = ((0, v.size),)
                data_compressed = bytes(v)
            data_pack[k] = (b'ndarray', DataShape, TypeCode, FragIdxs, data_compressed)
        else:
            if isinstance(v, np.ndarray):
                v = v.tolist()
            data_pack[k] = v
    with open(Output_Name, 'wb') as f:
        msgpack.dump(data_pack, f, use_bin_type=True)


def decompress_worker_thread(Data, X, DataType, FragIdx):
    X[FragIdx[0]:FragIdx[1]] = np.frombuffer(zlib.decompress(Data), dtype=DataType)


def load_fsz_thread(Input_Name, PoolSize=Job.CPUAvail, AlwaysReturnDict=False):
    """Read .fsz file"""
    with open(Input_Name, 'rb') as f:
        data_pack = msgpack.load(f, use_list=False, **msg_load_kwarg)
    DataDict = {}
    for k, v in data_pack.items():
        if isinstance(v, tuple) and v[0] == b'ndarray' and len(v) == 5:
            DataShape, TypeCode, FragIdxs, data_compressed = v[1:]
            DataType = typecode_to_dtype[TypeCode]
            if len(FragIdxs) > 1:
                X = np.zeros(DataShape, dtype=DataType)
                funcs = [(decompress_worker_thread, (data_compressed[i], X.ravel(), DataType, idx)) for i, idx in enumerate(FragIdxs)]
                run_functions_thread(funcs, PoolSize)
            else:
                X = np.frombuffer(data_compressed, dtype=DataType).reshape(DataShape)
        else:
            X = v
        DataDict[k] = X
    if len(DataDict) > 1 or AlwaysReturnDict:
        return DataDict
    else:
        return DataDict[list(DataDict.keys())[0]]


# def read_obj(Input_Name):
#     """Read a object from file. File format should be like .fsz"""
#     ObjDict = load_fsz_thread(Input_Name, PoolSize=1, AlwaysReturnDict=True)
#     return ObjDict

# def save_obj(Output_Name, Obj):
#     """Save a object as file. File format is like .fsz"""
#     OutDict={}
#     save_fsz_thread(Output_Name, OutDict, PoolSize=1)


def msg_dumps_and_zip(something):
    return zlib.compress(msgpack.dumps(something, use_bin_type=True))


def msg_unzip_and_loads(something):
    return msgpack.loads(zlib.decompress(something), use_list=False, **msg_load_kwarg)


def load_from_msgzip(filename):
    """Load something from msgzip file"""
    if not os.path.exists(filename):
        raise AssertionError('File not found: %s' % filename)
    with open(filename, 'rb') as f:
        something = f.read()
    return msg_unzip_and_loads(something)


def save_as_msgzip(filename, something):
    """Save something as msgzip into file"""
    with open(filename, 'wb') as f:
        f.write(msg_dumps_and_zip(something))

def read_file(Input_Name, DataOnly=False, *args, **kwargs):
    """Read .csv/txt/xyz/lammpstrj/npz/npy/fsz file, return (FileName,Data)"""
    assert os.path.isfile(Input_Name), 'ERROR: Fail to find file %s\n' % Input_Name
    temp = os.path.splitext(os.path.split(Input_Name)[1])
    FName = temp[0]
    ext = temp[1][1:]
    if ext == 'npy':
        Data = np.load(Input_Name, *args, **kwargs)
    elif ext == 'npz':
        Data = load_npz(Input_Name, *args, **kwargs)
    elif ext in ['fsz']:
        Data = load_fsz_thread(Input_Name, *args, **kwargs)
    elif ext == 'msgz':
        Data = load_from_msgzip(Input_Name)
    # elif ext in ['h5', 'hdf5']:
    #     Data = load_h5(Input_Name, *args, **kwargs)
    elif ext == 'csv':
        Data = np.loadtxt(open(Input_Name, "r"), delimiter=",")
    elif ext == 'txt':
        Data = np.loadtxt(open(Input_Name, "r"), delimiter="\t")
    elif ext == 'xyz':
        Data = load_xyz(Input_Name, *args, **kwargs)
    elif ext == 'lammpstrj':
        Data, Box = load_trj(Input_Name, *args, **kwargs)
        if DataOnly:
            return Data, Box
        else:
            return FName, Data, Box
    else:
        raise AssertionError('Unrecognized file format!')
    if DataOnly:
        return Data
    else:
        return FName, Data


def save_file(Output_Name, Output_Data, XYZPos=-1, XYZLen=3, Format='%.8g', Header='', Comments='', *args, **kwargs):
    """Save .csv/txt/xyz/npz/npy/fsz file"""
    temp = os.path.splitext(os.path.split(Output_Name)[1])
    ext = temp[1][1:]
    if ext == 'npy':
        return np.save(Output_Name, np.ascontiguousarray(Output_Data), *args, **kwargs)
    elif ext == 'npz':
        return save_npz(Output_Name, Output_Data, *args, **kwargs)
    elif ext == 'fsz':
        return save_fsz_thread(Output_Name, Output_Data, *args, **kwargs)
    elif ext == 'msgz':
        return save_as_msgzip(Output_Name, Output_Data)
    # elif ext in ['h5', 'hdf5']:
    #     save_h5(Output_Name, Output_Data, *args, **kwargs)
    elif ext == 'csv':
        if XYZPos >= 0:
            StrFmt = []
            for i in range(Output_Data.shape[1]):
                if i < XYZPos or i >= XYZPos + XYZLen:
                    StrFmt.append('%d')
                else:
                    StrFmt.append('%.12g')
            StrFmt = ','.join(StrFmt)
            return np.savetxt(Output_Name, Output_Data, fmt=StrFmt, header=Header, comments=Comments)
        else:
            return np.savetxt(Output_Name, Output_Data, fmt=Format, delimiter=",", header=Header, comments=Comments)
    elif ext == 'txt':
        if XYZPos >= 0:
            StrFmt = []
            for i in range(Output_Data.shape[1]):
                if i < XYZPos or i >= XYZPos + XYZLen:
                    StrFmt.append('%4d')
                else:
                    StrFmt.append('%.12g')
            StrFmt = '\t'.join(StrFmt)
            return np.savetxt(Output_Name, Output_Data, fmt=StrFmt, header=Header, comments=Comments)
        else:
            return np.savetxt(Output_Name, Output_Data, fmt=Format, delimiter="\t", header=Header, comments=Comments)
    elif ext == 'xyz':
        return save_xyz(Output_Name, Output_Data, *args, **kwargs)
    else:
        raise AssertionError('Unrecognized file format!')


def read_txt_lines(FileName, SPLIT=False, IncludeComment=False, INTACT=False, CommentMark = '#', encoding=''):
    """Read lines from txt file, skipping blank and comment lines. Return a list."""
    if not encoding:
        with open(FileName, 'rb') as f:
            cur_encoding = detect(f.read(10000))['encoding']
    else:
        cur_encoding = encoding
    with open(FileName, 'rt', encoding=cur_encoding) as f:
        lines = f.readlines()
    if INTACT:
        return lines
    TreatMethod = str.split if SPLIT else str.strip
    return list(TreatMethod(l) for l in lines if l.strip() and (IncludeComment or l.strip()[0] != CommentMark))


class IO_File(object):

    def __init__(self, FILENAME):
        self.rename(FILENAME)
        # self.Name = os.path.abspath(FILENAME)
        # temp = os.path.split(self.Name)
        # self.Path = temp[0]
        # temp = os.path.splitext(temp[1])
        # self.Proj = temp[0]
        # self.Ext = temp[1][1:]

    def rename(self, NewName):
        self.Name = os.path.abspath(NewName)
        temp = os.path.split(self.Name)
        self.Path = temp[0]
        temp = os.path.splitext(temp[1])
        self.Proj = temp[0]
        self.Ext = temp[1][1:]


class info_file(IO_File):

    def __init__(self, FILENAME, NPART=0, NACTIVE=0, TOPOTYPE='NULL'):
        super().__init__(FILENAME)
        self.NPart = NPART
        self.NActive = NACTIVE
        self.TopoType = TOPOTYPE.lower()
        if self.Proj[-5:] == '_Info':
            self.Proj = self.Proj[:-5]
        self.C_Name = self.Proj + '_Part'
        self.H_Name = self.Proj + '_H'
        self.S_Name = self.Proj + '_S'
        self.ProjSet = self.Proj
        self.Charge = []
        self.RelaxType = 'Z'
        self.PartKey = []
        self.FoundInfo = []
        self.Basis = '6-31g*'
        self.PartType = []

    def load(self):
        assert os.path.isfile(self.Name), "ERROR: Cannot find file %s!\n" % self.Name
        assert self.Ext.lower() == 'txt', "ERROR: Wrong file format %s!\n" % self.Ext
        LineList = [l.strip() for l in open(self.Name, 'rt').readlines()]
        LineList = [l for l in LineList if (l and not l.startswith('#'))]
        try:
            self.NPart, self.NActive, topo = (int(x) for x in LineList[0].split())
            self.TopoType = 'line' if topo == 0 else 'ring' if topo == 1 else 'null'
        except Exception as e:
            print('Exception %s: Fail to read line:\n\t%s' % (e, LineList[0]))
        for i in range(1, len(LineList)):
            ListTmp = LineList[i].split()
            if ListTmp[0] == 'C':
                self.C_Name = ListTmp[1]
            elif ListTmp[0] == 'H':
                self.H_Name = ListTmp[1]
            elif ListTmp[0] == 'S':
                self.S_Name = ListTmp[1]
            elif ListTmp[0] == 'P':
                self.ProjSet = ListTmp[1]
            elif ListTmp[0] == 'Charge':
                self.Charge = [int(x) for x in ListTmp[1:]]
                if len(self.Charge) < self.NPart:
                    self.Charge.extend([0 for _ in range(self.NPart - len(self.Charge))])
                    print('WARNING: Length of Charge list < NPart!\n\tSetting 0 charge for rest parts.')
                elif len(self.Charge) > self.NPart:
                    self.Charge = self.Charge[:self.NPart]
                    print('WARNING: Length of Charge list > NPart!\n\tIgnoring extra values.')
            elif ListTmp[0] == 'PartType':
                assert len(LineList) >= i + 1 + self.NPart
                self.PartType = [l.split()[1] for l in LineList[i + 1:i + 1 + self.NPart]]
            elif ListTmp[0] == 'Relax':
                self.RelaxType = ListTmp[1]
                if self.RelaxType in ['N', 'R']:
                    assert len(LineList) >= i + 1 + self.NPart
                    assert np.array([int(l.split()[0]) == n for n, l in enumerate(LineList[i + 1:i + 1 + self.NPart])]).all()
                    self.PartKey = [l.split()[1] for l in LineList[i + 1:i + 1 + self.NPart]]
            elif ListTmp[0] == 'DB':
                assert len(LineList) >= i + 1 + self.NPart
                self.FoundInfo = [l.split()[1:] for l in LineList[i + 1:i + 1 + self.NPart]]
            elif ListTmp[0] == 'Basis':
                self.Basis = ListTmp[1]
                assert self.Basis in Basis_List, 'ERROR: Unknown basis type %s!\n' % self.Basis
        return self

    def save(self):
        assert self.Ext.lower() == 'txt', "ERROR: Wrong file format %s!\n" % self.Ext
        print('%s info file!' % ('Updating' if os.path.isfile(self.Name) else 'Writing'))
        with open(self.Name, 'wt') as f:
            tp = 0 if self.TopoType == 'line' else 1 if self.TopoType == 'ring' else -1
            f.write('%d\t%d\t%d\n' % (self.NPart, self.NActive, tp))
            f.write('C\t%s\nH\t%s\nS\t%s\n' % (self.C_Name, self.H_Name, self.S_Name))
            if self.ProjSet != self.Proj:
                f.write('P\t%s\n' % (self.ProjSet))
            if self.Charge and np.array(self.Charge).any():
                Clist = '\t'.join([str(x) for x in self.Charge])
                f.write('Charge\t%s\n' % Clist)
            if self.Basis != '6-31g*':
                f.write('Basis\t%s\n' % (self.Basis))
            if self.PartType:
                f.write('PartType\n')
                for i, name in enumerate(self.PartType):
                    f.write('%d\t%s\n' % (i, name))
            f.write('Relax\t%s\n' % (self.RelaxType))
            if self.RelaxType in ['N', 'R']:
                for i, key in enumerate(self.PartKey):
                    f.write('%d\t%s\n' % (i, key))
            if self.FoundInfo:
                f.write('\nDB\n')
                for i, FI in enumerate(self.FoundInfo):
                    if len(FI) == 1:
                        f.write('%d\t%s\n' % (i, FI[0]))
                    else:
                        f.write('%d\t%s\t%s\n' % (i, FI[0], FI[1]))


class gjf_file(IO_File):
    '''.gjf file'''

    def __init__(self, FILENAME):
        super().__init__(FILENAME)
        self.ProjectName, self.CheckName = '', ''
        self.CPUNum = 1
        self.Memory = '40GB'
        self.Basis = 'b3lyp/6-31g'
        self.P_ON = True
        self.CtrlKWs = []
        self.Charge = 0
        self.SpinMul = 1
        self.Data = np.array([])
        self.RestLines = []

    def load(self):
        assert os.path.isfile(self.Name), "ERROR: Cannot find file %s!\n" % self.Name
        assert self.Ext.lower() == 'gjf', "ERROR: Wrong file format %s!\n" % self.Ext
        LineList = [l.strip() for l in open(self.Name, 'rt').readlines()]
        ParaLines = []
        PrintLine = ''
        StartLineNo = 0
        for i, l in enumerate(LineList):
            if l.startswith('%'):
                ParaLines.append(l)
            elif l.startswith('#'):
                PrintLine = l
                StartLineNo = i + 1
                break
        LCount = 0
        for i, l in enumerate(LineList[StartLineNo:]):
            if l:
                if LCount == 0:
                    self.ProjectName = l
                    LCount = 1
                elif LCount == 1:
                    self.Charge, self.SpinMul = [int(_) for _ in l.split()]
                    StartLineNo += i + 1
                    break
        Data = []
        for i, l in enumerate(LineList[StartLineNo:]):
            if l:
                ListTmp = l.split()
                if ListTmp[0] in Elements:
                    Data.append([atom_no(ListTmp[0])] + [float(_) for _ in ListTmp[1:]])
                else:
                    StartLineNo += i + 1
                    self.RestLines = [l for l in LineList[StartLineNo:] if l]
                    break
        else:
            self.RestLines = []
        for l in ParaLines:
            if l.startswith('%chk='):
                self.CheckName = l[5:l.rfind('.')].strip()
            elif l.startswith('%mem='):
                self.Memory = l[5:].strip()
            elif l.startswith('%nprocshared='):
                self.CPUNum = int(l[l.find('=') + 1:].strip())
        self.update_ctrlkw(PrintLine)
        # PrintInfos = PrintLine.split()
        # self.P_ON = True if PrintInfos[0] == '#p' else False
        # self.Basis = PrintInfos[1]
        # self.CtrlKWs = PrintInfos[2:]
        self.Data = np.array(Data)

        return self

    def save(self):
        with open(self.Name, 'w') as outfile:
            outfile.write('%%chk=%s\n' % (self.CheckName + '.chk'))
            if self.CPUNum > 1:
                outfile.write('%%nprocshared=%d\n' % (self.CPUNum))
            if self.Memory:
                outfile.write('%%mem=%s\n' % (self.Memory))
            outfile.write('#%s %s %s\n' % ('p' if self.P_ON else '', self.Basis, ' '.join(self.CtrlKWs)))
            outfile.write('\n%s\n\n' % (self.ProjectName))
            outfile.write('%d %d\n' % (self.Charge, self.SpinMul))
            for i in range(self.Data.shape[0]):
                outfile.write('%s\t%14.8f\t%14.8f\t%14.8f\n' % (atom_type(self.Data[i, 0]), self.Data[i, 1], self.Data[i, 2], self.Data[i, 3]))
            if self.RestLines:
                outfile.write('\n'.join(self.RestLines))
            outfile.write('\n\n\n')

    def update_ctrlkw(self, CtrlStr):
        PrintInfos = CtrlStr.split()
        self.P_ON = True if PrintInfos[0] == '#p' else False
        self.Basis = PrintInfos[1]
        self.CtrlKWs = PrintInfos[2:]
