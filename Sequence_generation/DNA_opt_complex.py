#!/usr/bin/env python
import numpy as np
import sys
import fs
import cvxpy as cp
import os

str_usage = """\
Usage: DNA_opt_complex.py (Mtx_File) (Options)
Generate DNA sequences for target DOS.
Options:
    -e start end            Ranges of energy, default: (-5.8,-5.1) eV
    -e0 start end           Ranges of zero-DOS energy, default: (-5.55,-5.45) eV
    -homo E_homo            Limit the HOMO energy
    -lumo E_lumo            Limit the LUMO energy
    -eV                     Use eV as unit of energy, default: Hartree(27.211eV)
    -m name_file            Name file convering shortname into full name, default: HOMO.txt
    -o output_file          Default: Fit-bms_(E0st,E0ed).png
    -l max_length           Max length of designed DNA, default: 50
    -cf CurveFile           File containing DNA DOS curves
    -tr DNA_trained         File containing trained DNAs, default: DNA_train.txt
    -curv u/d ratio         Parameter for target curve
    -freestyle t a b ...    Free style target curve, t:Type(u/d) (a,b):ERange
    -plotall                Plot all basemode curves into file
    -pdf                    Plot curves into pdf file
    -mix                    Use mix-type basemodes
    -ns                     No single base
    -notext                 Do not plot text
    -h                      Show this help page
    """

SinglePeakH = 0.025


def show_usage():
    print(str_usage)
    sys.exit()


def cutoff_func(x, Threshold=SinglePeakH):
    return cp.minimum(x, Threshold)


def vect_times(a, b):
    return cp.sum(cp.multiply(a, b))


def solve_mix(
    CurrentBase,
    MyTag,
    MixBMNames,
    MixBMCurves,
    HeadTailBase,
    Type=None,
    Verbose=False,
    IsDopant=True,
    TransFunc=None,
    TagLength=50,
):
    """Solve mix type problem"""
    with fs.timer('Fit Mix %s Type %s Dopants' % (CurrentBase, Type)):
        x = cp.Variable(len(MixBMNames), integer=True)
        l = np.array([len(_) - 1 for _ in MixBMNames])
        # obj = cp.Minimize(cp.norm(MixBMCurves.T @ x - MyTag, 2))
        if TransFunc is None:
            obj = cp.Minimize(cp.sum_squares(MixBMCurves.T @ x - MyTag))
        else:
            obj = cp.Minimize(cp.sum_squares(TransFunc(MixBMCurves.T @ x) - MyTag))
        # constriants = [cp.sum(cp.multiply(l, x)) <= 29, x >= 0]
        if IsDopant:
            constriants = [vect_times(l, x) <= int(0.7 * TagLength), x >= 0, cp.sum(x) >= 1]
        else:
            constriants = [vect_times(l, x) == TagLength - 1, x >= 0, cp.sum(x) >= 1]
        if 'C' in HeadTailBase.keys():
            v = HeadTailBase['C']
            vp = np.zeros_like(v)
            vm = np.zeros_like(v)
            vp[v > 0] = 1
            vm[v < 0] = 1
        if 'G' in HeadTailBase.keys():
            v = HeadTailBase['G']
            vp2 = np.zeros_like(v)
            vm2 = np.zeros_like(v)
            vp2[v > 0] = 1
            vm2[v < 0] = 1
        CG_CPLX = ('C' in HeadTailBase.keys()) and ('G' in HeadTailBase.keys())
        if Type is None:
            if IsDopant:
                for k, v in HeadTailBase.items():
                    if k == 'C':
                        if k == CurrentBase:
                            constriants.append(vect_times(vp, x) == 1)
                        else:
                            constriants.extend([vect_times(vp, x) <= 1, vect_times(vm, x) <= 1])
                    else:
                        if k == CurrentBase:
                            constriants.append(vect_times(v, x) == 1)
                        else:
                            constriants.extend([vect_times(v, x) >= -1, vect_times(v, x) <= 0])
            else:
                for k, v in HeadTailBase.items():
                    if k == 'C':
                        if k == CurrentBase:
                            constriants.append(vect_times(vp, x) <= 1)
                        else:
                            constriants.extend([vect_times(vp, x) <= 1, vect_times(vm, x) <= 1])
                    else:
                        if k == CurrentBase:
                            constriants.append(vect_times(v, x) <= 1)
                        else:
                            constriants.extend([vect_times(v, x) >= -1, vect_times(v, x) <= 0])
        elif Type in (0, 1, 2):
            if Type == 0:
                for k, v in HeadTailBase.items():
                    constriants.append(vect_times(v, x) == 0)
                    if k == 'C':
                        constriants.extend([vect_times(vp, x) <= 1, vect_times(vm, x) <= 1])
                    elif k == 'G':
                        constriants.extend([vect_times(vp2, x) <= 1, vect_times(vm2, x) <= 1])
                if CG_CPLX:
                    constriants.append(vect_times(vp + vp2, x) <= 1)
            else:
                t = 1
                for k, v in HeadTailBase.items():
                    if k == CurrentBase:
                        # constriants.append(vect_times(v, x) == 1)
                        if k == 'C':
                            constriants.append(vect_times(vp, x) == 1)
                        else:
                            constriants.append(vect_times(v, x) == 1)
                    else:
                        if t == Type:
                            # constriants.append(vect_times(v, x) == -1)
                            if k == 'C':
                                constriants.append(vect_times(vm, x) == 1)
                            else:
                                constriants.append(vect_times(v, x) == -1)
                        else:
                            # constriants.append(vect_times(v, x) == 0)
                            if k == 'C':
                                constriants.extend([vect_times(vp, x) <= 1, vect_times(vm, x) <= 1])
                            else:
                                constriants.append(vect_times(v, x) == 0)
                        t += 1
        else:
            assert isinstance(Type, str)
            assert len(Type) == 2
            BaseH, BaseL = Type
            for k, v in HeadTailBase.items():
                if k == 'C':
                    if k == BaseH:
                        constriants.append(vect_times(vp, x) == 1)
                    elif k == BaseL:
                        constriants.append(vect_times(vm, x) == 1)
                    else:
                        constriants.append(vect_times(vp + vm, x) == 0)
                elif k == 'G':
                    if k == BaseH:
                        constriants.append(vect_times(vp2, x) == 1)
                    elif k == BaseL:
                        constriants.append(vect_times(vm2, x) == 1)
                    else:
                        constriants.append(vect_times(vp2 + vm2, x) == 0)
                else:
                    if k == BaseH:
                        constriants.append(vect_times(v, x) == 1)
                    elif k == BaseL:
                        constriants.append(vect_times(v, x) == -1)
                    else:
                        constriants.append(vect_times(v, x) == 0)

        prob = cp.Problem(obj, constriants)
        try:
            # prob.solve(solver=cp.XPRESS, verbose=Verbose)
            prob.solve(solver=cp.MOSEK, verbose=Verbose)
            # prob.solve(solver=cp.ECOS_BB, verbose=Verbose)
        except Exception as e:
            print(e)
            print('Fail to solve problem for Base %s Type %s!' % (CurrentBase, Type))
            return None
    if prob.status == 'optimal' and x.value is not None:
        return np.rint(x.value)
    else:
        print("Status: ", prob.status)
        print('Fail to solve problem for Base %s Type %s!' % (CurrentBase, Type))
        return None


def solve_base(
    CurrentBase,
    MyTag,
    MyBMNames,
    MyBMCurves,
    Verbose=False,
    IsDopant=True,
    TransFunc=None,
    TagLength=50,
):
    """Solve single type problem"""
    with fs.timer('Fit Base %s Type Dopants' % CurrentBase):
        x = cp.Variable(len(MyBMNames), integer=True)
        l = np.array([len(_) - 1 for _ in MyBMNames])
        if TransFunc is None:
            obj = cp.Minimize(cp.sum_squares(MixBMCurves.T @ x - MyTag))
        else:
            obj = cp.Minimize(cp.sum_squares(TransFunc(MixBMCurves.T @ x) - MyTag))
        # obj = cp.Minimize(cp.norm(MyBMCurves.T @ x - MyTag, 2))
        # constriants = [cp.sum(cp.multipy(l, x)) <= 29, x >= 0]
        if IsDopant:
            constriants = [vect_times(l, x) <= int(0.6 * TagLength), x >= 0, cp.sum(x) >= 1]
        else:
            constriants = [vect_times(l, x) == TagLength - 1, x >= 0, cp.sum(x) >= 1]
        prob = cp.Problem(obj, constriants)
        prob.solve(solver=cp.MOSEK, verbose=Verbose)
        # prob.solve(solver=cp.ECOS_BB, verbose=Verbose)
        # prob.solve(solver=cp.SCS, verbose=Verbose)
    if prob.status == 'optimal' and x.value is not None:
        return np.rint(x.value)
    else:
        print("Status: ", prob.status)
        raise Exception('Fail to solve problem!')


if __name__ == "__main__":

    Num = {
        '.': (0, 1),
        'e': 2,
        'eV': 0,
        'm': 1,
        'o': 1,
        'l': 1,
        'plotall': 0,
        'e0': 2,
        'homo': 1,
        'lumo': 1,
        'mix': 0,
        'ns': 0,
        'pdf': 0,
        'tpl': (1, 2),
        'notext': 0,
        'cf': 1,
        'curv': (1, 2),
        'freestyle': (3, -1),
    }
    Para_Dict = fs.get_para(sys.argv[1:], Num_Dict=Num, help=show_usage)

    MtxFile = Para_Dict['.'][0] if len(Para_Dict['.'])>0 else os.path.abspath('../data/Mtx_trained.fsz')
    NameFile = Para_Dict['m'][0] if 'm' in Para_Dict else '../data/HOMO.txt'
    Unit = 1.0 / fs.hartree if 'eV' in Para_Dict else 1.0
    if 'e' in Para_Dict:
        Emin, Emax = (Unit * float(_) for _ in Para_Dict['e'])
    else:
        Emin, Emax = (-5.8 * Unit, -5.1 * Unit)
    FigHead = Para_Dict['o'][0] if 'o' in Para_Dict else 'Fit-bms'
    MaxLen = int(Para_Dict['l'][0]) if 'l' in Para_Dict else 50
    PLOTALL_ON = 'plotall' in Para_Dict
    if 'e0' in Para_Dict:
        E0_st, E0_ed = (Unit * float(_) for _ in Para_Dict['e0'])
    else:
        E0_st, E0_ed = (-5.55 / fs.hartree, -5.45 / fs.hartree)
    MIX_ON = 'mix' in Para_Dict
    NOSINGLE_ON = 'ns' in Para_Dict
    if 'homo' in Para_Dict:
        E0_st = Unit * float(Para_Dict['homo'][0])
        E0_ed = Emax
    if 'lumo' in Para_Dict:
        E0_ed = Unit * float(Para_Dict['lumo'][0])
        E0_st = Emin
    TPL_ON = 'tpl' in Para_Dict
    if TPL_ON:
        TplDNAName = Para_Dict['tpl'][0]
        TplScalse = float(Para_Dict['tpl'][1]) if len(Para_Dict['tpl']) > 1 else 1
    if 'curv' in Para_Dict:
        CurveType = Para_Dict['curv'][0].upper()
        CurvePara = Para_Dict['curv'][1] if len(Para_Dict['curv']) > 1 else ''
    else:
        CurveType, CurvePara = '', ''
    FREESTYLE_ON = 'freestyle' in Para_Dict
    if FREESTYLE_ON:
        i = 0
        Paras = Para_Dict['freestyle']
        NPara = len(Paras)
        FreePara = []
        while i < NPara:
            if i + 2 >= NPara:
                print('WARNING: Discarding unknown parameters "%s"!' % (' '.join(Paras[i:])))
                break
            prtmp = [Paras[i].upper(), float(Paras[i + 1]), float(Paras[i + 2])]
            if i + 3 < NPara and not Paras[i + 3].isalpha():
                prtmp.append(float(Paras[i + 3]))
                i += 4
            else:
                i += 3
            FreePara.append(tuple(prtmp))
        if not FreePara:
            raise RuntimeError('No information for target curve!')
    # LIMIT_HOMO = 'homo' in Para_Dict
    # E_homo = Unit * float(Para_Dict['homo']) if LIMIT_HOMO else 0
    # LIMIT_LUMO = 'lumo' in Para_Dict
    # E_lumo = Unit * float(Para_Dict['lumo']) if LIMIT_LUMO else 0

    Data = fs.read_file(MtxFile)[1]
    RefDNAName = Data['name'][0]
    BMNames = np.array(Data['name'][1:])
    BMx = Data['x']
    idx_energy = np.where((BMx >= Emin * fs.hartree) & (BMx <= Emax * fs.hartree))[0]
    BMx = BMx[idx_energy]
    BMCurves = Data['data'][:, idx_energy]
    if NOSINGLE_ON:
        BMNames = np.append(BMNames, [b for b in 'AGCT'])
        BMCurves = np.append(BMCurves, [np.zeros_like(BMx) for _ in 'AGCT'], axis=0)
    if MaxLen != 50:
        LengthScale = 50 / MaxLen
        BMCurves *= LengthScale

    if FREESTYLE_ON:
        idx_fit = []
        EScale = Unit * fs.hartree
        for pr in FreePara:
            idx_fit.extend(list(np.where((BMx >= pr[1] * EScale) & (BMx <= pr[2] * EScale))[0]))
        idx_fit = np.sort(list(set(idx_fit)))
    else:
        idx_fit = np.where((BMx >= E0_st * fs.hartree - 0.05) & (BMx <= E0_ed * fs.hartree + 0.05))[0]
    BMx_fit = BMx[idx_fit]
    BMCurves_fit = BMCurves[:, idx_fit]

    DNAs = fs.read_DNA_info(NameFile)

    X = fs.energy_axis((Emin, Emax), MultiplyHartree=True)
    TextLines = []
    if RefDNAName == 'None':
        REF_ON = False
        IsDopant = False
        Dopants = 'AGCT'
        RefCurve = np.zeros_like(X)
        NameSuffix = ''
    else:
        REF_ON = True
        IsDopant = True
        RefDNA = fs.find_DNA(RefDNAName, DNAs)
        MainBase = RefDNA.Base[0]
        Dopants = fs.dopant_of_DNA(MainBase)
        TotalLen = len(RefDNA.Base)
        NameSuffix = '-' + RefDNAName
    if TPL_ON or REF_ON:
        # Read DNA curves
        CurvDB = Para_Dict['cf'][0] if 'cf' in Para_Dict else os.path.abspath('../data/Stored_Curves_5.8-5.1.fsz')
        with fs.timer('Read DNA curves'):
            if os.path.isfile(CurvDB):
                CurvData = fs.read_file(CurvDB)[1]
                DNACurves = [CurvData[d.Base] for d in DNAs]
            else:
                raise RuntimeError('ERROR: File of DNA curves not found!')
        RefCurve = DNACurves[RefDNA.Base]

    TransFunc = None
    if FREESTYLE_ON:
        PrNames = []
        TagCurve = np.ones_like(BMx) * SinglePeakH
        for pr in FreePara:
            idx = np.where((BMx >= pr[1] * EScale) & (BMx <= pr[2] * EScale))[0]
            if pr[0] == 'D':
                TagCurve[idx] = 0
                PrNames.append('d(%g,%g)' % (pr[1] * EScale, pr[2] * EScale))
            elif pr[0] == 'U':
                height = float(pr[3]) if len(pr) > 3 else 0.2
                TagCurve[idx] = height
                PrNames.append('u%g(%g,%g)' % (height, pr[1] * EScale, pr[2] * EScale))
        TagCurve_fit = TagCurve[idx_fit]
        TestName = '-'.join(PrNames)
        FRONT_ON = False
    else:
        TestName = '(%g,%g)' % (E0_st * fs.hartree, E0_ed * fs.hartree)
        FRONT_ON = ('homo' in Para_Dict or 'lumo' in Para_Dict)
        if FRONT_ON:
            TagCurve = np.ones_like(X) * SinglePeakH
            idx_E0 = np.where((X > E0_st * fs.hartree) & (X < E0_ed * fs.hartree))[0]
            TagCurve[idx_E0] = 0
            TagCurve -= RefCurve
            # TransFunc = cutoff_func
        else:
            idx_E0 = np.where((X < E0_st * fs.hartree) | (X > E0_ed * fs.hartree))[0]
            if REF_ON:
                TagCurve = -RefCurve
                TagCurve[idx_E0] = 0
            elif TPL_ON:
                TplDNA = fs.find_DNA(TplDNAName, DNAs)
                TplCurve = DNACurves[TplDNA.Base]
                if CurveType == 'D':
                    TagCurve = -TplCurve
                    TagCurve[idx_E0] = 0
                    TagCurve += TplCurve
                    TagCurve *= TplScalse
                elif CurveType == 'U':
                    CurvePara = float(CurvePara) if CurvePara else 0.2
                    TagCurve = np.ones_like(X) * CurvePara - TplCurve * TplScalse
                    TagCurve[idx_E0] = 0
                    TagCurve += TplCurve * TplScalse
                else:
                    TagCurve = TplCurve * TplScalse
            else:
                TagCurve = np.zeros_like(X)
                TagCurve[idx_E0] = 0.2
        TagCurve_fit = np.interp(BMx_fit, X, TagCurve)
        TagCurve = np.interp(BMx, X, TagCurve)
    RefCurve = np.interp(BMx, X, RefCurve)

    if REF_ON:
        Title = 'Delta DOS (Base: %s)' % RefDNA.ShortName
        FitCurves = [TagCurve + RefCurve, RefCurve]
        Labels = ['Tag_' + TestName, 'Ref']
        Styles = [{'linestyle': '-'}, {'linestyle': '-.'}]
    elif TPL_ON:
        Title = 'DOS (Template: %s)' % TplDNAName
        FitCurves = [TagCurve]
        Labels = ['Tag_' + TestName]
        Styles = [{'linestyle': '-'}]
    else:
        Title = ' DOS '
        FitCurves = [TagCurve]
        Labels = ['Tag_' + TestName]
        Styles = [{'linestyle': '-'}]

    SingleBaseBM = {b: np.zeros_like(BMx_fit) for b in 'AGCT'}
    if MIX_ON:
        FigHead = FigHead + '-mix' + NameSuffix
        idx_Single = {}
        idx_MixBM = []
        for i, name in enumerate(BMNames):
            if len(name) == 1:
                idx_Single[name] = i
                SingleBaseBM[name] = BMCurves_fit[i]
            else:
                idx_MixBM.append(i)
        MixBMNames = BMNames[idx_MixBM]
        MixBMCurves_full = BMCurves[idx_MixBM, :]
        MixBMCurves = BMCurves_fit[idx_MixBM, :]
        HeadTailBase = {b: np.zeros_like(idx_MixBM).astype(int) for b in 'AGCT'}
        for i, name in enumerate(MixBMNames):
            HeadTailBase[name[0]][i] += 1
            HeadTailBase[name[-1]][i] -= 1
            MixBMCurves[i] -= SingleBaseBM[name[0]]

        funcs = []
        res = []
        MixTypes = [0] + [a + b for a in Dopants for b in Dopants if a != b]
        if FRONT_ON:
            BaseAndName = []
            for i in range(1, 11, 2):
                MyTag = TagCurve_fit * i
                # print(cp.installed_solvers())
                for t in MixTypes:
                    if t == 0:
                        b = Dopants[0]
                    else:
                        b = t[0]
                    BaseAndName.append([b, '%s-%d' % (t, i)])
                    funcs.append((solve_mix, (b, MyTag, MixBMNames, MixBMCurves, HeadTailBase, t), {'IsDopant': IsDopant, 'TransFunc': TransFunc, 'TagLength': MaxLen}))
            with fs.timer('Calculate best bms'):
                Solutions = fs.run_functions(funcs, PoolSize=min(3, fs.Job.CPUAvail))
            for (b, n), s in zip(BaseAndName, Solutions):
                if s is not None:
                    if (s == 0).all():
                        print('Skip all-0 result of Base %s Type %s!' % (b, n))
                    else:
                        res.append([b, n, s])
        else:
            BaseAndName = []
            MyTag = TagCurve_fit
            # print(cp.installed_solvers())
            for t in MixTypes:
                if t == 0:
                    b = Dopants[0]
                else:
                    b = t[0]
                BaseAndName.append([b, t])
                funcs.append((solve_mix, (b, MyTag, MixBMNames, MixBMCurves, HeadTailBase, t), {'IsDopant': IsDopant, 'TransFunc': TransFunc, 'TagLength': MaxLen}))
            with fs.timer('Calculate best bms'):
                Solutions = fs.run_functions(funcs, PoolSize=min(3, fs.Job.CPUAvail))
            for (b, n), s in zip(BaseAndName, Solutions):
                if s is not None:
                    if (s == 0).all():
                        print('Skip all-0 result of Base %s Type %s!' % (b, n))
                    else:
                        res.append([b, n, s])

        RMSEs = []
        Texts = []
        FitCs = []
        for b, n, MyCoef in res:
            print('%s:' % n)
            Tag_full = TagCurve - BMCurves[idx_Single[b]]
            MyFitCurve = np.zeros_like(Tag_full)
            CoefText = []
            for i, c in enumerate(MyCoef):
                if c:
                    CoefText.append('%3d * %s' % (c, MixBMNames[i]))
                    print(CoefText[-1])
                    MyFitCurve += c * MixBMCurves_full[i]
            coefs_name_count = [(MixBMNames[i], int(c)) for i, c in enumerate(MyCoef) if c]
            FOUND_SEQ = False
            try:
                my_seq = fs.sec2seq(coefs_name_count)
                print('Seq:', my_seq)
                FOUND_SEQ = True
            except Exception as e:
                print('Failed to generate sequence:', e)
                continue
            Error = Tag_full[idx_fit] - MyFitCurve[idx_fit]
            RMSE = np.sqrt((Error**2).mean())
            MAE = np.abs(Error).mean()
            RMSEs.append(RMSE)
            print("RMSE = ", RMSE)  # 均方根误差RMSE
            print("MAE = ", MAE, '\n')  # 平均绝对误差MAE
            if FOUND_SEQ:
                Texts.append("%s: %s\nRMSE = %10.4g   MAE = %10.4g" % (n, my_seq, RMSE, MAE))
            else:
                Texts.append("%s: %s\nRMSE = %10.4g   MAE = %10.4g" % (n, ' '.join(CoefText), RMSE, MAE))
            FitCs.append(MyFitCurve + BMCurves[idx_Single[b]] + RefCurve)
        Min5_Idx = np.argsort(RMSEs)[:5]
        for idx in Min5_Idx:
            b, n, MyCoef = res[idx]
            TextLines.append(Texts[idx])
            FitCurves.append(FitCs[idx])
            Labels.append('Fit_mix_%s' % n)
            Styles.append({'linestyle': '--'})

    else:
        FigHead = FigHead + NameSuffix
        BMSelectDict = {b: [] for b in Dopants + 'X'}
        idx_Single = {}
        for i, name in enumerate(BMNames):
            if len(name) == 1:
                idx_Single[name] = i
                SingleBaseBM[name] = BMCurves_fit[i]
            elif name[0] == name[-1]:
                BMSelectDict[name[0]].append(i)
            else:
                BMSelectDict['X'].append(i)

        funcs = []
        for b in Dopants:
            MyTag = TagCurve_fit - SingleBaseBM[b]
            MyBMNames = BMNames[BMSelectDict[b]]
            MyBMCurves = BMCurves_fit[BMSelectDict[b], :] - SingleBaseBM[b]
            funcs.append((solve_base, (b, MyTag, MyBMNames, MyBMCurves), {'IsDopant': IsDopant, 'TransFunc': TransFunc, 'TagLength': MaxLen}))

        res = {b: r for b, r in zip(Dopants, fs.run_functions(funcs))}

        for b in Dopants:
            MyCoef = res[b]
            print('%s:' % b)
            Tag_full = TagCurve - BMCurves[idx_Single[b]]
            MyFitCurve = np.zeros_like(Tag_full)
            CoefText = []
            for i, c in enumerate(MyCoef):
                if c:
                    CoefText.append('%3d * %s' % (c, BMNames[BMSelectDict[b][i]]))
                    print(CoefText[-1])
                    MyFitCurve += c * (BMCurves[BMSelectDict[b][i], :] - BMCurves[idx_Single[b]])
            coefs_name_count = [(MixBMNames[i], int(c)) for i, c in enumerate(MyCoef) if c]
            FOUND_SEQ = False
            try:
                my_seq = fs.sec2seq(coefs_name_count)
                print('Seq:', my_seq)
                FOUND_SEQ = True
            except Exception as e:
                print('Failed to generate sequence:', e)
                continue
            Error = Tag_full[idx_fit] - MyFitCurve[idx_fit]
            RMSE = np.sqrt((Error**2).mean())
            MAE = np.abs(Error).mean()
            print("RMSE = ", RMSE)  # 均方根误差RMSE
            print("MAE = ", MAE, '\n')  # 平均绝对误差MAE
            TextLines.append("%s: %s\nRMSE = %10.4g   MAE = %10.4g" % (b, my_seq, RMSE, MAE))
            # TextLines.append("%s: %s\nRMSE = %10.4g   MAE = %10.4g" % (b, ' '.join(CoefText), RMSE, MAE))
            FitCurves.append(MyFitCurve + BMCurves[idx_Single[b]] + RefCurve)
            Labels.append('Fit_%s' % b)
            Styles.append({'linestyle': '--'})
    if 'notext' not in Para_Dict:
        Text = (fs.linear_interp(BMx, 0.1), fs.linear_interp(FitCurves, 0.75), '\n'.join(TextLines))
    else:
        Text = []
    Ext_Fig = 'pdf' if 'pdf' in Para_Dict else 'png'
    fs.fig_dos(FigHead + '_%s.%s' % (TestName, Ext_Fig), np.array([BMx] + FitCurves).T, **{'Labels': Labels, 'Text': Text, 'Title': Title, 'Styles': Styles})

    exit()
