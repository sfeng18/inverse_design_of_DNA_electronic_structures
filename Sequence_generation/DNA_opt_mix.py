#!/usr/bin/env python
import numpy as np
import sys
import os
import fs
import signal
import cvxpy as cp
# from itertools import product

str_usage = """\
Usage: DNA_opt_mix.py (MtxFile) (Options)
Predict the sequence of DNA with target bit-like DOS automatically.
Options:
    -eb start end           Ranges of energy(eV) , default: -5.6 -5.3
    -e start end            Ranges of energy(eV) in plotting, default: -5.8 -5.1
    -o figure_head          Default: Opt-Tgt
    -edge len_ratio         Length of edges, Default: 0.25
    -cf CurveFile           File containing DNA DOS curves
    -pair No1 No2 ...       Target DNA pairs, default: 10 random pairs
    -reverse                Allow reserving DNA to generate sequence
    -ATfree                 Allow AT rich
    -long01                 Let all 0/1 to be full range low/high
    -h                      Show this help page
    """

SinglePeakH = 0.025
DOSHeight = 0.2


def show_usage():
    print(str_usage)
    sys.exit()


def cutoff_func(x, Threshold=SinglePeakH):
    return cp.minimum(x, Threshold)


def vect_times(a, b):
    return cp.sum(cp.multiply(a, b))


def solve_opt_tgt3(
    BMNames,
    BMCurves,
    TgtCurve,
    HeadTailBase,
    Len,
    Type=None,
    Verbose=False,
):
    """Solve mix type problem"""
    with fs.timer('Optimizing Type %s EndBases' % (Type)):
        x = cp.Variable(len(BMNames), integer=True)
        l = np.array([len(_) - 1 for _ in BMNames])
        LScale = Len / 50
        constraints = [x >= 0, cp.sum(x) >= 1]
        obj = cp.Minimize(cp.sum_squares(LScale * BMCurves.T @ x - TgtCurve))

        constraints.extend([vect_times(l, x) == Len - 1])
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
        if Type == 0:
            for k, v in HeadTailBase.items():
                constraints.append(vect_times(v, x) == 0)
                if k == 'C':
                    constraints.extend([vect_times(vp, x) <= 1, vect_times(vm, x) <= 1])
                elif k == 'G':
                    constraints.extend([vect_times(vp2, x) <= 1, vect_times(vm2, x) <= 1])
            if CG_CPLX:
                constraints.append(vect_times(vp + vp2, x) <= 1)
        else:
            assert isinstance(Type, str)
            assert len(Type) == 2
            BaseH, BaseL = Type
            # constraint_stat = {}
            if BaseH == BaseL:
                for k, v in HeadTailBase.items():
                    if k == 'C':
                        if k == BaseH:
                            constraints.append(vect_times(vp, x) == 1)
                            constraints.append(vect_times(vm, x) == 1)
                            # constraint_stat['C'] = (1, 1, 0)
                        else:
                            constraints.append(vect_times(vp + vm, x) == 0)
                            # constraint_stat['C'] = (0, 0, 0)
                    elif k == 'G':
                        if k == BaseH:
                            constraints.append(vect_times(vp2, x) == 1)
                            constraints.append(vect_times(vm2, x) == 1)
                            # constraint_stat['G'] = (1, 1, 0)
                        else:
                            constraints.append(vect_times(vp2 + vm2, x) == 0)
                            # constraint_stat['G'] = (0, 0, 0)
                    else:
                        constraints.append(vect_times(v, x) == 0)
                        # constraint_stat[k] = ('Inf', 'Inf', 0)
            else:
                for k, v in HeadTailBase.items():
                    if k == 'C':
                        if k == BaseH:
                            constraints.append(vect_times(vp, x) == 1)
                            constraints.append(vect_times(vm, x) == 0)
                            # constraint_stat['C'] = (1, 0, 1)
                        elif k == BaseL:
                            constraints.append(vect_times(vm, x) == 1)
                            constraints.append(vect_times(vp, x) == 0)
                            # constraint_stat['C'] = (0, 1, -1)
                        else:
                            constraints.append(vect_times(vp + vm, x) == 0)
                            # constraint_stat['C'] = (0, 0, 0)
                    elif k == 'G':
                        if k == BaseH:
                            constraints.append(vect_times(vp2, x) == 1)
                            constraints.append(vect_times(vm2, x) == 0)
                            # constraint_stat['G'] = (1, 0, 1)
                        elif k == BaseL:
                            constraints.append(vect_times(vm2, x) == 1)
                            constraints.append(vect_times(vp2, x) == 0)
                            # constraint_stat['G'] = (0, 1, -1)
                        else:
                            constraints.append(vect_times(vp2 + vm2, x) == 0)
                            # constraint_stat['G'] = (0, 0, 0)
                    else:
                        if k == BaseH:
                            constraints.append(vect_times(v, x) == 1)
                            # constraint_stat[k] = ('Inf', 'Inf', 1)
                        elif k == BaseL:
                            constraints.append(vect_times(v, x) == -1)
                            # constraint_stat[k] = ('Inf', 'Inf', -1)
                        else:
                            constraints.append(vect_times(v, x) == 0)
                            # constraint_stat[k] = ('Inf', 'Inf', 0)
        mosek_params = {
            # 'MSK_DPAR_MIO_REL_GAP_CONST': 0.05,
            # 'MSK_IPAR_MIO_MAX_NUM_ROOT_CUT_ROUNDS': 500,
            'MSK_IPAR_MIO_MAX_NUM_RELAXS': 100000,
        }
        prob = cp.Problem(obj, constraints)
        print(f'Type {Type} problem defined.')
        try:
            prob.solve(solver=cp.MOSEK, verbose=Verbose, mosek_params=mosek_params)
            print(f'Type {Type} problem solved.')
        except Exception as e:
            print(e)
            print('Fail to solve problem for Type %s!' % (Type))
            return None
    if prob.status == 'optimal' and x.value is not None:
        return np.rint(x.value)
    else:
        print("Status: ", prob.status)
        print('Fail to solve problem for Type %s!' % (Type))
        return None


def throw_error(e):
    print(e.__cause__)
    os.killpg(os.getpgid(os.getpid()), signal.SIGKILL)


if __name__ == "__main__":

    Num = {
        '.': (0, 1),
        'eb': 2,
        'e': 2,
        'o': 1,
        'edge': 1,
        'cf': 1,
        'pair': (2,-1),
        'tr': 1,
        'm': 1,
        'reverse': 0,
        'ATfree': 0,
        'long01': 0,
    }
    Para_Dict = fs.get_para(sys.argv[1:], Num_Dict=Num, help=show_usage)

    MtxFile = Para_Dict['.'][0] if len(Para_Dict['.'])>0 else os.path.abspath('../data/Mtx_trained.fsz')
    Unit = 1.0 / fs.hartree
    Emin, Emax = (float(_) for _ in Para_Dict['eb']) if 'eb' in Para_Dict else (-5.6, -5.3)
    Emin_plt, Emax_plt = (float(_) for _ in Para_Dict['e']) if 'e' in Para_Dict else (-5.8, -5.1)
    Emin_ht, Emax_ht = (Emin * Unit, Emax * Unit)
    FigHead = Para_Dict['o'][0] if 'o' in Para_Dict else 'Opt-Tgt'
    DNALen = 50

    Data = fs.read_file(MtxFile)[1]
    RefDNAName = Data['name'][0]
    BMNames = np.array(Data['name'][1:])
    BMx = Data['x']
    idx_energy = np.where((BMx >= Emin_plt) & (BMx <= Emax_plt))[0]
    BMx = BMx[idx_energy]
    BMCurves = Data['data'][:, idx_energy]

    idx_fit = np.where((BMx >= Emin) & (BMx <= Emax))[0][::10]
    BMx_fit = BMx[idx_fit]
    BMCurves_fit = BMCurves[:, idx_fit]
    BMCurves_fit[np.abs(BMCurves_fit) < 1e-8] = 1e-8

    NameSuffix = ''
    EndBases = 'ATCG'

    MixTypes = [a + b for a in EndBases for b in EndBases if a + b != 'TT']

    HeadTailBase = {b: np.zeros((len(BMNames),)).astype(int) for b in 'AGCT'}
    for i, name in enumerate(BMNames):
        HeadTailBase[name[0]][i] += 1
        HeadTailBase[name[-1]][i] -= 1

    DNAs_File = Para_Dict['m'][0] if 'm' in Para_Dict else '../data/HOMO.txt'
    DNAs = fs.read_DNA_info(DNAs_File)
    bmX = fs.energy_axis((Emin_plt * Unit, Emax_plt * Unit), MultiplyHartree=True)
    NbmX = len(bmX)
    print('bmX length: ', NbmX)

    # Read DNA curves
    CurvDB = Para_Dict['cf'][0] if 'cf' in Para_Dict else os.path.abspath('../data/Stored_Curves_5.8-5.1.fsz')
    with fs.timer('Read DNA curves'):
        if os.path.isfile(CurvDB):
            CurvData = fs.read_file(CurvDB)[1]
            DNACurves = [CurvData[d.Base] for d in DNAs]
        else:
            raise RuntimeError('ERROR: File of DNA curves not found!')

    test_DNAs = [(i, d) for i, d in enumerate(DNAs)]
    N_test = len(test_DNAs)
    print(f'{N_test} DNAs in test set.')

    if 'pair' in Para_Dict:
        input_DNA_Nos = [int(_) for _ in Para_Dict['pair']]
        pick_Nos = [(i1, i2) for i1, i2 in zip(input_DNA_Nos[:-1:2], input_DNA_Nos[1::2])]
    else:
        N_pairs = 10
        rand_Nos = np.random.randint(N_test, size=2 * N_pairs)
        pick_Nos = [(rand_Nos[2 * i], rand_Nos[2 * i + 1]) for i in range(N_pairs)]

    # pick_Nos = [(40, 291)]
    for i1, i2 in pick_Nos:
        test_No1, test_DNA1 = test_DNAs[i1]
        print(f'{i1} in test set, DNA1 No. {test_No1}\n{test_DNA1.Base}')
        OutputLines = [f'{i1} in test set, DNA No. {test_No1}\n{test_DNA1.Base}\n']
        my_curve1 = DNACurves[test_No1]
        my_curve_BMx1 = np.interp(BMx, bmX, my_curve1)
        test_No2, test_DNA2 = test_DNAs[i2]
        print(f'{i2} in test set, DNA2 No. {test_No2}\n{test_DNA2.Base}')
        OutputLines = [f'{i2} in test set, DNA No. {test_No2}\n{test_DNA2.Base}\n']
        my_curve2 = DNACurves[test_No2]
        my_curve_BMx2 = np.interp(BMx, bmX, my_curve2)
        for j in range(1):
            weight = j * 0.1
            TgtCurve_plot = weight * my_curve_BMx1 + (1 - weight) * my_curve_BMx2
            TgtCurve = TgtCurve_plot[idx_fit]
            print(f'weight: {weight}')
            OutputLines = [f'weight: {weight}\n']

            AllTypes = []
            funcs = []
            for t in MixTypes:
                AllTypes.append('%d-%s' % (DNALen, t))
                funcs.append((
                    solve_opt_tgt3,
                    (BMNames, BMCurves_fit, TgtCurve, HeadTailBase, DNALen, t),
                    {
                        # 'Verbose': True if t == 'CA' else False
                    },
                ))

            with fs.timer('Calculate best bms'):
                parallel_size = int(0.3 * fs.Job.CPUAvail) + 1
                task_num = len(funcs)
                if parallel_size > task_num:
                    parallel_size = task_num
                print(f'{task_num} tasks run in pool of size {parallel_size}.(CPUAvail={fs.Job.CPUAvail})')
                Err_Exit = {'error_callback': throw_error}
                Solutions = fs.run_functions(funcs, PoolSize=parallel_size, ApplyKwargs=Err_Exit)
            res = []
            for t, s in zip(AllTypes, Solutions):
                if s is not None:
                    if (s == 0).all():
                        print('Skip all-0 result of Type %s!' % (t))
                    else:
                        res.append([t, s])
            Errors = []
            Texts = []
            Out_Texts = []
            FitCs = []
            ExistCoef = []
            UniqueRes = []
            for t, MyCoef in res:
                if tuple(MyCoef) in ExistCoef:
                    continue
                ExistCoef.append(tuple(MyCoef))
                print('%s:' % t)
                MyCurve = np.zeros_like(BMx)
                CoefText = []
                coefs_name_count = []
                for i, c in enumerate(MyCoef):
                    if c:
                        CoefText.append('%3d * %s' % (c, BMNames[i]))
                        print(CoefText[-1])
                        MyCurve += c * BMCurves[i]
                        coefs_name_count.append((BMNames[i], int(c)))
                UniqueRes.append((t, coefs_name_count))
                MyCurve_TagRange = MyCurve[idx_fit]
                My_Error = np.sum(np.abs(MyCurve_TagRange - TgtCurve))
                Errors.append(My_Error)
                Texts.append(f"{t}: {My_Error}" + '\n' + ' '.join(CoefText))
                Out_Texts.append(f"{t}: {My_Error}" + "\n" + '\n'.join(CoefText))
                FitCs.append(MyCurve)
            Good_Idx = np.argsort(Errors)
            OutputLines.extend([Out_Texts[_] + '\n' for _ in Good_Idx])
            Best5_Idx = Good_Idx[:5]
            PlotCurves = [TgtCurve_plot]
            Labels = ['Tgt']
            TextLines = [f'Tgt:{weight}']
            FOUND_BEST = False
            for idx in Good_Idx:
                t, coefs_name_count = UniqueRes[idx]
                try:
                    my_seq = fs.sec2seq(coefs_name_count)
                    print('Best type:', t)
                    print('Seq:', my_seq)
                    OutputLines.append('Best type: %s\nSeq: %s\n' % (t, my_seq))
                    FOUND_BEST = True
                except Exception as e:
                    print('Fail at type:', t)
                    print(e)
                    continue
                if FOUND_BEST:
                    TextLines.append(f"{t}: {Errors[idx]}" + '\n' + my_seq)
                    PlotCurves.append(FitCs[idx])
                    Labels.append(t)
                    break

            MyName = 'MixNo%d-%d-%g_Best' % (test_No1, test_No2, weight)
            Title = MyName
            if 'notext' not in Para_Dict:
                Text = (fs.linear_interp(BMx, 0.1), fs.linear_interp(PlotCurves, 0.75), '\n'.join(TextLines))
            else:
                Text = []
            Ext_Fig = 'pdf' if 'pdf' in Para_Dict else 'png'
            fs.fig_shadow(FigHead + '-%dbp_%s.%s' % (DNALen, MyName, Ext_Fig), np.array([BMx] + PlotCurves).T, **{'Labels': Labels, 'Text': Text, 'Title': Title})

            output_name = f"Opt-Res_{MyName}.txt"
            with open(output_name, 'at') as f:
                f.writelines(OutputLines)

    exit()
