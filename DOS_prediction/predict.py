#!/usr/bin/env python
import numpy as np
import sys
import os
import fs

str_usage = """\
Usage: predict.py (MtxFile) (Options)
Predict the DOS of DNA.
Options:
    -eb start end           Ranges of energy(eV) , default: -5.6 -5.3
    -e start end            Ranges of energy(eV) in plotting, default: -5.8 -5.1
    -m name_file            Name file convering shortname into full name, default: ./data/HOMO.txt
    -o figure_head          Default: Pred_DOS
    -cf CurveFile           File containing DNA DOS curves
    -s seq1 (seq2...)       Sequence of target DNA
    -tgt No1 (No2 ...)      Target DNA numbers, default: 10 random numbers
    -h                      Show this help page
    """

SinglePeakH = 0.025
DOSHeight = 0.2


def show_usage():
    print(str_usage)
    sys.exit()


if __name__ == "__main__":

    Num = {
        '.': (0, 1),
        'eb': 2,
        'e': 2,
        'o': 1,
        'cf': 1,
        'm': 1,
        's': (1, -1),
        'tgt': (1, -1),
    }
    Para_Dict = fs.get_para(sys.argv[1:], Num_Dict=Num, help=show_usage)

    MtxFile = Para_Dict['.'][0] if len(Para_Dict['.']) > 0 else os.path.abspath('../data/Mtx_trained.fsz')
    Unit = 1.0 / fs.hartree
    Emin, Emax = (float(_) for _ in Para_Dict['eb']) if 'eb' in Para_Dict else (-5.6, -5.3)
    Emin_plt, Emax_plt = (float(_) for _ in Para_Dict['e']) if 'e' in Para_Dict else (-5.8, -5.1)
    Emin_ht, Emax_ht = (Emin * Unit, Emax * Unit)
    FigHead = Para_Dict['o'][0] if 'o' in Para_Dict else 'Pred_DOS'
    seqs = Para_Dict['s'] if 's' in Para_Dict else []

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

    CurvDB = Para_Dict['cf'][0] if 'cf' in Para_Dict else os.path.abspath('../data/Stored_Curves_5.8-5.1.fsz')
    DNAs_File = Para_Dict['m'][0] if 'm' in Para_Dict else os.path.abspath('../data/HOMO.txt')
    DNAs = fs.read_DNA_info(DNAs_File)
    bmX = fs.energy_axis((Emin_plt * Unit, Emax_plt * Unit), MultiplyHartree=True)
    NbmX = len(bmX)
    print('bmX length: ', NbmX)

    test_DNAs = [(i, d) for i, d in enumerate(DNAs)]
    N_test = len(test_DNAs)
    print(f'{N_test} DNAs in total.')

    with fs.timer('Read DNA curves'):
        if os.path.isfile(CurvDB):
            CurvData = fs.read_file(CurvDB)[1]
        else:
            CurvData = {}

    if not seqs:
        if 'tgt' in Para_Dict:
            pick_DNA_Nos = [int(_) for _ in Para_Dict['tgt']]
            pick_DNAs = [(f'pick-{i}', j, DNAs[j]) for i, j in enumerate(pick_DNA_Nos)]
        else:
            pick_Nos = np.random.randint(N_test, size=10)
            # pick_Nos = [np.random.randint(N_test)]
            pick_DNAs = [(f'pick-{i}', *test_DNAs[j]) for i, j in enumerate(pick_Nos)]
        seqs = [d.Base for i, j, d in pick_DNAs]

    BMNames_tuple = tuple(BMNames)
    plot_funcs = []
    for my_seq in seqs:
        bms = fs.cut_bms(my_seq)
        Unqbms = set()
        for bm in set(bms):
            if bm in BMNames:
                Unqbms.add(bm)
            else:
                bm_dis = [(_, fs.seq_dist(bm, _)) for _ in BMNames]
                closest_bm_dis = min(bm_dis, key=lambda x: x[1])
                anti_bm_dis = [(_, fs.seq_dist(fs.anti_seq(bm), _)) for _ in BMNames]
                closest_anti_bm_dis = min(anti_bm_dis, key=lambda x: x[1])
                if closest_bm_dis[1] < closest_anti_bm_dis[1]:
                    Unqbms.add(closest_bm_dis[0])
                else:
                    Unqbms.add(closest_anti_bm_dis[0])
                print(f'Warning: {bm} not found, use {closest_bm_dis[0]} instead.')
        if not Unqbms:
            print(f'Warning: No unique bms found for sequence {my_seq}.')
            continue
        my_curve_BMx = np.zeros_like(BMx)
        coefs, idx = [], []
        for bm in Unqbms:
            idx.append(BMNames_tuple.index(bm))
            MyCoef = bms.count(bm) if len(bm) > 1 else -bms.count(bm)
            coefs.append(MyCoef)
        my_curve_BMx = np.array(coefs).reshape((1, -1)).dot(BMCurves[idx, :]).reshape(-1)
        LengthScale = 50 / len(my_seq)
        my_curve_BMx *= LengthScale

        if len(CurvData) and my_seq in CurvData:
            plot_xy = np.array([bmX, my_curve_BMx, CurvData[my_seq]]).T
            Labels = ['Predicted', 'Calculated']
        else:
            plot_xy = np.array([bmX, my_curve_BMx]).T
            Labels = ['Predicted']
        fig_name = f"{FigHead}_{my_seq}.png"
        plot_funcs.append((
            fs.fig_XYMtx,
            (fig_name, plot_xy),
            {
                'Labels': Labels,
                'Title': my_seq,
                'AxisTitle': ('Energy (eV)', 'DOS')
            },
        ))

    fs.run_functions(plot_funcs)

    exit()
