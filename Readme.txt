Package for Inverse Design of DNA Electronic Structures
version 1.0

Description
    This package contains the essensial data and codes for generate DNA according to target density of states (DOS). To do this, you need to select the target from example DNAs (1372 in total).


Requirements
    To run the codes, you need a Python 3 environment (>=3.11) with the following packages installed:
    - numpy >= 1.26.4
    - scipy >= 1.12.0
    - matplotlib >= 3.8.4
    - ujson >= 5.4.0
    - zlib >= 1.2.13
    - chardet >= 4.0.0
    - msgpack-python >= 1.0.3
    - psutil >= 5.9.0
    - cvxpy >= 1.5.1
    - mosek >= 10.1.27

Installation
    The main python scripts can be run directly in Python 3 environment, without installation.

Documenation
    Content
        │  LICENSE
        │  Readme.txt
        │  requirements.txt
        │
        ├─data
        │      HOMO.txt
        │      Mtx_trained.fsz
        │      Stored_Curves_5.8-5.1.fsz
        │
        ├─DOS_prediction
        │      predict.py
        │
        ├─fs
        │      basis.py
        │      core.py
        │      data.py
        │      fileIO.py
        │      molecule.py
        │      molslc.py
        │      parallel.py
        │      plot.py
        │      requirements.txt
        │
        └─Sequence_generation
                DNA_opt_complex.py
                DNA_opt_mix.py
                DNA_opt_tgt_no.py
        
    Functionality of files in each folder

        requirements.txt        Installation requirement
        DNA_opt_tgt_no.py       Script for DNA generation

        Data
            This folder contains the database used by other scripts.
            Files:
                HOMO.txt                            List of example DNAs
                Mtx_trained.fsz                     Correlation between sequence and DOS 
                Stored_Curves_5.8-5.1.fsz            DOS of example DNAs
        fs
            This folder contains the functional tools used by other scripts.
        DOS_prediction
            This folder contains scripts for prediction of DNA's electronic density of states.
        Sequence_generation
            This folder contains scripts for generating the sequence of DNA for a target DOS.
    
    How to run the scripts?
        To predict the electronic density of states (DOS) of a target DNA, you can use the script predict.py in the folder DOS_prediction. Two possible ways to give the target DNAs are:
        1. Give the sequence of the target DNA using the option -s.
        2. Give the number of the target DNA using the option -tgt.
        The script will generate the DOS of the target DNA and plot the predicted DOS curve.

        Usage:
            predict.py (MtxFile) (Options)
            Options:
                -eb start end           Ranges of energy(eV) , default: -5.6 -5.3
                -e start end            Ranges of energy(eV) in plotting, default: -5.8 -5.1
                -m name_file            Name file convering shortname into full name, default: ./data/HOMO.txt
                -o figure_head          Default: Pred_DOS
                -cf CurveFile           File containing DNA DOS curves
                -s seq1 (seq2...)       Sequence of target DNA
                -tgt No1 (No2 ...)      Target DNA numbers, default: 10 random numbers
                -h                      Show this help page
        Example:
            # First, you need to change the directory to DOS_prediction
            >cd DOS_prediction
            # Then, you can run the script with the target DNA sequence
            >predict.py -s AGGCTAGATACGATCCT
            # Or, you can run the script with the target DNA number
            >predict.py -tgt 82 934

        To generate the sequence of DNA for a target DOS, you can use the scripts in the folder Sequence_generation. There are three scripts for different types of target DOS:
        1. DNA_opt_tgt_no.py: Mimic the DOS of a target DNA.
        2. DNA_opt_mix.py: Target DOS is a mixture of two DNAs.
        3. DNA_opt_complex.py: Target DOS is a complex bit-like curve.

        Usage:
            DNA_opt_tgt_no.py (MtxFile) (Options)
            Options:
                -eb start end           Ranges of energy(eV) , default: -5.6 -5.3
                -e start end            Ranges of energy(eV) in plotting, default: -5.8 -5.1
                -m name_file            Name file convering shortname into full name, default: HOMO.txt
                -o figure_head          Default: Opt-Tgt
                -cf CurveFile           File containing DNA DOS curves
                -l DNA_length           Default: 50
                -tgt No1 (No2 ...)      Target DNA numbers, default: 10 random numbers
                -ATend                  Only allow start by 'A' end by 'T'
                -h                      Show this help page
            
            DNA_opt_mix.py (MtxFile) (Options)
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

            DNA_opt_complex.py (Mtx_File) (Options)
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

        Examples:
        # First, you need to change the directory to Sequence_generation
        >cd Sequence_generation
        # Then, you can run the script with the target DNA number
        >DNA_opt_tgt_no.py -tgt 582 364
        # Or, you can run the script with the target DNA pairs
        >DNA_opt_mix.py -pair 40 291
        # Or, you can run the script with the target DNA complex curve, like a 010 3-bit curve
        >DNA_opt_complex.py -o Fit -eV -mix -ns -freestyle d -5.4 -5.3 u -5.5 -5.4 d -5.6 -5.5
        # The script will generate a DNA with length 50 by default. If you want to generate a DNA with a different length, you can use the -l option. For example:
        >DNA_opt_complex.py -o Fit-20bp -l 20 -eV -mix -ns -freestyle d -5.6 -5.55 u -5.5 -5.4
        # A much more complex curve can be designed by using the -freestyle option. For example:
        >DNA_opt_complex.py -o Fit-b6 -notext -eV -mix -ns -freestyle u -5.65 -5.583 u -5.583 -5.516 u -5.516 -5.45 d -5.45 -5.383 u -5.383 -5.316 d -5.316 -5.25




