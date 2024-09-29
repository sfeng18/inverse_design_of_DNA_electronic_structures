Package for Inverse Design of DNA Electronic Structures
version 1.0

Description
    This package contains the essensial data and codes for generate DNA according to target density of states (DOS). To do this, you need to select the target from example DNAs (1372 in total).


Requirements
    To run the codes, you need a Python 3 environment with the following packages installed:
    - numpy
    - scipy
    - matplotlib
    - ujson
    - zlib
    - chardet
    - msgpack-python
    - psutil
    - cvxpy
    - mosek

Installation
    The main python scripts can be run directly in Python 3 environment, without installation.

Documenation
    Content
        │  Readme.txt
        │  requirements.txt
        │  DNA_opt_tgt_no.py
        │  LICENSE
        │
        ├─Data
        │      HOMO.txt
        │      Mtx_trained.fsz
        │      Stored_Curves_5.8-5.1.fsz
        │

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
    
    How to run DNA_opt_tgt_no.py:
        Usage: DNA_opt_tgt_no.py (MtxFile) (Options)
        Predict the sequence of DNA to mimic the DOS of target DNA.
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


