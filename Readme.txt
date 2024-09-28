Package for Inverse Design of DNA Electronic Structures
version 1.0

Description
    This package contains the essensial data and codes for generate DNA according to target density of states (DOS). To do this, you need to select the target.


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
        │      dye_ABS.msgz
        │      dye_prop.msgz
        │      sdb.msgz
        │      train_idx.npy
        │      uvdb.msgz
        │      uvdb_trn.msgz
        │

    Functionality of files in each folder

        dye_ABS.msgz            Absorption spectra of dyes (in solution)
        dye_prop.msgz           Absorption spectra of dyes (in film)
        sdb.msgz                CD spectra of films
        train_idx.npy           The final training set used to train the forward prediction model
        uvdb.msgz               Transmission spectra of dyed films
        uvdb_trn.msgz           Transmission spectra of transparent films

        Data
            This folder contains the database used by other scripts.
            Files:
                dye_ABS.msgz            Absorption spectra of dyes (in solution)
                dye_prop.msgz           Absorption spectra of dyes (in film)
                sdb.msgz                CD spectra of films
                train_idx.npy           The final training set used to train the forward prediction model
                uvdb.msgz               Transmission spectra of dyed films
                uvdb_trn.msgz           Transmission spectra of transparent films
        fs
            This folder contains the functional tools used by other scripts.


