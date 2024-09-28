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
        │      HOMO.txt
        │      Mtx_trained.fsz
        │      Stored_Curvs_5.8-5.1.z01
        │      Stored_Curvs_5.8-5.1.zip
        │

    Functionality of files in each folder

        requirements.txt        Installation requirement
        DNA_opt_tgt_no.py       Script for DNA generation

        Data
            This folder contains the database used by other scripts.
            Files:
                HOMO.txt                            List of DNAs
                Mtx_trained.fsz                     Correlation between sequence and DOS 
                Stored_Curvs_5.8-5.1.z01            DOS of DNAs part II
                Stored_Curvs_5.8-5.1.zip            DOS of DNAs part I, need to unzip before use
        fs
            This folder contains the functional tools used by other scripts.


