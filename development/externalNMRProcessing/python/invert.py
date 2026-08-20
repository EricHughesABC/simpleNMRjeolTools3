#!python3

#A simple Python script to invert spectrum point-by-point
#How to use:
#In Jason, add "External command" to the processing list
#Set for external command parameters:
#Set "cmd" to "python" or to the full path to python.exe (e.g. C:\Users\<username>\AppData\Local\Programs\Python\Python38-32\python.exe)
#Set "arguments" to the path of the script (e.g. "C:\Users\<username>\Desktop\invert.py")
#Set "Data file" to "Spectrum as JJH5"

import sys
import h5py
import numpy as np

with h5py.File(sys.argv[1], "r+") as hf:
    print("Opening dataset:", sys.argv[1])
    dsetre = hf['/JasonDocument/DataPoints/0']
    dsetim = hf['/JasonDocument/DataPoints/1']
    print("dataset size:", len(dsetre))
    dsetre[...]=np.negative(dsetre)
    dsetim[...]=np.negative(dsetim)
    print("dataset changed")
    hf.close()
    