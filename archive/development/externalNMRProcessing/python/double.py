#!python3

#A simple Python script to double the processed spectrum by appending inverted spectrum at the end
#How to use:
#In Jason, add "External command" to the processing list at the end of the existing list
#Set for external command parameters:
#Set "cmd" to "python" or to the full path to python.exe (e.g. C:\Users\<username>\AppData\Local\Programs\Python\Python38-32\python.exe)
#Set "arguments" to the path of the script (e.g. "C:\Users\<username>\Desktop\double.py")
#Set "Data file" to "Spectrum as JJH5"
#Press "Apply"

import sys
import h5py
import numpy as np

with h5py.File(sys.argv[1], "r+") as hf:
    print("Opening dataset:", sys.argv[1])
    dsetre = hf['/JasonDocument/DataPoints/0']
    dsetim = hf['/JasonDocument/DataPoints/1']
    print("dataset initial size:", len(dsetre))
    dsetre=np.append(dsetre,np.negative(dsetre));
    dsetim=np.append(dsetim,np.negative(dsetim));
    length=hf['/JasonDocument'].attrs['Length']
    length[0]=length[0]*2
    hf['/JasonDocument'].attrs.modify('Length',length)
    print("dataset new size:", len(dsetre))
    
    del hf['/JasonDocument/DataPoints/0']
    hf.create_dataset('/JasonDocument/DataPoints/0', data=dsetre)
    del hf['/JasonDocument/DataPoints/1']
    hf.create_dataset('/JasonDocument/DataPoints/1', data=dsetim)
    
    print("dataset changed")
    hf.close()
    