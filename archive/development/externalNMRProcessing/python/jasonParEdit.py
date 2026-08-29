#! /usr/bin/env python

# ------------------------------------------------------------------------------- 
# --
# -- JEOL Ltd.
# -- 1-2 Musashino 3-Chome
# -- Akishima Tokyo 196-8558 Japan 
# -- Copyright 2024 
# -- 
#
# --++--------------------------------------------------------------------------- 
# -- 
# -- ModuleName : JASON
# -- ModuleType : Example external command script for JASON 
# -- Purpose : desktop software 
# -- Date : August 2024
# -- Author : Peter Kiraly
# -- Language : Python
# -- 
# --##---------------------------------------------------------------------------
#
# This script enables editing parameters ('SpectrometerFrequencies', 'SW', 'SpectrumRef') in the jason_parameters/specInfo group
#
# How to use:
#  In JASON, add "External command" to the processing list
# 
#  Set for external command parameters:
#
#   Set "cmd" to "python" or to the full path to python.exe 
#     (e.g. C:\Program Files\Python311\python.exe)
#
#   Set "arguments" to the path of the script, specifying the use of a temporary 
#   file
#     (e.g. "C:\Users\<username>\.jason\externalNMRProcessing\python\jasonParEdit.py -f $TMPFILE")
#     use the -p flag to set the name of the parameter to edit [SpectrometerFrequencies or SW or SpectrumRef] (default SpectrometerFrequencies)
#     use the -d flag to set the index of dimension when applicable [0-7 for 1st-8th dimensions] (default 0)
#     use the -v flag to set the new value [only numeric values supported] (default 500.0)
#     
#  Set "Data file" to "Spectrum as JJH5"
#
# Press "Apply"
#

import h5py
from argparse import ArgumentParser
import numpy as np

# Process commandline options
parser = ArgumentParser()
# temporary file name
parser.add_argument("-f", "--filename", action="store")
# parameter name
parser.add_argument("-p", "--par", action="store", default='SpectrometerFrequencies')
# dimension
parser.add_argument("-d", "--dim", action="store", type=int, default=0)
# new value
parser.add_argument("-v", "--vNum", action="store", type=float, default=500.0)

args = parser.parse_args()

# Open a file handle for the JASON datafile
f = h5py.File(args.filename, "r+")
print('Opening dataset: ', args.filename)

# check for supported parameter values
supportedPars = ['SpectrometerFrequencies', 'SW', 'SpectrumRef'] 
if args.par in supportedPars:
    # range check for dimension
    if args.dim>=0 and args.dim<8:  
        # checks if attribute exists in the file
        keys = f['JasonDocument/SpecInfo'].attrs.keys()
        if args.par in keys:
            # Modify parameter
            parameter = f['JasonDocument/SpecInfo'].attrs[args.par][:]
            parameter[args.dim] = args.vNum
            f['JasonDocument/SpecInfo'].attrs.modify(args.par, parameter)
            #report
            print('Parameter: ' + args.par +' updated')
        else:
            print('Error: parameter: ' + args.par +' does not exist.')            
    else:
        print('Error: dimension index must be in range of 0-7.')
else:
    print('Error: parameter: ' + args.par +' is not supported by this script.')         

f.close()
