#! /usr/bin/env python

# ------------------------------------------------------------------------------- 
# --
# -- JEOL Ltd.
# -- 1-2 Musashino 3-Chome
# -- Akishima Tokyo 196-8558 Japan 
# -- Copyright 2022 
# -- 
#
# --++--------------------------------------------------------------------------- 
# -- 
# -- ModuleName : noise_reduction.py
# -- ModuleType : Example external command script for JASON 
# -- Purpose : Reference deconvolution via external data processing in JASON 
# -- Date : August 2023
# -- Author : Iain J. Day / Peter Kiraly / Hugo Rocha
# -- Language : Python
# -- 
# --##---------------------------------------------------------------------------
#
# A script apply reference deconvolution to a pseudo2D spectrum. Based on the 1D MATLAB
# script written by Hugo Rocha and translated to Python by Iain J. Day
#
# How to use:
#  In Jason, add "External command" to the processing list
# 
#  Set for external command parameters:
#
#   Set "cmd" to "python" or to the full path to python.exe 
#     (e.g. C:\Program Files\Python311\python.exe)
#
#   Set "arguments" to the path of the script, specifying the use of a temporary 
#   file
#     (e.g. "C:\Users\<username>\.jason\externalNMRProcessing\python\reference_deconvolution_pseudo2D.py -f $TMPFILE")
#     use the -l flag to set the Lorentzian line width (default 0.0 Hz)
#     use the -g flag to set the Gaussian line width (default inf Hz)
#     use the -r flag to set the reference signal position in PPM (default 0.0 ppm)
#     use the -w flag to set the width of the reference signal in Hz (default 100.0 Hz)
#     use the -s flag to take into account satellites with the following two properties (default is false, no value given just use as -s)
#     use the -j flag to set the value of the heteronuclear J coupling for satellite of reference signal (default 6.6 Hz TMS)
#     use the -a flag to set the abundance of the heteronucleus for satellite of reference signal (default 4.67 Hz for 29Si)
#
#  Set "Data file" to "Spectrum as JJH5"
#
# Press "Apply"
#
# Note, use of peak reference processing function in JASON can help right before the external command to fix the reference position
#
# References:  Morris, GA, Barjat H, Horne TJ. 1997. Reference Deconvolution 
#   Methods. Progress in Nuclear Magnetic Resonance Spectroscopy. 31:197-257.
#
#   Castañar, L., Poggetto, G.D., Colbourne, A.A., Morris, G.A. and Nilsson, M., 
#   2018. The GNAT: A new tool for processing NMR data. Magnetic Resonance in 
#   Chemistry, 56(6), pp.546-558.

import h5py
from argparse import ArgumentParser
import numpy as np
from scipy.fft import fft, ifft, fftshift

# Process commandline options
parser = ArgumentParser()
# temporary file name
parser.add_argument("-f", "--filename", action="store")
# Lorentzian linewidth
parser.add_argument("-l", "--lhz", action="store", type=float, default=0.0)
# Gaussian linewidth
parser.add_argument("-g", "--ghz", action="store", type=float, default=np.inf)
# reference signal position in PPM
parser.add_argument("-r", "--refpos", action="store", type=float, default=0.0)
# reference region width to use in Hz
parser.add_argument("-w", "--width", action="store", type=float, default=100.0)
# use satellite doublet for reference signal
parser.add_argument("-s", "--sat", action="store_true")
# heteronuclear J-coupling of reference signal 
parser.add_argument("-j", "--jxy", action="store", type=float, default=6.6)
# abundance of reference satellite in % 
parser.add_argument("-a", "--abundance", action="store", type=float, default=4.67)

args = parser.parse_args()

# Open a file handle for the Jason datafile
f = h5py.File(args.filename, "r+")
print('Opening dataset: ', args.filename)

nptsF2 = f['JasonDocument'].attrs['Length'][0]
nptsF1 = f['JasonDocument'].attrs['Length'][1]

# Read in the spectrum
dsetreAll = f['JasonDocument/DataPoints/0'][()]

print('Dataset size: ', nptsF2, ' * ', nptsF1)
# Zero final output spectra
SPECTRA_ALL = np.array([], dtype=np.complex_)

for i in range(nptsF1):
    npts = nptsF2
#    print('i: ', i)
    dsetre = dsetreAll[i,:].copy()
    
    wholefid = ifft(np.flipud(fftshift(dsetre, 0)), npts)
    wholefid = wholefid[:npts//2]

# Zero SPECTRA
    SPECTRA = np.zeros(npts)

# Create the time shifted FID for the reference peak
    sw = f['JasonDocument/SpecInfo'].attrs['SW'][0]
    at = npts / sw
    sfrq = f['JasonDocument/SpecInfo'].attrs['SpectrometerFrequencies'][0]
    sw = sw / sfrq
    sref = f['JasonDocument/SpecInfo'].attrs['SpectrumRef'][0]
    x_offset = sref / sfrq
    sp = x_offset - sw / 2.0

    RDcentre = args.refpos
    widthP = args.width / sfrq
    speclim = np.array((RDcentre - 0.5* widthP, RDcentre + 0.5* widthP))

# Convert the integreation limits from ppm to number of points
    speclim = np.array(npts * (speclim - sp) / sw, dtype=int)

    temp = npts - np.round(npts * (RDcentre - sp) / sw)
    RDcentre = ((temp * sw) / npts) + sp

    if speclim[0] < 1:
        speclim[0] = 1

    if speclim[1] > npts:
        speclim[1] = npts

    speclim = npts - speclim
    exprefspec = dsetre.copy()

    exprefspec[:speclim[1]] = 0.0
    exprefspec[speclim[0]:] = 0.0

    expreffid = ifft(np.flipud(fftshift(exprefspec, 0)), npts)
    expreffid = expreffid[:npts//2]

# Create the perfect reference peak
    omega = 0.5 * sw * sfrq - (RDcentre - sp) * sfrq

    t = np.linspace(0, 0.5*at, npts//2)
    reffid = np.exp(1j * 2.0 * np.pi * omega * t)

    if args.sat:
        reffid = reffid + args.abundance/200.0*np.exp(1j * 2.0 * np.pi * (omega + args.jxy/2.0) * t) + args.abundance/200.0*np.exp(1j * 2.0 * np.pi * (omega - args.jxy/2.0) * t)

    reffid = reffid * np.exp(-t * np.pi * args.lhz - (t / args.ghz)**2)

    

# Create the correction fid
    corrfid = expreffid / reffid
    corrfid = corrfid / corrfid[0]
    endfid = wholefid / corrfid
    endfid[0] = 0.5 * endfid[0]

    SPECTRA = np.flipud(fftshift(fft(endfid, npts), 0))
        
# Store results per traces
    SPECTRA_ALL = np.append(SPECTRA_ALL, SPECTRA)

# Write out the changes to the original file

del f['/JasonDocument/DataPoints/0']
del f['/JasonDocument/DataPoints/1']
f.create_dataset('/JasonDocument/DataPoints/0', data=SPECTRA_ALL.real)
f.create_dataset('/JasonDocument/DataPoints/1', data=SPECTRA_ALL.imag)
print('Dataset changed')
f.close()
