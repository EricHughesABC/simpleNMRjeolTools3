#! /usr/bin/env python

# A script to read HDF5 files from Jason and plot a spectrum with Matplotlib.
# A t1-noise reduction algorithm is then applied, 
# taken from J. Biomol. NMR, 2 (1992) 485-494
# 
# NOTE: data must be baseline corrected prior to use
#
# IJD, 1/2/22

import sys
import h5py

import numpy as np


# Set some parameters (these could be read from a config file . . .)

eta = 3.0           # Real-peak threshold (peaks above this level are determined to be real)
itr_stop = 0.01     # Tolerance to stop the ANI iterations (not normally changed)
T = 10.0            # Smoothing threshold (low values give stronger t1 noise reduction)
smooth_itr = 10     # Number of smoothing iterations applied


# Open a file handle for the Jason datafile

f = h5py.File(sys.argv[1], "r+")
print('Opening dataset: ', sys.argv[1])


# Get the real part of the spectrum
# Need to swap the number of points around to get the right shape compared to
# how JASON stores the parameters

real_spec = f['JasonDocument/DataPoints/0'][()]

npts = real_spec.shape

t1red_spec = np.zeros(npts)


# Get some parameters associated with the spectrum from the HDF5 file

sw = f['JasonDocument/SpecInfo'].attrs['SW']
spec_freq = f['JasonDocument/SpecInfo'].attrs['SpectrometerFrequencies']
spec_ref = f['JasonDocument/SpecInfo'].attrs['SpectrumRef']


# Begin the ANI (Average NoIse) routine

print('Determining t1 noise profile . . .')

N = np.zeros(npts[1])

for i in range(npts[1]):
    # Loop over each t1 column of the spectrum
    # Find the locations of the zero crossings

    zero_crossings = np.where(np.diff(np.signbit(real_spec[:, i])))[0]


    # The peaks are the (abs) maxima between each zero crossing

    A_idx = np.zeros(len(zero_crossings)-1, dtype=int)

    for j in range(A_idx.shape[0]-1):
        A_idx[j] = j + np.abs(real_spec[zero_crossings[j:j+1], i]).argmax()


    # Build an array containing the peaks

    A = real_spec[A_idx, i]


    # Calculate the noise, and remove any real peaks, defined as peaks greater than
    # eta times the noise level

    while True:
        noise = np.abs(A).sum() / A.size

        A = A[A < eta * noise]

        new_noise = np.abs(A).sum() / A.size

        if noise - new_noise < itr_stop:
            N[i] = new_noise
            break


# Finally, calculate the error for each t1 point, and some derived quantities
# The elements of the smoothing array are set to a maximum value of 1.0

error = eta * N
L = T * error.min()
S = L / error
S[S>1.0] = 1.0


# Begin the RT1 (Reduce T1 noise) routine

print('Applying smoothing . . .')

for i in range(npts[0]):
    # Loop over each row in the spectrum, taking a copy of that row

    P = real_spec[i, :].copy()


    # Reduce the intensity of every point by the error

    P[P<error] = 0.0
    P = np.sign(P) * (np.abs(P) - error)


    # Apply t1 noise smoothing, checking that no point was adjusted by more than the error

    for j in range(smooth_itr):
        P = ((1.0 - S[1:-1]) * S[:-2] * P[:-2] + 2.0 * S[1:-1] * P[1:-1] + \
        (1.0 - S[1:-1]) * S[2:] * P[2:]) / \
        ((1.0 - S[1:-1]) * S[:-2] + 2.0 * S[1:-1] + (1.0 - S[1:-1]) * S[2:])

        P = np.append(np.append(0.0, P), 0.0)

        P[P > real_spec[i, :] + error] = real_spec[i, P > real_spec[i, :] + error] + error[P > real_spec[i, :] + error]
        P[P < real_spec[i, :] - error] = real_spec[i, P < real_spec[i, :] - error] - error[P < real_spec[i, :] - error]


    # Add the smoothed row to the output matrix

    t1red_spec[i, :] = P


# Write out the changes to the original file

del f['/JasonDocument/DataPoints/0']
f.create_dataset('/JasonDocument/DataPoints/0', data=t1red_spec)
print('Dataset denoised')
f.close()

