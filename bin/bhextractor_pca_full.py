#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2014-2015 James Clark <james.clark@ligo.org>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
"""

from __future__ import division

import os
import sys 

import numpy as np
import scipy.signal as signal
import scipy.io as sio

import lal
import lalsimulation as lalsim

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

def taper_start(input_data):
    """  
    Taper the start of the data
    """

    timeseries = lal.CreateREAL8TimeSeries('blah', 0.0, 0,
            1.0/16384, lal.StrainUnit, int(len(input_data)))
    timeseries.data.data = input_data

    lalsim.SimInspiralREAL8WaveTaper(timeseries.data,
        lalsim.SIM_INSPIRAL_TAPER_START)

    return timeseries.data.data

def window_wave(input_data):

    nonzero=np.argwhere(abs(input_data)>0)
    idx = range(nonzero[0],nonzero[-1])
    beta = 0.25
    win = lal.CreateTukeyREAL8Window(len(idx), beta)
    win.data.data[int(0.5*len(idx)):]=1.0

    input_data[idx] *= win.data.data

    return input_data

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 

catalogue_name=sys.argv[1]
theta=float(sys.argv[2])
phi=0.0

# Quick 'n' dirty sanity check (should use option parser really)
if catalogue_name not in ['Q','HR','RO3']:
    print >> sys.stderr, \
            "Error, catalogue '%s' not recognised "\
            "(must be Q, HR or RO3)"%catalogue_name

targetlen=5000
NRD_sampls=500 # number of samples to retain after peak
NINSP_sampls=2000 # number of samples to retain before peak

# Dictionary with the maximum waveform data lengths
maxlengths={'Q':10959, 'HR':31238, 'RO3':16493}

# Dictionary of NR sample rates. XXX Should probably get this from the data,
# really
NR_deltaT={'Q':0.155, 'HR':0.08, 'RO3':2./15}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct Waveforms

# --- Identify waveforms in catalogue
catalogue_path=os.environ['BHEX_PREFIX']+'/data/NR_data/'+catalogue_name+'-series'
waveforms = [ f for f in os.listdir(catalogue_path) if
        os.path.isdir(os.path.join(catalogue_path,f)) ]

# --- Load waveforms

# Preallocate a matrix for the catalogue
wflen=maxlengths[catalogue_name]
catalogue=np.zeros(shape=(wflen,len(waveforms))) + 0j

# Build waveforms in this catalogue
for w,waveform in enumerate(waveforms):
    print 'Building %s waveform'%waveform

    # Get files (modes) in this waveform dir
    waveform_dir = catalogue_path+'/'+waveform
    modes = [ f for f in os.listdir(waveform_dir) if
            os.path.isfile(os.path.join(waveform_dir,f)) ]

    # --- Sum modes
    # See line 95 of NRWaveInject.c
    hplus=np.zeros(wflen)
    hcross=np.zeros(wflen)
    for mode in modes:

        # Load mode data
        mode_data = np.loadtxt(waveform_dir+'/'+mode)

        # Extract l,m from filename
        spharm_degree = int(mode.split('_')[2].replace('l',''))
        spharm_order  = int(mode.split('_')[3].replace('m',''))

        # Compute spin -2 weighted spherical harmonic
        sYlm = lal.SpinWeightedSphericalHarmonic(theta, phi, 
                -2, spharm_degree, spharm_order)

        # Orient Waveforms
        hplus[0:len(mode_data)]  += mode_data[:,1]*np.real(sYlm) + mode_data[:,2]*np.imag(sYlm)
        hcross[0:len(mode_data)] += mode_data[:,2]*np.real(sYlm) - mode_data[:,1]*np.imag(sYlm)

    # Use complex waveform!
    catalogue[:,w] = hplus - 1j*hcross
    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Catalogue Conditioning

catalogue_real = np.real(catalogue)
catalogue_imag = np.imag(catalogue)

# Resampled catalogue
waveform_catalogue_real = np.zeros(shape=(targetlen,len(waveforms)))
waveform_catalogue_imag = np.zeros(shape=(targetlen,len(waveforms)))

print 'aligning peak times'

# Time axis
#time=mode_data[:,0]

# Find peak indices
peak_indices_real=np.argmax(abs(catalogue_real),axis=0)

# Align all waveform peaks to the latest-time peak
align_to_idx_real=max(peak_indices_real)

for w in xrange(len(waveforms)):
    print 'aligning & tapering %d of %d'%(w, len(waveforms))

    # ~~~ Align Peaks

    # Temp array to hold aligned waveform
    tmp_real = np.zeros(len(catalogue_real[:,w]))
    tmp_imag = np.zeros(len(catalogue_imag[:,w]))
    wf_real = np.copy(catalogue_real[:,w])
    wf_imag = np.copy(catalogue_imag[:,w])

    # Get the lengths of current waveform data to the left/right of the peak of
    # this waveform
    llen_real = len(wf_real[:peak_indices_real[w]])
    llen_imag = len(wf_imag[:peak_indices_real[w]])
    rlen_real = len(tmp_real[align_to_idx_real:])
    rlen_imag = len(tmp_imag[align_to_idx_real:])

    # populate left side of peak
    tmp_real[align_to_idx_real-llen_real:align_to_idx_real] = wf_real[:peak_indices_real[w]]
    tmp_imag[align_to_idx_real-llen_imag:align_to_idx_real] = wf_imag[:peak_indices_real[w]]

    # populate right side of peak
    tmp_real[align_to_idx_real:] = wf_real[peak_indices_real[w]:peak_indices_real[w]+rlen_real]
    tmp_imag[align_to_idx_real:] = wf_imag[peak_indices_real[w]:peak_indices_real[w]+rlen_real]

    # Re-populate catalogue with the aligned waveform
    catalogue_real[:,w] = np.copy(tmp_real)
    catalogue_imag[:,w] = np.copy(tmp_imag)

    # --- Resampling 
    # Resample to NR_deltaT = 0.01 (say)
    waveform_catalogue_real[:,w] = signal.resample(catalogue_real[:,w], targetlen)
    waveform_catalogue_imag[:,w] = signal.resample(catalogue_imag[:,w], targetlen)

    # --- Get rid of early-time numerical noise
    peak_idx = np.argmax(abs(waveform_catalogue_real[:,w]))
    zero_idx = max(0,peak_idx-NINSP_sampls)

    waveform_catalogue_real[:zero_idx,w] = 0.0
    waveform_catalogue_imag[:zero_idx,w] = 0.0

    # --- Get rid of late-time numerical noise
    zero_idx = min(peak_idx + NRD_sampls, targetlen)
    waveform_catalogue_real[zero_idx:,w] = 0.0
    waveform_catalogue_imag[zero_idx:,w] = 0.0

    # --- Windowing / tapering
    waveform_catalogue_real[:,w] = window_wave(waveform_catalogue_real[:,w])
    waveform_catalogue_imag[:,w] = window_wave(waveform_catalogue_imag[:,w])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA
#
# This part follows the procedure in http://arxiv.org/abs/0810.5707
#
# H = catalogue matrix (columns=waveforms)
# 
print 'Performing PCA'

# --- Combine real and imaginary parts of the catalogue
waveform_catalogue = waveform_catalogue_real #- 1j*waveform_catalogue_imag
#waveform_catalogue = waveform_catalogue_real - 1j*waveform_catalogue_imag

# --- Make catalogue a matrix for matrix arithmetic
H = np.matrix(waveform_catalogue)

# --- 1) Compute H^T * H - left side of equation (5) in 0810.5707
HTH = H.T * H

# --- 2) Compute eigenvectors (V) and eigenvalues (S) of H^T * H
#S, V = np.linalg.eigh(HTH)
S, V = np.linalg.eig(HTH)

# --- 3) Sort eigenvectors in descending order
idx = np.argsort(S)[::-1]
V = V[:,idx]
S = S[idx]

# --- 4) Sorted Eigenvectors of covariance matrix, C = PCs
# To get these, note:  H.(H^T.H) = M*C.H, since C = 1/M * H.H^T
# so H. eigenvectors of HTH are the eigenvectors of C
# i..e.,  C.U = s*U
#         (1/M)*H.H^T . U = s*U
#         U = H.V, V = eigenvectors of H^T.H
U = H*V 

for i in xrange(np.shape(U)[1]):
    U[:,i] /= np.linalg.norm(U[:,i])

# PC coefficients
Betas = H.T*U

# i.e to reconstruct H[:,0]:
h = np.zeros(len(H[:,0]),dtype=complex)
u = np.array(U)
for n in xrange(len(S)):
    h += Betas[0,n]*u[:,n]

# And to interpolate, we could do something like:

# mass ratios in Q-series:
Qs = [1.75, 1.60, 1.15, 1.30, 2.05, 1.45, 2.20, 1.90, 2.00, 1.5, 2.5, 2.35, 1.0]

interp_val = 2.1
idx=np.argsort(Qs)

q_axis = np.array(Qs)[idx]

beta0=np.concatenate(np.array(Betas[idx,0]))
beta1=np.concatenate(np.array(Betas[idx,1]))
beta2=np.concatenate(np.array(Betas[idx,2]))

beta0_target=np.interp(interp_val, np.array(Qs)[idx], beta0)
beta1_target=np.interp(interp_val, np.array(Qs)[idx], beta1)
beta2_target=np.interp(interp_val, np.array(Qs)[idx], beta2)

h_interp = beta0_target*u[:,0] + beta1_target*u[:,1] + beta2_target*u[:,2]

sys.exit()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save results

PCA_path=os.environ['BHEX_PREFIX']+'/data/'+'PCA_data'
if not os.path.exists(PCA_path): os.makedirs(PCA_path)
PCA_outname=PCA_path + '/' + catalogue_name + '_PCs_' + 'theta-%.0f'%theta
PC_dict={'PCs_final':np.array(U), 'EigVals':np.array(S)}
sio.savemat(PCA_outname, PC_dict)

catalogue_path=os.environ['BHEX_PREFIX']+'/data/'+'signal_data'
if not os.path.exists(catalogue_path): os.makedirs(catalogue_path)
catalogue_outname=catalogue_path + '/' + catalogue_name + '_catalogue_' + 'theta-%.0f'%theta
#waveform_dict={'MDC_final':waveform_catalogue}
waveform_dict={'MDC_final':waveform_catalogue}
sio.savemat(catalogue_outname, waveform_dict)

print "PCA complete"
print "PC matrix written to %s.mat"%PCA_outname
print "Catalogue matrix written to %s.mat"%catalogue_outname











