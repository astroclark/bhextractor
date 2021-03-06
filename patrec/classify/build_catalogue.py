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
#from os import os.environ,os.listdir,os.makedirs
#from os.path import os.isdir, os.isfile, join
import sys 

__author__ = "James Clark <james.clark@ligo.org>"

import numpy as np
import scipy.signal as signal
import scipy.io as sio

import lal
import lalsimulation as lalsim

import pycbc
import pycbc.filter


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 

catalogue_name=sys.argv[1]
theta=float(sys.argv[2])
phi=0.0

# Quick 'n' dirty sanity check (should use option parser really)
catalogue_names = ['Q','HR','RO3']
if catalogue_name not in ['Q','HR','RO3']:
    print >> sys.stderr, \
            "Error, catalogue '%s' not recognised "\
            "(must be Q, HR or RO3)"%catalogue_name

NRD_sampls=500 # number of samples to retain after peak

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
times=np.zeros(shape=(wflen,len(waveforms)))

# Get waveforms in this catalogue
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


    # XXX: for now we only consider optimally oriented signals so we'll only add
    # the plus polarisation to the catalogue
    catalogue[:,w] = hplus - 1j*hcross
    

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Align Peaks

catalogue_real = np.real(catalogue)
catalogue_imag = np.imag(catalogue)

print 'aligning peak times'

# Time axis
time=mode_data[:,0]

# Find peak indices
peak_indices_real=np.argmax(abs(catalogue_real),axis=0)
#peak_indices_imag=np.argmax(abs(catalogue_imag),axis=0)

# Align all waveform peaks to the latest-time peak
align_to_idx_real=max(peak_indices_real)
#align_to_idx_imag=max(peak_indices_imag)

resamp_catalogue_real = np.zeros(shape=(max(maxlengths.values()),
    len(waveforms)))
resamp_catalogue_imag = np.zeros(shape=(max(maxlengths.values()),
    len(waveforms)))

for w in range(len(waveforms)):

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

    # --- Tapering is done on lal TimeSeries objects
    hoft_real = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0,
            NR_deltaT[catalogue_name], lal.StrainUnit, maxlengths[catalogue_name])
    hoft_imag = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0,
            NR_deltaT[catalogue_name], lal.StrainUnit, maxlengths[catalogue_name])

    hoft_real.data.data = tmp_real
    hoft_imag.data.data = tmp_real

    lalsim.SimInspiralREAL8WaveTaper(hoft_real.data,
            lalsim.SIM_INSPIRAL_TAPER_STARTEND)
    lalsim.SimInspiralREAL8WaveTaper(hoft_imag.data,
            lalsim.SIM_INSPIRAL_TAPER_STARTEND)

    # Resample so that we have equal sampling for all catalogues
    real_part = signal.resample(hoft_real.data.data, max(maxlengths.values()))
    imag_part = signal.resample(hoft_real.data.data, max(maxlengths.values()))

    # Populate catalogue with the aligned waveform
    resamp_catalogue_real[:,w] = np.copy(real_part)
    resamp_catalogue_imag[:,w] = np.copy(imag_part)

    # --- Get rid of late-time numerical noise
    zero_tail_idx_real = np.argmax(abs(resamp_catalogue_real[:,w]))+NRD_sampls
    zero_tail_idx_imag = np.argmax(abs(resamp_catalogue_imag[:,w]))+NRD_sampls
    catalogue_real[zero_tail_idx_real:,w] = 0.0
    catalogue_imag[zero_tail_idx_real:,w] = 0.0

# --- Combine real and imaginary parts of the catalogue
waveform_catalogue = resamp_catalogue_real - 1j*resamp_catalogue_imag

#catalogue_path=os.environ['BHEX_PREFIX']+'/data/'+'signal_data'
catalogue_path='./'

#if not os.path.exists(catalogue_path): os.makedirs(catalogue_path)

catalogue_outname=catalogue_path + '/' + catalogue_name + '_catalogue_' + 'theta-%.0f'%theta

waveform_dict={'MDC_final':waveform_catalogue}
sio.savemat(catalogue_outname, waveform_dict)

sys.exit()


# ---------- MASS DEPENDENT CALCULATIONS -----------
# Determine current sample rate.  XXX: This assumes a mass scale
NR_deltaT = Mtot * lal.MTSUN_SI * NR_deltaT[catalogue_name]
Mscale = Mtot * lal.MRSUN_SI / (Dist * 1e9 * lal.PC_SI)

# Resampled catalogue
resamp_catalogue_real = np.zeros(shape=(wflen*NR_deltaT*fs,len(waveforms)))
resamp_catalogue_imag = np.zeros(shape=(wflen*NR_deltaT*fs,len(waveforms)))

# Resample and taper
for w in range(len(waveforms)):



        # --- Resampling (Note that LAL only handles downsampling by a factor of 2)
        # XXX apply mass scaling here.
        resampled_wf_real = Mscale*signal.resample(hoft_real.data.data, wflen*NR_deltaT * fs)
        resampled_wf_imag = Mscale*signal.resample(hoft_imag.data.data, wflen*NR_deltaT * fs)

        # High-pass it above 10 Hz
        tmp_real = pycbc.types.TimeSeries(initial_array=resampled_wf_real, delta_t=1.0/fs)
        tmp_real = pycbc.filter.highpass(tmp_real, frequency=10, filter_order=8, attenuation=0.1)
        resampled_wf_real = np.copy(tmp_real.data)
        tmp_imag = pycbc.types.TimeSeries(initial_array=resampled_wf_imag, delta_t=1.0/fs)
        tmp_imag = pycbc.filter.highpass(tmp_imag, frequency=10, filter_order=8, attenuation=0.1)
        resampled_wf_imag = np.copy(tmp_imag.data)

        # --- Compute SNR and rescale to SNR=1
        hoft_real = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0,
                1/float(fs), lal.StrainUnit, len(resampled_wf_real))
        hoft_real.data.data = resampled_wf_real
        rho, hrss, _, _, _ = optimal_snr(hoft_real,freqmin=10,freqmax=100)
        resampled_wf_real *= 1.0/rho
        hoft_imag = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0,
               1/float(fs), lal.StrainUnit, len(resampled_wf_imag))
        hoft_imag.data.data = resampled_wf_imag
        rho, hrss, _, _, _ = optimal_snr(hoft_imag,freqmin=10,freqmax=100)
        resampled_wf_imag *= 1.0/rho

        # Populate catalogue with scaled, resampled, tapered data. 
        resamp_catalogue_real[:,w] = resampled_wf_real
        resamp_catalogue_imag[:,w] = resampled_wf_imag

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create standardised catalogue 

targetlen=4*fs
waveform_catalogue_real = np.zeros(shape=(targetlen,len(waveforms)))
waveform_catalogue_imag = np.zeros(shape=(targetlen,len(waveforms)))
peak_indices_real=np.argmax(resamp_catalogue_real,axis=0)
#peak_indices_imag=np.argmax(resamp_catalogue_imag,axis=0)
align_to_idx_real=np.floor(0.75*targetlen)
#align_to_idx_imag=np.floor(0.75*targetlen)
for w in range(len(waveforms)):

    # current waveform
    wf_real=resamp_catalogue_real[:,w]
    wf_imag=resamp_catalogue_imag[:,w]

    # Get the lengths of current waveform data to the left/right of the peak of
    # this waveform
    llen_real = len(wf_imag[:peak_indices_real[w]])
    rlen_real = len(wf_imag[peak_indices_real[w]:])
    llen_imag = len(wf_imag[:peak_indices_real[w]])
    rlen_imag = len(wf_imag[peak_indices_real[w]:])


    # populate left side of peak
    waveform_catalogue_real[align_to_idx_real-llen_real:align_to_idx_real,w] = wf_real[:peak_indices_real[w]]
    waveform_catalogue_imag[align_to_idx_real-llen_imag:align_to_idx_real,w] = wf_imag[:peak_indices_real[w]]
    # populate right side of peak
    waveform_catalogue_real[align_to_idx_real:align_to_idx_real+rlen_real,w] = wf_real[peak_indices_real[w]:peak_indices_real[w]+rlen_real]
    waveform_catalogue_imag[align_to_idx_real:align_to_idx_real+rlen_imag,w] = wf_imag[peak_indices_real[w]:peak_indices_real[w]+rlen_imag]

del resamp_catalogue_real,resamp_catalogue_imag

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA
#
# This should be a straight copy of the matlab PCA scripts.  It looks like there
# are some differences in eig() between numpy and matlab but the results should
# come out the same.
#
# H = catalogue matrix (columns=waveforms)
# 
print 'Performing PCA'

# --- Combine real and imaginary parts of the catalogue
waveform_catalogue = waveform_catalogue_real - 1j*waveform_catalogue_imag

# --- Make catalogue a matrix for matrix arithmetic
#H = np.matrix(waveform_catalogue)
H = np.matrix(waveform_catalogue)

# --- 1) Compute catalogue covariance
C = H.T * H / np.shape(H)[0]

# --- 2) Compute eigenvectors (V) and eigenvalues (S) of C
S, V = np.linalg.eig(C)

# --- 3) Sort eigenvectors in descending order
idx = np.argsort(S)[::-1]
V = V[:,idx]
S = S[idx]

# --- 4) Compute the eigenvectors of the real covariance matrix U and normalise
U = H*V 
U /= np.linalg.norm(U,axis=0)

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











