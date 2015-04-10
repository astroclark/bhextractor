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
import scipy.linalg

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

def highpass(timeseries, delta_t=1./2048, knee=10., order=8, attn=0.1):

    tmp = pycbc.types.TimeSeries(initial_array=timeseries, delta_t=delta_t)

    return np.array(pycbc.filter.highpass(tmp, frequency=knee, filter_order=order,
            attenuation=attn).data)

def comp_norm(timeseries, delta_t=1./2048, snr=True, flow=10.):
    """
    If snr=True (default) returns the optimal SNR in aLIGO ZDHP. 
    If snr=False, returns the hrss
    """

    tmp = pycbc.types.TimeSeries(initial_array=timeseries, delta_t=delta_t)

    if snr:

        # make psd
        flen = len(tmp.to_frequencyseries())
        delta_f = np.diff(tmp.to_frequencyseries().sample_frequencies)[0]
        psd = aLIGOZeroDetHighPower(flen, delta_f, low_freq_cutoff=flow)

        return pycbc.filter.sigma(tmp, psd=psd, low_frequency_cutoff=flow)

    else:

        return pycbc.filter.sigma(tmp, low_frequency_cutoff=flow)

def freqseries(timeseries, delta_t=1./2048):
    tmp = pycbc.types.TimeSeries(initial_array=timeseries, delta_t=delta_t)
    return np.array(tmp.to_frequencyseries().data)

def pca(catalogue):

    U, S, Vt = scipy.linalg.svd(catalogue, full_matrices=False)
    V = Vt.T

    # sort the PCs by descending order of the singular values (i.e. by the
    # proportion of total variance they explain)
    ind = np.argsort(S)[::-1]
    U = U[:,ind]
    S = S[ind]
    V = V[:,ind]

    # See e.g.,:
    # http://en.wikipedia.org/wiki/Principal_component_analysis#Singular_value_decomposition

    # Score matrix:
    PCs = U * S 

    return PCs, V, S**2 #Betas



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

# Dictionary with the maximum waveform data lengths
maxlengths={'Q':10959, 'HR':31238, 'RO3':16493}

# Dictionary of NR sample rates. XXX Should probably get this from the data,
# really
NR_deltaT={'Q':0.155, 'HR':0.08, 'RO3':2./15}

# Samples to discard due to NR noise
NRD_sampls=500 # keep this many samples after the peak
NINSP_sampls=2000 # discard this many samples from the start


# We shall save the PCs as 4 second long time series sampled at 2048 Hz and
# scale the waveforms to 250 Msun
fs       = 2048
catalogue_len = 4 * fs
Mtot          = 250.
Dist          = 1.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct Waveforms

# --- Identify waveforms in catalogue
catalogue_path=os.environ['BHEX_PREFIX']+'/data/NR_data/'+catalogue_name+'-series'
waveforms = [ f for f in os.listdir(catalogue_path) if
        os.path.isdir(os.path.join(catalogue_path,f)) ]

# ---------- MASS-DEPENDENT CALCULATIONS -----------
# Determine current sample rate.  XXX: GIVEN A MASS 
SI_deltaT = Mtot * lal.MTSUN_SI * NR_deltaT[catalogue_name]
Mscale = Mtot * lal.MRSUN_SI / (Dist * 1e9 * lal.PC_SI)

# Get the length of the NR data
wflen=maxlengths[catalogue_name]

# Length to resample to: length of (e.g.) 250 Msun waveform in seconds x desired sample rate
resamp_len = np.ceil(wflen*SI_deltaT*fs)
catalogue=np.zeros(shape=(resamp_len,len(waveforms))) + 0j


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

    #
    # Clean up the waveform
    #

    # Get rid of junk
    peak_idx = np.argmax(abs(hplus))
    zero_idx = min(peak_idx + NRD_sampls, len(hplus))
    hplus[zero_idx:]  = 0.0
    hcross[zero_idx:] = 0.0

    #hplus[:NINSP_sampls] = 0.0
    #hcross[:NINSP_sampls] = 0.0
    hplus[:peak_idx-NINSP_sampls] = 0.0
    hcross[:peak_idx-NINSP_sampls] = 0.0
    

    # --- Windowing / tapering
    hplus=window_wave(hplus)
    hcross=window_wave(hcross)

    # --- Resampling 
    hplus = signal.resample(hplus, resamp_len)
    hcross = signal.resample(hcross, resamp_len)

    # --- Filtering
    hplus = highpass(hplus)
    hcross = highpass(hcross)

    # --- Set polarisations to unit norm
    # XXX: might need to be careful with relative +,x amplitudes...
    #hplus /= np.linalg.norm(hplus)
    #hcross /= np.linalg.norm(hcross)

    # Use complex waveform!
    catalogue[:,w] = hplus - 1j*hcross
    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Catalogue Conditioning

catalogue_real = np.real(catalogue)
catalogue_imag = np.imag(catalogue)

# Standardised catalogue (4 seconds long)
waveform_catalogue_real = np.zeros(shape=(catalogue_len,len(waveforms)))
waveform_catalogue_imag = np.zeros(shape=(catalogue_len,len(waveforms)))

print 'aligning peak times'

# Find peak indices
peak_idx=np.argmax(abs(catalogue_real),axis=0)

# Align all waveform peaks to the 3/4 of the way through the final catalogue
#align_idx=np.floor(0.75*catalogue_len)
align_idx=np.floor(0.5*catalogue_len)

for w in xrange(len(waveforms)):
    print 'aligning %d of %d'%(w, len(waveforms))

    non_zero_idx = \
            np.argwhere(abs(catalogue_real[:,w])>1e-2*max(abs(catalogue_real[:,w])))
    trunc_wav_real = \
            catalogue_real[non_zero_idx[0]:non_zero_idx[-1],w]
    trunc_wav_imag = \
            catalogue_imag[non_zero_idx[0]:non_zero_idx[-1],w]

    peak_idx=np.argmax(abs(trunc_wav_real))
    start_idx = align_idx - peak_idx

    waveform_catalogue_real[start_idx:start_idx+len(trunc_wav_real),w] = trunc_wav_real
    waveform_catalogue_imag[start_idx:start_idx+len(trunc_wav_real),w] = trunc_wav_imag

    # ~~~ Align Peaks
#   if catalogue_len >= resamp_len:
#       # populate the final catalogue with the full resampled waveform
#       waveform_catalogue_real[start_idx:resamp_len+start_idx,w] = catalogue_real[:,w]
#       waveform_catalogue_imag[start_idx:resamp_len+start_idx,w] = catalogue_imag[:,w]
#
#   else:
#       # populate the final catalogue with the truncated waveform
#       waveform_catalogue_real[start_idx:,w] = catalogue_real[:catalog_len-start_idx,w]
#       waveform_catalogue_imag[start_idx:,w] = catalogue_imag[:catalog_len-start_idx,w]

    # --- Normalisation (apply identical scaling to real, imag)
    N = comp_norm(waveform_catalogue_real[:,w], snr=False)
    waveform_catalogue_real[:,w] /= N
    waveform_catalogue_imag[:,w] /= N

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA
#

print 'Performing PCA'

Heng=False
if Heng==True:
    # This part follows the procedure in http://arxiv.org/abs/0810.5707
    #
    # H = catalogue matrix (columns=waveforms)
    # 
    # WARNING: I don't understand why the cross polarisation winds up looking so
    # small...


    # --- Combine real and imaginary parts of the catalogue
    waveform_catalogue = waveform_catalogue_real - 1j*waveform_catalogue_imag
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

    U = np.array(U)
    S = np.array(S)

    # normalise PCs
    for i in xrange(np.shape(U)[1]):
        U[:,i] /= np.linalg.norm(U[:,i])

else:
    # Use SVD

    waveform_catalogue = waveform_catalogue_real - 1j*waveform_catalogue_imag

    U, V, S = pca(waveform_catalogue)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save results

PCA_path=os.environ['BHEX_PREFIX']+'/data/'+'PCA_data'

if not os.path.exists(PCA_path): os.makedirs(PCA_path)

PCA_outname=PCA_path + '/' + catalogue_name + '_PCs_' + 'theta-%.0f'%theta


#
# Save FFTs of PCs to binary file 
#

# These are the template basis functions LIB will use to find the signal.

U_fdomain_plus = np.zeros(shape=(0.5*np.shape(U)[0]+1, np.shape(U)[1]), dtype=complex)
U_fdomain_cross = np.zeros(shape=(0.5*np.shape(U)[0]+1, np.shape(U)[1]), dtype=complex)
for i in xrange(np.shape(U)[1]):
    #U_fdomain_plus[:,i] = freqseries(np.real(U[:,i]))
    U_fdomain_plus[:,i] = freqseries(np.real(U[:,i]))
    U_fdomain_cross[:,i] = freqseries(np.imag(U[:,i]))

# Save to TEXT files
fp = open("%s_plus.dat"%PCA_outname, 'w')
fc = open("%s_cross.dat"%PCA_outname, 'w')

#U_fdomain_plus.tofile(fp)
#U_fdomain_cross.tofile(fc)

for i in xrange(np.shape(U_fdomain_plus)[0]):
    for j in xrange(np.shape(U_fdomain_plus)[1]):
        fp.write("%.15f %.10f\t"%(
            np.real(U_fdomain_plus[i,j]), np.imag(U_fdomain_plus[i,j])
            ))
        fc.write("%.15f %.10f\t"%(
            np.real(U_fdomain_cross[i,j]), np.imag(U_fdomain_cross[i,j])
            ))
    fp.write("\n")
    fc.write("\n")

fp.close()
fc.close()

#
# Save dictionary with PCs to mat file
#
PC_dict={'PCs_final':U, 'EigVals':S, 'Betas':V, \
        'pcs_plus':U_fdomain_plus, 'pcs_cross':U_fdomain_cross}
sio.savemat(PCA_outname, PC_dict)

catalogue_path=os.environ['BHEX_PREFIX']+'/data/'+'signal_data'
if not os.path.exists(catalogue_path): os.makedirs(catalogue_path)
catalogue_outname=catalogue_path + '/' + catalogue_name + '_catalogue_' + 'theta-%.0f'%theta
#waveform_dict={'MDC_final':waveform_catalogue}
waveform_dict={'MDC_final':waveform_catalogue}
sio.savemat(catalogue_outname, waveform_dict)

print "PCA complete"
print "PC matrix written to %s.mat and %s.dat"%(PCA_outname, PCA_outname)
print "Catalogue matrix written to %s.mat"%catalogue_outname











