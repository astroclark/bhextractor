#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2015-2016 James Clark <james.clark@ligo.org>
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
bw-cwb_compare.py

"""

import sys, os
import timeit
import pycbc.types
import pycbc.filter
import numpy as np
import scipy.interpolate
import scipy.signal
import triangle
from matplotlib import pyplot as pl

#def extract_wave(inwave, datalen=4.0, sample_rate = 2048):
def extract_wave(inwave, datalen=4.0, sample_rate = 4096):
    """
    Pull out datalen seconds around the peak of the timeseries contained in the
    numpy array inwave, sampled at sample_rate Hz
    """
    peakidx = np.argmax(inwave)
    nsamp = datalen * sample_rate
    return inwave[int(peakidx-0.5*nsamp): int(peakidx+0.5*nsamp)]

def resample(inwave, datalen=4.0, target_delta_t=1./1024):
    """
    Resample the timeseries in the numpy array inwave with sample spacing
    inwave_delta_t to a new sample spacing target_delta_t
    """

    return scipy.signal.resample(inwave, int(datalen/target_delta_t))

def compute_match(h1, h2, asd_data, delta_t=1./1024, f_min=30.0):
    """
    Compute the match between time series h1 and h2 using the ASD data contained
    in asd_data.
        h1: numpy array of h(t) samples
        h2: numpy array of h(t) samples
        asd_data: 2xN array of [frequency, ASD_value] samples

    NOTE:  h1 is assumed to be strain (e.g., CWB)
           h2 is assumed to be response (e.g., BW)
           h1 is FFT'd and whitened by the ASD.   The ASD is NOT therefore
           passed to pycbc.filter.match()
    """

    # Insert into pycbc time series
    h1_pycbc = pycbc.types.TimeSeries(h1, delta_t=delta_t)
    h2_pycbc = pycbc.types.TimeSeries(h2, delta_t=delta_t)

    # Retrieve the frequency samples for the full time series
    freq_axis = h1_pycbc.to_frequencyseries().sample_frequencies.data

    # Interpolate the asd to those frequencies (should already be the same,
    # really)
    asd = np.interp(freq_axis, asd_data[:,0], asd_data[:,1])

    # Whiten h1
    #H1_pycbc = h1_pycbc.to_frequencyseries()
    #H1_pycbc.data /= asd

    # Put this psd into a pycbc frequency series so we can compute match
    #psd = pycbc.types.FrequencySeries(asd**2, delta_f=np.diff(freq_axis)[0])
    psd = None

    #return pycbc.filter.match(h1_py, h2_py, psd=psd, low_frequency_cutoff=f_min)
    return pycbc.filter.match(h1_pycbc, h2_pycbc, psd=psd,
            low_frequency_cutoff=f_min)[0]
    #return pycbc.filter.match(H1_pycbc, h2_pycbc, psd=None,
    #        low_frequency_cutoff=f_min)[0]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

delta_t = 1./1024
datalen= 4.0
f_min = 50.0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Data


#
# --- Reconstruction data
#
print "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed")

#
# BayesWave
#
bw_times = np.loadtxt(os.path.join(event_file_dir, "bw/waveforms/timesamp.dat"))
h1_bw_samples = np.loadtxt(os.path.join(event_file_dir, 
    "bw/waveforms/signal_recovered_whitened_waveform.dat.0"))
l1_bw_samples = np.loadtxt(os.path.join(event_file_dir, 
    "bw/waveforms/signal_recovered_whitened_waveform.dat.1"))

h1_bw_asd = np.loadtxt(os.path.join(event_file_dir, "bw/IFO0_asd.dat"))
l1_bw_asd = np.loadtxt(os.path.join(event_file_dir, "bw/IFO1_asd.dat"))

#
# CWB
#
h1_cwb_estimate = np.loadtxt(os.path.join(event_file_dir, "cwb/H1_wf_signal.dat"))
l1_cwb_estimate = np.loadtxt(os.path.join(event_file_dir, "cwb/L1_wf_signal.dat"))

h1_cwb_asd = np.loadtxt(os.path.join(event_file_dir, "cwb/h1_asd.dat"))
l1_cwb_asd = np.loadtxt(os.path.join(event_file_dir, "cwb/l1_asd.dat"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unify sample rates & data lengths
#

# Extract the 4 second chunk in the middle
h1_cwb_estimate = extract_wave(h1_cwb_estimate, datalen=datalen)
l1_cwb_estimate = extract_wave(l1_cwb_estimate, datalen=datalen)

# Downsample the number of posterior samples (useful for testing)
# XXX: This is TEMPORARY
#nsampls = 50
#h1_bw_samples = h1_bw_samples[np.random.random_integers(low=0,high=nsampls,size=nsampls),:]
#l1_bw_samples = l1_bw_samples[np.random.random_integers(low=0,high=nsampls,size=nsampls),:]


# Resample to common rate; in practice, this is just downsampling cwb to bw
h1_cwb_estimate = resample(h1_cwb_estimate)
l1_cwb_estimate = resample(l1_cwb_estimate)



# Interpolate the ASD to the waveform frequencies
#h1_asd = np.interp(freq_axis, h1_bw_asd_data[:,0], h1_bw_asd_data[:,1])
#l1_asd = np.interp(freq_axis, l1_bw_asd_data[:,0], l1_bw_asd_data[:,1])
#
#mean_asd = np.sqrt(scipy.stats.hmean([h1_asd**2, l1_asd**2]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Match calculations
#

h1_matches = np.zeros(len(h1_bw_samples))
l1_matches = np.zeros(len(l1_bw_samples))

big_then = timeit.time.time()
print "...Evaluating matches..."
for s, (h1_bw_sample, l1_bw_sample) in enumerate(zip(h1_bw_samples, l1_bw_samples)):

#    print '-----------------------------'
#    print "Evaluating sample waveform %d of %d"%( s, len(h1_bw_samples) )

    then = timeit.time.time()

    h1_matches[s] = compute_match(h1_cwb_estimate, h1_bw_sample,
            delta_t=delta_t, asd_data=h1_bw_asd, f_min=f_min)

    l1_matches[s] = compute_match(l1_cwb_estimate, l1_bw_sample,
            delta_t=delta_t, asd_data=l1_bw_asd, f_min=f_min)

    now = timeit.time.time()
    
#    print "...match calculation took %.3f sec..."%(now-then)

big_now = timeit.time.time()
print "...all match calculations took %.3f sec..."%(big_now-big_then)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot match distributions
#

samples = np.array([h1_matches, l1_matches]).T

trifig = triangle.corner(samples, quantiles=[0.25, 0.50, 0.75], labels=['H1 Match', 
    'L1 Match'], plot_contours=True,
    plot_datapoints=True)
title = "H1, L1 Matches: \n ( BayesWave | CWB ), f$_{\mathrm{low}}$ = %d Hz"%f_min
trifig.suptitle(title, fontsize=16)
trifig.subplots_adjust(top=0.88)

