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
bhex_cwb_match.py

Compute matches between CWB reconstruction and NR waveforms
"""

import sys, os
import lal
import pycbc.types
import pycbc.filter
import numpy as np
import scipy.stats
import scipy.interpolate
import timeit
from matplotlib import pyplot as pl



def extract_wave(inwave, datalen=4.0, sample_rate = 4096):
    extract_len = 1
    peakidx = np.argmax(inwave)
    nsamp = extract_len * sample_rate

    extracted = inwave[int(peakidx-0.5*nsamp): int(peakidx+0.5*nsamp)]

    win = lal.CreateTukeyREAL8Window(len(extracted), 0.1)
    extracted *= win.data.data

    output = np.zeros(datalen*sample_rate)
    output[0.5*datalen*sample_rate-0.5*nsamp:
            0.5*datalen*sample_rate+0.5*nsamp] = np.copy(extracted)

    return output

def resample(inwave, current_delta_t = 1./4096, target_delta_t=1./2048):
    """
    Resample the timeseries in the numpy array inwave with sample spacing
    inwave_delta_t to a new sample spacing target_delta_t
    """

    num_new_samps = int(len(inwave) * current_delta_t / target_delta_t)

    return scipy.signal.resample(inwave, num_new_samps)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

strain_deltaT = 1./2048
response_deltaT = 1./4096

datalen= 4.0
f_min = 30.0

#
# --- Reconstruction data
#
print "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed")

#
# Astrophysical strain
#
h1_strain = np.loadtxt(os.path.join(event_file_dir, "cwb/H1_wf_strain.dat"))
l1_strain = np.loadtxt(os.path.join(event_file_dir, "cwb/L1_wf_strain.dat"))

# Extract the 4 second chunk in the middle
#h1_strain = extract_wave(h1_strain, datalen=datalen, sample_rate=1./strain_deltaT)
#l1_strain = extract_wave(l1_strain, datalen=datalen, sample_rate=1./strain_deltaT)

#
# Whitened response
#

h1_response = np.loadtxt(os.path.join(event_file_dir, "cwb/H1_wf_signal.dat"))
l1_response = np.loadtxt(os.path.join(event_file_dir, "cwb/L1_wf_signal.dat"))

# Extract the 4 second chunk in the middle
#   h1_response = extract_wave(h1_response, datalen=datalen,
#           sample_rate=1./response_deltaT)
#   l1_response = extract_wave(l1_response, datalen=datalen,
#           sample_rate=1./response_deltaT)

# Downsample
#   h1_response = resample(h1_response)
#   l1_response = resample(l1_response)

h1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "cwb/h1_asd.dat"))
l1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "cwb/l1_asd.dat"))

#h1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/IFO0_asd.dat"))
#l1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/IFO1_asd.dat"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conditioning
#

h1_strain   = pycbc.types.TimeSeries(h1_strain, delta_t=strain_deltaT)
l1_strain   = pycbc.types.TimeSeries(l1_strain, delta_t=strain_deltaT)

#h1_response = pycbc.types.TimeSeries(h1_response, delta_t=strain_deltaT)
#l1_response = pycbc.types.TimeSeries(l1_response, delta_t=strain_deltaT)
h1_response = pycbc.types.TimeSeries(h1_response, delta_t=response_deltaT)
l1_response = pycbc.types.TimeSeries(l1_response, delta_t=response_deltaT)


h1_strain.data /= pycbc.filter.sigma(h1_strain)
l1_strain.data /= pycbc.filter.sigma(l1_strain)
h1_response.data /= pycbc.filter.sigma(h1_response)
l1_response.data /= pycbc.filter.sigma(l1_response)

sample_freqs = h1_strain.to_frequencyseries().sample_frequencies.data

h1_asd = np.interp(sample_freqs, h1_cwb_asd_data[:,0], h1_cwb_asd_data[:,1])
l1_asd = np.interp(sample_freqs, l1_cwb_asd_data[:,0], l1_cwb_asd_data[:,1])

# Whiten the reconstructed strain
H1_strain_white = h1_strain.to_frequencyseries()
H1_strain_white.data /= h1_asd
h1_strain_white = H1_strain_white.to_timeseries()
h1_strain_white.data /= pycbc.filter.sigma(h1_strain_white)

L1_strain_white = l1_strain.to_frequencyseries()
L1_strain_white.data /= l1_asd
l1_strain_white = L1_strain_white.to_timeseries()
l1_strain_white.data /= pycbc.filter.sigma(l1_strain_white)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Matches
#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots
#

#time = np.arange(0, datalen, strain_deltaT) 
time = np.arange(0, len(h1_strain)*strain_deltaT, strain_deltaT) 

h1_strain_white_time = time - time[np.argmax(abs(h1_strain_white))]
l1_strain_white_time = time - time[np.argmax(abs(l1_strain_white))]

#time = np.arange(0, datalen, response_deltaT) 
time = np.arange(0, len(h1_response)*response_deltaT, response_deltaT) 

h1_response_time = time - time[np.argmax(abs(h1_response))]
l1_response_time = time - time[np.argmax(abs(l1_response))]


# --- Time Series
f, ax = pl.subplots(nrows=2, figsize=(10,8))

ax[0].plot(h1_response_time, h1_response, label='Response')
ax[0].plot(h1_strain_white_time, -1*h1_strain_white, 
        label='Rec. Strain / $\sqrt{S(f)}$')
ax[0].set_title('H1 Response')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Amplitude')
ax[0].minorticks_on()
ax[0].set_xlim(-0.25, .15)
ax[0].legend()

ax[1].plot(l1_response_time, l1_response, label='Response')
ax[1].plot(l1_strain_white_time, -1*l1_strain_white, 
        label='Rec. Strain / $\sqrt{S(f)}$')
ax[1].set_title('L1 Response')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Amplitude')
ax[1].minorticks_on()
ax[1].set_xlim(-0.25, .15)
ax[1].legend()

f.tight_layout()

# --- Frequency series
f, ax = pl.subplots(nrows=2, figsize=(10,8))
H1_response=h1_response.to_frequencyseries()
L1_response=l1_response.to_frequencyseries()

H1_strain_white=h1_strain_white.to_frequencyseries()
L1_strain_white=l1_strain_white.to_frequencyseries()

sample_freqs=H1_response.sample_frequencies
fnew = sample_freqs+16.5
test = np.interp(fnew, sample_freqs, abs(H1_response.data))

ax[0].plot(H1_response.sample_frequencies, abs(H1_response), label='Rec. IFO Response')
ax[0].plot(H1_strain_white.sample_frequencies, abs(H1_strain_white),
        label='Rec. Strain / $\sqrt{S(f)}$')
ax[0].set_title('H1 Response')
ax[0].set_xlabel('Frequency [Hz]')
ax[0].set_ylabel('Amplitude')
ax[0].minorticks_on()
ax[0].set_xlim(9, 512)
ax[0].legend(loc='upper right')
ax[0].axvline(30, color='r')

ax[1].plot(L1_response.sample_frequencies, abs(L1_response), label='Response')
ax[1].plot(L1_strain_white.sample_frequencies, abs(L1_strain_white),
        label='Rec. Strain / $\sqrt{S(f)}$')
ax[1].set_title('L1 Response')
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel('Amplitude')
ax[1].minorticks_on()
ax[1].set_xlim(9, 512)
ax[1].legend(loc='upper right')
ax[1].axvline(30, color='r')

f.tight_layout()


pl.show()

sys.exit()
