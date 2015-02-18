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

__author__ = "James Clark <james.clark@ligo.org>"

import numpy as np
import scipy.signal as signal
import scipy.io as sio

from matplotlib import pyplot as pl
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import lal
import lalsimulation as lalsim
import pycbc.types

import cwt

# -----------------------------
#
# --- Initial Setup
# 
# -----------------------------

#
# Input 
#
outputfilename='test'
mtot=250 # total mass in solar masses (don't go lower than 200)
inclination=0.0 # initial inclination of orbital plane
mass_ratio=10  # mass ratio q = mass1 / mass2
distance=1e9*lal.PC_SI # distance from Earth to source

# Data configuration
sample_rate = 512
deltaT = 1./sample_rate
fLow_si = 10.0 # waveform turns on at 10 Hz

# phase at coalescence
phiC=0.0


# ---------------------------------
#
# --- Generate EOBNRv2 waveform --- 
# 
# ---------------------------------

m2=mtot/(1.+mass_ratio)
m1=mtot-m2

# we don't use hx right now so just direct it to _
hplus, hcross = lalsim.SimIMREOBNRv2AllModes(phiC, deltaT,
        m1*lal.MSUN_SI, m2*lal.MSUN_SI, fLow_si, distance,
        inclination)

time = np.arange(0,hplus.data.length/float(sample_rate),1.0/sample_rate)


# ------------------------
#
# --- Wavelet Analysis --- 
# 
# ------------------------

# Trim data down to wherever the waveform is greater than 1% of the maximum
# value and create a pycbc TimeSeries object from this data
idx = abs(hplus.data.data)>0.01*abs(hplus.data.data)
signal = pycbc.types.TimeSeries(initial_array=hplus.data.data[idx], delta_t =
        deltaT)

#
# Construct CWT (continuous wavelet transform) using pyCWT module
#

# Range of wavelet scales we're interested in
scales = 1+np.arange(256)

# Construct the 'mother wavelet'; we'll begin using a Morlet wavelet
# (sine-Gaussian) but we should/could investigate others
# Could also experiment with values of f0.  See cwt.Morlet for info
mother_wavelet = cwt.Morlet(len_signal = len(signal), scales = scales,
        sampf=sample_rate, f0=1)

# Compute the CWT
wavelet = cwt.cwt(signal.data, mother_wavelet)

# ------------------------
#
# --- Plot The Results --- 
# 
# ------------------------

# Compute the center frequencies corresponding to each scale
freqs = 0.5 * sample_rate * wavelet.motherwavelet.fc / wavelet.motherwavelet.scales

# Choose the number of colour levels to plo
collevs=np.linspace(0, max(map(max,abs(wavelet.coefs)**2)), 100)

# Open the figure
fig, ax_cont = pl.subplots(figsize=(10,5))

# plot a contour map of the wavelet coefficient magnitudes
ax_cont.contourf(signal.sample_times, freqs, np.abs(wavelet.coefs)**2,
        levels=collevs, cmap=cm.gnuplot2)

ax_cont.set_yscale('log')
ax_cont.set_xlabel('Time [s]')
ax_cont.set_ylabel('Frequency [Hz]')

ax_cont.set_xlim(min(signal.sample_times),max(signal.sample_times))

divider = make_axes_locatable(ax_cont)

# time-series
ax_ts = divider.append_axes("top", 1.5, sharex=ax_cont)
ax_ts.plot(signal.sample_times, signal)

# fourier spectrum
signal_frequency_spectrum = signal.to_frequencyseries()

ax_fs = divider.append_axes("right", 1.5, sharey=ax_cont)
ax_fs.semilogx(abs(signal_frequency_spectrum)**2,
        signal_frequency_spectrum.sample_frequencies)


# Adjust axis limits
ax_ts.set_xlim(min(signal.sample_times),max(signal.sample_times))
ax_cont.set_xlim(min(signal.sample_times),max(signal.sample_times))
ax_cont.set_ylim(min(signal_frequency_spectrum.sample_frequencies),
        max(signal_frequency_spectrum.sample_frequencies))
ax_fs.set_ylim(min(signal_frequency_spectrum.sample_frequencies),
        max(signal_frequency_spectrum.sample_frequencies))

pl.setp(ax_ts.get_xticklabels()+ax_fs.get_yticklabels(),visible=False)

#ax_fs.tick_params(axis='both', which='major', labelsize=8)
ax_fs.tick_params(axis='x', which='major', labelsize=8)

pl.show()

#pl.tight_layout()
pl.savefig(outputfilename+'.png')

