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
mass_scaling.py

Generates a few different mass BBH waveforms to investigate/check scaling
behaviour.
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

from matplotlib import pyplot as pl

import lal
import lalsimulation as lalsim

def startFreqHz(startFreq,m1,m2):

    mtotal=(m1+m2)*lal.LAL_MTSUN_SI
    startFreqHz=startFreq/(lal.LAL_TWOPI*mtotal)

    return startFreqHz

# -----------------------------------------
#
# --- Signal 1: low mass
# 
# -----------------------------------------
sample_rate = 1024
deltaT = 1./sample_rate

# Stellar params
phiC=0.0
distance=1e9*lal.LAL_PC_SI
mtot=100
inclination=0.0
mass_ratio=0.5

m2=mtot/(1.+mass_ratio)
m1=mtot-m2
fLow_si = 10.

# we don't use hx right now so just direct it to _
hlowmass, _ = lalsim.SimIMREOBNRv2AllModes(phiC, deltaT,
        m1*lal.LAL_MSUN_SI, m2*lal.LAL_MSUN_SI, fLow_si, distance,
        inclination)
# Taper start of waveform
lalsim.SimInspiralREAL8WaveTaper(hlowmass.data,
        lalsim.LAL_SIM_INSPIRAL_TAPER_STARTEND)
time_lowmass = np.arange(0, hlowmass.data.length*deltaT, deltaT)

idx_peak = np.argmax(abs(hlowmass.data.data))
time_lowmass -= time_lowmass[idx_peak]

#lal.ResizeREAL8TimeSeries(hlowmass,0,Nlow)

# -----------------------------------------
#
# --- Signal 2: high mass
# 
# -----------------------------------------
# Stellar params
mtot=200
m2=mtot/(1.+mass_ratio)
m1=mtot-m2

# we don't use hx right now so just direct it to _
hhighmass, _ = lalsim.SimIMREOBNRv2AllModes(phiC, deltaT,
        m1*lal.LAL_MSUN_SI, m2*lal.LAL_MSUN_SI, fLow_si, distance,
        inclination)

# Taper start of waveform
lalsim.SimInspiralREAL8WaveTaper(hhighmass.data,
        lalsim.LAL_SIM_INSPIRAL_TAPER_STARTEND)
lal.ResizeREAL8TimeSeries(hhighmass,0,hlowmass.data.length)

time_highmass = np.arange(0, hhighmass.data.length*deltaT, deltaT)

#   from pycbc import filter
#   # make pycbc time series object:
#   h_tmp = filter.matchedfilter.TimeSeries(\
#           initial_array=hhighmass.data.data, delta_t=hhighmass.deltaT)
#   h_tmp = filter.highpass(h_tmp,10,filter_order=12,attenuation=0.01)
#   hhighmass.data.data = h_tmp.data

idx_peak = np.argmax(abs(hhighmass.data.data))
time_highmass -= time_highmass[idx_peak]

# -----------------------------------------
#
# --- Signal 3: geometric/NR waveform
# 
# -----------------------------------------

Mcurrent = 100
# we don't use hx right now so just direct it to _
hnumrel = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0, deltaT,
        lal.lalStrainUnit, hlowmass.data.length)

# Copy 100 Msun waveform and rescale to geometric units
hnumrel.data.data = distance * np.copy(hlowmass.data.data) / (Mcurrent *
        lal.LAL_MRSUN_SI)

NRdeltaT = deltaT / (Mcurrent * lal.LAL_MTSUN_SI)

time_numrel = np.arange(0, hnumrel.data.length*NRdeltaT, NRdeltaT)

idx_peak = np.argmax(abs(hnumrel.data.data))
time_numrel -= time_numrel[idx_peak]

# --------------------------------------------------------
#
# --- Signal 4: rescaled geometric/NR waveform to 100 Msun
# 
# -------------------------------------------------------

Mtarget = 200

deltaT_new = NRdeltaT * (Mtarget * lal.LAL_MTSUN_SI)

# we don't use hx right now so just direct it to _
hrescaled = lal.CreateREAL8TimeSeries('hoft', lal.LIGOTimeGPS(), 0.0,
        deltaT_new, lal.lalStrainUnit, hnumrel.data.length)

# Copy 100 Msun waveform and rescale to geometric units
hrescaled.data.data = np.copy(hnumrel.data.data) \
        * Mtarget * lal.LAL_MRSUN_SI / distance

time_rescaled = np.arange(0, hrescaled.data.length*deltaT_new, deltaT_new)

# --- Finally high-pass the rescaled data to eliminate the unnecessary
#     low-frequencies (unnecessary really, but good for sanity checking)
import pycbc
from pycbc import filter
# make pycbc time series object:

hrescaled_tmp = pycbc.types.TimeSeries(\
        initial_array=hrescaled.data.data, delta_t=hrescaled.deltaT)
hhighmass_tmp = pycbc.types.TimeSeries(\
        initial_array=hhighmass.data.data, delta_t=hhighmass.deltaT)
hhighmass_tmp.resize(len(hrescaled_tmp.data))
sys.exit()

#hrescaled_tmp = filter.highpass(hrescaled_tmp,10,filter_order=12,attenuation=0.9)
hrescaled.data.data = hrescaled_tmp.data

idx_peak = np.argmax(abs(hrescaled.data.data))
time_rescaled -= time_rescaled[idx_peak]

# -----------------------------------------
#
# --- PLOTS
# 
# -----------------------------------------
fig1,ax1=pl.subplots(nrows=3,ncols=2,figsize=(8,8))

#hhighmass.data.data[time_highmass<-0.6] = 0.0
#hrescaled.data.data[time_rescaled<-0.6] = 0.0

#hhighmass_tmp = hhighmass.data.data[time_highmass>-0.6]
#hrescaled_tmp = hrescaled.data.data[time_highmass>-0.6]

# --- Time-domain waveform plot

ax1[0][0].plot(time_lowmass,hlowmass.data.data,label=r'M$_{\mathrm{tot}}=%d$'%100)
ax1[0][0].set_xlim(-2,.1)
ax1[0][0].set_ylim(-3e-21,3e-21)
ax1[0][0].set_xlabel('Time [s]')
ax1[0][0].legend(loc='upper left')

ax1[1][0].plot(time_numrel,hnumrel.data.data,label=r'NR')
ax1[1][0].set_xlabel('Time / M$_{\odot}$')
ax1[1][0].legend(loc='upper left')

ax1[2][0].plot(time_highmass,hhighmass.data.data,label=r'M$_{\mathrm{tot}}=%d$'%Mtarget)
ax1[2][0].set_xlim(-2,.1)
ax1[2][0].set_ylim(-3e-21,3e-21)
ax1[2][0].set_xlabel('Time [s]')
ax1[2][0].legend(loc='upper left')

ax1[2][0].plot(time_rescaled,hrescaled.data.data,label=r'NR to: %d'%Mtarget)
ax1[2][0].set_xlabel('Time [s]')
ax1[2][0].set_xlim(-2,.1)
ax1[2][0].set_ylim(-3e-21,3e-21)
ax1[2][0].legend(loc='upper left')

# --- PSDs of waveforms

freq, Pxx_den = signal.periodogram(hlowmass.data.data, 1/deltaT,
        window=np.ones(hlowmass.data.length))
ax1[0][1].semilogy(freq,Pxx_den,label=r'M$_{\mathrm{tot}}=%d$'%100)
ax1[0][1].set_xlim(0,200)
ax1[0][1].set_ylim(1e-48,1e-44)
ax1[0][1].axvline(10,color='k')
ax1[0][1].set_xlabel('Frequency [Hz]')

freq, Pxx_den = signal.periodogram(hnumrel.data.data, 1/NRdeltaT,
        window=np.ones(hlowmass.data.length))
ax1[1][1].semilogy(freq,Pxx_den,label=r'NR')
ax1[1][1].set_xlim(0,.1)
ax1[1][1].set_ylim(1e-5,1)
ax1[1][1].set_xlabel('Frequency [M$_{\odot}$]')

freq1, Pxx_den1 = signal.periodogram(hhighmass.data.data, 1/deltaT,
        scaling='density')
ax1[2][1].semilogy(freq1,Pxx_den1,label=r'M$_{\mathrm{tot}}=%d$'%200)
ax1[2][1].set_xlim(0,200)
ax1[2][1].set_ylim(1e-48,1e-43)
ax1[2][1].axvline(10,color='k')
ax1[2][1].set_xlabel('Frequency [Hz]')

freq2, Pxx_den2 = signal.periodogram(hrescaled.data.data, 1/deltaT_new,
        scaling='density')
ax1[2][1].semilogy(freq2,Pxx_den2,label=r'NR to: %d'%Mtarget)
ax1[2][1].set_xlim(0,200)
ax1[2][1].set_ylim(1e-48,1e-43)
ax1[2][1].axvline(10,color='k')
ax1[2][1].set_xlabel('Frequency [Hz]')

print max(Pxx_den2) / max(Pxx_den1)

fig1.tight_layout()
pl.show()




