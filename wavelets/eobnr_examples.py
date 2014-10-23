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


# -----------------------------
#
# --- Initial Setup
# 
# -----------------------------
sample_rate = 512
deltaT = 1./sample_rate

# Stellar params
phiC=0.0
distance=1e9*lal.PC_SI
mtot=500
inclination=45.0
mass_ratio=100


# -----------------------------
#
# --- EOBNRv2
# 
# -----------------------------

m2=mtot/(1.+mass_ratio)
m1=mtot-m2
fLow_si = 10.0

# we don't use hx right now so just direct it to _
hplus, hcross = lalsim.SimIMREOBNRv2AllModes(phiC, deltaT,
        m1*lal.MSUN_SI, m2*lal.MSUN_SI, fLow_si, distance,
        inclination)

time = np.arange(0,hplus.data.length/float(sample_rate),1.0/sample_rate)


######################
# Wavelets !

#idx=np.concatenate(np.argwhere(abs(hplus.data.data)>0.01*abs(hplus.data.data)))
T_max=time[np.argwhere(abs(hplus.data.data)>0.01*abs(hplus.data.data))[-1]]
idx=range(len(time))
data = hplus.data.data[idx]
time = time[idx]

#
# CWT from pyCWT
#

import cwt

scales = 1+np.arange(256)
mother_wavelet = cwt.SDG(len_signal = len(data), scales = scales,
        normalize = True, fc = 'center')
wavelet = cwt.cwt(data, mother_wavelet)

# --- Plotting
freqs = 0.5 * sample_rate * wavelet.motherwavelet.fc / wavelet.motherwavelet.scales

collevs=np.linspace(0, max(map(max,abs(wavelet.coefs)**2)), 100)
fig, ax_cont = pl.subplots(figsize=(10,5))
ax_cont.contourf(time,freqs,np.abs(wavelet.coefs)**2, levels=collevs,
        cmap=cm.gnuplot2)
ax_cont.set_yscale('log')
ax_cont.set_xlim(min(time),T_max)
ax_cont.set_ylim(min(freqs),max(freqs))
ax_cont.set_xlabel('Time [s]')
ax_cont.set_ylabel('Frequency [Hz]')

divider = make_axes_locatable(ax_cont)

# time-series
ax_ts = divider.append_axes("top", 1.2, sharex=ax_cont)
ax_ts.plot(time, data)
ax_cont.set_xlim(min(time),T_max)
ax_ts.set_ylim(-1.1*max(abs(data)), 1.1*max(abs(data)))

# fourier spectrum
freq_fourier, Pxx = signal.periodogram(data, fs=sample_rate)
ax_fs = divider.append_axes("right", 1.2, sharey=ax_cont)
ax_fs.semilogx(Pxx,freq_fourier)
ax_fs.set_ylim(min(freqs),max(freqs))
ax_fs.set_xlim(0.01*max(Pxx),1.1*max(Pxx))

pl.setp(ax_ts.get_xticklabels()+ax_fs.get_yticklabels(),visible=False)

pl.tight_layout()
pl.show()






#   #
#   # CWT in scipy
#   #
#   scales  = np.arange(1,256+1)
#   wavelet = signal.ricker
#   cwtcmatr = signal.cwt(data, wavelet, widths)
#
#
#   #
#   # DWT in pyrange(0,2*numpy.pi,numpy.pi/8.)
#   #data = numpy.sin(x**2)
#   scales = numpy.arange(10)
#
#   mother_wavelet = SDG(len_signal = len(data), scales = np.arange(10), normalize = True, fc ='center')
#
#   wavelet = cwt(data, mother_wavelet)
#   wavelet.scalogram(origin = 'bottom')
#
#
#
#   # --- pywt stuff
#   import pywt
#
#   wavelet = 'db4'
#   level = 4
#   order = "freq"  # "normal"
#   interpolation = 'nearest'
#   cmap = cm.hot
#
#   # --- perform wavelet decomposition at levels 1->maxlevel
#   #wp = pywt.WaveletPacket(data, wavelet, 'sym', maxlevel=level)
#   wd = pywt.wavedec(data, wavelet, 'sym', level=level)
#   sys.exit()
#
#   # --- Get the nodes of the decomposition at the level you specified
#   #nodes = wp.get_level(level, order=order)
#
#   # Get the time & frequency resolution at this level
#   Fres = sample_rate / (2**level)
#   Tres = 1.0/Fres
#   time_axis = np.arange(min(time), max(time)+Tres, Tres)
#
#   #labels = [n.path for n in nodes]
#   scales=len(nodes)
#   freqs = pywt.scal2frq(wavelet, (1+np.arange(0,scales))[::-1], 1./sample_rate)
#   labels = ['%.2f'%f for f in freqs]
#   values = np.array([n.data for n in nodes], 'd')
#   values = abs(values)**2
#
#
#   # plotting:
#
#   f = pl.figure()
#   f.subplots_adjust(hspace=0.2, bottom=.03, left=.07, right=.97, top=.92)
#   pl.subplot(2, 1, 1)
#   pl.title("linchirp signal")
#   pl.plot(time, data, 'b')
#
#   ax = pl.subplot(2, 1, 2)
#   pl.title("Wavelet packet coefficients at level %d" % level)
#   pl.imshow(values, interpolation=interpolation, cmap=cmap, aspect="auto",
#       origin="lower", extent=[0, max(time), 0, max(freqs)])
#
#   pl.show()

