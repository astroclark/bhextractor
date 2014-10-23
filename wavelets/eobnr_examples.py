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

import lal
import lalsimulation as lalsim


# -----------------------------
#
# --- Initial Setup
# 
# -----------------------------
sample_rate = 1024
deltaT = 1./sample_rate

# Stellar params
phiC=0.0
distance=1e9*lal.PC_SI
mtot=100
inclination=0.0

mass_ratios=np.arange(1,110,10)
Nwaveforms=len(mass_ratios)

catalogue_name='eobnr-%.0d-%.0d-%d'%(min(mass_ratios), max(mass_ratios),
        len(mass_ratios))

# -----------------------------
#
# --- EOBNRv2
# 
# -----------------------------

m2=mtot/(1.+max(mass_ratios))
m1=mtot-m2
fLow_si = 10.0

# we don't use hx right now so just direct it to _
hplus, hcross = lalsim.SimIMREOBNRv2AllModes(phiC, deltaT,
        m1*lal.MSUN_SI, m2*lal.MSUN_SI, fLow_si, distance,
        inclination)

time = np.arange(0,hplus.data.length/float(sample_rate),1.0/sample_rate)


######################
# Wavelets !

data = np.copy(hplus.data.data)

#
# CWT from pyCWT
#
#x = numpy.arange(0,2*numpy.pi,numpy.pi/8.)
#data = numpy.sin(x**2)

import cwt

scales = 1+np.arange(100)
mother_wavelet = cwt.SDG(len_signal = len(data), scales = scales,
        normalize = True, fc = 'center')
wavelet = cwt.cwt(data, mother_wavelet)

# --- Plotting
freqs = 0.5 * sample_rate * wavelet.motherwavelet.fc / wavelet.motherwavelet.scales

pl.figure()
collevs=np.linspace(min(abs(data)**2), max(abs(data)**2), 1000)
pl.contourf(time,freqs,np.abs(wavelet.coefs)**2)#, levels=collevs)
pl.show()

sys.exit()

#
# CWT in scipy
#
scales  = np.arange(1,256+1)
wavelet = signal.ricker
cwtcmatr = signal.cwt(data, wavelet, widths)


#
# DWT in pyrange(0,2*numpy.pi,numpy.pi/8.)
#data = numpy.sin(x**2)
scales = numpy.arange(10)

mother_wavelet = SDG(len_signal = len(data), scales = np.arange(10), normalize = True, fc ='center')

wavelet = cwt(data, mother_wavelet)
wavelet.scalogram(origin = 'bottom')



# --- pywt stuff
import pywt

wavelet = 'db4'
level = 4
order = "freq"  # "normal"
interpolation = 'nearest'
cmap = cm.hot

# --- perform wavelet decomposition at levels 1->maxlevel
#wp = pywt.WaveletPacket(data, wavelet, 'sym', maxlevel=level)
wd = pywt.wavedec(data, wavelet, 'sym', level=level)
sys.exit()

# --- Get the nodes of the decomposition at the level you specified
#nodes = wp.get_level(level, order=order)

# Get the time & frequency resolution at this level
Fres = sample_rate / (2**level)
Tres = 1.0/Fres
time_axis = np.arange(min(time), max(time)+Tres, Tres)

#labels = [n.path for n in nodes]
scales=len(nodes)
freqs = pywt.scal2frq(wavelet, (1+np.arange(0,scales))[::-1], 1./sample_rate)
labels = ['%.2f'%f for f in freqs]
values = np.array([n.data for n in nodes], 'd')
values = abs(values)**2


# plotting:

f = pl.figure()
f.subplots_adjust(hspace=0.2, bottom=.03, left=.07, right=.97, top=.92)
pl.subplot(2, 1, 1)
pl.title("linchirp signal")
pl.plot(time, data, 'b')

ax = pl.subplot(2, 1, 2)
pl.title("Wavelet packet coefficients at level %d" % level)
pl.imshow(values, interpolation=interpolation, cmap=cmap, aspect="auto",
    origin="lower", extent=[0, max(time), 0, max(freqs)])

pl.show()

