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
inclination=0.0

fLow_si = 10.0

# -----------------------------
#
# --- EOBNRv2 Catalogue(s)
# 
# -----------------------------

mass_ratio=100

m2=mtot/(1.+mass_ratio)
m1=mtot-m2


# non-spinning
hplus_nospin, hcross_nospin = lalsim.SimIMREOBNRv2AllModes(phiC, deltaT,
        m1*lal.MSUN_SI, m2*lal.MSUN_SI, fLow_si, distance, inclination)

## spin-aligned
#s1z=0.5
#s2z=0.5
#hplus_spin, hcross_spin = lalsim.SimIMRSpinAlignedEOBWaveform(phiC, deltaT,
#        m1*lal.MSUN_SI, m2*lal.MSUN_SI, fLow_si, distance, inclination, s1z,
#        s2z, 2)

time_nospin = np.arange(0,hplus_nospin.data.length/float(sample_rate),1.0/sample_rate)
#time_spin = np.arange(0,hplus_spin.data.length/float(sample_rate),1.0/sample_rate)

pl.figure()
pl.plot(time_nospin,hplus_nospin.data.data)
#pl.plot(time_spin,hplus_spin.data.data)
pl.show()


