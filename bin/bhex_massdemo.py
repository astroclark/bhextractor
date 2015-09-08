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
bhextractor_plotpca.py

Construct waveform catalogues and PCA for plotting and diagnostics
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as pl
import bhextractor_pca as bhex



# -------------------------------
# USER INPUT

catalogue_name='Q'

# END USER INPUT
# -------------------------------

# -------------------------------
# ANALYSIS

#
# Setup and then build the catalogue
#
lowmass_catalogue = bhex.waveform_catalogue(catalogue_name=catalogue_name, fs=512,
        catalogue_len=4, mtotal_ref=250, Dist=1.)

highmass_catalogue = bhex.waveform_catalogue(catalogue_name=catalogue_name, fs=512,
        catalogue_len=4, mtotal_ref=500, Dist=1.)

lowmass_times  = lowmass_catalogue.sample_times.data
highmass_times = highmass_catalogue.sample_times.data

# -------------------------------
# Interpolation example

lowmass_hplus = np.real(lowmass_catalogue.aligned_catalogue[0,:])
lowmass_amp = lowmass_catalogue.aligned_amplitudes[0,:]
lowmass_phase = lowmass_catalogue.aligned_phases[0,:]

lowmass_rec = lowmass_amp*np.exp(1j*lowmass_phase)

highmass_hplus = np.real(highmass_catalogue.aligned_catalogue[0,:])
highmass_amp = highmass_catalogue.aligned_amplitudes[0,:]
highmass_phase = highmass_catalogue.aligned_phases[0,:]

highmass_rec = highmass_amp*np.exp(1j*highmass_phase)

#
# Interpolate
#
amp_ratio = 500./250.
new_times = amp_ratio*lowmass_times
highmass_hplus_interp = np.interp(lowmass_times, new_times, lowmass_hplus)
highmass_amp_interp   = np.interp(lowmass_times, new_times, lowmass_amp)
highmass_phase_interp = np.interp(lowmass_times, new_times, lowmass_phase)


highmass_hplus /= np.linalg.norm(highmass_hplus)
highmass_hplus_interp /= np.linalg.norm(highmass_hplus_interp)
#highmass_amp_interp   /= np.linalg.norm(highmass_amp_interp) 
#highmass_phase_interp /= np.linalg.norm(highmass_phase_interp) 

highmass_rec_interp = np.real(highmass_amp_interp *
        np.exp(1j*highmass_phase_interp))
highmass_rec_interp /= np.linalg.norm(highmass_rec_interp)


f, ax = pl.subplots(figsize=(10,10))
#f, ax = pl.subplots(nrows=3, figsize=(10,10))
#ax[0].plot(lowmass_times, highmass_hplus, label='highmass')
#ax[0].plot(lowmass_times, highmass_hplus_interp, label='lowmass rescaled', linewidth=2)
#ax[0].plot(lowmass_times, highmass_rec_interp, label='A,$\phi$ interp rec')

ax.plot(lowmass_times, lowmass_hplus, label='lowmass', marker='+')
ax.plot(lowmass_times, highmass_hplus, label='highmass', marker='+')
ax.plot(lowmass_times, highmass_hplus_interp, label='lowmass rescaled',
        linewidth=2, marker='x')
ax.plot(lowmass_times, highmass_rec_interp, label='A,$\phi$ interp rec', marker='o')

#ax[1].plot(highmass_amp, label='highmass')
#ax[1].plot(highmass_amp_interp, label='lowmass rescaled')

#ax[2].plot(highmass_phase, label='highmass')
#ax[2].plot(highmass_phase_interp, label='lowmass rescaled')

pl.show()

import sys
sys.exit()
# -------------------------------
# Plot

f, ax = pl.subplots(nrows=3, figsize=(10,10))
ax[0].plot(lowmass_times, lowmass_hplus, label='low')
ax[0].plot(highmass_times, highmass_hplus, label='high')
ax[0].set_ylabel('h$_+$(t)')

ax[1].plot(lowmass_times, lowmass_amp, label='low')
ax[1].plot(highmass_times, highmass_amp, label='high')
ax[1].set_ylabel('|h(t)|')

ax[2].plot(lowmass_times, lowmass_phase, label='low')
ax[2].plot(highmass_times, highmass_phase, label='high')
ax[2].set_ylabel('arg[h(t)]')

pl.show()




