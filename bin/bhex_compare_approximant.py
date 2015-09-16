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
gatech_EOBNR_comparison.py

Builds waveforms from the 2,2 modes of the GATech NR catalogue(s) and computes
matches with approximants (particularly EOBNR) available to pycbc.
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Imports & local definitions

import sys
import numpy as np
import bhex_wavedata as bwave

import lal
import pycbc.types
from pycbc.waveform import get_td_waveform
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

from matplotlib import pyplot as pl

def component_masses(total_mass, mass_ratio):
    """
    Return m1 and m2, given total mass and mass ratio (m1/m2)

    m1, m2 = component_masses(total_mass, mass_ratio)
    """

    m1 = mass_ratio * total_mass / (1.0 + mass_ratio)
    m2 = total_mass - m1

    return m1, m2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User Options (mass configurations, noise curves etc)

f_low = 10. # Min freq for waveform generation
f_min = 30. # Min freq for match calculation

sample_rate = 1024
datalen= 4.0

total_mass = 150.
#mass_ratio = 1.

distance=1.

#series_names = ['Eq-series']#, 'HRq-series']
series_names = ['HR-series']

#bounds=None
# Change to e.g.:
#
bounds=dict()
#bounds['q'] = [1, 2] 
bounds['a1'] = [0, 0]
bounds['a2'] = [0, 0]
bounds['q'] = [15, np.inf] 
#
# to only use simulations with q<=2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Select & Generate NR waveforms

# Pick out the simulations we want:

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''

simulations_list = bwave.simulation_details(series_names=series_names,
        param_bounds=bounds, Mmin30Hz=total_mass)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
NR_catalogue = bwave.waveform_catalogue(simulations_list, ref_mass=total_mass,
        sample_rate=sample_rate, datalen=datalen, distance=distance)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate EOBNR

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'Beginning Match Calculations'
print ''

matches = np.zeros(simulations_list.nsimulations)
mass_ratios = np.zeros(simulations_list.nsimulations)

# Retrieve the mass ratio:
f1, ax1 = pl.subplots()
f2, ax2 = pl.subplots()
for s,simulation in enumerate(simulations_list.simulations):
    
    print 'Matching waveform %d of %d'%(s, simulations_list.nsimulations)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # --- Condition the NR waveform ---

    # Get the NR (plus) wave and put it in a pycbc TimeSeries object
    hplus_NR = \
            pycbc.types.TimeSeries(np.real(NR_catalogue.SIComplexTimeSeries[s,:]),
                    delta_t=1./sample_rate)
    hplus_NR.data = bwave.window_wave(hplus_NR.data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # --- Generate Approximant ---
    #
    mass_ratios[s] = simulation['q']

    mass1, mass2 = component_masses(total_mass, mass_ratios[s])

    hplus_EOBNR, _ = get_td_waveform(approximant="SEOBNRv2",
            distance=distance,
            mass1=mass1,
            mass2=mass2,
            spin1z=0.99,
            spin2z=0.99,
            f_lower=f_low,
            delta_t=1.0/sample_rate)
            #f_lower=fpeak_NR,
    # divide out the spherical harmonic (2,2) amplitude (this is just for nice
    # plots / sanity - it does not affect match)
    sY22 = lal.SpinWeightedSphericalHarmonic(0,0,-2,2,2)
    hplus_EOBNR.data /= np.real(sY22)

    # (1-sided) Planck window for smooth FFT
    hplus_EOBNR.data = bwave.window_wave(hplus_EOBNR.data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Match Calculation
    #

    # Make the timeseries consistent lengths
    tlen = max(len(hplus_NR), len(hplus_EOBNR))
    hplus_EOBNR.resize(tlen)
    hplus_NR.resize(tlen)

    # Generate a noise curve:
    Hplus_NR = hplus_NR.to_frequencyseries()
    delta_f = Hplus_NR.delta_f
    flen = len(Hplus_NR)
    psd = aLIGOZeroDetHighPower(flen, delta_f, f_low) 

    m, snr_max_loc = pycbc.filter.match(hplus_EOBNR, hplus_NR,
            low_frequency_cutoff=f_min, psd=psd)
    matches[s] = np.copy(m)
    # match() returns the index for peak snr as well as the match

    NR_max_loc = np.argmax(hplus_NR)

    ax1.plot(hplus_EOBNR.sample_times, hplus_EOBNR, linestyle='-')
    ax1.plot(hplus_NR.sample_times - hplus_NR.sample_times[NR_max_loc],
            -1*hplus_NR, linestyle='--')
    ax1.minorticks_on()

    ax2.loglog(hplus_EOBNR.to_frequencyseries().sample_frequencies,
            abs(hplus_EOBNR.to_frequencyseries()))
    ax2.loglog(hplus_NR.to_frequencyseries().sample_frequencies,
            abs(hplus_NR.to_frequencyseries()), linestyle='--')
    ax2.axvline(30, color='k')
    ax2.minorticks_on()

    pl.show()
#    sys.exit()
