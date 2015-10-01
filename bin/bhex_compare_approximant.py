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
from bhex_utils import bhex_wavedata as bwave

import lal
import pycbc.types
from pycbc.waveform import get_td_waveform
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower
from pylal import spawaveform

from matplotlib import pyplot as pl

def component_masses(total_mass, mass_ratio):
    """
    Return m1 and m2, given total mass and mass ratio (m1/m2)

    m1, m2 = component_masses(total_mass, mass_ratio)
    """

    m1 = mass_ratio * total_mass / (1.0 + mass_ratio)
    m2 = total_mass - m1

    return m1, m2


def cartesian_spins(spin_magnitude, spin_theta):
    """
    Compute cartesian spin components.  Only does z-component for now
    """

    if np.isnan(spin_magnitude) or np.isnan(spin_theta):
        return 0.0
    else:
        spin_magnitude * np.cos(spin_theta * np.pi / 180.0)
    return spin_z

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User Options (mass configurations, noise curves etc)

f_low = 10. # Min freq for waveform generation
f_min = 30. # Min freq for match calculation

SI_deltaT = 1./1024
SI_datalen= 4.0

max_minMass = 65 # The maximum value of the smallest mass we want to generate
total_mass = 72. # The mass at which we generate the catalogue

distance=1. # Distance in Mpc to the source


# Make catalogue
#bounds=dict()
#bounds['a1'] = [0, 0]
#bounds['a2'] = [0, 0]
#bounds['q'] = [1, 1] 

# Alternatively (easier):
bounds = bwave.bounds_dict("NonSpinning")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Select & Generate NR waveforms

# Pick out the simulations we want:

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''

simulations = bwave.simulation_details( param_bounds=bounds,
        Mmin30Hz=max_minMass)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
NR_catalogue = bwave.waveform_catalogue(simulations, ref_mass=total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
        trunc_time=False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate EOBNR

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'Beginning Match Calculations'
print ''

matches = np.zeros(simulations.nsimulations)
mass_ratios = np.zeros(simulations.nsimulations)
effective_spins = np.zeros(simulations.nsimulations)

spin_magnitudes=[]
spin_angles=[]

# Retrieve the mass ratio:
#f1, ax1 = pl.subplots()
#f2, ax2 = pl.subplots()
for s,simulation in enumerate(simulations.simulations):
    
    print 'Matching waveform %d of %d'%(s, simulations.nsimulations)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # --- Condition the NR waveform ---

    # Get the NR (plus) wave and put it in a pycbc TimeSeries object
    hplus_NR = pycbc.types.TimeSeries(np.real(
        NR_catalogue.SIComplexTimeSeries[s,:]),
        delta_t=SI_deltaT)
    hplus_NR.data = bwave.window_wave(hplus_NR.data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # --- Generate approximant with this NR waveform's parameters ---
    #
    mass_ratios[s] = simulation['q']
    spin_magnitudes.append((simulation['a1'], simulation['a2']))
    spin_angles.append((simulation['th1L'], simulation['th2L']))

    mass1, mass2 = component_masses(total_mass, mass_ratios[s])

    spin1z = cartesian_spins(simulation['a1'], simulation['th1L'])
    spin2z = cartesian_spins(simulation['a2'], simulation['th2L'])

    effective_spins[s] = spawaveform.computechi(mass1, mass2, spin1z, spin2z)


    print "Generating EOBNR"
    hplus_EOBNR, _ = get_td_waveform(approximant="SEOBNRv2",
            distance=distance,
            mass1=mass1,
            mass2=mass2,
            spin1z=spin1z,
            spin2z=spin2z,
            f_lower=f_low,
            delta_t=SI_deltaT)

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

    matches[s], snr_max_loc = pycbc.filter.match(hplus_EOBNR, hplus_NR,
            low_frequency_cutoff=f_min, psd=psd)
    # match() returns the index for peak snr as well as the match

#   NR_max_loc = np.argmax(hplus_NR)
#
#   ax1.plot(hplus_EOBNR.sample_times, hplus_EOBNR, linestyle='-')
#   ax1.plot(hplus_NR.sample_times - hplus_NR.sample_times[NR_max_loc],
#           -1*hplus_NR, linestyle='--')
#   ax1.minorticks_on()
#
#   ax2.loglog(hplus_EOBNR.to_frequencyseries().sample_frequencies,
#           abs(hplus_EOBNR.to_frequencyseries()), label='SEOBNR')
#   ax2.loglog(hplus_NR.to_frequencyseries().sample_frequencies,
#           abs(hplus_NR.to_frequencyseries()), linestyle='--', label='NR')
#   ax2.axvline(f_min, color='r')
#   ax2.minorticks_on()
#
#   ax2.legend()
#
#   ax2.set_title('Match=%.2f'%matches[s])

#pl.show()
#sys.exit()

