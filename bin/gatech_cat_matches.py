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
burst_nr_match.py

Compute matches between burst reconstruction and NR waveforms
"""

import sys, os
from optparse import OptionParser
import ConfigParser
import subprocess
import cPickle as pickle

import numpy as np
import scipy.optimize
import timeit

import lal
import pycbc.types
from pycbc.waveform import get_td_waveform
import pycbc.filter

import burst_nr_utils as bnru
from bhex_utils import bhex_wavedata as bwave

from matplotlib import pyplot as pl

__author__ = "James Clark <james.clark@ligo.org>"
git_version_id = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
__version__ = "git id %s" % git_version_id

gpsnow = subprocess.check_output(['lalapps_tconvert', 'now']).strip()
__date__ = subprocess.check_output(['lalapps_tconvert', gpsnow]).strip()

def scale_approx(wave, target_total_mass, init_total_mass):
    """ 
    Scale the waveform to total_mass.  Assumes the waveform is initially
    generated at init_total_mass defined in this script.
    """
    amp = abs(wave.data[:])
    phase = bnru.phase_of(wave.data[:])

    scale_ratio = target_total_mass / init_total_mass
    amp *= scale_ratio

    peakidx = np.argmax(amp)

    interp_times = scale_ratio * wave.sample_times.data[:]

    resamp_amp = np.interp(wave.sample_times.data[:], interp_times, amp)
    resamp_phase = np.interp(wave.sample_times.data[:], interp_times, phase)
 
    return resamp_amp, resamp_phase


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse input

#
# --- Catalogue Definition
#
bounds = dict()
bounds['q'] = [1,1]
bounds['a1'] = [0.0, 0.0]
bounds['a2'] = [0.0, 0.0]

#
# --- Plotting options
#
nMassPoints = 5
maxMass = 500.0

#
# --- Time Series Config
#
deltaT = 1./8192
datalen = 4.0

#
# --- Noise Spectrum
#
asd_file = \
        "/home/jclark/Projects/bhextractor/data/noise_curves/early_aligo.dat"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

init_total_mass = 100   # Generate a catalogue at this mass; shouldn't matter,
                        # we rescale anyway

distance=100. # Mpc (doesn't really matter)


#
# --- Generate initial catalogue
#
print >> sys.stdout,  '~~~~~~~~~~~~~~~~~~~~~'
print >> sys.stdout,  'Selecting Simulations'
print >> sys.stdout,  ''
then = timeit.time.time()
simulations = \
        bwave.simulation_details(param_bounds=bounds)

sys.exit()

print >> sys.stdout,  '~~~~~~~~~~~~~~~~~~~~~'
print >> sys.stdout,  'Building NR catalogue'
print >> sys.stdout,  ''
catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
        SI_deltaT=deltaT, SI_datalen=datalen, distance=distance,
        trunc_time=False)
now = timeit.time.time()
print >> sys.stdout,  "...catalogue construction took %.1f..."%(now-then)

asd_data = np.loadtxt(asd_file)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Match Calculations
#

# For each waveform in the catalogue:
#   1) GeIn [204]: nerate approximant(s) @ 5 mass scales between the min allowed mass / max
#      mass
#   2) Compute matches at the 5 mass scales
#   3) Append the 5 match values (per approximant) to the
#      simulations.simulations
#
#   By the end then, simulations.simulations is a list of dictionaries where
#   each dictionary is 1 GAtech waveform with all physical attributes, as well
#   as matches at 5 mass scales with some selection of approximants


# Loop over waves in NR catalogue
for w, wave in enumerate(catalogue.SIComplexTimeSeries):

    print >> sys.stdout,  "________________________________"
    print >> sys.stdout,  "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w+1, simulations.nsimulations)

    # Set up the Masses we're going to study
    masses = np.linspace(simulations.simulations[w]['Mmin30Hz'], maxMass,
            nMassPoints)

    # matches is going to be a list of tuples: (mass, match)
    matches = []

    # Get the NR (plus) wave and put it in a pycbc TimeSeries object
    hplus_NR = pycbc.types.TimeSeries(np.real(wave), delta_t=deltaT)
    #hplus_NR.data = bwave.window_wave(hplus_NR.data)


    # --- Generate SEOBNR at init_total_mass and rescale later

    # Extract physical parameters
    mass1, mass2 = bwave.component_masses(init_total_mass, simulations.simulations[w]['q'])
    spin1z = simulations.simulations[w]['a1z']
    spin2z = simulations.simulations[w]['a2z']

    # ~~~~~~~~~~~~~~~~
    # SEOBNR
    hplus_SEOBNR, hcross_SEOBNR = get_td_waveform(approximant="SEOBNRv2",
            distance=distance,
            mass1=mass1,
            mass2=mass2,
            spin1z=spin1z,
            spin2z=spin2z,
            f_lower=10.0,
            delta_t=deltaT)

    # (1-sided) Planck window for smooth FFT
    hplus_SEOBNR.data = bwave.window_wave(hplus_SEOBNR.data)
    hcross_SEOBNR.data = bwave.window_wave(hcross_SEOBNR.data)

    # divide out the spherical harmonic (2,2) amplitude (this is just for nice
    # plots / sanity - it does not affect match)
    sY22 = lal.SpinWeightedSphericalHarmonic(0,0,-2,2,2)
    hplus_SEOBNR.data /= np.real(sY22)
#   hcross_SEOBNR.data /= np.imag(sY22)
#
    # Make the timeseries consistent lengths
    tlen = max(len(hplus_NR), len(hplus_SEOBNR))
    hplus_SEOBNR.resize(tlen)
    hcross_SEOBNR.resize(tlen)
    hplus_NR.resize(tlen)

    # ~~~~~~~~~~~~~~~~


    # Interpolate the ASD to the waveform frequencies (this is convenient so that we
    # end up with a PSD which overs all frequencies for use in the match calculation
    # later
    asd = np.interp(hplus_NR.to_frequencyseries().sample_frequencies,
            asd_data[:,0], asd_data[:,1])

    # Now insert into a pycbc frequency series so we can use pycbc.filter.match()
    # later
    noise_psd = pycbc.types.FrequencySeries(asd**2, delta_f =
            hplus_NR.to_frequencyseries().delta_f)

    for m,mass in enumerate(masses):

        # --- Generate the NR waveform at this mass

        # Scale the NR waveform to the mass we want
        amp, phase = bnru.scale_wave(hplus_NR, mass, init_total_mass)

        hplus_NR_new = pycbc.types.TimeSeries(np.real(amp*np.exp(1j*phase)),
                delta_t=deltaT)

        # Scale the SEOBNR waveform to this mass
        amp, phase = scale_approx(hplus_SEOBNR, mass, init_total_mass)

        hplus_SEOBNR_new = pycbc.types.TimeSeries(np.real(amp*np.exp(1j*phase)),
                delta_t=deltaT)
        hplus_SEOBNR_new.data = bwave.window_wave(hplus_SEOBNR_new.data)

        match, _ = pycbc.filter.match(hplus_SEOBNR_new, hplus_NR_new,
                low_frequency_cutoff=30.0, psd=noise_psd)

        print "~~~~~~~~~~~~~~~~~~~~~~~"
        print "Mass, mismatch (%)"
        print mass, 100*(1-match)

#
#       pl.figure()
#       pl.plot(hplus_NR.sample_times-hplus_NR.sample_times[np.argmax(hplus_NR)],
#               hplus_NR, label='NR')
#       pl.plot(hplus_SEOBNR.sample_times-hplus_SEOBNR.sample_times[np.argmax(hplus_SEOBNR)],
#               hplus_SEOBNR, label='SEOBNR')
#       pl.xlim(-0.1, 0.05)
 
        pl.figure() 
        pl.loglog(hplus_NR_new.to_frequencyseries().sample_frequencies,
                abs(hplus_NR_new.to_frequencyseries()), label='NR')
 
        pl.loglog(hplus_SEOBNR_new.to_frequencyseries().sample_frequencies,
                abs(hplus_SEOBNR_new.to_frequencyseries()), label='SEOBNR')
 
        pl.axvline(30, color='r')
 
    pl.show()
#
    sys.exit()






