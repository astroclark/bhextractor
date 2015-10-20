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

    scaling_data = np.copy(wave.data[:])

    amp = abs(scaling_data)

    scale_ratio = target_total_mass / init_total_mass
    scaling_data *= scale_ratio

    peakidx = np.argmax(amp)

    interp_times = scale_ratio * wave.sample_times.data[:]

    resamp_wave = np.interp(wave.sample_times.data[:], interp_times,
            scaling_data)

    return resamp_wave


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse input

#
# --- Catalogue Definition
#
bounds = dict()
bounds['q'] = [2,2]
bounds['a1'] = [0.1, 0.2]
bounds['a2'] = [0.5, 0.7]
bounds['Mmin30Hz'] = [-np.inf,  60]

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

plot_snr = 8

#
# --- Generate initial catalogue
#
print >> sys.stdout,  '~~~~~~~~~~~~~~~~~~~~~'
print >> sys.stdout,  'Selecting Simulations'
print >> sys.stdout,  ''
then = timeit.time.time()
simulations = \
        bwave.simulation_details(param_bounds=bounds)

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
    hplus_NR.data[:] = bnru.taper(hplus_NR.data[:], delta_t=hplus_NR.delta_t)


    # Extract physical parameters
    mass1, mass2 = bwave.component_masses(init_total_mass, simulations.simulations[w]['q'])
    spin1z = simulations.simulations[w]['a1z']
    spin2z = simulations.simulations[w]['a2z']



    f, ax = pl.subplots(nrows = len(masses), figsize=(8,15))
    for m,mass in enumerate(masses):
        print 'mass: ', mass

        # --- Scale the NR waveform at this mass

        # Scale the NR waveform to the mass we want
        hplus_NR_new = pycbc.types.TimeSeries(bnru.scale_wave(hplus_NR, mass,
            init_total_mass), delta_t=deltaT)

        # --- Generate the SEOBNR waveform to this mass

        mass1, mass2 = bwave.component_masses(mass, simulations.simulations[w]['q'])
        hplus_SEOBNR, _ = get_td_waveform(approximant="SEOBNRv2",
                distance=distance,
                mass1=mass1,
                mass2=mass2,
                spin1z=spin1z,
                spin2z=spin2z,
                f_lower=10.0,
                delta_t=deltaT)

        hplus_SEOBNR.data = bnru.taper(hplus_SEOBNR.data,
                delta_t=hplus_SEOBNR.delta_t)

     
        # Make the timeseries consistent lengths
        tlen = max(len(hplus_NR_new), len(hplus_SEOBNR))
        hplus_SEOBNR.resize(tlen)
        hplus_NR_new.resize(tlen)

        # Interpolate the ASD to the waveform frequencies (this is convenient so that we
        # end up with a PSD which overs all frequencies for use in the match calculation
        # later
        asd = np.interp(hplus_NR_new.to_frequencyseries().sample_frequencies,
                asd_data[:,0], asd_data[:,1])


        # Now insert ASD into a pycbc frequency series so we can use
        # pycbc.filter.match() later
        noise_psd = pycbc.types.FrequencySeries(asd**2, delta_f =
                hplus_NR_new.to_frequencyseries().delta_f)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        Hf = abs(hplus_SEOBNR.to_frequencyseries())
        inband = noise_psd.sample_frequencies.data>30
        upp_bound = \
                noise_psd.sample_frequencies[inband][np.argwhere(Hf.data[inband]<1e-2*Hf[inband].max())[0]]
#        upp_bound = 0.5*1./deltaT

        match, _ = pycbc.filter.match(hplus_SEOBNR, hplus_NR_new,
                low_frequency_cutoff=30.0, psd=noise_psd,
                high_frequency_cutoff=upp_bound)

        print "~~~~~~~~~~~~~~~~~~~~~~~"
        print "Mass, mismatch (%)"
        print mass, 100*(1-match)

        # Normalise to unit SNR
        hplus_SEOBNR.data[:] /= pycbc.filter.sigma(hplus_SEOBNR,
                psd=noise_psd, low_frequency_cutoff=30, high_frequency_cutoff=upp_bound)

        hplus_NR_new.data[:] /= pycbc.filter.sigma(hplus_NR_new,
                psd=noise_psd, low_frequency_cutoff=30, high_frequency_cutoff=upp_bound)

        Hplus_SEOBNR = hplus_SEOBNR.to_frequencyseries()
        Hplus_NR_new = hplus_NR_new.to_frequencyseries()
 
        ax[m].loglog(Hplus_NR_new.sample_frequencies,
                   abs(Hplus_NR_new), label='NR')
 
        ax[m].loglog(Hplus_SEOBNR.sample_frequencies,
                   abs(Hplus_SEOBNR), label='SEOBNR')

        ax[m].loglog(noise_psd.sample_frequencies, np.sqrt(noise_psd),
                label='noise psd', color='k', linestyle='--')

        ax[m].set_title('M$_{\mathrm{tot}}$=%.2f M$_{\odot}$, mismatch=%.2f %%'%(
            mass, 100*(1-match)))

        ax[m].legend(loc='lower left')
 
        ax[m].axvline(30, color='r')
        ax[m].axvline(upp_bound, color='r')

        ax[m].set_xlabel('Frequency [Hz]')
        ax[m].set_ylabel('|H(f)| [arb units]')
        ax[m].set_ylim(1e-27, 1e-19)
 
    f.tight_layout()
    pl.show()
#
    sys.exit()






