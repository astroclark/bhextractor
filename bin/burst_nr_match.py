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

from bhex_utils import bhex_wavedata as bwave
import burst_nr_utils as bnru

__author__ = "James Clark <james.clark@ligo.org>"
git_version_id = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
__version__ = "git id %s" % git_version_id

gpsnow = subprocess.check_output(['lalapps_tconvert', 'now']).strip()
__date__ = subprocess.check_output(['lalapps_tconvert', gpsnow]).strip()

def parser():

    # --- Command line input
    parser = OptionParser()
    parser.add_option("-t", "--user-tag", default="TEST", type=str)
    parser.add_option("-o", "--output-dir", type=str, default=None)
    parser.add_option("-a", "--algorithm", type=str, default=None)

    (opts,args) = parser.parse_args()

    if len(args)==0:
        print >> sys.stderr, "ERROR: require config file"
        sys.exit()

    algorithms=["BW", "CWB"]
    if opts.algorithm is not None and opts.algorithm not in algorithms:
        print >> sys.stderr, "ERROR: algorithm %s not recognised"%opts.algorithm
        print >> sys.stderr, "must be in ", algorithms
        sys.exit(-1)


    # --- Read config file
    configparser = ConfigParser.ConfigParser()
    configparser.read(args[0])

    # --- Where did the reconstruction come from?
    if opts.algorithm is not None:
        # override from the commandline
        configparser.set('analysis','algorithm',opts.algorithm)

    # Check algorithm is defined (might have been in the ini file)
    if configparser.has_option('analysis', 'algorithm'):
        alg = configparser.get('analysis','algorithm')
        if alg not in algorithms:
            print >> sys.stderr, "ERROR: algorithm %s not recognised"%alg
            print >> sys.stderr, "must be in ", algorithms
            sys.exit(-1)
    else:
        print >> sys.stderr, "ERROR: algorithm not defined"
        print >> sys.stderr, "must be in ", algorithms
        print >> sys.stderr, "and defined in [analysis] of ini or with --algorithm"
        sys.exit(-1)


    return opts, args, configparser

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse input

opts, args, cp = bnru.parser()
config = bnru.configuration(cp)

#
# --- Catalogue Definition
#
bounds = dict()
bounds['Mchirpmin30Hz'] = [-np.inf, config.min_chirp_mass]


#
# --- Reconstruction data
#
print >> sys.stdout,  "Loading data"
reconstruction_data = np.loadtxt(config.reconstruction)
asd_data = np.loadtxt(config.spectral_estimate)

# If BayesWave, select the user-specified number of samples for which we will
# compute matches (useful for speed / development work)
if config.algorithm=='BW':
    print 'reducing sample size'
    idx = np.random.random_integers(low=0, high=len(reconstruction_data)-1,
            size=config.nsampls)
    reconstruction_data = reconstruction_data[idx]
elif config.algorith=='CWB':
    # Make it iterable so that the BW/CWB codes can be consistent
    reconstruction_data = [reconstruction_data]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

init_total_mass = 100   # Generate a catalogue at this mass; shouldn't matter,
                        # we rescale anyway
distance=1. # Mpc

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
        SI_deltaT=config.deltaT, SI_datalen=config.datalen, distance=distance,
        trunc_time=False)
now = timeit.time.time()
print >> sys.stdout,  "...catalogue construction took %.1f..."%(now-then)


# Useful time/freq samples
time_axis = np.arange(0, config.datalen, config.deltaT)
freq_axis = np.arange(0, catalogue.SI_flen*catalogue.SI_deltaF,
        catalogue.SI_deltaF)

# Interpolate the ASD to the waveform frequencies (this is convenient so that we
# end up with a PSD which overs all frequencies for use in the match calculation
# later - In practice, this will really just pad out the spectrum at low
# frequencies)
asd = np.interp(freq_axis, asd_data[:,0], asd_data[:,1])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter Estimation
#

#
# --- Minimise mismatch over waveforms and total mass
#

# Preallocate
matches = np.zeros(shape=(simulations.nsimulations, len(reconstruction_data)))
masses  = np.zeros(shape=(simulations.nsimulations, len(reconstruction_data)))

# Loop over waves in NR catalogue
for w, wave in enumerate(catalogue.SIComplexTimeSeries):

    print >> sys.stdout,  "________________________________"
    print >> sys.stdout,  "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w+1, simulations.nsimulations)


    # Find best-fitting mass (in terms of match)
    print >> sys.stdout,  "Optimising for total mass for each sampled waveform..."

    # Find min/max allowable mass to which we can scale the waveform
    min_mass = bnru.mtot_from_mchirp(config.min_chirp_mass,
            simulations.simulations[w]['q'])
    max_mass = bnru.mtot_from_mchirp(config.max_chirp_mass,
            simulations.simulations[w]['q'])

    for s, sample in enumerate(reconstruction_data):


        print >> sys.stdout,  '-----------------------------'
        print >> sys.stdout,  "Evaluating sample waveform %d of %d"%( s,
                len(reconstruction_data) )
        print >> sys.stdout,  " (NR: %s [%d/%d])"%( simulations.simulations[w]['wavename'], w+1,
                simulations.nsimulations)

        #
        # Optimise match over total mass
        #

        then = timeit.time.time()
        ifo_response=True
        result = scipy.optimize.fmin(bnru.mismatch, x0=config.mass_guess, args=(
            init_total_mass, [min_mass, max_mass], wave, sample, asd, config.deltaT,
            catalogue.SI_deltaF, ifo_response), full_output=True,
            retall=True, disp=False)
        now = timeit.time.time()
        print >> sys.stdout,  "...mass optimisation took %.3f sec..."%(now-then)

        matches[w, s] = 1-result[1]
        masses[w, s] = result[0][0]

        print >> sys.stdout,  "Best matching mass [match]: %.2f [%.2f]"%(
                masses[w,s], matches[w,s])

    bestidx=np.argmax(matches[w, :])

    print >> sys.stdout,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print >> sys.stdout,  " Best matching mass [match]: %.2f [%.2f]"%(
            masses[w,bestidx], max(matches[w,:]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dump data

filename=opts.user_tag+'_'+config.algorithm+'_'+gpsnow+'.pickle'

#np.savez(filename, matches=matches, masses=masses)

# Dump results and configuration to pickle
pickle.dump([matches, masses, config, simulations, __author__, __version__,
    __date__], open(filename, "wb"))



