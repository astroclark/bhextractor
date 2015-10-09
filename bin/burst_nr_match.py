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
bhex_bw_match.py

Compute matches between BW reconstruction and NR waveforms
"""

import sys, os
from optparse import OptionParser
import ConfigParser

import numpy as np
import scipy.optimize
import scipy.stats
import scipy.interpolate
import timeit
from matplotlib import pyplot as pl

from bhex_utils import bhex_wavedata as bwave
import lal
import pycbc.filter
from pylal import spawaveform
import pycbc.types

def scale_wave(wave, total_mass):
    """
    Scale the waveform to total_mass.  Assumes the waveform is initially
    generated at init_total_mass defined in this script.
    """
    amp = abs(wave)
    phase = bwave.phase_of(wave)

    scale_ratio = total_mass / init_total_mass
    amp *= scale_ratio

    peakidx = np.argmax(amp)
    interp_times = scale_ratio * time_axis - \
            peakidx*config.deltaT*(scale_ratio-1)

    resamp_amp = np.interp(time_axis, interp_times, amp)
    resamp_phase = np.interp(time_axis, interp_times, phase)
 
    return resamp_amp*np.exp(1j*resamp_phase)


def mismatch(total_mass, mass_bounds, tmplt_wave_data, event_wave_data, asd=None,
        delta_t=1./512, delta_f=0.25, ifo_response=False):
    """
    Compute mismatch (1-match) between the tmplt wave and the event wave, given
    the total mass.  Uses event_wave and psd which are defined globally in the
    script.
    """
    min_mass, max_mass = mass_bounds

    if (total_mass >= min_mass) and (total_mass <= max_mass):

        # Rescale the template to this total mass
        tmplt = scale_wave(tmplt_wave_data, total_mass)

        # Convert the wave to a pycbc timeseries object
        hp = pycbc.types.TimeSeries(np.real(tmplt[:]), delta_t=config.deltaT)

        if ifo_response and asd is not None:
            # Whiten the template
            Hp = hp.to_frequencyseries()
            Hp.data /= asd

            # IFFT (just simplifies the code below) 
            hp = Hp.to_timeseries()


        if asd is not None and not ifo_response:
            psd = pycbc.types.FrequencySeries(asd**2, delta_f=delta_f)
        else:
            psd = None

        # Put the reconstruction data in a TimeSeries
        event_wave = pycbc.types.TimeSeries(event_wave_data, delta_t=config.deltaT)

        try:
            match, _ = pycbc.filter.match(hp, event_wave, psd=psd,
                    low_frequency_cutoff=config.f_min)
        except ZeroDivisionError:
            match = np.nan

        return 1-match

    else:

        return 1.


def mtot_from_mchirp(mc, q):
    eta = q/(1+q)**2.0
    return mc * eta**(-3./5)


def parser():

    # --- Command line input
    parser = OptionParser()
    parser.add_option("-t", "--user-tag", default="TEST", type=str)
    parser.add_option("-o", "--output-dir", type=str, default=None)
    parser.add_option("-a", "--algorithm", type=str, default="BW")

    (opts,args) = parser.parse_args()

    if len(args)==0:
        print >> sys.stderr, "ERROR: require config file"
        sys.exit()

    algorithms=["BW", "CWB"]
    if opts.algorithm not in algorithms:
        print >> sys.stderr, "ERROR: algorithm %s not recognised"%opts.algorithm
        print >> sys.stderr, "must be in ", algorithms
        sys.exit()


    # --- Read config file
    configparser = ConfigParser.ConfigParser()
    configparser.read(args[0])

    return opts, args, configparser

class configuration:

    def __init__(self, configparser)

    if configparser.has_section('analysis'):

        deltaT=configparser.get_float('analysis', 'deltaT')
        datalen=configparser.get_float('analysis', 'datalen')
        f_min=configparser.get_float('analysis', 'f_min')
        nsampls=configparser.get_int('analysis', 'nsampls')

        mass_guess=configparser.get_float('analysis', 'mass_guess')
        min_chirp_mass=configparser.get_float('analysis', 'min_chirp_mass')
        max_chirp_mass=configparser.get_float('analysis', 'max_chirp_mass')

    else:

        print >> sys.stderr, "ERROR: ini file requires [analysis] section"
        sys.exit()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT


opts, args, cp = parser()
config = 




use_waves=['geo', 'h1', 'l1']
#use_waves=['geo'] # XXX: can only do GEO for software injections right now

#
# --- Catalogue Definition
#
bounds = bwave.bounds_dict(usertag)
bounds = dict()
bounds['Mchirpmin30Hz'] = [-np.inf, min_chirp_mass]


sys.exit()

#
# --- Reconstruction data
#
print >> sys.stdout,  "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/injections")
       # "data/observed")

# Geocentric waveform:
if 'geo' in use_waves:
    geo_wave_samples = np.loadtxt(os.path.join(event_file_dir,
        "bw/G187068/G187068_geo_waveforms.dat"))
        #"bw/job_45906/waveforms/waveform_geo_1000.dat"))
        #"bw/geo_waveforms/waveform_geo_1000.dat"))
    # Downsample the number of posterior samples (useful for testing)
    geo_wave_samples=geo_wave_samples[np.random.random_integers(low=0,high=config.nsampls,size=config.nsampls),:]
else:
    geo_wave_samples=np.zeros(shape=(config.nsampls, config.datalen/config.deltaT))

# Detector responses:
if 'h1' in use_waves:
    h1_wave_samples = np.loadtxt(os.path.join(event_file_dir,
        "bw/G187068/G187068_H1.dat"))
        #"bw/waveforms/signal_recovered_whitened_waveform.dat.0"))
    h1_wave_samples=h1_wave_samples[np.random.random_integers(low=0,high=config.nsampls,size=config.nsampls),:]
else:
    h1_wave_samples=np.zeros(shape=(1000, config.datalen/config.deltaT))

if 'l1' in use_waves:
    l1_wave_samples = np.loadtxt(os.path.join(event_file_dir,
        "bw/G187068/G187068_L1.dat"))
        #"bw/waveforms/signal_recovered_whitened_waveform.dat.1"))
    l1_wave_samples=l1_wave_samples[np.random.random_integers(low=0,high=config.nsampls,size=config.nsampls),:]
else:
    l1_wave_samples=np.zeros(shape=(1000, config.datalen/config.deltaT))


# PSD estimates
# --- The event

#h1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/IFO0_asd.dat"))
#l1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/IFO1_asd.dat"))

# --- Software injection
#h1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/job_45906/IFO0_asd.dat"))
#l1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/job_45906/IFO1_asd.dat"))

# --- Hardware injection
h1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/job_45906/IFO0_asd.dat"))
l1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/job_45906/IFO1_asd.dat"))

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
# later)
h1_asd = np.interp(freq_axis, h1_bw_asd_data[:,0], h1_bw_asd_data[:,1])
l1_asd = np.interp(freq_axis, l1_bw_asd_data[:,0], l1_bw_asd_data[:,1])

mean_asd = np.sqrt(scipy.stats.hmean([h1_asd**2, l1_asd**2]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter Estimation
#

#
# --- Minimise mismatch over waveforms and total mass
#

# Preallocate
geo_matches = np.zeros(shape=(simulations.nsimulations, len(geo_wave_samples)))
geo_masses  = np.zeros(shape=(simulations.nsimulations, len(geo_wave_samples)))

h1_matches = np.zeros(shape=(simulations.nsimulations, len(h1_wave_samples)))
h1_masses  = np.zeros(shape=(simulations.nsimulations, len(h1_wave_samples)))

l1_matches = np.zeros(shape=(simulations.nsimulations, len(l1_wave_samples)))
l1_masses  = np.zeros(shape=(simulations.nsimulations, len(l1_wave_samples)))

# Loop over waves in NR catalogue
for w, wave in enumerate(catalogue.SIComplexTimeSeries):

    print >> sys.stdout,  "________________________________"
    print >> sys.stdout,  "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w+1, simulations.nsimulations)


    # Find best-fitting mass (in terms of match)
    print >> sys.stdout,  "Optimising for total mass for each sampled waveform..."

    # Find min/max allowable mass
    min_mass = mtot_from_mchirp(min_chirp_mass, simulations.simulations[w]['q'])
    max_mass = mtot_from_mchirp(max_chirp_mass, simulations.simulations[w]['q'])

    for s, (geo_sample, h1_sample, l1_sample) in enumerate(zip(geo_wave_samples,
        h1_wave_samples, l1_wave_samples)):


        print >> sys.stdout,  '-----------------------------'
        print >> sys.stdout,  "Evaluating sample waveform %d of %d"%( s, len(geo_wave_samples) )
        print >> sys.stdout,  " (NR: %s [%d/%d])"%( simulations.simulations[w]['wavename'], w+1,
                simulations.nsimulations)


        if 'geo' in use_waves:
            #
            # GEOCENTRIC WAVEFORM
            #

            then = timeit.time.time()
            ifo_response=False
            geo_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
                [min_mass, max_mass], wave, geo_sample, mean_asd, config.deltaT,
                catalogue.SI_deltaF, ifo_response), full_output=True,
                retall=True, disp=False)
            now = timeit.time.time()
            print >> sys.stdout,  "...mass optimisation took %.3f sec..."%(now-then)

            geo_matches[w, s] = 1-geo_result[1]
            geo_masses[w, s] = geo_result[0][0]

            print >> sys.stdout,  "geo: Best matching mass [match]: %.2f [%.2f]"%(
                    geo_masses[w,s], geo_matches[w,s])

        if 'h1' in use_waves:
            #
            # H1 RESPONSE
            #
            ifo_response=True
            then = timeit.time.time()
            h1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
                [min_mass, max_mass], wave, h1_sample, h1_asd, config.deltaT,
                catalogue.SI_deltaF, ifo_response), full_output=True,
                retall=True, disp=False)
            now = timeit.time.time()
            print >> sys.stdout,  "...mass optimisation took %.3f sec..."%(now-then)

            h1_matches[w, s] = 1-h1_result[1]
            h1_masses[w, s] = h1_result[0][0]

            print >> sys.stdout,  "h1: Best matching mass [match]: %.2f [%.2f]"%(
                    h1_masses[w,s], h1_matches[w,s])

        if 'l1' in use_waves:
            #
            # L1 RESPONSE
            #

            ifo_response=True
            then = timeit.time.time()
            l1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
                [min_mass, max_mass], wave, l1_sample, l1_asd, config.deltaT,
                catalogue.SI_deltaF, ifo_response), full_output=True,
                retall=True, disp=False)
            now = timeit.time.time()
            print >> sys.stdout,  "...mass optimisation took %.3f sec..."%(now-then)

            l1_matches[w, s] = 1-l1_result[1]
            l1_masses[w, s] = l1_result[0][0]

            print >> sys.stdout,  "l1: Best matching mass [match]: %.2f [%.2f]"%(
                    l1_masses[w,s], l1_matches[w,s])

    geo_bestidx=np.argmax(geo_matches[w, :])
    h1_bestidx=np.argmax(h1_matches[w, :])
    l1_bestidx=np.argmax(l1_matches[w, :])

    print >> sys.stdout,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print >> sys.stdout,  "Maximising over all waveforms:"
    if 'geo' in use_waves:
        print >> sys.stdout,  "geo: Best matching mass [match]: %.2f [%.2f]"%(
                geo_masses[w,geo_bestidx], max(geo_matches[w,:]))
    if 'h1' in use_waves:
        print >> sys.stdout,  "h1: Best matching mass [match]: %.2f [%.2f]"%(
                h1_masses[w,h1_bestidx], max(h1_matches[w,:]))
    if 'l1' in use_waves:
        print >> sys.stdout,  "l1: Best matching mass [match]: %.2f [%.2f]"%(
                l1_masses[w,l1_bestidx], max(l1_matches[w,:]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dump data
import uuid
filename=usertag+'_'+str(uuid.uuid4())
np.savez(filename, geo_matches=geo_matches, geo_masses=geo_masses,
        h1_matches=h1_matches, h1_masses=h1_masses,
        l1_matches=l1_matches, l1_masses=l1_masses)



