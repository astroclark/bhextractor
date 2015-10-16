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
burst_nr_utils.py

Utilities for matching burst reconstructions to NR waveforms
"""

import sys
from optparse import OptionParser
import ConfigParser
import subprocess

import numpy as np

import lal
import pycbc.filter
import pycbc.types

__author__ = "James Clark <james.clark@ligo.org>"
git_version_id = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
__version__ = "git id %s" % git_version_id

gpsnow = subprocess.check_output(['lalapps_tconvert', 'now']).strip()
__date__ = subprocess.check_output(['lalapps_tconvert', gpsnow]).strip()

def phase_of(z):
    return np.unwrap(np.angle(z))

def scale_wave(wave, target_total_mass, init_total_mass):
    """
    Scale the waveform to total_mass.  Assumes the waveform is initially
    generated at init_total_mass defined in this script.
    """
    amp = abs(wave.data[:])
    phase = phase_of(wave.data[:])

    scale_ratio = target_total_mass / init_total_mass
    amp *= scale_ratio

    peakidx = np.argmax(amp)

    interp_times = scale_ratio * wave.sample_times.data[:] - \
            peakidx*wave.delta_t*(scale_ratio-1)

    resamp_amp = np.interp(wave.sample_times.data[:], interp_times, amp)
    resamp_phase = np.interp(wave.sample_times.data[:], interp_times, phase)
 
    return resamp_amp, resamp_phase

def extract_wave(inwave, datalen=4.0, sample_rate = 4096):
    extract_len = 0.5 # retain this many seconds of reconstruction
    delta = 0.15 # center the retained data on the peak of the waveform minus
                 # this many seconds
    peakidx = np.argmax(abs(inwave)) - delta*sample_rate
    nsamp = extract_len * sample_rate

    extracted = inwave[int(peakidx-0.5*nsamp): int(peakidx+0.5*nsamp)]

    win = lal.CreateTukeyREAL8Window(len(extracted), 0.1)
    extracted *= win.data.data

    output = np.zeros(datalen*sample_rate)
    output[0.5*datalen*sample_rate-0.5*nsamp:
            0.5*datalen*sample_rate+0.5*nsamp] = np.copy(extracted)

    return output

def mismatch(target_total_mass, init_total_mass, mass_bounds, tmplt_wave_data,
        event_wave_data, asd=None, delta_t=1./512, delta_f=0.25,
        ifo_response=False, f_min=30.0):
    """
    Compute mismatch (1-match) between the tmplt wave and the event wave, given
    the total mass.  Uses event_wave and psd which are defined globally in the
    script.
    """
    min_mass, max_mass = mass_bounds

    if (target_total_mass >= min_mass) and (target_total_mass <= max_mass):

        # Convert the real part of the wave to a pycbc timeseries object
        init_tmplt = pycbc.types.TimeSeries(np.real(tmplt_wave_data[:]),
                delta_t=delta_t)

        # Rescale the template to this total mass
        tmplt_amp, tmplt_phase = scale_wave(init_tmplt, target_total_mass, init_total_mass)

        tmplt = pycbc.types.TimeSeries(np.real(tmplt_amp*np.exp(1j*tmplt_phase)),
                delta_t=delta_t)

        if ifo_response and asd is not None:
            # Whiten the template
            Tmplt = tmplt.to_frequencyseries()
            Tmplt.data /= asd

            # IFFT (just simplifies the code below) 
            tmplt = Tmplt.to_timeseries()


        if asd is not None and not ifo_response:
            psd = pycbc.types.FrequencySeries(asd**2, delta_f=delta_f)
        else:
            psd = None

        # Put the reconstruction data in a TimeSeries
        event_wave = pycbc.types.TimeSeries(event_wave_data, delta_t=delta_t)

        try:
            match, _ = pycbc.filter.match(tmplt, event_wave, psd=psd,
                    low_frequency_cutoff=f_min)
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

class configuration:

    def __init__(self, configparser):

        self.sample_rate=configparser.getint('analysis', 'sample_rate')
        self.deltaT=1./self.sample_rate
        self.datalen=configparser.getfloat('analysis', 'datalen')
        self.f_min=configparser.getfloat('analysis', 'f_min')
        self.algorithm=configparser.get('analysis', 'algorithm')

        self.nsampls=configparser.getint('parameters', 'nsampls')
        self.mass_guess=configparser.getfloat('parameters', 'mass_guess')
        self.min_chirp_mass=configparser.getfloat('parameters', 'min_chirp_mass')
        self.max_chirp_mass=configparser.getfloat('parameters', 'max_chirp_mass')

        self.reconstruction=configparser.get('paths', 'reconstruction')
        self.spectral_estimate=configparser.get('paths', 'spectral-estimate')

################################################################################
# MAIN

def main():
    print >> sys.stdout, sys.argv[0]
    print >> sys.stdout, __version__
    print >> sys.stdout, __date__
    return 0

#
# End definitions
#

if __name__ == "__main__":
    result = main()
