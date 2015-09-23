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
bhextractor_plot_catalogue.py

Construct catalogues and principal component analysis for NR BBH waveforms

This script simply builds and plots a catalogue (defined by the user in this
script)

"""

import sys, os
from bhex_utils import bhex_wavedata as bwave
from bhex_utils import bhex_pca as bpca
#from pycbc.psd import aLIGOZeroDetHighPower
import pycbc.types
import numpy as np
import scipy.optimize
from matplotlib import pyplot as pl

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
            peakidx*SI_deltaT*(scale_ratio-1)

    resamp_amp = np.interp(time_axis, interp_times, amp)
    resamp_phase = np.interp(time_axis, interp_times, phase)

    return resamp_amp*np.exp(1j*resamp_phase)


def mismatch(total_mass, tmplt_wave_data, event_wave_data, delta_t=1./512):
    """
    Compute mismatch (1-match) between the tmplt wave and the event wave, given
    the total mass.  Uses event_wave and psd which are defined globally in the
    script.
    """

    # Rescale the template to this total mass
    tmplt = scale_wave(tmplt_wave_data, total_mass)

    # Convert the wave to a pycbc timeseries object
    hp = pycbc.types.TimeSeries(np.real(tmplt[:]), delta_t=SI_deltaT)
    hc = pycbc.types.TimeSeries(-1*np.imag(tmplt[:]), delta_t=SI_deltaT)

    event_wave = pycbc.types.TimeSeries(event_wave_data, delta_t=SI_deltaT)

    match, _ = pycbc.filter.match(hp, event_wave, psd=None,
            low_frequency_cutoff=f_min)

    return 1-match

def whiten_wave(wave, psd):
    """
    Whiten the complex waveform wave and return a wave with unit-hrss in each
    polarisation
    """

    hp = pycbc.types.TimeSeries(np.real(wave[:]), delta_t=SI_deltaT)
    hc = pycbc.types.TimeSeries(-1*np.imag(wave[:]), delta_t=SI_deltaT)

    # XXX Should have some checks here on sample rate / datalen but everything uses
    # 4s, 1024 Hz so should be ok...for now

    Hp = hp.to_frequencyseries()
    Hp.data /= np.sqrt(psd[:,1])
    Hc = hc.to_frequencyseries()
    Hc.data /= np.sqrt(psd[:,1])

    hp_white = Hp.to_timeseries()
    hp_white /= pycbc.filter.sigma(hp_white)
    hc_white = Hc.to_timeseries()
    hc_white /= pycbc.filter.sigma(hc_white)

    return hp_white - 1j*hc_white

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Some useful info:
#

#valid_series = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
#        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series"
#        "TP-series"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

SI_deltaT = 1./1024
SI_datalen= 4.0

f_min = 30.0

# Initial guess at the mass
mass_guess = 100 + 100*np.random.random()

init_total_mass = 100.  # Select waveforms which can go down to at least this mass
                        # (and generate SI catalogue at this mass)
distance=1. # Mpc

#
# --- Catalogue Definition
#
series_names = ['HR-series']#, 'HRq-series', 'RO3-series'] # (see above for valid choices)

bounds=None
#bounds=dict()
#bounds['q'] = [1, 1]
#bounds['a1'] = [0,0]
#bounds['a2'] = [0,1]

#
# --- Reconstruction data
#

event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed")
h1_wave_samples = np.loadtxt(os.path.join(event_file_dir, "waveforms/waveform_H1_1000.dat"))
l1_wave_samples = np.loadtxt(os.path.join(event_file_dir, "waveforms/waveform_L1_1000.dat"))

geo_wave_samples = np.loadtxt(os.path.join(event_file_dir, "geo_waveforms/waveform_geo_1000.dat"))

h1_psd = np.loadtxt(os.path.join(event_file_dir, \
        "lalinferencenest-0-H1-1126259462.39-14.datH1-PSD.dat"))
l1_psd = np.loadtxt(os.path.join(event_file_dir, \
        "lalinferencenest-0-L1-1126259462.39-26.datL1-PSD.dat"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

#
# --- Generate initial catalogue
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
simulations = \
        bwave.simulation_details(series_names=series_names, param_bounds=bounds,
                Mmin30Hz=init_total_mass)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
        trunc_time=True)

time_axis = np.arange(0, SI_datalen, SI_deltaT)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter Estimation
#

#
# --- Minimise mismatch over waveforms and total mass
#

# Loop over waves

geo_matches = np.zeros(shape=(simulations.nsimulations, len(geo_wave_samples)))
geo_masses  = np.zeros(shape=(simulations.nsimulations, len(geo_wave_samples)))

h1_matches = np.zeros(shape=(simulations.nsimulations, len(h1_wave_samples)))
h1_masses  = np.zeros(shape=(simulations.nsimulations, len(h1_wave_samples)))
l1_matches = np.zeros(shape=(simulations.nsimulations, len(l1_wave_samples)))
l1_masses  = np.zeros(shape=(simulations.nsimulations, len(l1_wave_samples)))

# XXX: Hack to use median, rather than sampled waveforms
#geo_wave_samples = [np.median(geo_wave_samples, axis=0)]

#h1_wave_samples = [np.median(h1_wave_samples, axis=0)]
#l1_wave_samples = [np.median(l1_wave_samples, axis=0)]


for w, wave in enumerate(catalogue.SIComplexTimeSeries):

    print "________________________________"
    print "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w, simulations.nsimulations)

    #
    # --- Whiten the catalogue waveform
    #
    h1_wave = whiten_wave(wave, h1_psd)
    l1_wave = whiten_wave(wave, l1_psd)


    # Find best-fitting mass (in terms of match)
    print "Optimising for total mass for each sampled waveform..."

    for s, (geo_sample, h1_sample, l1_sample) in enumerate(zip(geo_wave_samples,
        h1_wave_samples, l1_wave_samples)):


        print "Evaluating waveform %d of %d"%( s, len(geo_wave_samples) )

        geo_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(wave,
            geo_sample, SI_deltaT), full_output=True, retall=True, disp=True)

#        h1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(h1_wave,
#            h1_sample, SI_deltaT), full_output=True, retall=True, disp=True)

#        l1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(l1_wave,
#            l1_sample, SI_deltaT), full_output=True, retall=True, disp=True)

        geo_matches[w, s] = 1-geo_result[1]
        geo_masses[w, s] = geo_result[0][0]

#        h1_matches[w, s] = 1-h1_result[1]
#        h1_masses[w, s] = h1_result[0][0]

#        l1_matches[w, s] = 1-l1_result[1]
#        l1_masses[w, s] = l1_result[0][0]

    geo_bestidx=np.argmax(geo_matches[w, :])
    h1_bestidx=np.argmax(h1_matches[w, :])
    l1_bestidx=np.argmax(l1_matches[w, :])

    print "geo: Best matching mass [match]: %.2f [%.2f]"%(
            geo_masses[w,geo_bestidx], max(geo_matches[w,:]))
    print "H1: Best matching mass [match]: %.2f [%.2f]"%(
            h1_masses[w,h1_bestidx], max(h1_matches[w,:]))
    print "L1: Best matching mass [match]: %.2f [%.2f]"%(
            l1_masses[w,l1_bestidx], max(l1_matches[w,:]))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting



