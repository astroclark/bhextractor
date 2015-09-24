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
import scipy.stats
import timeit
from matplotlib import pyplot as pl

def scale_wave(wave, total_mass):
    """
    Scale the waveform to total_mass.  Assumes the waveform is initially
    generated at init_total_mass defined in this script.
    """
    amp = abs(wave)
    phase = bwave.phase_of(wave)

    scale_ratio = total_mass / init_total_mass
    #scale_ratio = init_total_mass / total_mass
    amp *= scale_ratio

    peakidx = np.argmax(amp)
    interp_times = scale_ratio * time_axis - \
            peakidx*SI_deltaT*(scale_ratio-1)

#   resamp_hplus = np.interp(time_axis, interp_times, np.real(wave))
#   resamp_hcross = np.interp(time_axis, interp_times, np.imag(wave))
#
#   return resamp_hplus + 1j*resamp_hcross

    #resamp_amp = np.interp(time_axis, interp_times, amp)
    #resamp_phase = np.interp(time_axis, interp_times, phase)

    resamp_amp = np.interp(time_axis, interp_times, amp)
    resamp_phase = np.interp(time_axis, interp_times, phase)

    return resamp_amp*np.exp(1j*resamp_phase)


def mismatch(total_mass, tmplt_wave_data, event_wave_data, asd=None,
        delta_t=1./512, delta_f=0.25):
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

    if asd is not None:
        psd = pycbc.types.FrequencySeries(asd**2, delta_f=delta_f)
    else:
        psd = None

    try:
        match, _ = pycbc.filter.match(hp, event_wave, psd=psd,
                low_frequency_cutoff=f_min)
    except ZeroDivisionError:
        match = np.nan

    return 1-match

def whiten_wave(wave, asd):
    """
    Whiten the complex waveform wave and return a wave with unit-hrss in each
    polarisation
    """

    hp = pycbc.types.TimeSeries(np.real(wave[:]), delta_t=SI_deltaT)
    hc = pycbc.types.TimeSeries(-1*np.imag(wave[:]), delta_t=SI_deltaT)

    # XXX Should have some checks here on sample rate / datalen but everything uses
    # 4s, 1024 Hz so should be ok...for now

    Hp = hp.to_frequencyseries()
    Hp.data /= asd
    Hc = hc.to_frequencyseries()
    Hc.data /= asd

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

f_min = 50.0

# Initial guess at the mass
mass_guess = 74.0#100 + 100*np.random.random()

init_total_mass = 100.  # Select waveforms which can go down to at least this mass
                        # (and generate SI catalogue at this mass)
distance=1. # Mpc

#
# --- Catalogue Definition
#
series_names = ['HRq-series']#, 'HRq-series', 'RO3-series'] # (see above for valid choices)

bounds=None
bounds=dict()
bounds['q'] = [1, 1]
#bounds['a1'] = [0,0]
#bounds['a2'] = [0,1]

# XXX SANITY CHECKING SCALING AND BEHAVIOUR OF NR WAVE wrt SEOBNR

#   series_names = ['HRq-series']#, 'HRq-series', 'RO3-series'] # (see above for valid choices)
#   bounds=dict()
#   bounds['q'] = [1, 1]
#   from pycbc.waveform import get_td_waveform
#   hplus_EOBNR, _ = get_td_waveform(approximant="SEOBNRv2",
#           distance=distance,
#           mass1=38,
#           mass2=38,
#           spin1z=0.0,
#           spin2z=0.0,
#           f_lower=10,
#           delta_t=1.0/1024)
#   emax=np.argmax(abs(hplus_EOBNR))
#
#   import lal
#   sY22 = lal.SpinWeightedSphericalHarmonic(0,0,-2,2,2)
#   hplus_EOBNR.data /= np.real(sY22)
#   hplus_EOBNR.data = bwave.window_wave(hplus_EOBNR.data)
#
#
#   simulations = \
#           bwave.simulation_details(series_names=series_names, param_bounds=bounds,
#                   Mmin30Hz=100.0)
#
#   catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
#           SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
#           trunc_time=False)
#
#   # Useful time/freq samples
#   wave = np.real(catalogue.SIComplexTimeSeries[0,:])
#
#   nmax = np.argmax(abs(wave))
#
#   time_axis = np.arange(0, SI_datalen, SI_deltaT)  
#
#   # rescale 
#   wave = scale_wave(wave, 76)
#
#   pl.figure()
#   pl.plot(hplus_EOBNR.sample_times - hplus_EOBNR.sample_times[emax], hplus_EOBNR,
#           label='SEOBNR')
#   pl.plot(time_axis-time_axis[nmax],wave, label='NR')
#
#   hplus_NR = pycbc.types.TimeSeries(np.real(wave), delta_t=1./1024)
#   tlen = max(len(hplus_NR), len(hplus_EOBNR))
#   hplus_EOBNR.resize(tlen)
#   hplus_NR.resize(tlen)
#
#   match = pycbc.filter.match(hplus_EOBNR,
#           hplus_NR,low_frequency_cutoff=f_min)
#   pl.title('match=%.2f'%match[0])
#
#   pl.xlim(-0.15,0.1)
#
#
#   pl.figure()
#   pl.loglog(hplus_EOBNR.to_frequencyseries().sample_frequencies,
#           abs(hplus_EOBNR.to_frequencyseries()))
#
#   pl.loglog(hplus_NR.to_frequencyseries().sample_frequencies,
#           abs(hplus_NR.to_frequencyseries()))
#
#   pl.axvline(f_min, color='r')
#
#   pl.show()
#
#   sys.exit()
#
#   # --XXX

#
# --- Reconstruction data
#
print "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed")

h1_wave_samples = np.loadtxt(os.path.join(event_file_dir, "waveforms/waveform_H1_1000.dat"))
l1_wave_samples = np.loadtxt(os.path.join(event_file_dir, "waveforms/waveform_L1_1000.dat"))
geo_wave_samples = np.loadtxt(os.path.join(event_file_dir, "geo_waveforms/waveform_geo_1000.dat"))

h1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "IFO0_asd.dat"))
l1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "IFO1_asd.dat"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

#
# --- Generate initial catalogue
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
then = timeit.time.time()
simulations = \
        bwave.simulation_details(series_names=series_names, param_bounds=bounds,
                Mmin30Hz=init_total_mass)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
        trunc_time=False)
now = timeit.time.time()
print "...catalogue construction took %.1f..."%(now-then)

#   # XXX Zero padding data
#   geo_wave_samples_new = np.zeros(shape=(len(geo_wave_samples), SI_datalen /
#       SI_deltaT))
#   for s,sample_wave in enumerate(geo_wave_samples):
#       geo_wave_samples_new[s, :len(sample_wave)] = np.copy(sample_wave)
#   geo_wave_samples = np.copy(geo_wave_samples_new)
#
# Useful time/freq samples
time_axis = np.arange(0, SI_datalen, SI_deltaT)
freq_axis = np.arange(0, catalogue.SI_flen*catalogue.SI_deltaF,
        catalogue.SI_deltaF)

# Interpolate the ASD to the waveform frequencies
h1_asd = np.interp(freq_axis, h1_bw_asd_data[:,0], h1_bw_asd_data[:,1])
l1_asd = np.interp(freq_axis, l1_bw_asd_data[:,0], l1_bw_asd_data[:,1])

mean_asd = np.sqrt(scipy.stats.hmean([h1_asd**2, l1_asd**2]))

#sys.exit()
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

    print "________________________________"
    print "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w, simulations.nsimulations)

    #
    # --- Whiten the catalogue waveform
    #
    #h1_wave = whiten_wave(wave, h1_asd)
    #l1_wave = whiten_wave(wave, l1_asd)


    # Find best-fitting mass (in terms of match)
    print "Optimising for total mass for each sampled waveform..."

    for s, (geo_sample, h1_sample, l1_sample) in enumerate(zip(geo_wave_samples,
        h1_wave_samples, l1_wave_samples)):


        print '-----------------------------'
        print "Evaluating waveform %d of %d"%( s, len(geo_wave_samples) )

        then = timeit.time.time()
        geo_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(wave,
            geo_sample, mean_asd, SI_deltaT), full_output=True, retall=True,
            disp=False)
        now = timeit.time.time()
        print "...mass optimisation took %.1f sec..."%(now-then)

        geo_matches[w, s] = 1-geo_result[1]
        geo_masses[w, s] = geo_result[0][0]

#       h1_result = scipy.optimize.fmin(mismatch, x0=geo_masses[w,s], args=(h1_wave,
#           h1_sample, None, SI_deltaT), full_output=True, retall=True,
#           disp=False)
#
#       h1_matches[w, s] = 1-h1_result[1]
#       h1_masses[w, s] = h1_result[0][0]
#
#       l1_result = scipy.optimize.fmin(mismatch, x0=geo_masses[w,s], args=(l1_wave,
#           l1_sample, None, SI_deltaT), full_output=True, retall=True,
#           disp=False)
#
#       l1_matches[w, s] = 1-l1_result[1]
#       l1_masses[w, s] = l1_result[0][0]

        print "geo: Best matching mass [match]: %.2f [%.2f]"%(
                geo_masses[w,s], geo_matches[w,s])
        print "H1: Best matching mass [match]: %.2f [%.2f]"%(
                h1_masses[w,s], h1_matches[w,s])
        print "L1: Best matching mass [match]: %.2f [%.2f]"%(
                l1_masses[w,s], l1_matches[w,s])


    geo_bestidx=np.argmax(geo_matches[w, :])
    h1_bestidx=np.argmax(h1_matches[w, :])
    l1_bestidx=np.argmax(l1_matches[w, :])

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Maximising over all waveforms:"
    print "geo: Best matching mass [match]: %.2f [%.2f]"%(
            geo_masses[w,geo_bestidx], max(geo_matches[w,:]))
    print "H1: Best matching mass [match]: %.2f [%.2f]"%(
            h1_masses[w,h1_bestidx], max(h1_matches[w,:]))
    print "L1: Best matching mass [match]: %.2f [%.2f]"%(
            l1_masses[w,l1_bestidx], max(l1_matches[w,:]))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting



