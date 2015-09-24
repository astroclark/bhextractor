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
import scipy.interpolate
import timeit
from matplotlib import pyplot as pl

import lalmetaio
import lalinspiral
from glue.ligolw import lsctables, table, utils, ligolw, ilwd
table.use_in(ligolw.LIGOLWContentHandler)
lsctables.use_in(ligolw.LIGOLWContentHandler)

def swig_row(wavefile, total_mass):
    """
    Make a sim_inspiral row which points numrel_data to the wavefile and sets up
    m1 and m2 to correspond to the total mass of the system
    """

    swigrow = lalmetaio.SimInspiralTable()

    swigrow.mass1=0.5*total_mass
    swigrow.mass2=0.5*total_mass
    swigrow.numrel_data="file://localhost"+os.path.realpath(wavefile)
    swigrow.numrel_mode_max=2
    swigrow.numrel_mode_min=2
    swigrow.f_lower=10.0
    swigrow.distance=1.0

    return swigrow

def scale_wave(wave_file, total_mass):
    """
    Scale the waveform to total_mass.  Assumes the waveform is initially
    generated at init_total_mass defined in this script.

    """

    row = swig_row(wave_file, total_mass)
    Hp, Hc = lalinspiral.NRInjectionFromSimInspiral(row, SI_deltaT)

    hp = pycbc.types.TimeSeries(Hp.data.data[:], delta_t=Hp.deltaT, epoch=Hp.epoch)
    hc = pycbc.types.TimeSeries(Hc.data.data[:], delta_t=Hp.deltaT, epoch=Hp.epoch)

    return hp, hc


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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Some useful info:
#

#valid_series = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
#        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series"
#        "TP-series"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

#SI_deltaT = 1./1024
SI_deltaT = 1./4096
SI_datalen= 4.0
SI_deltaF=1./SI_datalen
SI_flen=len(np.fft.fftfreq(int(SI_datalen/SI_deltaT), SI_deltaT))
f_min = 40.0

nsampls=500

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
#bounds=dict()
#bounds['q'] = [1, 1]
#bounds['a1'] = [0,0]
#bounds['a2'] = [0,1]

# XXX SANITY CHECKING SCALING AND BEHAVIOUR OF NR WAVE wrt SEOBNR

if 0:
    mtot=70
    series_names = ['HRq-series']#, 'HRq-series', 'RO3-series'] # (see above for valid choices)
    bounds=dict()
    bounds['q'] = [1, 1]
    from pycbc.waveform import get_td_waveform
    hplus_EOBNR, _ = get_td_waveform(approximant="SEOBNRv2",
            distance=distance,
            mass1=0.5*mtot,
            mass2=0.5*mtot,
            spin1z=0.0,
            spin2z=0.0,
            f_lower=10,
            delta_t=SI_deltaT)
    emax=np.argmax(abs(hplus_EOBNR))
   
    #import lal
    #sY22 = lal.SpinWeightedSphericalHarmonic(0,0,-2,2,2)
    #hplus_EOBNR.data /= np.real(sY22)
    hplus_EOBNR.data = bwave.window_wave(hplus_EOBNR.data)
   
   
    simulations =  bwave.simulation_details(series_names=series_names, param_bounds=bounds,
            Mmin30Hz=100.0)

    # rescale 
    then = timeit.time.time()
    waveframe = simulations.simulations[0]['wavefile'].replace('asc','gwf')
    hplus_NR, _ = scale_wave(waveframe, mtot)
    now = timeit.time.time()
    print "Waveform scaling took %.3f sec"%(now-then)
    #hplus_NR.data = bwave.window_wave(hplus_NR.data)
   
    tlen = max(len(hplus_NR), len(hplus_EOBNR))
    hplus_EOBNR.resize(tlen)
    hplus_NR.resize(tlen)
   
    match = pycbc.filter.match(hplus_EOBNR,
            hplus_NR,low_frequency_cutoff=f_min)
   
    pl.figure()
    pl.loglog(hplus_EOBNR.to_frequencyseries().sample_frequencies,
            abs(hplus_EOBNR.to_frequencyseries()), label='SEOBNRv2')
   
    pl.loglog(hplus_NR.to_frequencyseries().sample_frequencies,
            abs(hplus_NR.to_frequencyseries()), label='NR')
    pl.title('M$_{\mathrm{tot}}$=%d M$_{\odot}$: match=%.2f'%(mtot,match[0]))
    pl.xlabel('Frequency [Hz]')
    pl.ylabel('|H(f)|')
    pl.legend()
   
    pl.axvline(f_min, color='r')
   
    pl.show()
   
    sys.exit()
#
#   # --XXX

#
# --- Reconstruction data
#
print "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed")

geo_wave_samples_in = np.loadtxt(os.path.join(event_file_dir, "geo_waveforms/waveform_geo_1000.dat"))

# Downsample the number of posterior samples (useful for testing)
geo_wave_samples_in=geo_wave_samples_in[np.random.random_integers(low=0,high=nsampls,size=nsampls),:]


# XXX Zero padding data
#   geo_wave_samples = np.zeros(shape=(len(geo_wave_samples_in),
#           SI_datalen/SI_deltaT))
#   for s, sample in enumerate(geo_wave_samples_in):
#       geo_wave_samples[s,:len(sample)] = sample

# XXX Resampling data
geo_wave_samples = np.zeros(shape=(len(geo_wave_samples_in),
        int(SI_datalen/SI_deltaT)))
for s, sample in enumerate(geo_wave_samples_in):
    geo_wave_samples[s,:] = scipy.signal.resample(sample,
            int(SI_datalen/SI_deltaT))

h1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "IFO0_asd.dat"))
l1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "IFO1_asd.dat"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
then = timeit.time.time()
simulations = \
        bwave.simulation_details(series_names=series_names, param_bounds=bounds,
                Mmin30Hz=init_total_mass)

now = timeit.time.time()
print "...catalogue construction took %.1f..."%(now-then)


# Useful time/freq samples
time_axis = np.arange(0, SI_datalen, SI_deltaT)
freq_axis = np.arange(0, SI_flen * SI_deltaF, SI_deltaF)

# Interpolate the ASD to the waveform frequencies
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

# Loop over waves in NR catalogue
for w in xrange(simulations.nsimulations):

    print "________________________________"
    print "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w+1, simulations.nsimulations)

    waveframe = simulations.simulations[w]['wavefile'].replace('asc','gwf')

    print "Optimising for total mass for each sampled waveform..."

    for s, geo_sample in enumerate(geo_wave_samples):


        print '-----------------------------'
        print "Evaluating sample waveform %d of %d"%( s, len(geo_wave_samples) )
        print " (NR: %s [%d/%d])"%( simulations.simulations[w]['wavename'], w+1,
                simulations.nsimulations)

        # Find best-fitting mass (in terms of match)
        then = timeit.time.time()
        geo_result = scipy.optimize.fmin(mismatch, x0=mass_guess,
                args=(waveframe, geo_sample, mean_asd, SI_deltaT,
                    SI_deltaF),
            full_output=True, retall=True, disp=False)
        now = timeit.time.time()
        print "...mass optimisation took %.3f sec..."%(now-then)

        geo_matches[w, s] = 1-geo_result[1]
        geo_masses[w, s] = geo_result[0][0]

        print "geo: Best matching mass [match]: %.2f [%.2f]"%(
                geo_masses[w,s], geo_matches[w,s])


    geo_bestidx=np.argmax(geo_matches[w, :])

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Maximising over all waveforms:"
    print "geo: Best matching mass [match]: %.2f [%.2f]"%(
            geo_masses[w,geo_bestidx], max(geo_matches[w,:]))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dump data
filename=series_names[0] + 'v2'
np.savez(filename, geo_matches=geo_matches, geo_masses=geo_masses)



