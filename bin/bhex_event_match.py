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

#    resamp_hplus = np.interp(time_axis, interp_times, np.real(wave))
#    resamp_hcross = np.interp(time_axis, interp_times, np.imag(wave))
#    return resamp_hplus #+ 1j*resamp_hcross

    resamp_amp = np.interp(time_axis, interp_times, amp)
    resamp_phase = np.interp(time_axis, interp_times, phase)
 
#    f = scipy.interpolate.interp1d(interp_times, amp, kind='cubic')
#    resamp_amp = f(time_axis)
#    resamp_phase = f(time_axis)
    

    return resamp_amp*np.exp(1j*resamp_phase)


def mismatch(total_mass, mass_bounds, tmplt_wave_data, event_wave_data, asd=None,
        delta_t=1./512, delta_f=0.25):
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

    else:

        return 1.


def mtot_from_mchirp(mc, q):
    eta = q/(1+q)**2.0
    return mc * eta**(-3./5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

#SI_deltaT = 1./1024
SI_deltaT = 1./4096
SI_datalen= 4.0
f_min = 40.0

#nsampls=500

# Initial guess at the mass
mass_guess = 74.0#100 + 100*np.random.random()

init_total_mass = 100.  # Select waveforms which can go down to at least this mass
                        # (and generate SI catalogue at this mass)
distance=1. # Mpc

usertag=sys.argv[1]

#
# --- Catalogue Definition
#
bounds=None
#bounds=dict()
#bounds['q'] = [4, 6]
#bounds['a1'] = [0,0]
#bounds['a2'] = [0,1]

# XXX SANITY CHECKING SCALING AND BEHAVIOUR OF NR WAVE wrt SEOBNR

if 0:
    mtot=70
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
   
    import lal
    sY22 = lal.SpinWeightedSphericalHarmonic(0,0,-2,2,2)
    hplus_EOBNR.data /= np.real(sY22)
    hplus_EOBNR.data = bwave.window_wave(hplus_EOBNR.data)
   
   
    simulations = \
            bwave.simulation_details(param_bounds=bounds, Mmin30Hz=100.0)
   
    catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
            SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
            trunc_time=False)
   
    # Useful time/freq samples
    wave = np.real(catalogue.SIComplexTimeSeries[0,:])
   
    nmax = np.argmax(abs(wave))
   
    time_axis = np.arange(0, SI_datalen, SI_deltaT)  
   
    # rescale 
    then = timeit.time.time()
    wave = scale_wave(wave, mtot)
    now = timeit.time.time()

    print "Waveform scaling took %.3f sec"%(now-then)
   
   
    hplus_NR = pycbc.types.TimeSeries(np.real(wave), delta_t=SI_deltaT)
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
#geo_wave_samples_in=geo_wave_samples_in[np.random.random_integers(low=0,high=nsampls,size=nsampls),:]


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

#
# --- Generate initial catalogue
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
then = timeit.time.time()
simulations = \
        bwave.simulation_details(param_bounds=bounds, Mmin30Hz=init_total_mass)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
        trunc_time=False)
now = timeit.time.time()
print "...catalogue construction took %.1f..."%(now-then)


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

# Loop over waves in NR catalogue
for w, wave in enumerate(catalogue.SIComplexTimeSeries):

    print "________________________________"
    print "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w+1, simulations.nsimulations)


    # Find best-fitting mass (in terms of match)
    print "Optimising for total mass for each sampled waveform..."

    # Find minimum mass
    min_mass = mtot_from_mchirp(27, simulations.simulations[w]['q'])
    max_mass = mtot_from_mchirp(34, simulations.simulations[w]['q'])

    for s, geo_sample in enumerate(geo_wave_samples):


        print '-----------------------------'
        print "Evaluating sample waveform %d of %d"%( s, len(geo_wave_samples) )
        print " (NR: %s [%d/%d])"%( simulations.simulations[w]['wavename'], w+1,
                simulations.nsimulations)

        then = timeit.time.time()
        geo_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
            [min_mass, max_mass], wave, geo_sample, mean_asd, SI_deltaT,
            catalogue.SI_deltaF), full_output=True, retall=True, disp=False)
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
import uuid
filename=usertag+'_'+str(uuid.uuid4())
np.savez(filename, geo_matches=geo_matches, geo_masses=geo_masses)



