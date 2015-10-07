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
bhex_cwb_match.py

Compute matches between CWB reconstruction and NR waveforms
"""

import sys, os
from bhex_utils import bhex_wavedata as bwave
import pycbc.filter
import pycbc.types
import lal
from pylal import spawaveform
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
            peakidx*deltaT*(scale_ratio-1)

    resamp_amp = np.interp(time_axis, interp_times, amp)
    resamp_phase = np.interp(time_axis, interp_times, phase)
    

    return resamp_amp*np.exp(1j*resamp_phase)


def mismatch(total_mass, mass_bounds, tmplt_wave_data, event_wave_data, asd=None,
        delta_t=1./512, delta_f=0.25, whitened_response=False):
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
        hp = pycbc.types.TimeSeries(np.real(tmplt[:]), delta_t=deltaT)

        if whitened_response and asd is not None:
            # Whiten the template
            Hp = hp.to_frequencyseries()
            Hp.data /= asd

            # IFFT (just simplifies the code below) 
            hp = Hp.to_timeseries()


        if asd is not None and not whitened_response:
            psd = pycbc.types.FrequencySeries(asd**2, delta_f=delta_f)
        else:
            psd = None

        # Put the reconstruction data in a TimeSeries
        event_wave = pycbc.types.TimeSeries(event_wave_data, delta_t=deltaT)

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

def extract_wave(inwave, datalen=4.0, sample_rate = 4096):

    extract_len = 0.5
    peakidx = np.argmax(inwave) - 500
    nsamp = extract_len * sample_rate

    extracted = inwave[int(peakidx-0.5*nsamp): int(peakidx+0.5*nsamp)]

    win = lal.CreateTukeyREAL8Window(len(extracted), 0.5)
    extracted *= win.data.data

    output = np.zeros(int(datalen*sample_rate))

    output[int(0.5*datalen*sample_rate-0.5*nsamp):
            int(0.5*datalen*sample_rate+0.5*nsamp)] = np.copy(extracted)

    return output

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

deltaT = 1./4096
datalen= 4.0
f_min = 30.0

#Initial guess at the mass
mass_guess = 72.0
min_chirp_mass = 27.0
max_chirp_mass = 34.0
init_total_mass = 100   # Generate a catagg

distance=1. # Mpc

usertag=sys.argv[1]

# --- Catalogue Definition
#
bounds = bwave.bounds_dict(usertag)
bounds = dict()
bounds['Mchirpmin30Hz'] = [-np.inf, min_chirp_mass]


#
# --- Reconstruction data
#
print "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed")

h1_wave_in = np.loadtxt(os.path.join(event_file_dir, "cwb/H1_wf_signal.dat"))
l1_wave_in = np.loadtxt(os.path.join(event_file_dir, "cwb/L1_wf_signal.dat"))


# Extract the 4 second chunk in the middle
h1_wave = extract_wave(h1_wave_in, datalen=datalen, sample_rate=1./deltaT)
l1_wave = extract_wave(l1_wave_in, datalen=datalen, sample_rate=1./deltaT)

h1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "cwb/h1_asd.dat"))
l1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "cwb/l1_asd.dat"))

#h1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/IFO0_asd.dat"))
#l1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "bw/IFO1_asd.dat"))

#sys.exit()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

#
# --- Generate initial catalogue
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
then = timeit.time.time()
simulations = bwave.simulation_details(param_bounds=bounds)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
        SI_deltaT=deltaT, SI_datalen=datalen, distance=distance,
        trunc_time=False)
now = timeit.time.time()
print "...catalogue construction took %.1f..."%(now-then)


# Useful time/freq samples
time_axis = np.arange(0, datalen, deltaT)
freq_axis = np.arange(0, catalogue.SI_flen*catalogue.SI_deltaF,
        catalogue.SI_deltaF)

# Interpolate the ASD to the waveform frequencies
h1_asd = np.interp(freq_axis, h1_cwb_asd_data[:,0], h1_cwb_asd_data[:,1])
l1_asd = np.interp(freq_axis, l1_cwb_asd_data[:,0], l1_cwb_asd_data[:,1])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter Estimation
#

#
# --- Minimise mismatch over waveforms and total mass
#

# Preallocate
h1_matches = np.zeros(shape=(simulations.nsimulations))
h1_masses  = np.zeros(shape=(simulations.nsimulations))

l1_matches = np.zeros(shape=(simulations.nsimulations))
l1_masses  = np.zeros(shape=(simulations.nsimulations))

# Loop over waves in NR catalogue
for w, wave in enumerate(catalogue.SIComplexTimeSeries):

    print "________________________________"
    print "Computing match for %s (%d/%d)"%(simulations.simulations[w]['wavename'],
            w+1, simulations.nsimulations)


    # Find best-fitting mass (in terms of match)
    print "Optimising for total mass for CWB reconstructions..."

    # Find minimum mass
    min_mass = mtot_from_mchirp(min_chirp_mass, simulations.simulations[w]['q'])
    max_mass = mtot_from_mchirp(max_chirp_mass, simulations.simulations[w]['q'])


    # --- H1
    then = timeit.time.time()

    ifo_response=True
    then = timeit.time.time()
    h1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
        [min_mass, max_mass], wave, h1_wave, h1_asd, deltaT,
        catalogue.SI_deltaF, ifo_response), full_output=True,
        retall=True, disp=False)

    now = timeit.time.time()
    print "...H1 mass optimisation took %.3f sec..."%(now-then)

    h1_matches[w] = 1-h1_result[1]
    h1_masses[w] = h1_result[0][0]

    print "H1: Best matching mass [match]: %.2f [%.2f]"%(
            h1_masses[w], h1_matches[w])

    # --- L1
    then = timeit.time.time()

    ifo_response=True
    then = timeit.time.time()
    l1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
        [min_mass, max_mass], wave, l1_wave, l1_asd, deltaT,
        catalogue.SI_deltaF, ifo_response), full_output=True,
        retall=True, disp=False)

    now = timeit.time.time()
    print "...L1 mass optimisation took %.3f sec..."%(now-then)

    l1_matches[w] = 1-l1_result[1]
    l1_masses[w] = l1_result[0][0]

    print "L1: Best matching mass [match]: %.2f [%.2f]"%(
            l1_masses[w], l1_matches[w])



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dump data
import uuid
filename=usertag+'_CWB_'+str(uuid.uuid4())
np.savez(filename, h1_matches=h1_matches, h1_masses=h1_masses,
        l1_matches=l1_matches, l1_masses=l1_masses)



