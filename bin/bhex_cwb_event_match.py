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

    resamp_amp = np.interp(time_axis, interp_times, amp)
    resamp_phase = np.interp(time_axis, interp_times, phase)
    

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

def extract_wave(inwave, datalen=4.0, sample_rate = 2048):
    peakidx = np.argmax(h1_wave_in)
    nsamp = datalen * sample_rate
    return inwave[int(peakidx-0.5*nsamp): int(peakidx+0.5*nsamp)]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

SI_deltaT = 1./2048
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
if usertag=="NoConstraint":
    bounds=None

elif usertag=="NonSpinning":
    # NonSpinning
    bounds=dict()
    bounds['a1'] = [0,0]
    bounds['a2'] = [0,0]

elif usertag=="AlignedSpinUp":
    # AlignedSpinUp
    bounds=dict()
    bounds['a1'] = [0.001,np.inf]
    bounds['a2'] = [0.001,np.inf]
    bounds['th1L'] = [0,0]
    bounds['th2L'] = [0,0]

elif  usertag=="AlignedSpinDown":
    # AlignedSpinDown
    bounds=dict()
    bounds['a1'] = [0.001,np.inf]
    bounds['a2'] = [0.001,np.inf]
    bounds['th1L'] = [180,180]
    bounds['th2L'] = [180,180]

elif usertag=="BigBHSpinUp":
    # BigBHSpinUp
    bounds=dict()
    bounds['a1'] = [0.001,np.inf]
    bounds['a2'] = [0, 0]
    bounds['th1L'] = [0,0]

elif usertag=="BigBHSpinDown":
    # BigBHSpinDown
    bounds=dict()
    bounds['a1'] = [0.001,np.inf]
    bounds['a2'] = [0,0]
    bounds['th1L'] = [180,180]

elif usertag=="SmallBHSpinUp":
    # SmallBHSpinUp
    bounds=dict()
    bounds['a1'] = [0,0]
    bounds['a2'] = [0.001,np.inf]
    bounds['th2L'] = [0,0]

elif usertag=="SmallBHSpinDown":
    # SmallBHSpinDown
    bounds=dict()
    bounds['a1'] = [0, 0]
    bounds['a2'] = [0.001,np.inf]
    bounds['th1L'] = [180,180]

else:
    print >> sys.stderr, "Configuration not recognised"
    sys.exit(-1)



#
# --- Reconstruction data
#
print "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed/cwb")

h1_wave_in = np.loadtxt(os.path.join(event_file_dir, "H1_wf_strain.dat"))
l1_wave_in = np.loadtxt(os.path.join(event_file_dir, "L1_wf_strain.dat"))


# Extract the 4 second chunk in the middle
h1_wave = extract_wave(h1_wave_in, datalen=SI_datalen)
l1_wave = extract_wave(l1_wave_in, datalen=SI_datalen)

# resample to 4096 
#h1_wave = scipy.signal.resample(h1_wave, 2*len(h1_wave))
#l1_wave = scipy.signal.resample(l1_wave, 2*len(l1_wave))

h1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "h1_asd.dat"))
l1_cwb_asd_data = np.loadtxt(os.path.join(event_file_dir, "l1_asd.dat"))

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
    min_mass = mtot_from_mchirp(27, simulations.simulations[w]['q'])
    max_mass = mtot_from_mchirp(34, simulations.simulations[w]['q'])


    # --- H1
    then = timeit.time.time()

    h1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
        [min_mass, max_mass], wave, h1_wave, h1_asd, SI_deltaT,
        catalogue.SI_deltaF), full_output=True, retall=True, disp=False)

    now = timeit.time.time()
    print "...H1 mass optimisation took %.3f sec..."%(now-then)

    h1_matches[w] = 1-h1_result[1]
    h1_masses[w] = h1_result[0][0]

    print "H1: Best matching mass [match]: %.2f [%.2f]"%(
            h1_masses[w], h1_matches[w])

    # --- L1
    then = timeit.time.time()

    l1_result = scipy.optimize.fmin(mismatch, x0=mass_guess, args=(
        [min_mass, max_mass], wave, l1_wave, l1_asd, SI_deltaT,
        catalogue.SI_deltaF), full_output=True, retall=True, disp=False)

    now = timeit.time.time()
    print "...L1 mass optimisation took %.3f sec..."%(now-then)

    l1_matches[w] = 1-l1_result[1]
    l1_masses[w] = l1_result[0][0]

    print "L1: Best matching mass [match]: %.2f [%.2f]"%(
            h1_masses[w], h1_matches[w])



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dump data
import uuid
filename=usertag+'_CWB_'+str(uuid.uuid4())
np.savez(filename, h1_matches=h1_matches, h1_masses=h1_masses,
        l1_matches=l1_matches, l1_masses=l1_masses)



