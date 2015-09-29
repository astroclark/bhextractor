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
import pycbc.types
import numpy as np
import timeit
from matplotlib import pyplot as pl
import triangle

pl.rcParams.update({'axes.labelsize': 16})
pl.rcParams.update({'xtick.labelsize':16})
pl.rcParams.update({'ytick.labelsize':16})
pl.rcParams.update({'legend.fontsize':16})

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

def make_labels(simulations):
    """
    Return a list of strings with suitable labels for e.g., box plots
    """

    labels=[]
    for sim in simulations:
        labelstr = \
                r"$q=%.1f$, $a_1=%.1f$, $a_2=%.1f$, $\theta_1=%.1f$, $\theta_2=%.1f$, $\phi_{12}=%.1f$"%(
                        sim['q'], sim['a1'], sim['a2'], sim['th1L'], sim['th2L'],
                        sim['th12'])
        labels.append(labelstr)

    return labels

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
f_min = 40.0


init_total_mass = 100.  # Select waveforms which can go down to at least this mass
                        # (and generate SI catalogue at this mass)
distance=1. # Mpc

#
# --- Catalogue Definition
#
series_names = [sys.argv[1].split('_')[0]]
match_file = sys.argv[1]
match_results = np.load(match_file)
geo_matches = match_results['geo_matches']
geo_masses = match_results['geo_masses']

bounds=None
#bounds=dict()
#bounds['q'] = [1, 1]
#bounds['a1'] = [0,0]
#bounds['a2'] = [0,1]

#
#   # --XXX

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
#catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
#        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
#        trunc_time=False)
now = timeit.time.time()
print "...catalogue construction took %.1f..."%(now-then)


# Useful time/freq samples
time_axis = np.arange(0, SI_datalen, SI_deltaT)
#freq_axis = np.arange(0, catalogue.SI_flen*catalogue.SI_deltaF,
#        catalogue.SI_deltaF)
#sys.exit()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Manipulation and derived FOMs
#

mean_matches = np.mean(geo_matches, axis=1)
median_matches = np.median(geo_matches, axis=1)
match_sort = np.argsort(median_matches)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting

# Plots to make:
# box plots for distributions of mass / match
#   sort these plots in descending match
#   label by parsing parameters from the simulations object

fmatchbox, axmatchbox = pl.subplots(figsize=(12,8))
match_box = axmatchbox.boxplot(geo_matches[match_sort].T, whis='range', showcaps=True,
        showmeans=True, showfliers=False,
        vert=False)
axmatchbox.set_xlabel('Mass-optimised Match')
axmatchbox.set_ylabel('Waveform Parameters')
axmatchbox.grid(linestyle='-', color='grey')
axmatchbox.minorticks_on()

ylabels=make_labels(np.array(simulations.simulations)[match_sort])
axmatchbox.set_yticklabels(ylabels)#, rotation=90)

fmatchbox.tight_layout()

# --- masses
fmassbox, axmassbox = pl.subplots(figsize=(12,8))
mass_box = axmassbox.boxplot(geo_masses[match_sort].T, whis='range', showcaps=True,
        showmeans=True, showfliers=False,
        vert=False)
axmassbox.set_xlabel('Match-optimised mass')
axmassbox.set_ylabel('Waveform Parameters')
axmassbox.grid(linestyle='-', color='grey')
axmassbox.minorticks_on()

ylabels=make_labels(np.array(simulations.simulations)[match_sort])
axmassbox.set_yticklabels(ylabels)#, rotation=90)

fmassbox.tight_layout()

# 1- and 2-D Histograms of mass, match (do with a triangle plot) for the
# waveform with the highest median match

samples = np.array([geo_matches[match_sort[-1],:], geo_masses[match_sort[-1],:]]).T
trifig = triangle.corner(samples, quantiles=[0.25, 0.50, 0.75], labels=['Match', 
    'M$_{\mathrm{tot}}$ [M$_{\odot}$]'], plot_contours=True,
    plot_datapoints=True)
title = make_labels([simulations.simulations[match_sort[-1]]])
trifig.suptitle(title[0], fontsize=16)
trifig.subplots_adjust(top=0.9)


pl.show()

sys.exit()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load reconstruction data

#
# --- Reconstruction data
#
print "Loading data"
event_file_dir = os.path.join(os.environ.get('BHEX_PREFIX'),
        "data/observed")

geo_wave_samples_in = np.loadtxt(os.path.join(event_file_dir, "geo_waveforms/waveform_geo_1000.dat"))

# XXX Resampling data
geo_wave_samples = np.zeros(shape=(len(geo_wave_samples_in),
        int(SI_datalen/SI_deltaT)))
for s, sample in enumerate(geo_wave_samples_in):
    geo_wave_samples[s,:] = scipy.signal.resample(sample,
            int(SI_datalen/SI_deltaT))

h1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "IFO0_asd.dat"))
l1_bw_asd_data = np.loadtxt(os.path.join(event_file_dir, "IFO1_asd.dat"))


# Interpolate the ASD to the waveform frequencies
h1_asd = np.interp(freq_axis, h1_bw_asd_data[:,0], h1_bw_asd_data[:,1])
l1_asd = np.interp(freq_axis, l1_bw_asd_data[:,0], l1_bw_asd_data[:,1])

mean_asd = np.sqrt(scipy.stats.hmean([h1_asd**2, l1_asd**2]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting



