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

def make_labels(simulations):
    """
    Return a list of strings with suitable labels for e.g., box plots
    """

    labels=[]
    for sim in simulations:

        # check nans
        vals = []
        for val in [sim['q'], sim['a1'], sim['a2'], sim['th1L'], sim['th2L']]:
            if np.isnan(val):
                val = 0.0
            vals.append(val)

        labelstr = \
                r"$q=%.2f$, $a_1=%.2f$, $a_2=%.2f$, $\theta_1=%.2f$, $\theta_2=%.2f$"%(
                        vals[0], vals[1], vals[2], vals[3], vals[4])
        labels.append(labelstr)

    return labels



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT


init_total_mass = 65.  # Select waveforms which can go down to at least this mass
                        # (and generate SI catalogue at this mass)
distance=1. # Mpc

#
# --- Catalogue Definition
#
match_file = sys.argv[1]
match_results = np.load(match_file)
geo_matches = match_results['geo_matches']
geo_masses = match_results['geo_masses']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

usertag=sys.argv[2]

#
# --- Catalogue Definition
#
bounds = bwave.bounds_dict(usertag)


#
# --- Generate initial catalogue
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
then = timeit.time.time()
simulations = \
        bwave.simulation_details(param_bounds=bounds, Mmin30Hz=init_total_mass)



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

if sum(mean_matches==0):
    axmatchbox.set_ylim((len(mean_matches) -
        np.where(mean_matches==0)[0])[0]+0.5,len(mean_matches)+0.5)

    axmatchbox.set_xlim(0.8,1)

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

if sum(mean_matches==0):
    axmassbox.set_ylim((len(mean_matches) -
        np.where(mean_matches==0)[0])[0]+0.5,len(mean_matches)+0.5)

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

fmatchmax, axmatchmass = pl.subplots(nrows=2, sharex=True)
axmatchmass[0].plot(median_masses, median_matches, 'ms', label='Median value')
axmatchmass[0].set_ylim(0.8, 0.95)
axmatchmass[0].legend()
axmatchmass[0].minorticks_on()
axmatchmass[1].set_ylim(0.8, 0.95)
axmatchmass[1].legend()
axmatchmass[1].set_xlabel('Total Mass [M$_{\odot}$]')
axmatchmass[1].minorticks_on()
axmatchmass[0].set_ylabel('Match')
axmatchmass[1].set_ylabel('Match')
pl.subplots_adjust(hspace=0)

pl.show()


