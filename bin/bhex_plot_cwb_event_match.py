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
match_file = sys.argv[1]
match_results = np.load(match_file)
h1_matches = match_results['h1_matches']
h1_masses = match_results['h1_masses']
l1_matches = match_results['l1_matches']
l1_masses = match_results['l1_masses']

#bounds=None
bounds=dict()

#bounds['th1L'] = [0,0]
#bounds['th2L'] = [0,0]

bounds['a1'] = [0,0]
bounds['a2'] = [0,0]

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Manipulation
h1_match_sort = np.argsort(h1_matches)
l1_match_sort = np.argsort(l1_matches)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting


fmatchplot, axmatchplot = pl.subplots(figsize=(12,8))
h1_match_plot = axmatchplot.plot(h1_matches[h1_match_sort],
        range(1,len(h1_matches)+1), 'rs', label='H1 response')
l1_match_plot = axmatchplot.plot(l1_matches[h1_match_sort],
        range(1,len(h1_matches)+1), 'go', label='L1 response')
axmatchplot.set_xlabel('Mass-optimised Match')
axmatchplot.set_ylabel('Waveform Parameters')
axmatchplot.grid(linestyle='-', color='grey')
axmatchplot.minorticks_on()

axmatchplot.set_yticks(range(1,len(h1_matches)+1))

if sum(h1_matches==0):
    axmatchplot.set_ylim((len(h1_matches) -
        np.where(h1_matches==0)[0])[0]+0.5,len(h1_matches)+0.5)

    axmatchplot.set_xlim(0.8,1)

ylabels=make_labels(np.array(simulations.simulations)[h1_match_sort])
axmatchplot.set_yticklabels(ylabels)#, rotation=90)

fmatchplot.tight_layout()

pl.show()
sys.exit()


