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

min_chirp_mass = 27.0
max_chirp_mass = 34.0

#
# --- Catalogue Definition
#
match_file = sys.argv[1]
match_results = np.load(match_file)
matches = match_results['matches']
total_masses = match_results['total_masses']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

usertag=sys.argv[2]

#
# --- Catalogue Definition
#
bounds = bwave.bounds_dict(usertag)
bounds = dict()
bounds['Mchirpmin30Hz'] = [-np.inf, min_chirp_mass]


#
# --- Generate initial catalogue
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
then = timeit.time.time()
simulations = bwave.simulation_details(param_bounds=bounds)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Manipulation
h1_match_sort = np.argsort(matches)
l1_match_sort = np.argsort(l1_matches)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCATTER PLOTS


# --- Mass vs Mass-ratio Scatter plot
f, ax = pl.subplots()

scat = ax.scatter(mass_ratios, total_masses, c=matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()
#ax.set_xlim(mass_ratios[matches>0.5][0]-0.1,
#        mass_ratios[matches>0.5][-1]+0.1)
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Mass ratio (q=m$_1$/m$_2$)')
ax.set_ylabel('Total Mass [M$_{\odot}$]')

f.tight_layout()

f.savefig("CWB_%s_massratio-totalmass.png"%ifo_label)

# --- Chirp Mass vs Mass-ratio Scatter plot
f, ax = pl.subplots()

scat = ax.scatter(mass_ratios, mchirps, c=matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper right')

ax.minorticks_on()
#ax.set_xlim(mass_ratios[matches>0.5][0]-0.1,
#        mass_ratios[matches>0.5][-1]+0.1)
ax.grid()
#ax.set_ylim(27,34)

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Mass ratio (q=m$_1$/m$_2$)')
ax.set_ylabel('$\mathcal{M}_{\mathrm{chirp}}$ [M$_{\odot}$]')

f.tight_layout()

f.savefig("CWB_%s_massratio-chirpmass.png"%ifo_label)

# --- Chi vs Mass-ratio Scatter plot
f, ax = pl.subplots()

scat = ax.scatter(mass_ratios, chis, c=matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()
#ax.set_xlim(mass_ratios[matches>0.5][0]-0.1,
#        mass_ratios[matches>0.5][-1]+0.1)
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Mass ratio (q=m$_1$/m$_2$)')
ax.set_ylabel('Effective Spin ($\chi$)')

f.tight_layout()

f.savefig("CWB_%s_massratio-chi.png"%ifo_label)

# --- Chi vs Mass Scatter plot
f, ax = pl.subplots()

scat = ax.scatter(total_masses, chis, c=matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()

idx=np.argwhere(matches>0.5)

#ax.set_xlim(total_masses[idx[0]]-1, total_masses[idx[-1]]+1)
#ax.set_ylim(chis[idx[0]], chis[idx[-1]])
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Total Mass [M$_{\odot}$]')
ax.set_ylabel('Effective Spin ($\chi$)')

f.tight_layout()

f.savefig("CWB_%s_totalmass-chi.png"%ifo_label)

# --- Chi vs Chirp Mass Scatter plot
f, ax = pl.subplots()

scat = ax.scatter(mchirps, chis, c=matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()

idx=np.argwhere(matches>0.5)

#ax.set_xlim(27,34)
#ax.set_ylim(chis[idx[0]], chis[idx[-1]])
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_ylabel('Effective Spin ($\chi$)')
ax.set_xlabel('$\mathcal{M}_{\mathrm{chirp}}$ [M$_{\odot}$]')

f.tight_layout()

f.savefig("CWB_%s_chirpmass-chi.png"%ifo_label)

# --- Mass vs Chirp Scatter plot
f, ax = pl.subplots()

scat = ax.scatter(mchirps, total_masses, c=matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper right')

ax.minorticks_on()
#ax.set_xlim(27,34)
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('$\mathcal{M}_{\mathrm{chirp}}$ [M$_{\odot}$]')
ax.set_ylabel('Total Mass [M$_{\odot}$]')

f.tight_layout()

f.savefig("CWB_%s_chirpmass-totalmass.png"%ifo_label)

pl.show()
sys.exit()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BoxPlots

fmatchplot, axmatchplot = pl.subplots(figsize=(12,8))
h1_match_plot = axmatchplot.plot(matches[h1_match_sort],
        range(1,len(matches)+1), 'rs', label='H1 response')
axmatchplot.set_xlabel('Mass-optimised Match')
axmatchplot.set_ylabel('Waveform Parameters')
axmatchplot.grid(linestyle='-', color='grey')
axmatchplot.minorticks_on()

axmatchplot.set_yticks(range(1,len(matches)+1))

if sum(matches==0):
    axmatchplot.set_ylim((len(matches) -
        np.where(matches==0)[0])[0]+0.5,len(matches)+0.5)

#    axmatchplot.set_xlim(0.8,1)

ylabels=make_labels(np.array(simulations.simulations)[h1_match_sort])
axmatchplot.set_yticklabels(ylabels)#, rotation=90)

fmatchplot.tight_layout()

pl.show()

