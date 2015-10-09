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
import lal
from pylal import spawaveform
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


#
# --- Catalogue Definition
#
match_file = sys.argv[1]
match_results = np.load(match_file)

matches = match_results['h1_matches']
total_masses = match_results['h1_masses']

ifo_label='H1'

nsamples = np.shape(matches)[1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate The Catalogue

usertag=sys.argv[2]

#
# --- Catalogue Definition
#
min_chirp_mass = 27.0
max_chirp_mass = 34.0
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
# Manipulation and derived FOMs
#

# Remove zero-match waveforms:
mean_matches = np.mean(matches, axis=1)

nonzero_match = mean_matches>0
matches = matches[nonzero_match]
total_masses = total_masses[nonzero_match]

# XXX: bit hacky..
simulations.simulations = np.array(simulations.simulations)[nonzero_match]
simulations.nsimulations = len(simulations.simulations)


# Continue
mean_matches = np.mean(matches, axis=1)

median_matches = np.median(matches, axis=1)
match_sort = np.argsort(median_matches)

median_total_masses = np.median(total_masses, axis=1)
std_total_masses = np.std(total_masses, axis=1)

mass_ratios = np.zeros(simulations.nsimulations)
chis = np.zeros(shape=(simulations.nsimulations, nsamples))
chirp_masses = np.zeros(shape=(simulations.nsimulations, nsamples))

m1 = np.zeros(shape=(simulations.nsimulations, nsamples))
m2 = np.zeros(shape=(simulations.nsimulations, nsamples))

for s, sim in enumerate(simulations.simulations):

    mass_ratios[s] = sim['q']

    spin1z = bwave.cartesian_spins(sim['a1'], sim['th1L'])
    spin2z = bwave.cartesian_spins(sim['a2'], sim['th2L'])

    for n in xrange(nsamples):

        mass1, mass2 = bwave.component_masses(total_masses[s, n], mass_ratios[s])

        chirp_masses[s, n] = spawaveform.chirpmass(mass1, mass2) \
                / lal.MTSUN_SI
        chis[s, n] = spawaveform.computechi(mass1, mass2, spin1z, spin2z)

        print mass1, mass2, spin1z, spin2z, chis[s,n]

median_mchirps = np.median(chirp_masses, axis=1)
std_mchirps = np.std(chirp_masses, axis=1)

median_chis = np.median(chis, axis=1)
std_chis = np.std(chis, axis=1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCATTER PLOTS


# --- Mass vs Mass-ratio Scatter plot
f, ax = pl.subplots()

err = ax.errorbar(mass_ratios, median_total_masses, std_total_masses, color='k',
        linestyle='None', label='1$\sigma$', ecolor='grey', zorder=-1)

scat = ax.scatter(mass_ratios, median_total_masses, c=median_matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()
#ax.set_xlim(mass_ratios[median_matches>0.5][0]-0.1,
#        mass_ratios[median_matches>0.5][-1]+0.1)
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Mass ratio (q=m$_1$/m$_2$)')
ax.set_ylabel('Total Mass [M$_{\odot}$]')

f.tight_layout()

f.savefig("BW_%s_massratio-totalmass.png"%ifo_label)

# --- Chirp Mass vs Mass-ratio Scatter plot
f, ax = pl.subplots()

err = ax.errorbar(mass_ratios, median_mchirps, std_mchirps, color='k',
        linestyle='None', label='1$\sigma$', ecolor='grey', zorder=-1)

scat = ax.scatter(mass_ratios, median_mchirps, c=median_matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper right')

ax.minorticks_on()
#ax.set_xlim(mass_ratios[median_matches>0.5][0]-0.1,
#        mass_ratios[median_matches>0.5][-1]+0.1)
ax.grid()
#ax.set_ylim(27,34)

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Mass ratio (q=m$_1$/m$_2$)')
ax.set_ylabel('$\mathcal{M}_{\mathrm{chirp}}$ [M$_{\odot}$]')

f.tight_layout()

f.savefig("BW_%s_massratio-chirpmass.png"%ifo_label)

# --- Chi vs Mass-ratio Scatter plot
f, ax = pl.subplots()

err = ax.errorbar(mass_ratios, median_chis, std_chis, color='k',
        linestyle='None', label='1$\sigma$', ecolor='grey', zorder=-1)

scat = ax.scatter(mass_ratios, median_chis, c=median_matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()
#ax.set_xlim(mass_ratios[median_matches>0.5][0]-0.1,
#        mass_ratios[median_matches>0.5][-1]+0.1)
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Mass ratio (q=m$_1$/m$_2$)')
ax.set_ylabel('Effective Spin ($\chi$)')

f.tight_layout()

f.savefig("BW_%s_massratio-chi.png"%ifo_label)

# --- Chi vs Mass Scatter plot
f, ax = pl.subplots()

err = ax.errorbar(median_total_masses, median_chis, xerr=std_total_masses,
        yerr=std_chis, color='k', linestyle='None', label='1$\sigma$', ecolor='grey', zorder=-1)

scat = ax.scatter(median_total_masses, median_chis, c=median_matches, s=50,
        label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()

idx=np.argwhere(median_matches>0.5)

#ax.set_xlim(median_total_masses[idx[0]]-1, median_total_masses[idx[-1]]+1)
#ax.set_ylim(median_chis[idx[0]], median_chis[idx[-1]])
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_xlabel('Total Mass [M$_{\odot}$]')
ax.set_ylabel('Effective Spin ($\chi$)')

f.tight_layout()

f.savefig("BW_%s_totalmass-chi.png"%ifo_label)

# --- Chi vs Chirp Mass Scatter plot
f, ax = pl.subplots()

err = ax.errorbar(median_mchirps, median_chis, xerr=std_mchirps,
        yerr=std_chis, color='k', linestyle='None', label='1$\sigma$', ecolor='grey', zorder=-1)

scat = ax.scatter(median_mchirps, median_chis, c=median_matches, s=50,
        zorder=1)
        #label='Median', zorder=1)

scat.set_clim(0.8,1)

ax.legend(loc='upper left')

ax.minorticks_on()

idx=np.argwhere(median_matches>0.5)

#ax.set_xlim(27,34)
#ax.set_ylim(median_chis[idx[0]], median_chis[idx[-1]])
ax.grid()

colbar = f.colorbar(scat) 
colbar.set_label('FF')

ax.set_ylabel('Effective Spin ($\chi$)')
ax.set_xlabel('$\mathcal{M}_{\mathrm{chirp}}$ [M$_{\odot}$]')

f.tight_layout()

f.savefig("BW_%s_chirpmass-chi.png"%ifo_label)

# --- Mass vs Chirp Scatter plot
f, ax = pl.subplots()

err = ax.errorbar(median_mchirps, median_total_masses, xerr=std_mchirps,
        yerr=std_total_masses, color='k', linestyle='None', label='1$\sigma$', ecolor='grey', zorder=-1)

scat = ax.scatter(median_mchirps, median_total_masses, c=median_matches, s=50,
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

f.savefig("BW_%s_chirpmass-totalmass.png"%ifo_label)


#sys.exit()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BOX PLOTS


# --- Match vs Waveform boxes
fmatchbox, axmatchbox = pl.subplots(figsize=(12,8))
match_box = axmatchbox.boxplot(matches[match_sort].T, whis='range', showcaps=True,
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


# --- Mass vs Waveform
#   fmassbox, axmassbox = pl.subplots(figsize=(12,8))
#   mass_box = axmassbox.boxplot(total_masses[match_sort].T, whis='range', showcaps=True,
#           showmeans=True, showfliers=False,
#           vert=False)
#   axmassbox.set_xlabel('Match-optimised mass')
#   axmassbox.set_ylabel('Waveform Parameters')
#   axmassbox.grid(linestyle='-', color='grey')
#   axmassbox.minorticks_on()
#
#   if sum(mean_matches==0):
#       axmassbox.set_ylim((len(mean_matches) -
#           np.where(mean_matches==0)[0])[0]+0.5,len(mean_matches)+0.5)
#
#   ylabels=make_labels(np.array(simulations.simulations)[match_sort])
#   axmassbox.set_yticklabels(ylabels)#, rotation=90)
#
#   fmassbox.tight_layout()

#
# Summary of best waveform
#

# 1- and 2-D Histograms of mass, match (do with a triangle plot) for the
# waveform with the highest median match

samples = np.array([matches[match_sort[-1],:], total_masses[match_sort[-1],:]]).T
trifig = triangle.corner(samples, quantiles=[0.25, 0.50, 0.75], labels=['Match', 
    'M$_{\mathrm{tot}}$ [M$_{\odot}$]'], plot_contours=True,
    plot_datapoints=True)
title = make_labels([simulations.simulations[match_sort[-1]]])
trifig.suptitle(title[0], fontsize=16)
trifig.subplots_adjust(top=0.9)

pl.show()

