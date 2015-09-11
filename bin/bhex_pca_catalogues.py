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
bhex_pca_catalogues.py

Construct catalogues and principal component analysis for NR BBH waveforms

This script performs a PCA decomposition of two training catalogues.  A test
catalogue is then projected onto each new basis and we study the distribution of
the test waveforms in the PC spaces.

"""

import sys
from bhex_utils import bhex_wavedata as bwave
from bhex_utils import bhex_pca as bpca
from pycbc.psd import aLIGOZeroDetHighPower
import pycbc.types
import numpy as np
from matplotlib import pyplot as pl

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Some useful info:
#

#valid_series = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
#        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series"
#        "TP-series"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT (do not edit beyond this section)

# Options common to train & test data:
SI_deltaT = 1./512
SI_datalen= 4.0
total_mass = 150. 
distance=1. # Mpc

#
# Define the First train catalogue
#
train1_series_names = ['HRq-series']
train1_bounds=dict()
train1_bounds['a1'] = [0, 0]
train1_bounds['a2'] = [0, 0]
train1_bounds['q'] = [-np.inf, 3] 

#
# Define the Second train catalogue
#
train2_series_names = ['RO3-series']
train2_bounds=dict()
train2_bounds['q'] = [-np.inf, 3] 

#
# Define the test catalogue
#
# Remember: intuitive results <=> train == test
test_series_names = ['RO3-series']
test_bounds=dict()
test_bounds['q'] = [-np.inf, 3] 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build Catalogues

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'Building Training Catalogues'
print ''
print 'The principal component bases are constructed from these waveforms'
print ''
train1_simulations = \
        bwave.simulation_details(series_names=train1_series_names,
                param_bounds=train1_bounds, Mmin30Hz=total_mass)

train1_catalogue = bwave.waveform_catalogue(train1_simulations,
        ref_mass=total_mass, SI_deltaT=SI_deltaT, SI_datalen=SI_datalen,
        distance=distance)

# Second training set

train2_simulations = \
        bwave.simulation_details(series_names=train2_series_names,
                param_bounds=train2_bounds, Mmin30Hz=total_mass)

train2_catalogue = bwave.waveform_catalogue(train2_simulations,
        ref_mass=total_mass, SI_deltaT=SI_deltaT, SI_datalen=SI_datalen,
        distance=distance)

print ''
print 'Building Testing Catalogue'
print ''
print 'The principal component basis is constructed from these waveforms'
print ''
test_simulations = \
        bwave.simulation_details(series_names=test_series_names,
                param_bounds=test_bounds, Mmin30Hz=total_mass)

test_catalogue = bwave.waveform_catalogue(test_simulations, ref_mass=total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Do the PCA
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Performing PCA'
print ''
pca1 = bpca.waveform_pca(train1_catalogue, test_catalogue)
pca2 = bpca.waveform_pca(train2_catalogue, test_catalogue)

# Characterise reconstructions
# XXX build a PSD and call projection_fidelity()
psd = aLIGOZeroDetHighPower(pca1.SI_flen, pca1.SI_deltaF, pca1.fmin)
euclidean_distances1, projections1, matches1 = pca1.compute_projection_fidelity(psd=psd)
euclidean_distances2, projections2, matches2 = pca2.compute_projection_fidelity(psd=psd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Plots
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Plotting results'
print ''



#
# Matches
#
npcs = np.arange(1,pca.ntrain+1)

# Adopt the same approach as for explained variance here
match_types = ['matches_hplus','matches_ampphase']

# Dictionary to store the number of PCs required for 95% match for each
# decomposition
MatchThresh=0.97
npcMatchThresh = dict()

xmin=np.inf
xmax=-np.inf
for m, mat in enumerate(match_types):

    matches = getattr(pca,mat)
    med_matches = np.median(matches, axis=0)

    if npcs[np.argwhere(med_matches>MatchThresh)[0][0]]>xmax:
        xmax = npcs[np.argwhere(med_matches>MatchThresh)[0][0]]

    npcMatchThresh[mat] = npcs[np.argwhere(med_matches>MatchThresh)[0][0]]

    f, ax = bpca.plot_fidelity_by_npc(matches)
    ax.set_xlim(0.5,pca.ntrain+0.5)#xmax+0.5)
    ax.set_xticks(npcs)

    ax.set_ylabel('Match with %s reconstruction'%(mat.replace('matches_','')))

# --- hplus

# --- A(t) & phi(t)

# --- A(f) & phi(f)


