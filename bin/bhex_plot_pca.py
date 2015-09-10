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
bhextractor_plot_pca.py

Construct catalogues and principal component analysis for NR BBH waveforms

This script simply builds a catalogue (defined by the user in this
script) and produces plots of the PC basis and related characteristics.

See bhex_pca_projections.py for projections of test catalogues onto a basis and
reconstruction fidelity measures.

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
# Define the train catalogue
#
train_series_names = ['HRq-series']
train_bounds=dict()
train_bounds['a1'] = [0, 0]
train_bounds['a2'] = [0, 0]
train_bounds['q'] = [-np.inf, 3] 

#
# Define the test catalogue
#
# Remember: intuitive results <=> train == test
#   test_series_names = ['HR-series']
#   test_bounds=dict()
#   test_bounds['a1'] = [0, 0]
#   test_bounds['a2'] = [0, 0]
#   test_bounds['q'] = [-np.inf, 3] 
test_series_names = ['HRq-series']
test_bounds=dict()
test_bounds['a1'] = [0, 0]
test_bounds['a2'] = [0, 0]
test_bounds['q'] = [-np.inf, 3] 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build Catalogues

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'Building Training Catalogue'
print ''
print 'The principal component basis is constructed from these waveforms'
print ''
train_simulations = \
        bwave.simulation_details(series_names=train_series_names,
                param_bounds=train_bounds, Mmin30Hz=total_mass)

train_catalogue = bwave.waveform_catalogue(train_simulations,
        ref_mass=total_mass, SI_deltaT=SI_deltaT, SI_datalen=SI_datalen,
        distance=distance, trunc_time=True)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Do the PCA
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Performing PCA'
print ''
pca = bpca.waveform_pca(train_catalogue, train_catalogue)

# Characterise reconstructions
# XXX build a PSD and call projection_fidelity()
psd = aLIGOZeroDetHighPower(pca.SI_flen, pca.SI_deltaF, pca.fmin)
pca.compute_projection_fidelity(psd=psd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Plots
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Plotting results'
print ''


#
# Explained variance: Time Series
#

# PCA decompositions whose variances we will plot:
pc_types = ['NRhplusTimeSeriesPCA', 'NRAmpTimeSeriesPCA', \
        'NRPhaseTimeSeriesPCA']

# Dictionary to store the number of PCs required for 95% variance for each
# decomposition
VarThresh=0.99
npcVarThresh = dict()

fev, axev = pl.subplots()
xmax=-np.inf
npcs = np.arange(1,pca.ntrain+1)
for p, pcs in enumerate(pc_types):

    ev = np.cumsum(getattr(pca,pcs).explained_variance_ratio_)

    s=axev.step(npcs-0.5, ev, label=pcs.replace('NR','').replace('PCA',''))

    if npcs[np.argwhere(ev>VarThresh)[0][0]] > xmax:
        xmax=npcs[np.argwhere(ev>VarThresh)[0][0]]

    npcVarThresh[pcs] = npcs[np.argwhere(ev>VarThresh)[0][0]]
    axev.axvline(npcVarThresh[pcs],color=s[0].get_color(), linestyle='--')

axev.set_xlabel('Number of PCs')
axev.set_ylabel('explained variance [%]')
axev.legend(loc='lower right')
axev.minorticks_on()
axev.set_xticks(npcs)
axev.set_ylim(0.9,1)
axev.set_xlim(0.5,pca.ntrain+0.5)#xmax+0.5)

#
# Matches
#

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

#
# The principal components themselves
#

# --- hplus
NRtime = np.arange(0, train_catalogue.NR_datalen, train_catalogue.NR_deltaT)
NRtime -= NRtime[np.argmax(train_catalogue.NRAmpTimeSeries[0])]

# zoom on significant amplitude regime
xmin=-np.inf
xmax=np.inf
for w in xrange(pca.ntrain):
    idx = np.argwhere(train_catalogue.NRAmpTimeSeries[w,:] >
            1e-3*max(train_catalogue.NRAmpTimeSeries[w,:]))
    if NRtime[idx][0]>xmin:
        xmin=NRtime[idx][0]
    if NRtime[idx][-1]<xmax:
        xmax=NRtime[idx][-1]

for p,pcs in enumerate(pc_types):

    f, ax = pl.subplots()
    for p in xrange(npcVarThresh[pcs]):
        ax.plot(NRtime, getattr(pca, pcs).components_[p,:], label='PC %d'%(p+1))
        ax.set_xlim(xmin-100,xmax+50)
        #ax.set_xlim(-1000,500)
    ax.set_xlabel('Time [M$_{\odot}$]')
    ax.set_ylabel('PC amplitude for %s'%pcs.replace('NR','').replace('PCA',''))
    ax.minorticks_on()
    ax.legend(loc='upper left')
        

pl.show()


