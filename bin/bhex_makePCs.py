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
bhextractor_makePCs.py

Construct catalogues and principal component analysis for NR BBH waveforms; save
the PCA results to file using the file dump method on the PCA result object

This is a pretty simple script but it might be preferable to set up the
configuration in a config parser & ini file.
"""

import sys
from bhex_utils import bhex_wavedata as bwave
from bhex_utils import bhex_pca as bpca

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Some useful info:
#

#valid_series = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
#        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series"
#        "TP-series"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

# Options common to train & test data:
SI_deltaT = 1./512
SI_datalen= 4.0
total_mass = 150. 
distance=1. # Mpc

#train_series_names = ['HRq-series'] # (see above for valid choices)
train_series_names = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series",
        "TP-series"]
#train_series_names = ['HRq-series'] # (see above for valid choices)

#
# Modify for imposing parameter bounds on the catalogue:
#
#train_bounds=None
train_bounds=dict()
train_bounds['a1'] = [0, 0]
train_bounds['a2'] = [0, 0]
#train_bounds['q'] = [-np.inf, 3] 

catalogue_name = 'NonSpinning'

save_pcs = ['NRhplusTimeSeriesPCA', 'NRhcrossTimeSeriesPCA',
'NRAmpTimeSeriesPCA', 'NRPhaseTimeSeriesPCA',  'SIhplusTimeSeriesPCA',
'SIhcrossTimeSeriesPCA', 'SIAmpTimeSeriesPCA', 'SIPhaseTimeSeriesPCA']
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do the work:

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
train_simulations = \
        bwave.simulation_details(series_names=train_series_names,
                param_bounds=train_bounds, Mmin30Hz=total_mass)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
train_catalogue = bwave.waveform_catalogue(train_simulations, ref_mass=total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance)

# Dump catalogue to pickle
train_catalogue.file_dump(catalogue_name=catalogue_name)

#
# Do the PCA
#
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'Performing PCA (training only)'
print ''
pca = bpca.waveform_pca(train_catalogue)
pca.file_dump(pca_attrs=save_pcs, pcs_filename=catalogue_name)

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'DONE.'



