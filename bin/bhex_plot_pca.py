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
sample_rate = 512
SI_datalen= 4.0
total_mass = 150. 
distance=1. # Mpc

# Define the train catalogue
train_series_names = ['HR-series']
train_bounds=dict()
train_bounds['a1'] = [0, 0]
train_bounds['a2'] = [0, 0]
train_bounds['q'] = [-np.inf, 3] 

# Define the test catalogue
test_series_names = ['HR-series']
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
        ref_mass=total_mass, sample_rate=sample_rate, SI_datalen=SI_datalen,
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
        sample_rate=sample_rate, SI_datalen=SI_datalen, distance=distance)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Do the PCA
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Performing PCA'
print ''
pca = bpca.waveform_pca(train_catalogue, test_catalogue)

# Characterise reconstructions
# XXX build a PSD and call projection_fidelity()
psd = aLIGOZeroDetHighPower(





