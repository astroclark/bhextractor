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
bhextractor_plotpca.py

Construct waveform catalogues and PCA for plotting and diagnostics
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as pl
import bhextractor_pca as bhex


# -------------------------------
# USER INPUT

train_catalogue_name='Q'
test_catalogue_name='Q'

# END USER INPUT
# -------------------------------

# -------------------------------
# ANALYSIS

#
# Setup and then build the catalogue
#
train_catalogue = bhex.waveform_catalogue(catalogue_name=train_catalogue_name, fs=512,
        catalogue_len=4, mtotal_ref=250, Dist=1.)


test_catalogue = bhex.waveform_catalogue(catalogue_name=test_catalogue_name, fs=512,
        catalogue_len=4, mtotal_ref=250, Dist=1.)

#
# Do the PCA
#
pca = bhex.waveform_pca(training_catalogue=train_catalogue,
        testing_catalogue=test_catalogue)


#
# Find the coefficients for projecting each test waveform onto the PC basis
# formed from the training data
#
betas_hplus = np.zeros(shape=(test_catalogue.nwaves, train_catalogue.nwaves))
betas_amp = np.zeros(shape=(test_catalogue.nwaves, train_catalogue.nwaves))
betas_phase = np.zeros(shape=(test_catalogue.nwaves, train_catalogue.nwaves))
for n, name in enumerate(test_catalogue.waveform_names):
    betas_hplus[n,:] = pca.projection_plus[name]
    betas_amp[n,:] = pca.projection_amp[name]
    betas_phase[n,:] = pca.projection_phase[name]

# -------------------------------
# Plotting

#
#
#



