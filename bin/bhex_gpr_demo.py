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
"""

import sys
from bhex_utils import bhex_wavedata as bwave
from bhex_utils import bhex_pca as bpca
from pycbc.psd import aLIGOZeroDetHighPower
import pycbc.types
import numpy as np
from matplotlib import pyplot as pl

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Training data input
# XXX: probably going to want a config parser

SI_deltaT = 1.0/512
SI_datalen= 4.0

total_mass = 150. 
distance=1. # Mpc

train_series_names = ['HR-series']

train_bounds=dict()
train_bounds['a1'] = [0, 0]
train_bounds['a2'] = [0, 0]
#train_bounds['q'] = [-np.inf, 3] 
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
train_simulations = \
        bwave.simulation_details(series_names=train_series_names,
                param_bounds=train_bounds, Mmin30Hz=total_mass)


print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
train_catalogue = bwave.waveform_catalogue(train_simulations,
        ref_mass=total_mass, SI_deltaT=SI_deltaT, SI_datalen=SI_datalen,
        distance=distance)

#
# Do the PCA
#
print '~~~~~~~~~~~~~~~~~~~~~'
print 'Performing PCA'
print ''
pca = bpca.waveform_pca(train_catalogue, train_catalogue)


#
# Extract the PC coefficients for each waveform
#
mass_ratios = np.array([pca.test_catalogue_data[w]['q'] for w in
    xrange(pca.ntest)])

beta1_measured = np.array([pca.test_catalogue_data[w]['NRAmpTimeSeriesBetas'][0]
    for w in xrange(pca.ntest)])
beta2_measured = np.array([pca.test_catalogue_data[w]['NRAmpTimeSeriesBetas'][1] for
    w in xrange(pca.ntest)])
beta3_measured = np.array([pca.test_catalogue_data[w]['NRAmpTimeSeriesBetas'][3] for
    w in xrange(pca.ntest)])

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# GP fitting experiments
#

# Get initial polynomial fit
pfit_beta1 = np.polyfit(mass_ratios, beta1_measured, deg=3)

def poly_beta(x,pfit):
    """
    Polynomial function, given parameters pfit
    """
    deg = len(pfit)-1

    y = pfit[-1]
    for n in xrange(deg):
        power = deg - n
        #print n, power
        y += pfit[n]*x**power

    return y
    
mass_ratios_fit = np.arange(mass_ratios.min(), mass_ratios.max(), 0.01)
beta1_fit = poly_beta(mass_ratios_fit, pfit_beta1)

from sklearn import gaussian_process

# jiggery-pokery
X = np.atleast_2d(mass_ratios).T
x = np.atleast_2d(mass_ratios_fit).T
y = np.atleast_2d(beta1_measured).T


#gp = gaussian_process.GaussianProcess(theta0=1e-2, thetaL=1e-4, thetaU=1e-1)
gp = gaussian_process.GaussianProcess(corr='cubic', theta0=1e-2, thetaL=1e-4, thetaU=1e-1,
                             random_start=100)
gp.fit(X, y)  
y_pred, sigma2_pred = gp.predict(x, eval_MSE=True)
sigma = np.sqrt(sigma2_pred)


fig = pl.figure()
pl.plot(mass_ratios_fit, beta1_fit, 'r:', label=u'$f(x) = ax^3 + bx^2 +cx + d$')
pl.plot(mass_ratios, beta1_measured, 'r.', markersize=10, label=u'Observations')
pl.plot(x, y_pred, 'b-', label=u'Prediction')

low = np.concatenate(y_pred) - 1.9600 * sigma
upp = np.concatenate(y_pred) + 1.9600 * sigma

pl.fill_between(mass_ratios_fit, y1=low, y2=upp,
        alpha=.5,label='95% confidence interval')

pl.legend(loc='upper left')

pl.xlabel('Mass Ratio (q)')
pl.ylabel('First PC coefficient')

pl.minorticks_on()

pl.show()



