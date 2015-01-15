#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2014-2015 James Clark <james.clark@ligo.org>
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
Generate a catalogue with 3 waveforms:
    1) Gaussian pulse
    2) Chirp
    3) IMR
"""

from __future__ import division

import os
import sys 

__author__ = "James Clark <james.clark@ligo.org>"

import numpy as np
import scipy.signal as signal
import scipy.stats as stats
import scipy.io as sio 

import lal 
import lalsimulation as lalsim

import pycbc
import pycbc.filter

def truncparms(low,upp,mu,sigma):
    a = (low - mu) / sigma
    b = (upp - mu) / sigma
    return a, b

#
# Common Variables
#
Fs=1024
datalen=1.
time_axis = np.arange(-0.5*datalen, 0.5*datalen, 1./Fs)
Nwaveforms = 1000
Ntypes = 2

catalogue = np.zeros(shape=(len(time_axis), Ntypes*Nwaveforms))

#
# Sine-Gaussians ('Gaussian Pulse')
#

# center frequency
fc_mu=500
fc_sigma=100
fc_low=10
fc_upp=1000
a, b = truncparms(fc_low, fc_upp, fc_mu, fc_sigma)
fcs = stats.truncnorm.rvs(a, b, loc=fc_mu, scale=fc_sigma, size=Nwaveforms)

# bandwidth
bw_mu=0.005
bw_sigma=0.001
bw_low=0.0001
bw_upp=0.1
a, b = truncparms(bw_low, bw_upp, bw_mu, bw_sigma)
bws = stats.truncnorm.rvs(a, b, loc=bw_mu, scale=bw_sigma, size=Nwaveforms)

# populate catalogue
win = lal.CreateTukeyREAL8Window(len(time_axis), 0.1)
for i in xrange(Nwaveforms):
    catalogue[:,i] = win.data.data*signal.gausspulse(time_axis, fcs[i], bws[i])
    catalogue[:,i] /= np.sqrt(np.dot(catalogue[:,i], catalogue[:,i]))


#
# Chirps
#
chirp_times = np.arange(0, datalen, 1.0/Fs)
idx = chirp_times<=0.5
win = lal.CreateTukeyREAL8Window(int(sum(idx)), 0.1)

tau=0.5
f0=100
t1=0.5
f1_mu = 500
f1_sigma = 100
f1_low = 10
f1_upp = 1000
a, b = truncparms(f1_low, f1_upp, f1_mu, f1_sigma)
f1s = stats.truncnorm.rvs(a, b, loc=f1_mu, scale=f1_sigma, size=Nwaveforms)
phis = 180*np.random.random(Nwaveforms)

for i in xrange(Nwaveforms):
    catalogue[:,i+Nwaveforms] = signal.chirp(chirp_times, f0 = f0, t1 = t1, f1 = f1s[i],
            phi=phis[i])
    catalogue[:,i+Nwaveforms] *= np.exp(chirp_times / tau)
    catalogue[:,i+Nwaveforms][chirp_times>0.5] = 0.0
    catalogue[:,i+Nwaveforms][idx]*=win.data.data
    catalogue[:,i] /= np.sqrt(np.dot(catalogue[:,i], catalogue[:,i]))

#catalogue += 0.1*np.random.randn(Fs*datalen, Ntypes*Nwaveforms)


if 0:
    #
    # PCAT magic: Lifting the following from GMM.py in PCAT
    #
    import PCA, GMM

    score_matrix, principal_components, means, stds, eigenvalues = \
            PCA.PCA(catalogue, components_number=10)

    principal_components_number=10

    reduced_score_matrix = score_matrix[:,:principal_components_number]

    mat, tmp, tmp1 = PCA.matrix_whiten(reduced_score_matrix, std=True)

    #labels = GMM.gaussian_mixture(mat,upper_bound=5)
    labels = GMM.gaussian_mixture(reduced_score_matrix,upper_bound=5)

    colored_clusters = GMM.color_clusters( score_matrix, labels )

    GMM.print_cluster_info(colored_clusters)


#sys.exit()
#
# PCA
#
#H = np.matrix(waveform_catalogue)
H = np.matrix(catalogue)

# --- 1) Compute catalogue covariance
C = H.T * H / np.shape(H)[0]

# --- 2) Compute eigenvectors (V) and eigenvalues (S) of C
S, V = np.linalg.eigh(C)

# --- 3) Sort eigenvectors in descending order
idx = np.argsort(S)[::-1]
V = V[:,idx]
S = S[idx]

# --- 4) Compute the eigenvectors of the real covariance matrix U and normalise
U = H*V 
for i in xrange(Ntypes*Nwaveforms):
    U /= np.linalg.norm(U[:,i])

# Coeffs
coeffs = H.T*U

# 'scores'
Z = np.dot(catalogue, coeffs)

# whiten
for i in xrange(np.shape(catalogue)[0]):
    Z[i,:] -= np.mean(Z[i,:])
    Z[i,:] /= np.std(Z[i,:])

# reduce
Z = Z[:50,:]

#
# GMM
#
from sklearn import mixture

niterations=10
print ''
for i in xrange(niterations):
    bic=np.zeros(5)
    min_bic=np.inf
    for c in xrange(5):
        #print c
        g = mixture.GMM(n_components=c+1)
        g.fit(Z)
        bic[c] = g.bic(Z)
        
        if bic[c] < min_bic:
            min_bic=bic[c]
            best_g = g

    best_g.fit(Z)
    labels=best_g.predict(Z)

    nclusters=np.argmin(bic)+1
    print 'favoured # of clusters:', nclusters

    print 'Cluster membership:'
    for label in xrange(nclusters):
        print 'cluster %d: %.2f percent'%(label,
                float(sum(labels==label))/len(labels) )



