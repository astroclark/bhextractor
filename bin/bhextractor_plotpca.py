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
bhextractor_plotpca.py

Construct waveform catalogues and PCA for plotting and diagnostics
"""
import numpy as np
from matplotlib import pyplot as pl
import bhextractor_pca as bhex

def image_matches(match_matrix, waveform_names, title=None, mismatch=False):
    """
    Make a Nice plot.
    """

    # XXX: THIS SHOULD NOT LIVE HERE!!! Move it to a plotting script which
    # builds the classes in this module!

    from matplotlib import pyplot as pl

    if mismatch:
        match_matrix = 1-match_matrix
        text_thresh = 0.1
        clims = (0,0.2)
        bar_label = 'mismatch'
    else:
        text_thresh = 0.85
        clims = (0.75,1.0)
        bar_label = 'match'

    #fig, ax = pl.subplots(figsize=(15,8))
    #fig, ax = pl.subplots(figsize=(8,4))
    fig, ax = pl.subplots()
    nwaves = np.shape(match_matrix)[0]
    npcs = np.shape(match_matrix)[1]

    im = ax.imshow(match_matrix, interpolation='nearest', origin='lower',
            aspect='auto')

    for x in xrange(nwaves):
        for y in xrange(npcs):
            if match_matrix[x,y]<text_thresh:
                ax.text(y, x, '%.2f'%(match_matrix[x,y]), \
                    va='center', ha='center', color='w')
            else:
                ax.text(y, x, '%.2f'%(match_matrix[x,y]), \
                    va='center', ha='center', color='k')

    ax.set_xticks(range(0,npcs))
    ax.set_yticks(range(0,nwaves))

    xlabels=range(1,npcs+1)
    ax.set_xticklabels(xlabels)

    ax.set_yticklabels(waveform_names)

    im.set_clim(clims)
    im.set_cmap('gnuplot2')

    ax.set_xlabel('Number of PCs')
    ax.set_ylabel('Waveform type')

    if title is not None:
        ax.set_title(title)

    #c=pl.colorbar(im, ticks=np.arange(clims[0],clims[1]+0.05,0.05),
    #        orientation='horizontal')
    #c.set_label(bar_label)

    fig.tight_layout()

    return fig, ax

def image_euclidean(euclidean_matrix, waveform_names, title=None, clims=None):

    #clims = (0.0,0.10)
    if clims is None:
        clims = (0.0, euclidean_matrix.max())
    text_thresh = 0.5*max(clims)
    bar_label = '$||\Phi - \Phi_r||$'

    #fig, ax = pl.subplots(figsize=(15,8))
    #fig, ax = pl.subplots(figsize=(7,7))
    fig, ax = pl.subplots()
    nwaves = np.shape(euclidean_matrix)[0]
    npcs = np.shape(euclidean_matrix)[1]

    im = ax.imshow(euclidean_matrix, interpolation='nearest', origin='lower',
            aspect='auto')

    for x in xrange(nwaves):
        for y in xrange(npcs):
            if euclidean_matrix[x,y]<text_thresh:
                ax.text(y, x, '%.2f'%(euclidean_matrix[x,y]), \
                    va='center', ha='center', color='k')
            else:
                ax.text(y, x, '%.2f'%(euclidean_matrix[x,y]), \
                    va='center', ha='center', color='w')

    ax.set_xticks(range(0,npcs))
    ax.set_yticks(range(0,nwaves))

    xlabels=range(1,npcs+1)
    ax.set_xticklabels(xlabels)

    ylabels=[name.replace('_lessvisc','') for name in waveform_names]
    ax.set_yticklabels(ylabels)

    im.set_clim(clims)
    im.set_cmap('gnuplot2_r')

    ax.set_xlabel('Number of PCs')
    ax.set_ylabel('Waveform type')

    if title is not None:
        ax.set_title(title)

    #c=pl.colorbar(im, ticks=np.arange(clims[0],clims[1]+0.05,0.05),
    #        orientation='horizontal')
    #c.set_label(bar_label)

    fig.tight_layout()

    return fig, ax


# -------------------------------
# USER INPUT

catalogue_name='Q'
theta=90.0

# END USER INPUT
# -------------------------------

# -------------------------------
# ANALYSIS

#
# Setup and then build the catalogue
#
catalogue = bhex.waveform_catalogue(catalogue_name=catalogue_name, fs=2048,
        catalogue_len=4, mtotal_ref=250, Dist=1., theta=theta)

#
# Do the PCA
#
pca = bhex.waveform_pca(catalogue)

#
# Compute matches
#
pca.compute_matches()


# -------------------------------
# PLOTTING

#
# Waveform Catalogue
#

print "Plotting Catalogue"

fig, ax = pl.subplots(figsize=(10,10),
        nrows=np.shape(catalogue.aligned_catalogue)[0], ncols=3)

for w in xrange(len(catalogue.waveform_names)):

    ax[w][0].plot(catalogue.sample_times,
            np.real(catalogue.aligned_catalogue[w,:]), label='+')
    ax[w][0].plot(catalogue.sample_times,
            np.imag(catalogue.aligned_catalogue[w,:]), label='x')

    ax[w][0].set_xlim(-1.0,0.1)

    ax[w][1].plot(catalogue.sample_frequencies, catalogue.ampSpectraPlus[w,:],
            label='+')
    ax[w][1].plot(catalogue.sample_frequencies, catalogue.ampSpectraCross[w,:],
            label='x')
    ax[w][1].set_xlim(9,128)
    ax[w][1].set_ylim(1e-5,1e-2)
    ax[w][1].axvline(10,color='k',linestyle='--')
    ax[w][1].set_yscale('log')

    ax[w][2].plot(catalogue.sample_frequencies, catalogue.phaseSpectraPlus[w,:],
            label='+')
    ax[w][2].plot(catalogue.sample_frequencies, catalogue.phaseSpectraCross[w,:],
            label='x')
    ax[w][2].set_xlim(9,128)
    ax[w][2].set_ylim(-512,512)
    ax[w][2].axvline(10,color='k',linestyle='--')

fig.tight_layout()

#
# Principal Components
#

print "Plotting PCs"

fig, ax = pl.subplots(figsize=(10,10),
        nrows=np.shape(catalogue.aligned_catalogue)[0]+1, ncols=1)

ax[0].plot(catalogue.sample_times, pca.pca_plus.mean_, label='Mean')
ax[0].set_xlim(-1.0,0.1)
ax[0].legend(loc='upper left')

for w in xrange(len(catalogue.waveform_names)):

    ax[w+1].plot(catalogue.sample_times,
         pca.pca_plus.components_[w,:], label=r'$\beta_{%d}$'%(w+1))

    ax[w+1].legend(loc='upper left')
    ax[w+1].set_xlim(-1.0,0.1)


fig.tight_layout()


#
# Explained Variance
#
fe, ax = pl.subplots()
ax.plot(range(1,len(catalogue.waveform_names)+1),
        1-pca.pca_plus.explained_variance_ratio_)
ax.set_xlabel('Number of PCs')
ax.set_ylabel('Explained Variance')
ax.set_ylim(0,1)
ax.set_xlim(1,len(catalogue.waveform_names)+1)

#
# Matches for 250 Msun
#
imfig, imax = image_matches(pca.matches, catalogue.waveform_names)
imfig.tight_layout()

#
# Euclidean distances for 250 Msun
#
imfig, imax = image_euclidean(pca.euclidean_distances, catalogue.waveform_names)
imfig.tight_layout()

pl.show()




