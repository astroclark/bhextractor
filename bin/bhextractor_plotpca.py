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
# Compute matches
#
pca.compute_matches()


# -------------------------------
# PLOTTING

#
# Waveform Catalogue
#

print "Plotting Training Catalogue"

fig1, ax1 = pl.subplots(figsize=(10,10),
        nrows=np.shape(train_catalogue.aligned_catalogue)[0], ncols=1, sharex='col')

fig2, ax2 = pl.subplots(figsize=(10,10),
        nrows=np.shape(train_catalogue.aligned_catalogue)[0], ncols=1, sharex='col')

fig3, ax3 = pl.subplots(figsize=(10,10),
        nrows=np.shape(train_catalogue.aligned_catalogue)[0], ncols=1, sharex='col')

#fig4, ax4 = pl.subplots(figsize=(10,10),
#        nrows=np.shape(train_catalogue.aligned_catalogue)[0], ncols=1, sharex='col')

for w in xrange(len(train_catalogue.waveform_names)):

    ax1[w].plot(train_catalogue.sample_times,
            np.real(train_catalogue.aligned_catalogue[w,:]), label='+')
    ax1[w].plot(train_catalogue.sample_times,
            np.imag(train_catalogue.aligned_catalogue[w,:]), label='x')

    ax1[w].set_xlim(-1.5,0.25)

    ax2[w].plot(train_catalogue.sample_frequencies, train_catalogue.ampSpectraPlus[w,:],
            label='+')
    ax2[w].plot(train_catalogue.sample_frequencies, train_catalogue.ampSpectraCross[w,:],
            label='x')
    ax2[w].set_xlim(0,128)
    ax2[w].set_ylim(1e-5,1e-2)
    ax2[w].axvline(10,color='k',linestyle='--')
    ax2[w].set_yscale('log')

    ax3[w].plot(train_catalogue.sample_times_ampalign,
            train_catalogue.aligned_amplitudes[w,:], label='|h(t)|')
    ax3[w].set_xlim(-1.5,0.25)

#    ax4[w].plot(train_catalogue.sample_times_ampalign, train_catalogue.aligned_phases[w,:],
#            label='arg[h(t)]')
#    ax4[w].set_xlim(-1.5,0.25)

fig1.suptitle('h$_{+,\\times}$(t)')
fig1.subplots_adjust(hspace=0)
fig1.savefig('train_catalogue_timeseries.png')

fig2.suptitle('$|\\tilde{H}_{+,\\times}|$')
fig2.subplots_adjust(hspace=0)
fig2.savefig('train_catalogue_amplitudespectra.png')

fig3.suptitle('|h(t)|')
fig3.subplots_adjust(hspace=0)
fig3.savefig('train_catalogue_timeseries_amp.png')

#   fig4.suptitle('arg[h(t)]')
#   fig4.subplots_adjust(hspace=0)
#   fig4.savefig('train_catalogue_timeseries_phase.png')

pl.show()

import sys
sys.exit()
sys.exit()

#
# Principal Components
#

print "Plotting PCs"

fig, ax = pl.subplots(figsize=(10,10),
        nrows=np.shape(train_catalogue.aligned_catalogue)[0]+1, ncols=1, sharex='col')

fig2, ax2 = pl.subplots(figsize=(10,10),
        nrows=np.shape(train_catalogue.aligned_catalogue)[0]+1, ncols=1, sharex='col')

fig3, ax3 = pl.subplots(figsize=(10,10),
        nrows=np.shape(train_catalogue.aligned_catalogue)[0]+1, ncols=1, sharex='col')

ax[0].plot(train_catalogue.sample_times, pca.pca_plus.mean_, label='Mean')
ax[0].set_xlim(-1.5,0.1)
ax[0].legend(loc='upper left')

ax2[0].plot(train_catalogue.sample_times_ampalign, pca.pca_amp.mean_, label='Mean')
ax2[0].set_xlim(-1.5,0.1)
ax2[0].legend(loc='upper left')

ax3[0].plot(train_catalogue.sample_times_ampalign, pca.pca_phase.mean_, label='Mean')
ax3[0].set_xlim(-1.5,0.1)
ax3[0].legend(loc='upper left')

for w in xrange(len(train_catalogue.waveform_names)):

    ax[w+1].plot(train_catalogue.sample_times,
         pca.pca_plus.components_[w,:], label=r'$\beta_{%d}$'%(w+1))

    ax2[w+1].plot(train_catalogue.sample_times_ampalign,
         pca.pca_amp.components_[w,:], label=r'$\beta_{%d}$'%(w+1))

    ax3[w+1].plot(train_catalogue.sample_times_ampalign,
         pca.pca_phase.components_[w,:], label=r'$\beta_{%d}$'%(w+1))

    ax[w+1].legend(loc='upper left')
    ax[w+1].set_xlim(-1.5,0.1)

    ax2[w+1].legend(loc='upper left')
    ax2[w+1].set_xlim(-1.5,0.1)

    ax3[w+1].legend(loc='upper left')
    ax3[w+1].set_xlim(-1.5,0.1)

fig.suptitle(r'h${_+,\times}$(t)')
fig.subplots_adjust(hspace=0)
fig.savefig('hplus_pcs.png')

fig2.suptitle(r'|h(t)|')
fig2.subplots_adjust(hspace=0)
fig2.savefig('amplitude_pcs.png')

fig3.suptitle(r'arg[h(t)]')
fig3.subplots_adjust(hspace=0)
fig3.savefig('phase_pcs.png')


#
# Explained Variance
#
#   fe, ax = pl.subplots()
#   ax.step(np.arange(1,len(train_catalogue.waveform_names)+1)-0.5,
#           np.cumsum(pca.pca_plus.explained_variance_ratio_), color='k',
#           label='h$_+$')
#   ax.step(np.arange(1,len(train_catalogue.waveform_names)+1)-0.5,
#           np.cumsum(pca.pca_plus.explained_variance_ratio_), color='grey',
#           label=r'h$_{\times}$', linestyle='--')
#   ax.set_xlabel('Number of PCs')
#   ax.set_ylabel('Explained Variance')
#   ax.set_ylim(0,1)
#   ax.set_xlim(1,len(train_catalogue.waveform_names)+1)
#   ax.grid()
#   ax.minorticks_on()
#   ax.legend(loc='lower right')
#   fe.savefig('explained_variance_pluscross.png')

#
# Explained Variance
#
fe, ax = pl.subplots()
ax.step(np.arange(1,len(train_catalogue.waveform_names)+1)-0.5,
        np.cumsum(pca.pca_amp.explained_variance_ratio_), color='k',
        label='|h(t)|')
ax.step(np.arange(1,len(train_catalogue.waveform_names)+1)-0.5,
        np.cumsum(pca.pca_phase.explained_variance_ratio_), color='grey',
        linestyle='--', label='arg[h(t)]')
ax.set_xlabel('Number of PCs')
ax.set_ylabel('Explained Variance')
ax.set_ylim(0,1)
ax.set_xlim(1,len(train_catalogue.waveform_names)+1)
ax.grid()
ax.minorticks_on()
ax.legend(loc='lower right')
fe.savefig('explained_variance_ampphase.png')

#
# Matches for 250 Msun
#
#   imfig, imax = image_matches(pca.matches_plus, test_catalogue.waveform_names)
#   imax.set_title(r"Matches for +/$\times$ PCA (%s waveforms trained with %s)"%(\
#           test_catalogue_name, train_catalogue_name))
#   imfig.tight_layout()
#   imfig.savefig('plus_matches.png')

imfig, imax = image_matches(pca.matches_ampphase, test_catalogue.waveform_names)
imax.set_title(r"Matches for A,$\phi$ PCA (%s waveforms trained with %s)"%(\
        test_catalogue_name, train_catalogue_name))
imfig.tight_layout()
imfig.savefig('ampphase_matches.png')

#
# Euclidean distances for 250 Msun
#
#   imfig, imax = image_euclidean(pca.euclidean_distances_plus, test_catalogue.waveform_names)
#   imax.set_title(r"Euclidean Distances for +/$\times$ PCA (%s waveforms trained with%s)"%(\
#           test_catalogue_name, train_catalogue_name))
#   imfig.tight_layout()
#   imfig.savefig('plus_euclidean_distances.png')

pl.show()




