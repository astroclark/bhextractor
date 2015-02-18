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
"""

from __future__ import division

import sys

import numpy as np
import scipy.io as sio

import matplotlib
from matplotlib import pyplot as pl

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower



fig_width_pt = 800  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
fig_height_pt = fig_width_pt*golden_mean
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_height_pt*inches_per_pt
fig_size =  [fig_width,fig_height]

matplotlib.rcParams.update(
        {'axes.labelsize': 12,
        'text.fontsize':   12,  
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'text.usetex': True,
        'figure.figsize': fig_size,
        'font.family': "serif",
        'font.serif': ["Times"]
        })  

matplotlib.rcParams.update(
        {'savefig1.dpi': 200,
        'xtick.major.size':6,
        'xtick.minor.size':4,
        'ytick.major.size':6,
        'ytick.minor.size':4
        })
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

def compute_betas(waveforms, PCs):
    """
    """
    betas=np.array(np.matrix(waveforms).T * PCs)
    return betas


def comp_match(timeseries1, timeseries2, delta_t=1./2048, flow=10.,
        weighted=False):
    """
    """

    #tmp1 = pycbc.types.TimeSeries(initial_array=timeseries1, delta_t=delta_t)
    #tmp2 = pycbc.types.TimeSeries(initial_array=timeseries2, delta_t=delta_t)

    tmp1 = pycbc.types.FrequencySeries(initial_array=np.fft.fft(timeseries1),
            delta_f = 1.0 / len(timeseries1) / (1./delta_t))
    tmp2 = pycbc.types.FrequencySeries(initial_array=np.fft.fft(timeseries2),
            delta_f = 1.0 / len(timeseries2) / (1./delta_t))

    if weighted:

        # make psd
        flen = len(tmp1.to_frequencyseries())
        delta_f = np.diff(tmp1.to_frequencyseries().sample_frequencies)[0]
        psd = aLIGOZeroDetHighPower(flen, delta_f, low_freq_cutoff=flow)


        return pycbc.filter.match(tmp1, tmp2, psd=psd,
                low_frequency_cutoff=flow)[0]

    else:

        #return pycbc.filter.match(tmp1, tmp2, low_frequency_cutoff=flow)[0]
        return pycbc.filter.match(tmp1, tmp2, low_frequency_cutoff=0.0)[0]

def minimal_match_by_npc(waveforms,betas,PCs):
    """
    """
    match=np.zeros(shape=(np.shape(waveforms)[1],np.shape(waveforms)[1]))
    npcs=np.arange(np.shape(waveforms)[1])+1

    w=0
    for w in range(np.shape(waveforms)[1]):

        # --- Normalise the test waveform to unit norm
        target_wf = waveforms[:,w] / np.linalg.norm(waveforms[:,w])

        for npcidx, npc in enumerate(npcs):

            # --- Reconstruct the wf-th waveform using n PCs
            rec_wf = np.zeros(np.shape(waveforms)[0])
            for n in range(npc):
                rec_wf+=betas[w,n] * PCs[:,n]

            # --- Normalise the reconstructed waveform to unit norm
            #rec_wf /= np.linalg.norm(rec_wf)

            match[w,npcidx] = comp_match(rec_wf,target_wf)

    # --- Find minimal match as a function of nPC
    minimal_match = np.zeros(len(npcs))
    average_match = np.zeros(len(npcs))
    maximum_match = np.zeros(len(npcs))
    for npc, _ in enumerate(npcs):
        minimal_match[npc] = min(match[:,npc])
        average_match[npc] = np.mean(match[:,npc])
        maximum_match[npc] = max(match[:,npc])

    return minimal_match, average_match, maximum_match, npcs

def eigenenergy(eigenvalues):
    """
    """
    eigenvalues=abs(np.concatenate(eigenvalues))
    # See:
    # http://en.wikipedia.org/wiki/Principal_component_analysis#Compute_the_cumulative_energy_content_for_each_eigenvector
    gp = sum(eigenvalues)
    gj=np.zeros(len(eigenvalues))
    for e,i in enumerate(eigenvalues):
        for a in range(e+1):
            gj[e]+=eigenvalues[a]
    gj/=gp

    return gj

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 

labels=['Q','HR','RO3']
cols=['k','gray','k']
styles=['-','-','--']
# comma-separated lists of files
PC_files=sys.argv[1].split(',')
WF_files=sys.argv[2].split(',')

f1,ax1 = pl.subplots()
f2,ax2 = pl.subplots()
fcounter=0
for pc_file, wf_file in zip(PC_files,WF_files):

    PC_dict = sio.loadmat(pc_file)
    waveform_dict = sio.loadmat(wf_file)

    PCs = PC_dict['PCs_final']
    waveforms = waveform_dict['MDC_final']
    betas = PC_dict['Betas']

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the beta values by projecting data onto PCs
    #betas=compute_betas(waveforms,PCs)

    # XXX: NEED TO HANDLE COMPLEX MATCH???

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Compute match as a function of nPC
    minimal_match, average_match, maximum_match, npcs =\
            minimal_match_by_npc(waveforms,betas,PCs)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Compute cumulative energy content for each eigenvector
    #
    cum_energy = eigenenergy(PC_dict['EigVals'])
    npcs = range(1,len(cum_energy)+1)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Plot Results

    # --- Minimal Match
    ax1.plot(npcs, minimal_match, color=cols[fcounter],
            linestyle=styles[fcounter], marker='o', 
            label=labels[fcounter]+': min match')

    #ax1.plot(npcs, average_match, color=cols[fcounter],
    #        linestyle=styles[fcounter], marker='x',
    #        label=labels[fcounter]+': mean match')
    #ax1.plot(npcs, maximum_match, color=cols[fcounter],
    #        linestyle=styles[fcounter], marker='s',
    #        label=labels[fcounter]+': max match')

    # --- Eigenenergy
    ax2.plot(npcs, cum_energy, color=cols[fcounter], linestyle=styles[fcounter],
            marker='o', label=labels[fcounter])

    fcounter+=1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make Plots nice
ax1.set_xlabel('Number of PCs')
ax1.set_ylabel('Minimal Match')
#ax1.set_title('Match for worst-reconstruction in catalogue')
ax1.minorticks_on()
ax1.grid(which='major',color='grey',linestyle='-')
#ax1.grid(which='minor',color='grey',linestyle=':')
ax1.set_xlim(0.5,10)
ax1.set_ylim(0,1)
ax1.legend(loc='lower right')

f1.savefig('match.png')
#pl.close(f1)

ax2.set_xlabel('Number of PCs')
#ax2.set_ylabel('Cumulative Eigenvector Energy')
#ax2.set_title('Distribution of Energy Among Eigenvectors')
ax2.set_ylabel('$E(k)$')
ax2.minorticks_on()
ax2.grid(which='major',color='grey',linestyle='-')
#ax2.grid(which='minor',color='grey',linestyle=':')
ax2.set_xlim(0.5,10)
ax2.legend(loc='lower right')
#f2.tight_layout()
pl.subplots_adjust(bottom=0.1,top=0.95,left=0.1)


f2.tight_layout()
f2.savefig('eigenenergy.png')
f2.savefig('eigenenergy.pdf')

#pl.show()
