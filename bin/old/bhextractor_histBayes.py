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
import os

from glob import glob

import numpy as np
import scipy.io as sio
import matplotlib
from matplotlib import pyplot as pl
import itertools

#font = {'family' : 'Helvetica',
#        'weight' : 'bold',
#        'size'   : 12}
#matplotlib.rc('font', **font)

fig_width_pt = 332  # Get this from LaTeX using \showthe\columnwidth
#fig_height_pt = 300 
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
#fig_height_pt = fig_width_pt*golden_mean
fig_height_pt = fig_width_pt*0.4
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_height_pt*inches_per_pt
fig_size =  [fig_width,fig_height]

matplotlib.rcParams.update(
        {'axes.labelsize': 8,
        'text.fontsize':   8,  
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# defs

def makebins(lowest,highest,increment):
    """
    """
    bins=[]
    edges=np.arange(lowest,highest+increment,increment)
    print edges
    for edge in edges:
        if edge<0:
            bins.append((-np.inf,edge))
        elif edge>0:
            bins.append((edge,np.inf))

    # insert a bin around zero
    bins.insert(np.argwhere(edges>0)[0],(edges[edges<0][-1],edges[edges>0][0]))

    labels=[]
    for b in bins:
        if abs(b[0])==abs(b[1]):
            labels.append('$[%d,%d]$'%(b[0],b[1]))
        elif b[0]==-np.inf:
            labels.append('$<%d$'%b[1])
        elif b[1]==np.inf:
            labels.append('$>%d$'%b[0])

    return bins,labels

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 
catalogues=['Q', 'HR', 'RO3']

# --- Injection run results
resultsfile=sys.argv[1]

num_pcs={'Q':2, 'HR':4, 'RO3':5}
num_wfs={'Q':10, 'HR':10, 'RO3':10}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data extraction & Plotting
# 

# The (unique) permuations for pairs of catalogues
catcombos=list(itertools.combinations(catalogues,2))

# --- Loading data from injection runs
results_data=np.load(resultsfile)
logB  = results_data['logB']

# Note log B is indexed by:
# [injection, injection #, recovery catalogue, noise realisation]

fig, ax = pl.subplots(nrows=1,ncols=3,sharey='row')
hatchstyles=['none','/','X']
colors={'Q':'b','HR':'g','RO3':'r'}

# Construct bins
bins,labels=makebins(-25,25,10)
width=0.5

# --- Loop through model selection combinations
for c,combo in enumerate(catcombos):
    print "Extracting Bayes factors for %s vs %s: "%(combo[0],combo[1])

    # find index for numerator evidence:
    numerator_index = [i for i, cat in enumerate(catalogues) if cat==combo[0]]
    denominator_index = [i for i, cat in enumerate(catalogues) if cat==combo[1]]

    # loop through injections
    for idx,inj in enumerate(combo):
        print "...%s injections"%inj
        injection_index = [i for i, cat in enumerate(catalogues) if cat==inj]

        # Extract the Bayes factors for all of these injection types
        current_logB = logB[injection_index,:,numerator_index,:] \
                - logB[injection_index,:,denominator_index,:]
        current_logB = np.concatenate(current_logB[0])

        # --- Get counts for SMEE fig 7 histograms
        ntotal=len(current_logB)
        counts=np.zeros(len(bins))
        for b, current_bin in enumerate(bins):
            counts[b] = sum( (current_logB>current_bin[0]) * 
                    (current_logB<current_bin[1]) )

        # --- Plotting
        ax[c].bar(np.arange(len(bins))-0.5*width,counts/ntotal, color=colors[inj],
                hatch=hatchstyles[idx], width=width, label='%s'%inj)
        ax[c].set_xticks(range(len(bins)))

    # --- Plot labels etc
    #ax[c].legend(frameon=False,loc='best')
    ax[c].legend(loc='upper center', bbox_to_anchor=(0.55, 1.25),
                      fancybox=False, frameon=False, ncol=2)
                      #fancybox=True, ncol=1)
    ax[c].set_ylim([0,1.0])
    #ax[c].minorticks_on()
    ax[c].set_xlabel('log B$_{\mathrm{%s,%s}}$'%(combo[0],combo[1]))
    xtickNames = pl.setp(ax[c], xticklabels=labels)
    pl.setp(xtickNames, rotation=45)
    #pl.show()
    #sys.exit()
ax[0].set_ylabel('Fractional Count')



#pl.subplots_adjust(hspace=0.2,wspace=0.2,bottom=0.2)
fig.tight_layout()
pl.subplots_adjust(hspace=0,wspace=0.05,bottom=0.3,top=0.85)
pl.savefig('bayesfactors.pdf')
pl.show()



