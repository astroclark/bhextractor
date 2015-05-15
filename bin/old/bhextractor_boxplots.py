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
import matplotlib.ticker as mtick
from matplotlib.patches import Polygon

font = {'family' : 'Helvetica',
        'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)

import lal

def model_wf(pcs, *betas):
    """
    Construct a template from the sum of PCs weighted by the coefficients in
    betas
    """
    model_wf = np.zeros(np.shape(pcs)[0])
    for n in range(len(betas)):
        model_wf+=betas[n] * pcs[:,n] 
    return model_wf 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 
catalogues=['Q', 'HR', 'RO3']

num_pcs={'Q':2, 'HR':4, 'RO3':5}
num_wfs={'Q':10, 'HR':10, 'RO3':10}

bhexpath=os.getenv('BHEX_PREFIX')

catalogue_files={'Q':'%s/data/signal_data/Q_catalogue_theta-0.mat'%(bhexpath),
        'HR':'%s/data/signal_data/HR_catalogue_theta-0.mat'%(bhexpath),
        'RO3':'%s/data/signal_data/RO3_catalogue_theta-0.mat'%(bhexpath)}

pc_files={'Q':'%s/data/PCA_data/Q_PCs_theta-0.mat'%(bhexpath),
        'HR':'%s/data/PCA_data/HR_PCs_theta-0.mat'%(bhexpath),
        'RO3':'%s/data/PCA_data/RO3_PCs_theta-0.mat'%(bhexpath)}


# --- Injection run results
resultsfile=sys.argv[1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sanity
if not bhexpath:
    print >> sys.stderr, "ERROR, must set BHEX_PREFIX"
    sys.exit()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data extraction & Plotting

# --- Labels for the different waveforms
Q_labels=['q=1.0','q=1.15','q=1.3','q=1.45','q=1.5','q=1.6',
        'q=1.75', 'q=1.9', 'q=2.0', 'q=2.05']

HR_labels=['q=1.0,\na=0.0','q=1.0,\na=0.1','q=1.0,\n a=0.2','q=1.0,\n a=0.3', \
        'q=1.0,\n a=0.4', 'q=1.0,\n a=0.5','q=1.0,\n a=0.6','q=1.0,\n a=0.7',\
        'q=1.0,\n a=0.8','q=1.5,\n a=0.0']

RO3_labels=['q=1.5,\n a=0.4,\n $\Theta$=60$^{\circ}$', \
        'q=1.5,\n a=0.6,\n $\Theta$=45$^{\circ}$',
        'q=1.5,\n a=0.6,\n $\Theta$=60$^{\circ}$',
        'q=1.5,\n a=0.6,\n $\Theta$=90$^{\circ}$',
        'q=2.0,\n a=0.4,\n $\Theta$=60$^{\circ}$',
        'q=2.0,\n a=0.6,\n $\Theta$=45$^{\circ}$',
        'q=2.0,\n a=0.6,\n $\Theta$=60$^{\circ}$',
        'q=2.0,\n a=0.6,\n $\Theta$=90$^{\circ}$',
        'q=2.0,\n a=0.6,\n $\Theta$=135$^{\circ}$',
        'q=2.0,\n a=0.6,\n $\Theta$=180$^{\circ}$']
wf_labels={'Q':Q_labels, 'HR':HR_labels, 'RO3':RO3_labels}

# --- Loading data from injection runs
results_data=np.load(resultsfile)
logB = results_data['logB']
match = results_data['match']

# --- Load targets and PCs first to save overheads
target_waveforms=[]
princ_comps=[]
for wf in catalogues:
    target_waveforms.append(sio.loadmat(catalogue_files[wf])['MDC_final'])
    princ_comps.append(sio.loadmat(pc_files[wf])['PCs_final'])

# --- Loop through injected catalogues (subplots) and plot each recovered
# catalogue's responses against the injected waveform in that subplot

cols=np.array(['black','grey'])

for i,inj in enumerate(catalogues):

    # --- Seperate figures for each injection catalogue
    f1,ax1=pl.subplots(figsize=(10,5)) # for bar plots of means
    f2,ax2=pl.subplots(figsize=(10,5)) # for box/whisker plots

    # x-locations for boxplot
    xlocs=np.arange(len(catalogues)-1)+1
    xticklocs_boxplot=np.zeros(num_wfs[inj])

    xlocs2=np.arange(len(catalogues))+1
    xticklocs_boxplot2=np.zeros(num_wfs[inj])

    for w in range(num_wfs[inj]):

        # For this injection catalogue and waveform #, pre-allocate array for
        # odds of this injection catalogue vs the other catalogues
        logB_sigvsig = np.zeros(shape=(len(catalogues),len(logB[0][0][0][:])))
        match_tmp = np.zeros(shape=(len(catalogues),len(logB[0][0][0][:])))

        # Loop through recovery catalogues
        for c,rec in enumerate(catalogues):

            # --- Bayes factor for injected catalogue vs recovered catalogues
            logB_sigvsig[c] = (logB[i][w][c][:]-logB[i][w][i][:])/np.log(10)
            match_tmp[c]    = match[i][w][c][:]

        # --- don't use the results of injected model
        logB_sigvsig = np.delete(logB_sigvsig,np.s_[i],0)

        # --- Plot the mean values as a bar on the injection subplot
        #mean_logB_sigvsig = np.median(logB_sigvsig,axis=1)
        #delta_logB_sigvsig = \
        #        np.std(logB_sigvsig,axis=1)/np.sqrt(len(logB_sigvsig))
        #sortidx=np.argsort(abs(mean_logB_sigvsig))[::-1]

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Bars
        #ax1.bar(w*np.ones(len(catalogues)-1), height=mean_logB_sigvsig[sortidx],
        #        color=cols[sortidx], width=0.5,
        #        yerr=delta_logB_sigvsig[sortidx], ecolor='r')

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Bayes Factor Boxes
        xlocs+=4
        boxdat=logB_sigvsig.T
        bp=ax1.boxplot(boxdat, positions=xlocs, sym='', widths=.75)

        pl.setp(bp['boxes'], color='black')
        pl.setp(bp['whiskers'], color='black')
        pl.setp(bp['medians'], color='red')

        xticklocs_boxplot[w]=np.mean(xlocs)

        # Color the box faces
        # -------------------
        # Lifted from: http://matplotlib.org/examples/pylab_examples/boxplot_demo2.html
        boxColors = ['black','grey']
        numBoxes = 2
        medians = range(numBoxes)
        for box_index in range(numBoxes):
            box = bp['boxes'][box_index]
            boxX = []
            boxY = []

            for j in range(5):
                boxX.append(box.get_xdata()[j])
                boxY.append(box.get_ydata()[j])
            boxCoords = zip(boxX,boxY)
            # Alternate between Dark Khaki and Royal Blue
            k = box_index % 2
            boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
            ax1.add_patch(boxPolygon)
            # Now draw the median lines back over what we just filled in
            med = bp['medians'][box_index]
            medianX = []
            medianY = []

            for j in range(2):
                medianX.append(med.get_xdata()[j])
                medianY.append(med.get_ydata()[j])
            ax1.plot(medianX, medianY, 'r', linewidth=2)
            medians[box_index] = medianY[0]


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Waveform-match boxes
        # --- don't use the results of injected model
        xlocs2+=5
        boxdat=match_tmp.T
        bp2=ax2.boxplot(boxdat,
                positions=xlocs2, sym='', widths=0.75)#,
                #notch=False, sym='', bootstrap=5000)
        
        pl.setp(bp2['boxes'], color='black')
        pl.setp(bp2['whiskers'], color='black')
        pl.setp(bp2['medians'], color='red')
        
        xticklocs_boxplot2[w]=np.mean(xlocs2)

        # Color the box faces
        # -------------------
        # Lifted from: http://matplotlib.org/examples/pylab_examples/boxplot_demo2.html
        boxColors = ['white','grey','black']
        numBoxes = 3
        medians = range(numBoxes)
        for box_index in range(numBoxes):
            box = bp2['boxes'][box_index]
            boxX = []
            boxY = []

            for j in range(5):
                boxX.append(box.get_xdata()[j])
                boxY.append(box.get_ydata()[j])
            boxCoords = zip(boxX,boxY)
            # Alternate between Dark Khaki and Royal Blue
            k = box_index #% 2
            boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
            ax2.add_patch(boxPolygon)
            # Now draw the median lines back over what we just filled in
            med = bp2['medians'][box_index]
            medianX = []
            medianY = []

            for j in range(2):
                medianX.append(med.get_xdata()[j])
                medianY.append(med.get_ydata()[j])
            ax2.plot(medianX, medianY, 'r', linewidth=2)
            medians[box_index] = medianY[0]

    # --- Label Box Plots
    for c,cat in enumerate(catalogues): 
        if cat==inj: rmidx=c
    these_cats=np.delete(catalogues,np.s_[rmidx])
    # Dummy feature for legend
    cols=['black','grey']
    for c in range(len(these_cats)):
        ax1.bar(-100+i,0,color=cols[c],label='X=%s'%these_cats[c])

    # ---------------
    # - logB plots - 
    # ---------------

    ax1.minorticks_on()
    ax1.set_xlim(0, xlocs[-1]+1)
    ax1.set_xticks(xticklocs_boxplot)
    for x in xticklocs_boxplot: 
        #ax1.axvline(x,color='black',linestyle='--')
        x1=x-1.0
        x2=x+1.0
        ylow=ax1.get_ylim()[0]
        yupp=ax1.get_ylim()[1]
        ax1.fill_between([x1,x2],y1=ylow,y2=yupp,color='grey',alpha=0.2)
    #sys.exit()

    # Grab the labels to use from the dictionary
    xlabels=wf_labels[inj][:len(xticklocs_boxplot)]
    ax1.set_xticklabels(xlabels,rotation=0,weight='bold')

    ax1.yaxis.grid(True, linestyle='-', which='major', color='grey')
    ax1.yaxis.grid(True, linestyle=':', which='minor', color='grey')
    ax1.axhline(0,color='m')
    ax1.set_axisbelow(True)
    ax1.set_ylabel('log$_{\\mathbf{10}}$ [prob(X) / prob(%s)]'%inj, weight='bold')
    ax1.set_xlabel('%s Injection Parameters'%inj, weight='bold')
    
    ax1.set_xlim(min(xticklocs_boxplot)-1,max(xticklocs_boxplot)+1)
    ax1.legend(loc='lower left')

    ax1.set_title('Injected: %s, Recovered: X (see legend)'%inj)
    #pl.subplots_adjust(bottom=0.2)

    # ---------------
    # - match plots - 
    # ---------------
    mcols=['w','grey','k']
    for c in range(len(catalogues)):
        ax2.bar(-100+i,0,color=mcols[c],label='X=%s'%catalogues[c])
    ax2.minorticks_on()
    ax2.set_xlim(0, xlocs2[-1]+1)
    ax2.set_xticks(xticklocs_boxplot2)
    ax2.set_ylim(-1,1)
    for x in xticklocs_boxplot2: 
        #ax1.axvline(x,color='black',linestyle='--')
        x1=x-1.4
        x2=x+1.4
        ylow=ax2.get_ylim()[0]
        yupp=ax2.get_ylim()[1]
        ax2.fill_between([x1,x2],y1=ylow,y2=yupp,color='grey',alpha=0.2)

    # Grab the labels to use from the dictionary
    xlabels=wf_labels[inj][:len(xticklocs_boxplot2)]
    ax2.set_xticklabels(xlabels,rotation=0,weight='bold')
    ax2.set_xlabel('%s Injection Parameters'%inj, weight='bold')

    ax2.yaxis.grid(True, linestyle='-', which='major', color='grey')
    ax2.yaxis.grid(True, linestyle=':', which='minor', color='grey')
    ax2.set_axisbelow(True)
    ax2.set_ylabel('match', weight='bold')
    
    ax2.set_xlim(min(xticklocs_boxplot2)-2,max(xticklocs_boxplot2)+2)
    ax2.legend(loc='lower right')


    ax2.set_title('Injected: %s, Recovered: X (see legend)'%inj)

    # --- Finalise and save
    f1.tight_layout()
    f2.tight_layout()

    f1.savefig('Inj-%s-logB_sigvsig-%s.png'%(inj,resultsfile.replace('.npz','')))
    f2.savefig('Inj-%s-match-%s.png'%(inj,resultsfile.replace('.npz','')))


pl.show()



