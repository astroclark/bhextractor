# -*- coding: utf-8 -*-
#!/usr/bin/env python
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
#import scipy.optimize as optimize
#import scipy.signal as signal

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
#num_wfs={'Q':13, 'HR':10, 'RO3':10}

smeepath=os.getenv('SMEEBBH_PREFIX')

catalogue_files={'Q':'%s/PCA/bhextractor/Q_catalogue_theta-0.mat'%(smeepath),
        'HR':'%s/PCA/bhextractor/HR_catalogue_theta-0.mat'%(smeepath),
        'RO3':'%s/PCA/bhextractor/RO3_catalogue_theta-0.mat'%(smeepath)}

pc_files={'Q':'%s/PCA/bhextractor/Q_PCs_theta-0.mat'%(smeepath),
        'HR':'%s/PCA/bhextractor/HR_PCs_theta-0.mat'%(smeepath),
        'RO3':'%s/PCA/bhextractor/RO3_PCs_theta-0.mat'%(smeepath)}


Mtot=250.
Dist=1.0
Mscale= Mtot * lal.LAL_MRSUN_SI / (Dist * 1e9 * lal.LAL_PC_SI)

if not smeepath:
    print >> sys.stderr, "ERROR, must set SMEEBBH_PREFIX"
    sys.exit()

resultsdir=sys.argv[1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop through results files

# --- Preallocate main results arrays
logB_out   = []
match_out  = []

# --- Load targets and PCs first to save overheads
target_waveforms=[]
princ_comps=[]
for wf in catalogues:
    target_waveforms.append(sio.loadmat(catalogue_files[wf])['MDC_final'])
    princ_comps.append(sio.loadmat(pc_files[wf])['PCs_final'])

# --- Begin loading all data
# Loop through injected waveform
#tmp=['HR']
#for i,inj in enumerate(tmp):
for i,inj in enumerate(catalogues):

    # number of injections performed from this catalogue
    nwf=num_wfs[inj]

    print >> sys.stdout, "Loading {0} Injections".format(inj)

    # Loop through injected waveform catalogue number 
    a=0
    allmat=glob(resultsdir+'/'+'*smee_output_%s*mat'%inj)

    logB_by_inj   = []
    match_by_inj  = []

    for wf in range(nwf):
        
        # pick out the target waveform
        target_wf = target_waveforms[i][:,wf]

        print >> sys.stdout, "Frac Complete: {0:.2f}".format(
                float(a)/len(allmat))

        logB_by_rec   = []
        match_by_rec  = []

        # Loop through waveform used for recovery
        for r,rec in enumerate(catalogues):

            # --- Get the PCs for this catalogue
            npc=num_pcs[rec]
            pcs = princ_comps[r][:,:npc]

            # --- construct file pattern

            # --- Glob for mat files from this injection / recovery combo
            filepat='smee_output_{inj}{wf}_model{rec}_PCs{npc}_*'.format(
                    inj=inj, rec=rec, wf=int(wf+1), npc=npc)
            matching_files=glob(resultsdir+'/'+filepat)
            #for i in range(5): print matching_files[i]


            #print 'Inj: %s-%d, Rec: %s, Nsim: %d'%(
            #        inj, wf, rec, len(matching_files))

            # --- Preallocate arrays for these files
            Zsig=np.zeros(len(matching_files))
            logB=np.zeros(len(matching_files))
            match=np.zeros(len(matching_files))

            # Loop through noise realisations
            for f,mf in enumerate(matching_files):
                #print '%d %d'%(a,len(allmat))
                a+=1

                data=sio.loadmat(mf)

                # populate signal evidence, log B
                Zsig[f]=data['Z_end']
                logB[f]=data['Bayes']

                # --- Compute match between reconstructed and target waveforms
                # Reconstruct the waveform from the max-likelihood betas
                betas=data['betas']
                Lw=data['Lw']
                beta_maxL=betas[np.argmax(Lw),:]
                #mean_betas=np.mean(betas,axis=0)

                # The reconstructed waveform:
                rec_wf = model_wf(pcs, *beta_maxL)
                #rec_wf = model_wf(pcs, *mean_betas)

                # Compute match between unit-norm waveforms
                match[f] = np.dot(rec_wf/np.linalg.norm(rec_wf),
                        target_wf/np.linalg.norm(target_wf))

            logB_by_rec.append(logB)
            match_by_rec.append(match)

        logB_by_inj.append(logB_by_rec)
        match_by_inj.append(match_by_rec)

    logB_out.append(logB_by_inj)
    match_out.append(match_by_inj)
    

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save for further processing

# FIXME: input-check 
outname=sys.argv[2]

np.savez(outname, logB   = logB_out,
                  match  = match_out
                  )
