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

import sys 
import numpy as np
import scipy.io as sio 
from scipy import signal
import matplotlib
#matplotlib.use("Agg")
from matplotlib import pyplot as pl
import glob

import lal
import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

def model_wf(pcs, results_data):
    """
    Construct a template from the sum of PCs weighted by the coefficients in
    betas
    """

    likelihoods=results_data['Lw']
    betas=results_data['betas']
    beta_maxL=betas[np.argmax(likelihoods),:]

    model_wf = np.zeros(np.shape(pcs)[0])
    for n in range(len(beta_maxL)):
        model_wf+=beta_maxL[n] * pcs[:,n]
    return model_wf

def get_injected(injectionfile, injnum):
    injection_data = sio.loadmat(injectionfile)

    target_wf = injection_data['MDC_final'][:,wfnum-1]

    # Put into pycbc TimeSeries for further manipulation
    target_ts = pycbc.types.TimeSeries(target_wf,delta_t=1.0/16384)

    return target_ts

def get_recovered(pcfile, resultsfile):
    results_data = sio.loadmat(resultsfile)
    npcs = results_data['numPCs'][0][0]
    pcs = Mscale*sio.loadmat(pcfile)['PCs_final'][:,0:npcs]

    # recovered:
    rec_wf = model_wf(pcs, results_data)

    # Put into pycbc TimeSeries for further manipulation
    recon_ts  = pycbc.types.TimeSeries(rec_wf,delta_t=1.0/16384)

    return recon_ts, results_data['Bayes']

# data files
Qresultsfile='spring_results/Inj-RO3_Wf-9_Seed-9561/smee_output_RO39_modelQ_PCs2_detno1_SNR50_seed9561.mat'
HRresultsfile='spring_results/Inj-RO3_Wf-9_Seed-9561/smee_output_RO39_modelHR_PCs4_detno1_SNR50_seed9561.mat'
RO3resultsfile='spring_results/Inj-RO3_Wf-9_Seed-9561/smee_output_RO39_modelRO3_PCs5_detno1_SNR50_seed9561.mat'

wfnum=9
Qnpcs=2
HRnpcs=4
RO3npcs=5

QPCfile='spring_results/Q_PCs_theta-0.mat'
HRPCfile='spring_results/HR_PCs_theta-0.mat'
RO3PCfile='spring_results/RO3_PCs_theta-0.mat'

Qinjectionfile='spring_results/Q_catalogue_theta-0.mat'
HRinjectionfile='spring_results/HR_catalogue_theta-0.mat'
RO3injectionfile='spring_results/RO3_catalogue_theta-0.mat'

Mtot=250.
Dist=1.0
Mscale= Mtot * lal.MRSUN_SI / (Dist * 1e9 * lal.PC_SI)


# load data
target = get_injected(RO3injectionfile, wfnum-1)
Qrecon, QBayes = get_recovered(QPCfile, Qresultsfile) 
HRrecon, HRBayes = get_recovered(HRPCfile, HRresultsfile) 
RO3recon, RO3Bayes = get_recovered(RO3PCfile, RO3resultsfile) 

flen = 8192
psd = aLIGOZeroDetHighPower(flen, 1., 0.0)

Qmatch, Qmaxidx =  pycbc.filter.match(target, Qrecon, psd=psd)
HRmatch, HRmaxidx =  pycbc.filter.match(target, HRrecon, psd=psd)
RO3match, RO3maxidx =  pycbc.filter.match(target, RO3recon, psd=psd)

print Qmatch
print HRmatch
print RO3match

# apply time shit for easy plotting XXX: this doesn't look right.
QmaxL_recon=pycbc.types.TimeSeries(initial_array=np.copy(Qrecon.data),
              delta_t=Qrecon.delta_t, epoch=Qrecon.start_time+Qmaxidx/16384.)
HRmaxL_recon=pycbc.types.TimeSeries(initial_array=np.copy(HRrecon.data),
              delta_t=HRrecon.delta_t, epoch=HRrecon.start_time+HRmaxidx/16384.)
RO3maxL_recon=pycbc.types.TimeSeries(initial_array=np.copy(RO3recon.data),
              delta_t=RO3recon.delta_t, epoch=RO3recon.start_time+RO3maxidx/16384.)

target_peak, target_peak_idx = target.max_loc()
target/=np.linalg.norm(target)

Q_peak, Q_peak_idx = Qrecon.max_loc()
#QmaxL_recon *= target_peak/Q_peak
QmaxL_recon /= np.linalg.norm(QmaxL_recon)

HR_peak, HR_peak_idx = HRrecon.max_loc()
#HRmaxL_recon *= target_peak/HR_peak
HRmaxL_recon /= np.linalg.norm(HRmaxL_recon)

RO3_peak, RO3_peak_idx = RO3recon.max_loc()
#RO3maxL_recon *= target_peak/RO3_peak
RO3maxL_recon /= np.linalg.norm(RO3maxL_recon)


f,ax=pl.subplots(nrows=1,ncols=3,figsize=(12,4),sharey=True)
ax[0].plot(target.sample_times, target, label='Target (RO3)', color='k')
#ax[0].plot(Qrecon.sample_times, Qrecon, label='Q-series', color='r', linestyle='--')
ax[0].plot(Qrecon.sample_times, -1*QmaxL_recon,
        label='Q recon.', color='r', linestyle='-')
ax[0].legend()
ax[0].set_title('log B=%.2f\nmax-L match=%.2f'%(QBayes, Qmatch))

ax[1].plot(target.sample_times, target, label='Target (RO3)', color='k')
ax[1].plot(HRrecon.sample_times, HRmaxL_recon,
        label='HR recon.', color='r', linestyle='-')
ax[1].legend()
ax[1].set_title('log B=%.2f\nmax-L match=%.2f'%(HRBayes, HRmatch))

ax[2].plot(target.sample_times, target, label='Target (RO3)', color='k')
ax[2].plot(RO3recon.sample_times, RO3maxL_recon,
        label='RO3 recon.'%RO3Bayes, color='r', linestyle='-')
ax[2].legend()
ax[2].set_title('log B=%.2f\nmax-L match=%.2f'%(RO3Bayes, RO3match))

ax[0].set_xlim(0.5,0.8)
ax[1].set_xlim(0.5,0.8)
ax[2].set_xlim(0.5,0.8)

ax[0].minorticks_on()
ax[1].minorticks_on()
ax[2].minorticks_on()

#   ax[0].set_ylim(-5e-23,6e-23)
#   ax[1].set_ylim(-5e-23,6e-23)
#   ax[2].set_ylim(-5e-23,6e-23)
ax[0].set_ylim(-0.05, 0.06)
ax[1].set_ylim(-0.05, 0.06)
ax[2].set_ylim(-0.05, 0.06)

ax[0].set_ylabel('h(t)')

ax[0].set_xlabel('Time [s]')
ax[1].set_xlabel('Time [s]')
ax[2].set_xlabel('Time [s]')

f.tight_layout()


pl.show()

