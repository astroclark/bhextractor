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

import numpy as np
import scipy.io as sio
import scipy.optimize as optimize
import scipy.signal as signal

import lal

def model_wf(pcs, *parms):
    """
    Construct a template from the sum of PCs weighted by the coefficients in
    betas
    """
    betas=parms[:-1]
    phi=parms[-1]
    model_wf = np.zeros(np.shape(pcs)[0])
    for n in range(len(betas)):
        model_wf+=betas[n] * pcs[:,n] 
    model_wf = signal.filtfilt(b, a, model_wf)
    #model_wf*=np.exp(-1j*phi)
    return model_wf 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 
PC_dict = sio.loadmat(sys.argv[1])
waveform_dict = sio.loadmat(sys.argv[2])
Mtot=250.0
Dist=1 # Gpc

catalogue=sys.argv[1].split('_')[0]

pcs = PC_dict['PCs_final']
waveforms = waveform_dict['MDC_final']

# Dictionary of npcs to use for each catalogue
# TODO Integrate this with bhextractor_match.py
npcs={'Q':2,'HR':4,'RO3':5}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find max likelihood betas

# --- High pass
b, a = signal.butter(N=4, Wn=10./8192., btype='highpass')

overlap=np.zeros(np.shape(waveforms)[1])
for w in range(np.shape(waveforms)[1]):

    # --- unit norm target
    target_wf=waveforms[:,w] / np.linalg.norm(waveforms[:,w])

    target_wf = signal.filtfilt(b, a, target_wf)

    # --- find best fitting reconstruction parameters
    p0 = np.ones(npcs[catalogue]+1)
    p0[-1] = np.random.randn()
    popt, pcov = optimize.curve_fit(model_wf, xdata=pcs, ydata=target_wf, p0=p0)

    # --- best fitting waveform
    rec_wf = model_wf(pcs,*popt)

    # --- overlap
    overlap[w] = np.dot(rec_wf,target_wf)


