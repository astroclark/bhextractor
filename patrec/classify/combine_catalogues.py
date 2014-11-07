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

import os
#from os import os.environ,os.listdir,os.makedirs
#from os.path import os.isdir, os.isfile, join
import sys 

__author__ = "James Clark <james.clark@ligo.org>"

import numpy as np
import scipy.signal as signal
import scipy.io as sio

import lal
import lalsimulation as lalsim

import pycbc
import pycbc.filter


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUT 

# comma-separated list of catalogue files
catalogue_names=sys.argv[1]
theta=float(sys.argv[2])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load waveforms
catalogues = []
ncols=0
for c,cat in enumerate(catalogue_names.split(',')):

    catalogue = sio.loadmat(cat)['MDC_final']
    catalogues.append(catalogue)

    # columns = waveforms (should increment)
    ncols+=np.shape(catalogue)[1]

    # rows = samples (should be constant)
    nrows=np.shape(catalogue)[0]

# Matrix containing ALL initial waveform sets
allWaves = np.concatenate(catalogues,1)
allWaves_aligned = np.zeros(np.shape(allWaves), dtype=complex)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Align Peaks

catalogue_real = np.real(allWaves)
catalogue_imag = np.imag(allWaves)

print 'aligning peak times'

# Find peak indices
peak_indices_real=np.argmax(abs(catalogue_real),axis=0)
#peak_indices_imag=np.argmax(abs(catalogue_imag),axis=0)

# Align all waveform peaks to the latest-time peak
align_to_idx_real=max(peak_indices_real)
#align_to_idx_imag=max(peak_indices_imag)

for w in xrange(np.shape(allWaves)[1]):

    # Temp array to hold aligned waveform
    tmp_real = np.zeros(len(catalogue_real[:,w]))
    tmp_imag = np.zeros(len(catalogue_imag[:,w]))
    wf_real = np.copy(catalogue_real[:,w])
    wf_imag = np.copy(catalogue_imag[:,w])

    # Get the lengths of current waveform data to the left/right of the peak of
    # this waveform
    llen_real = len(wf_real[:peak_indices_real[w]])
    llen_imag = len(wf_imag[:peak_indices_real[w]])
    rlen_real = len(tmp_real[align_to_idx_real:])
    rlen_imag = len(tmp_imag[align_to_idx_real:])

    # populate left side of peak
    tmp_real[align_to_idx_real-llen_real:align_to_idx_real] = wf_real[:peak_indices_real[w]]
    tmp_imag[align_to_idx_real-llen_imag:align_to_idx_real] = wf_imag[:peak_indices_real[w]]

    # populate right side of peak
    tmp_real[align_to_idx_real:] = wf_real[peak_indices_real[w]:peak_indices_real[w]+rlen_real]
    tmp_imag[align_to_idx_real:] = wf_imag[peak_indices_real[w]:peak_indices_real[w]+rlen_real]

    allWaves_aligned[:,w] = tmp_real + 1j*tmp_imag



#catalogue_path=os.environ['BHEX_PREFIX']+'/data/'+'signal_data'
catalogue_path='./'

#if not os.path.exists(catalogue_path): os.makedirs(catalogue_path)

filename='allWaves' + 'theta-%.0f'%theta

datadict={'complex_strain':allWaves_aligned}
sio.savemat(filename, datadict)

sys.exit()

