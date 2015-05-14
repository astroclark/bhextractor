
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

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower


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

oriwave250 = np.copy(catalogue.aligned_catalogue[0,:])

#
# Do the PCA
#
pca = bhex.waveform_pca(catalogue)

#
# Build a 350 solar mass waveform from the 250 Msun PCs#

# Just use the first waveform
betas = pca.projection[catalogue.waveform_names[0]]

times = np.arange(0,len(catalogue.aligned_catalogue[0,:])/2048.,1./2048)

recwave350 = bhex.reconstruct_waveform(pca.pca, betas, len(catalogue.waveform_names),
        mtotal_target=350.0)


#
# Now make a catalogue at 350 solar masses and then compute the overlap
#

catalogue350 = bhex.waveform_catalogue(catalogue_name=catalogue_name, fs=2048,
        catalogue_len=4, mtotal_ref=350, Dist=1., theta=theta)

oriwave350 = np.copy(catalogue350.aligned_catalogue[0,:])


# Finally, compute the match between the reconstructed 350 Msun system and the
# system we generated at that mass in the first place

recwave350_pycbc = pycbc.types.TimeSeries(np.real(recwave350), delta_t=1./2048)
oriwave250_pycbc = pycbc.types.TimeSeries(np.real(oriwave250), delta_t=1./2048)
oriwave350_pycbc = pycbc.types.TimeSeries(np.real(oriwave350), delta_t=1./2048)
psd = aLIGOZeroDetHighPower(len(recwave350_pycbc.to_frequencyseries()),
        recwave350_pycbc.to_frequencyseries().delta_f, low_freq_cutoff=10.0)

print 'Match between 250 and 350 Msun catalogue waves: ',\
        pycbc.filter.match(oriwave250_pycbc.to_frequencyseries(),
                oriwave350_pycbc.to_frequencyseries(), psd=psd,
                low_frequency_cutoff=10)[0]

print 'Match between 350 reconstruction and 350 catalogue wave: ',\
        pycbc.filter.match(recwave350_pycbc.to_frequencyseries(),
                oriwave350_pycbc.to_frequencyseries(), psd=psd,
                low_frequency_cutoff=10)[0]



