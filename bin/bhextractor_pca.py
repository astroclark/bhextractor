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
import sys 

import numpy as np
import scipy.signal as signal
import scipy.io as sio
import scipy.linalg

from sklearn.decomposition import PCA 
from sklearn.decomposition import TruncatedSVD

import lal
import lalsimulation as lalsim

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

# *****************************************************************************
# Function Definitions 

def window_wave(input_data):

    nonzero=np.argwhere(abs(input_data)>0)
    idx = range(nonzero[0],nonzero[-1])
    beta = 0.25
    win = lal.CreateTukeyREAL8Window(len(idx), beta)
    win.data.data[int(0.5*len(idx)):]=1.0

    input_data[idx] *= win.data.data

    return input_data

def highpass(timeseries, delta_t=1./2048, knee=10., order=8, attn=0.1):

    tmp = pycbc.types.TimeSeries(initial_array=timeseries, delta_t=delta_t)

    return np.array(pycbc.filter.highpass(tmp, frequency=knee, filter_order=order,
            attenuation=attn).data)

def perform_pca(aligned_td_catalogue):
    """
    Use the scikits-learn PCA tools to do PCA on a peak-aligned time-domain
    waveform catalogue
    """

    pca = PCA()
    pca.fit(aligned_td_catalogue)

    return pca

def reconstruct_waveform(pca, betas, npcs, waveform_name, mtotal=250.0):
    """
    Reconstruct the specified waveform (waveform_name) from the PC basis
    contatined in pca and the projection coefficients betas

    pca is the object returned by sklearn.decomposition.PCA()
    """

    #betas=projection[waveform_name]

    reconstruction = np.zeros(len(pca.components_[0,:]), dtype=complex)
    for b in xrange(npcs):
        reconstruction += pca.components_[b,:] * betas[b]
    reconstruction += pca.mean_

    # TODO: mass scaling

    return reconstruction

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


# *****************************************************************************
# Class Definitions 

class waveform_pca:
    """
    The full PCA result, including projection coefficients etc. Takes an aligned
    waveform catalogue object as input.
    """

    def __init__(self, catalogue):

        self.catalogue = catalogue

        # Do the PCA
        self.pca = perform_pca(self.catalogue.aligned_catalogue)

        # Project the catalogue onto the new basis
        self.project_catalogue()

    def project_catalogue(self):
        """
        Project catalogue waveforms onto the new basis to find the
        coefficients for each waveform
        """

        print 'Projecting catalogue onto new basis'

        # Store the projection coefficients in a dictionary, keyed by the
        # waveform names
        self.projection = dict()

        for w in xrange(np.shape(self.catalogue.aligned_catalogue)[0]):

            # Center data
            waveform_centered = self.catalogue.aligned_catalogue[w,:] - \
                    self.pca.mean_

            self.projection[self.catalogue.waveform_names[w]] = \
                    np.concatenate(self.pca.transform(waveform_centered))

    def compute_matches(self, mtotal=250.0):
        """
        Compute the match as a function of nPCs, for each waveform in the
        catalogue.  Uses the aLIGO ZDHP noise curve.
        """

        if not hasattr(self, 'projection'):
            self.project_catalogue()

        self.matches = np.zeros(shape=(len(self.catalogue.waveform_names),
            len(self.catalogue.waveform_names)))

        # build psd for match calculation
        example_td = pycbc.types.TimeSeries(
                np.real(self.catalogue.aligned_catalogue[0,:]),
                delta_t=1.0/self.catalogue.fs)

        flen = len(example_td.to_frequencyseries())

        psd = aLIGOZeroDetHighPower(flen,
                example_td.to_frequencyseries().delta_f,
                low_freq_cutoff=10.0)#self.low_frequency_cutoff)

        print 'Computing match as a function of nPC for each waveform'
        for w in xrange(len(self.catalogue.waveform_names)):

            # retrieve projection coefficients
            betas = self.projection[self.catalogue.waveform_names[w]]

            for n in xrange(len(self.catalogue.waveform_names)):

                reconstructed_waveform = reconstruct_waveform(self.pca, betas,
                        n, self.catalogue.waveform_names[n], mtotal=mtotal)

                recwav = pycbc.types.TimeSeries(np.real(reconstructed_waveform),
                        delta_t=1./self.catalogue.fs)

                targetwav = pycbc.types.TimeSeries(
                        np.real(self.catalogue.aligned_catalogue[w,:]),
                        delta_t=1./self.catalogue.fs)

                self.matches[w,n] = pycbc.filter.match(recwav, targetwav,
                        psd=psd, low_frequency_cutoff=10.0)[0]


        # TODO:
        # Introduce mass scaling
        # Produce 'standard' explained variance / match plots
        # Save PCs

                     
class waveform_catalogue:
    """
    Object with the full waveform catalogue and all conditioning information
    """

    def __init__(self,catalogue_name='Q', theta=0.0, phi=0.0, fs=2048,
            catalogue_len=4, mtotal_ref=250, Dist=1.):

        # ######################################################
        # Other initialisation
        self.catalogue_name = catalogue_name
        self.theta = theta
        self.phi = phi
        self.fs = fs
        self.catalogue_len = catalogue_len * fs # XXX
        self.mtotal_ref = mtotal_ref
        self.Dist = Dist

        #
        # Build the catalogue
        #
        self.build_catalogue()

    def build_catalogue(self):
        """
        Return a complex array with the time-domain waveforms for the requested
        catalogue
        """

        # ######################################################
        # Some catalogue-specific stuff

        # Dictionary with the maximum waveform data lengths
        maxlengths={'Q':10959, 'HR':31238, 'RO3':16493}

        # Dictionary of NR sample rates. XXX Should probably get this from the data,
        # really
        NR_deltaT={'Q':0.155, 'HR':0.08, 'RO3':2./15}

        # Samples to discard due to NR noise
        NRD_sampls=500 # keep this many samples after the peak
        NINSP_sampls=2000 # discard this many samples from the start

        # --- Identify waveforms in catalogue
        catalogue_path=os.environ['BHEX_PREFIX']+'/data/NR_data/'+self.catalogue_name+'-series'
        waveforms = [ f for f in os.listdir(catalogue_path) if
                os.path.isdir(os.path.join(catalogue_path,f)) ]

        # ---------- MASS-DEPENDENT CALCULATIONS -----------
        # Determine current sample rate.  XXX: GIVEN A MASS 
        SI_deltaT = self.mtotal_ref * lal.MTSUN_SI * NR_deltaT[self.catalogue_name]
        Mscale = self.mtotal_ref * lal.MRSUN_SI / (self.Dist * 1e9 * lal.PC_SI)

        # Get the length of the NR data
        wflen=maxlengths[self.catalogue_name]

        # Length to resample to: length of (e.g.) 250 Msun waveform in seconds x desired sample rate
        resamp_len = np.ceil(wflen*SI_deltaT*self.fs)
        catalogue=np.zeros(shape=(len(waveforms), resamp_len), dtype=complex)

        self.waveform_names = []

        # Build waveforms in this catalogue
        for w, waveform in enumerate(waveforms):
            print 'Building %s waveform'%waveform

            self.waveform_names.append(waveform)

            # Get files (modes) in this waveform dir
            waveform_dir = catalogue_path+'/'+waveform
            modes = [ f for f in os.listdir(waveform_dir) if
                    os.path.isfile(os.path.join(waveform_dir,f)) ]

            # --- Sum modes
            # See line 95 of NRWaveInject.c
            hplus=np.zeros(wflen)
            hcross=np.zeros(wflen)
            for mode in modes:

                # Load mode data
                mode_data = np.loadtxt(waveform_dir+'/'+mode)

                # Extract l,m from filename
                spharm_degree = int(mode.split('_')[2].replace('l',''))
                spharm_order  = int(mode.split('_')[3].replace('m',''))

                # Compute spin -2 weighted spherical harmonic
                sYlm = lal.SpinWeightedSphericalHarmonic(self.theta, self.phi, 
                        -2, spharm_degree, spharm_order)

                # Orient Waveforms
                hplus[0:len(mode_data)]  += mode_data[:,1]*np.real(sYlm) + \
                        mode_data[:,2]*np.imag(sYlm)
                hcross[0:len(mode_data)] += mode_data[:,2]*np.real(sYlm) - \
                        mode_data[:,1]*np.imag(sYlm)

            #
            # Clean up the waveform
            #

            # Get rid of junk
            peak_idx = np.argmax(abs(hplus))
            zero_idx = min(peak_idx + NRD_sampls, len(hplus))
            hplus[zero_idx:]  = 0.0
            hcross[zero_idx:] = 0.0

            hplus[:peak_idx-NINSP_sampls] = 0.0
            hcross[:peak_idx-NINSP_sampls] = 0.0

            # --- Windowing / tapering
            hplus=window_wave(hplus)
            hcross=window_wave(hcross)

            # --- Resampling 
            hplus  = signal.resample(hplus, resamp_len)
            hcross = signal.resample(hcross, resamp_len)

            # --- Filtering
            hplus  = highpass(hplus)
            hcross = highpass(hcross)

            # Use complex waveform
            catalogue[w,:] = hplus - 1j*hcross
            
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Catalogue Conditioning

        # Standardised catalogue (4 seconds long)
        self.aligned_catalogue = \
                np.zeros(shape=(len(waveforms),self.catalogue_len),
                        dtype=complex)

        print 'aligning peak times'

        # Find peak indices
        peak_idx=np.argmax(abs(catalogue),axis=0)

        # Align all waveform peaks to the 0.5 of the way through the final catalogue
        align_idx=np.floor(0.5*self.catalogue_len)

        for w in xrange(len(waveforms)):
            print 'aligning %d of %d'%(w, len(waveforms))

            non_zero_idx = \
                    np.argwhere(abs(catalogue[w,:])>1e-2*max(abs(catalogue[w,:])))
            trunc_wav = \
                    catalogue[w,non_zero_idx[0]:non_zero_idx[-1]]

            peak_idx=np.argmax(abs(trunc_wav))
            start_idx = align_idx - peak_idx

            self.aligned_catalogue[w,start_idx:start_idx+len(trunc_wav)] = trunc_wav

            # --- Normalisation
            self.aligned_catalogue[w,:] /= \
                    np.linalg.norm(self.aligned_catalogue[w,:])

# *******************************************************************************
def main():

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # INPUT 

    catalogue_name=sys.argv[1]
    theta=float(sys.argv[2])
    phi=0.0

    # Quick 'n' dirty sanity check (should use option parser really)
    if catalogue_name not in ['Q','HR','RO3']:
        print >> sys.stderr, \
                "Error, catalogue '%s' not recognised "\
                "(must be Q, HR or RO3)"%catalogue_name

    #
    # Setup and then build the catalogue
    #
    catalogue = waveform_catalogue(catalogue_name=catalogue_name, fs=2048,
            catalogue_len=4, mtotal_ref=250, Dist=1., theta=theta)

    #
    # Do the PCA
    #
    pca = waveform_pca(catalogue)

    return pca


#
# End definitions
#
if __name__ == "__main__":
    pca_results = main()

    



