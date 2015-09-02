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
bhextractor_pca.py

The main driver module for constructing NR BBH waveform catalogues from
NINJA-style data and for performing PCA of that catalogue.
"""

from __future__ import division

import os
import sys 
import cPickle as pickle

import numpy as np
from scipy.spatial.distance import euclidean as euclidean_distance

from sklearn.decomposition import PCA 

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

# *****************************************************************************
# Function Definitions 


def highpass(timeseries, delta_t=1./512, knee=9., order=12, attn=0.1):

    tmp = pycbc.types.TimeSeries(initial_array=timeseries, delta_t=delta_t)

    return np.array(pycbc.filter.highpass(tmp, frequency=knee, filter_order=order,
            attenuation=attn).data)



def reconstruct_waveform(pca, betas, npcs, mtotal_ref=250.0,
        mtotal_target=250.0, fs=512, fd_interpolation=False):
    """
    Reconstruct the specified waveform (waveform_name) from the PC basis
    contatined in pca and the projection coefficients betas

    pca is the object returned by sklearn.decomposition.PCA()
    """

    #betas=projection[waveform_name]

    reconstruction = np.zeros(shape=(np.shape(pca.components_[0,:])))

    for b in xrange(npcs):
        reconstruction += betas[b]*pca.components_[b,:]
    reconstruction += pca.mean_

    #
    # Mass Scaling
    #
    if mtotal_ref==mtotal_target:
        return reconstruction
    else:
        amp_ratio = mtotal_target / mtotal_ref

        #
        # TD scaling
        #

        # Procedure: increasing the mass shortens the waveform; so interpolate
        # the original to time stamps which are mtotal_target / mtotal_ref x
        # closer together.  We do this by creating a new time axis which is
        # denser and shorter and then interpolate the *original* waveform FROM
        # those new timestamps TO the *original* time stamps.
        #from matplotlib import pyplot as pl

        delta_t = 1./fs
        times = np.arange(0,len(reconstruction)/float(fs), delta_t)

        # XXX 
        peakidx = 0.5*len(times)#np.argmax(abs(reconstruction))

        interp_times = amp_ratio*times - peakidx*delta_t*(amp_ratio-1.0)

        reconstruction = np.interp(times, interp_times, reconstruction)

    return reconstruction

def compute_projection(waveform, pca_result):
    """
    Compute the projection coefficients (betas) for the given waveform onto the
    PCA result
    """

    return np.concatenate(pca_result.transform(waveform))

def perform_pca(data_matrix):
    """
    Use the scikits-learn PCA tools to do PCA on a peak-aligned time-domain
    waveform catalogue
    """

    print 'performing PCA'

    #pca = PCA(whiten=True)
    pca = PCA(whiten=False)
    pca.fit(data_matrix)

    return pca


# *****************************************************************************
# Class Definitions 

class waveform_pca:
    """
    The full PCA result, including projection coefficients etc. Takes a
    waveform catalogue object as input.
    """

    def __init__(self, training_catalogue, testing_catalogue):

        self.training_catalogue = training_catalogue
        self.testing_catalogue  = testing_catalogue

        #
        # PCA of the NR data
        #
        print "Performing PCA of NR waveforms"
        self.NRwave_PCA = perform_pca(training_catalogue.NRdata)

        # Project the testing catalogue onto the new basis formed from the
        # training catalogue
        self.project_catalogue(testing_catalogue.NRdata)

        #
        # PCA of the SI data (if any)
        #
        if hasattr(training_catalogue, 'SIdata'):
            print "Performing PCA of SI waveforms"
            self.SIwave_PCA = perform_pca(training_catalogue.SIdata)

        if hasattr(training_catalogue, 'SIdata') and hasattr(testing_catalogue, 'SIdata'):
            self.project_catalogue(testing_catalogue.SIdata)
        else:
            print "No SI projections computed; one of training/testing catalogues
does not have SIdata"


    def add_pca(self, training_data):
        """
        Append PCA objects for the plus, cross, amplitude and phase time series

        """
        # TODO: 1) return pca objects for plus, cross, amplitude and phase

        hplus = np.real(training_data)
        hcross = -1*np.imag(training_data)
        amplitude = abs(training_data)
        phase = np.zeros(shape=np.shape(training_data))
        for idx in xrange(np.shape(training_data)[0]):
            phase[idx,:] = phase_of(training_data[idx,:])



        return pca


    def file_dump(self, identifier=None):
        """
        Dump the all the results to a python-friendly pickle and stick the PCs
        in ascii files; first row of the ascii file is the mean.

        (Pickling is pretty useful here, since I want to retain the object
        itself)
        """

        if identifier is None:
            identifier = "PCA"

        # Output path
        try:
            pcs_path=os.environ['BHEX_PREFIX']+"/data/PCA_data"
        except KeyError:
            print >> sys.stderr, "BHEX_PREFIX environment variable appears to be un-set"
            sys.exit()

        # Make output directory
        if not os.path.exists(pcs_path):
            os.makedirs(pcs_path)


        # Build filename
        filename = pcs_path + '/' +\
                self.training_catalogue.catalogue_name + "-" + \
                str(identifier)

        print "Performing data dump to %s.pickle"%filename
        
        # Pickle me!
        pickle.dump(self, open(filename+".pickle",'wb'))

        #
        # Ascii Dump
        #
        fplus  = filename+"_hplusPCs.asc"
        fcross = filename+"_hcrossPCs.asc"
        fcomplex_strain = filename+"_hPCs.dat"
        famp = filename+"_AmpPCs.dat"
        fphase = filename+"_PhasePCs.dat"

        # First row contains the mean waveform
        dims = np.shape(self.pca_plus.components_)
        output_array_plus  = np.zeros(shape=(dims[0]+1,dims[1]))
        output_array_cross = np.zeros(shape=(dims[0]+1,dims[1]))

        output_array = np.zeros(shape=(dims[0]+1,dims[1]), dtype=complex)
        output_array_amp = np.zeros(shape=(dims[0]+1,dims[1]))
        output_array_phase = np.zeros(shape=(dims[0]+1,dims[1]))

        output_array_plus[0,:] = self.pca_plus.mean_
        output_array_plus[1:,:] = self.pca_plus.components_
        output_array_cross[0,:] = self.pca_cross.mean_
        output_array_cross[1:,:] = self.pca_cross.components_

        output_array[0,:] = self.pca_plus.mean_ - 1j*self.pca_cross.mean_
        output_array[1:,:] = self.pca_plus.components_ - 1j*self.pca_cross.components_

        output_array_amp[0,:] = self.pca_amp.mean_
        output_array_amp[1:,:] = self.pca_amp.components_

        output_array_phase[0,:] = self.pca_phase.mean_
        output_array_phase[1:,:] = self.pca_phase.components_

        print "Performing data dump to %s"%fplus
        np.savetxt(fplus, output_array_plus)
        print "Performing data dump to %s"%fcross
        np.savetxt(fcross, output_array_cross)

        print "Performing data dump to %s"%fcomplex_strain
        fpcomplex_strain = open(fcomplex_strain,"wb")
        output_array.tofile(fpcomplex_strain)

        print "Performing data dump to %s"%famp
        fpamp = open(famp,"wb")
        output_array_amp.tofile(fpamp)

        print "Performing data dump to %s"%fphase
        fpphase = open(fphase,"wb")
        output_array_phase.tofile(fpphase)


    def project_catalogue(self):
        """
        Project catalogue waveforms onto the new basis to find the
        coefficients for each waveform
        """

        print 'Projecting time-domain testing catalogue onto new basis'

        # Store the projection coefficients in a dictionary, keyed by the
        # waveform names
        self.projection_plus  = dict()
        self.projection_cross = dict()
        self.projection_amp   = dict()
        self.projection_phase = dict()

        for w in xrange(np.shape(self.testing_catalogue.aligned_catalogue)[0]):

            self.projection_plus[self.testing_catalogue.waveform_names[w]] = \
                    compute_projection(np.real(self.testing_catalogue.aligned_catalogue[w,:]),
                            self.pca_plus)

            self.projection_cross[self.testing_catalogue.waveform_names[w]] = \
                    compute_projection(-1*np.imag(self.testing_catalogue.aligned_catalogue[w,:]),
                            self.pca_cross)

            self.projection_amp[self.testing_catalogue.waveform_names[w]] = \
                    compute_projection(self.testing_catalogue.aligned_amplitudes[w,:],
                            self.pca_amp)

            self.projection_phase[self.testing_catalogue.waveform_names[w]] = \
                    compute_projection(self.testing_catalogue.aligned_phases[w,:],
                            self.pca_phase)

    def compute_matches(self, mtotal_ref=250.0):
        """
        Compute the match as a function of nPCs, for each waveform in the
        catalogue.  Uses the aLIGO ZDHP noise curve.

        Since hx is just a phase shifted copy of h+, we only use h+ here.
        """

        if not hasattr(self, 'projection'):
            self.project_catalogue()

        self.matches_plus = np.zeros(shape=(self.testing_catalogue.nwaves,
            self.training_catalogue.nwaves))

        self.euclidean_distances_plus = np.zeros(shape=(self.testing_catalogue.nwaves,
            self.training_catalogue.nwaves)) 

        self.matches_ampphase = np.zeros(shape=(self.testing_catalogue.nwaves,
            self.training_catalogue.nwaves))

        self.euclidean_distances_amp = np.zeros(shape=(self.testing_catalogue.nwaves,
            self.training_catalogue.nwaves)) 

        self.euclidean_distances_phase = np.zeros(shape=(self.testing_catalogue.nwaves,
            self.training_catalogue.nwaves)) 

        # build psd for match calculation
        example_td = pycbc.types.TimeSeries(
                np.real(self.training_catalogue.aligned_catalogue[0,:]),
                delta_t=1.0/self.training_catalogue.fs)

        flen = len(example_td.to_frequencyseries())

        psd = aLIGOZeroDetHighPower(flen,
                example_td.to_frequencyseries().delta_f,
                low_freq_cutoff=10.0)

        print 'Computing match as a function of nPC for each waveform'
        for w in xrange(self.testing_catalogue.nwaves):

            # retrieve projection coefficients
            hplus_betas  = self.projection_plus[self.testing_catalogue.waveform_names[w]]
            amp_betas    = self.projection_amp[self.testing_catalogue.waveform_names[w]]
            phase_betas  = self.projection_phase[self.testing_catalogue.waveform_names[w]]
            
            target_hplus = pycbc.types.TimeSeries(
                    np.real(self.testing_catalogue.aligned_catalogue[w,:]),
                    delta_t=1./self.testing_catalogue.fs)

            target_amp   = self.testing_catalogue.aligned_amplitudes[w,:]
            target_phase = self.testing_catalogue.aligned_phases[w,:]

            for n in xrange(self.training_catalogue.nwaves):

                #
                # Plus reconstruction
                #
                reconstructed_hplus = np.real(reconstruct_waveform(self.pca_plus,
                        hplus_betas, n+1, mtotal_ref=mtotal_ref))
                recwav = pycbc.types.TimeSeries(reconstructed_hplus,
                        delta_t=1./self.training_catalogue.fs)

                self.matches_plus[w,n] = pycbc.filter.match(recwav,
                        target_hplus, psd=psd, low_frequency_cutoff=10.0)[0]

                self.euclidean_distances_plus[w,n] = euclidean_distance(
                        reconstructed_hplus, target_hplus.data)
                #
                # Amp/phase reconstruction
                #
                reconstructed_amp = reconstruct_waveform(self.pca_amp,
                        amp_betas, n+1, mtotal_ref=mtotal_ref)
                reconstructed_phase = reconstruct_waveform(self.pca_phase,
                        phase_betas, n+1, mtotal_ref=mtotal_ref)

                reconstructed_h = reconstructed_amp * \
                        np.exp(1j*reconstructed_phase)

                recwav = pycbc.types.TimeSeries(np.real(reconstructed_h),
                        delta_t=1./self.training_catalogue.fs)

                self.matches_ampphase[w,n] = pycbc.filter.match(recwav,
                        target_hplus, psd=psd, low_frequency_cutoff=10.0)[0]

                self.euclidean_distances_amp[w,n] =\
                        euclidean_distance(reconstructed_amp, target_amp)

                self.euclidean_distances_phase[w,n] =\
                        euclidean_distance(reconstructed_amp, target_amp)

                     


# *******************************************************************************
def main():

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # INPUT 

    catalogue_name=sys.argv[1]

    # Quick 'n' dirty sanity check (should use option parser really)
    if catalogue_name not in ['Q','HR','RO3']:
        print >> sys.stderr, \
                "Error, catalogue '%s' not recognised "\
                "(must be Q, HR or RO3)"%catalogue_name

    #
    # Setup and then build the catalogue
    #
    catalogue = waveform_catalogue(catalogue_name=catalogue_name, fs=512,
            catalogue_len=8, mtotal_ref=250, Dist=1.)

    #
    # Do the PCA
    #
    pca = waveform_pca(catalogue)
    #pca.file_dump()

    return pca


#
# End definitions
#
if __name__ == "__main__":
    pca_results = main()
    pca_results.file_dump()

    



