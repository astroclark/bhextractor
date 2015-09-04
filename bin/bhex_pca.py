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

import bhex_wavedata as bwave

# *****************************************************************************
# Function Definitions 


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

def perform_pca(data_matrix):
    """
    Use the scikits-learn PCA tools to do PCA on a peak-aligned time-domain
    waveform catalogue
    """

    #pca = PCA(whiten=True)
    pca = PCA(whiten=False)
    pca.fit(data_matrix)

    return pca




#   def compute_matches(self, mtotal_ref=250.0):
#       """
#       Compute the match as a function of nPCs, for each waveform in the
#       catalogue.  Uses the aLIGO ZDHP noise curve.
#
#       Since hx is just a phase shifted copy of h+, we only use h+ here.
#       """
#
#       if not hasattr(self, 'projection'):
#           self.project_catalogue()
#
#       self.matches_plus = np.zeros(shape=(self.test_catalogue.nwaves,
#           self.training_catalogue.nwaves))
#
#       self.euclidean_distances_plus = np.zeros(shape=(self.test_catalogue.nwaves,
#           self.training_catalogue.nwaves)) 
#
#       self.matches_ampphase = np.zeros(shape=(self.test_catalogue.nwaves,
#           self.training_catalogue.nwaves))
#
#       self.euclidean_distances_amp = np.zeros(shape=(self.test_catalogue.nwaves,
#           self.training_catalogue.nwaves)) 
#
#       self.euclidean_distances_phase = np.zeros(shape=(self.test_catalogue.nwaves,
#           self.training_catalogue.nwaves)) 
#
#       # build psd for match calculation
#       example_td = pycbc.types.TimeSeries(
#               np.real(self.training_catalogue.aligned_catalogue[0,:]),
#               delta_t=1.0/self.training_catalogue.fs)
#
#       flen = len(example_td.to_frequencyseries())
#
#       psd = aLIGOZeroDetHighPower(flen,
#               example_td.to_frequencyseries().delta_f,
#               low_freq_cutoff=10.0)
#
#       print 'Computing match as a function of nPC for each waveform'
#       for w in xrange(self.test_catalogue.nwaves):
#
#           # retrieve projection coefficients
#           hplus_betas  = self.projection_plus[self.test_catalogue.waveform_names[w]]
#           amp_betas    = self.projection_amp[self.test_catalogue.waveform_names[w]]
#           phase_betas  = self.projection_phase[self.test_catalogue.waveform_names[w]]
#           
#           target_hplus = pycbc.types.TimeSeries(
#                   np.real(self.test_catalogue.aligned_catalogue[w,:]),
#                   delta_t=1./self.test_catalogue.fs)
#
#           target_amp   = self.test_catalogue.aligned_amplitudes[w,:]
#           target_phase = self.test_catalogue.aligned_phases[w,:]
#
#           for n in xrange(self.training_catalogue.nwaves):
#
#               #
#               # Plus reconstruction
#               #
#               reconstructed_hplus = np.real(reconstruct_waveform(self.pca_plus,
#                       hplus_betas, n+1, mtotal_ref=mtotal_ref))
#               recwav = pycbc.types.TimeSeries(reconstructed_hplus,
#                       delta_t=1./self.training_catalogue.fs)
#
#               self.matches_plus[w,n] = pycbc.filter.match(recwav,
#                       target_hplus, psd=psd, low_frequency_cutoff=10.0)[0]
#
#               self.euclidean_distances_plus[w,n] = euclidean_distance(
#                       reconstructed_hplus, target_hplus.data)
#               #
#               # Amp/phase reconstruction
#               #
#               reconstructed_amp = reconstruct_waveform(self.pca_amp,
#                       amp_betas, n+1, mtotal_ref=mtotal_ref)
#               reconstructed_phase = reconstruct_waveform(self.pca_phase,
#                       phase_betas, n+1, mtotal_ref=mtotal_ref)
#
#               reconstructed_h = reconstructed_amp * \
#                       np.exp(1j*reconstructed_phase)
#
#               recwav = pycbc.types.TimeSeries(np.real(reconstructed_h),
#                       delta_t=1./self.training_catalogue.fs)
#
#               self.matches_ampphase[w,n] = pycbc.filter.match(recwav,
#                       target_hplus, psd=psd, low_frequency_cutoff=10.0)[0]
#
#               self.euclidean_distances_amp[w,n] =\
#                       euclidean_distance(reconstructed_amp, target_amp)
#
#               self.euclidean_distances_phase[w,n] =\
#                       euclidean_distance(reconstructed_amp, target_amp)

# *****************************************************************************
# Class Definitions 

class waveform_pca:
    """
    The full PCA result, including projection coefficients etc. Takes a
    waveform catalogue object as input.
    """

    def __init__(self, training_catalogue, test_catalogue=None):

        # XXX: don't need to store these, right?
#       self.training_catalogue = training_catalogue
#       if test_catalogue is not None:
#           self.test_catalogue  = test_catalogue

        #
        # PCA of the NR data
        #
        print "PCA of NR waveforms"

        # XXX: Probably the smart place to truncate the catalogue so we have
        # equal length waveforms

        self.NRhplusTimeSeriesPCA = \
                perform_pca(np.real(training_catalogue.NRComplexTimeSeries))
        self.NRhcrossTimeSeriesPCA = \
                perform_pca(-1*np.imag(training_catalogue.NRComplexTimeSeries))
        self.NRAmpTimeSeriesPCA = \
                perform_pca(training_catalogue.NRAmpTimeSeries)
        self.NRPhaseTimeSeriesPCA = \
                perform_pca(training_catalogue.NRPhaseTimeSeries)




        #
        # PCA of the SI data (if any)
        #
        if hasattr(training_catalogue, 'ref_mass'):
            print "PCA of SI waveforms"
            
            # XXX: Can add any extra conditioning (e.g., filters) here

            self.SIhplusTimeSeriesPCA = \
                    perform_pca(np.real(training_catalogue.SIComplexTimeSeries))
            self.SIhcrossTimeSeriesPCA = \
                    perform_pca(-1*np.imag(training_catalogue.SIComplexTimeSeries))
            self.SIAmpTimeSeriesPCA = \
                    perform_pca(training_catalogue.SIAmpTimeSeries)
            self.SIPhaseTimeSeriesPCA = \
                    perform_pca(training_catalogue.SIPhaseTimeSeries)
 

        # Project the test catalogue onto the new basis formed from the
        # training catalogue
        if test_catalogue is not None:
            print "Evaluating projection of test catalogue onto new basis"
            self.project_test_catalogue(test_catalogue)

    def project_test_catalogue(self, test_catalogue):
        """
        Project catalogue waveforms onto the new basis to find the
        coefficients for each waveform
        
        Updates the simulation_details object in the test catalogue with the
        projection onto this basis

        """

        print 'Projecting time-domain test catalogue onto new basis'

        self.test_catalogue_data = \
                list(test_catalogue.simulation_details.simulations)

        for w in xrange(len(self.test_catalogue_data)):

            self.test_catalogue_data[w]\
                    ['NRhplusTimeSeriesBetas'] = \
                    np.concatenate( 
                            self.NRhplusTimeSeriesPCA.transform(
                                np.real(test_catalogue.NRComplexTimeSeries[w,:])
                                ) 
                            )

            self.test_catalogue_data[w]\
                    ['NRhcrossTimeSeriesBetas'] = \
                    np.concatenate( 
                            self.NRhcrossTimeSeriesPCA.transform(
                                -1*np.imag(test_catalogue.NRComplexTimeSeries[w,:])
                                ) 
                            )

            self.test_catalogue_data[w]\
                    ['NRAmpTimeSeriesBetas'] = \
                    np.concatenate( 
                            self.NRAmpTimeSeriesPCA.transform(
                                test_catalogue.NRAmpTimeSeries[w,:]
                                ) 
                            )

            self.test_catalogue_data[w]\
                    ['NRPhaseTimeSeriesBetas'] = \
                    np.concatenate( 
                            self.NRPhaseTimeSeriesPCA.transform(
                                test_catalogue.NRPhaseTimeSeries[w,:]
                                ) 
                            )

            if hasattr(test_catalogue, 'ref_mass'):

                self.test_catalogue_data[w]\
                        ['SIhplusTimeSeriesBetas'] = \
                        np.concatenate( 
                                self.SIhplusTimeSeriesPCA.transform(
                                    np.real(test_catalogue.SIComplexTimeSeries[w,:])
                                    ) 
                                )

                self.test_catalogue_data[w]\
                        ['SIhcrossTimeSeriesBetas'] = \
                        np.concatenate( 
                                self.SIhcrossTimeSeriesPCA.transform(
                                    -1*np.imag(test_catalogue.SIComplexTimeSeries[w,:])
                                    ) 
                                )

                self.test_catalogue_data[w]\
                        ['SIAmpTimeSeriesBetas'] = \
                        np.concatenate( 
                                self.SIAmpTimeSeriesPCA.transform(
                                    test_catalogue.SIAmpTimeSeries[w,:]
                                    ) 
                                )

                self.test_catalogue_data[w]\
                        ['SIPhaseTimeSeriesBetas'] = \
                        np.concatenate( 
                                self.SIPhaseTimeSeriesPCA.transform(
                                    test_catalogue.SIPhaseTimeSeries[w,:]
                                    ) 
                                )


        return 0


    def file_dump(self, pca_attrs=['NRhplusTimeSeriesPCA'] , pcs_filename=None):
        """
        Dump the all the results to a python-friendly pickle and stick the PCs
        in ascii files; first row of the ascii file is the mean.

        (Pickling is pretty useful here, since I want to retain the object
        itself)

        pca_attr is a list of the PCA attributes you wish to dump
        """

        if pcs_filename is None:
            print >> sys.stderr, "ERROR: you must provide a name for catalogue \
file dumps"
            sys.exit(-1)

        # Output path
        try:
            pcs_path=os.path.join(os.environ['BHEX_PREFIX'],"data/PCA_data")
        except KeyError:
            print >> sys.stderr, "BHEX_PREFIX environment variable appears to be un-set"
            sys.exit()

        # Make output directory
        if not os.path.exists(pcs_path):
            os.makedirs(pcs_path)

        # Build filename for pickle with everything
        pcs_filename = os.path.join(pcs_path, pcs_filename)

        print "Performing data dump to %s.pickle"%pcs_filename
        
        # Pickle me!
        pickle.dump(self, open(pcs_filename+".pickle",'wb'))

        #
        # Ascii Dump
        #
        for pca_attr in pca_attrs:

            this_name  = os.path.join(pcs_path, pcs_filename + "_" + pca_attr + ".asc")
            print "Dumping to %s"%this_name
            
            pcaObj = getattr(self, pca_attr)


            # First row contains the mean waveform
            dims = np.shape(pcaObj.components_)
            output_array  = np.zeros(shape=(dims[0]+1,dims[1]))

            output_array[0,:]  = pcaObj.mean_
            output_array[1:,:] = pcaObj.components_

            np.savetxt(this_name, output_array)


        return 0

                     


# *******************************************************************************
def main():

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Training data input
    # XXX: probably going to want a config parser

    sample_rate = 512
    datalen= 4.0

    total_mass = 150. 
    distance=1. # Mpc

    train_series_names = ['HR-series']

    train_bounds=dict()
    train_bounds['a1'] = [0, 0]
    train_bounds['a2'] = [0, 0]
    train_bounds['q'] = [-np.inf, 3] 

    print '~~~~~~~~~~~~~~~~~~~~~'
    print 'Selecting Simulations'
    print ''
    train_simulations = \
            bwave.simulation_details(series_names=train_series_names,
                    param_bounds=train_bounds, Mmin30Hz=total_mass)

    print '~~~~~~~~~~~~~~~~~~~~~'
    print 'Building NR catalogue'
    print ''
    train_catalogue = bwave.waveform_catalogue(train_simulations,
            ref_mass=total_mass, sample_rate=sample_rate, datalen=datalen,
            distance=distance)


    #
    # Do the PCA
    #
    print '~~~~~~~~~~~~~~~~~~~~~'
    print 'Performing PCA'
    print ''
    pca = waveform_pca(train_catalogue, train_catalogue)

    return pca


#
# End definitions
#
if __name__ == "__main__":
    pca_results = main()

    



