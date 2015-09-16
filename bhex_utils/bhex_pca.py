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
import copy

import numpy as np
from scipy.spatial.distance import euclidean as euclidean_distance

from sklearn.decomposition import PCA 

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

from bhex_utils import bhex_wavedata as bwave

# *****************************************************************************
# Function Definitions 

def plot_fidelity_by_npc(fidelity_matrix, title=None,
        figsize=None, ylabel = '$\mathcal{M}$', legloc=None):

    from matplotlib import pyplot as pl

    # Calculations
    min_fidelity = fidelity_matrix.min(axis=0)
    max_fidelity = fidelity_matrix.max(axis=0)
    mean_fidelity = fidelity_matrix.mean(axis=0)
    median_fidelity = np.median(fidelity_matrix, axis=0)
    rms_fidelity = np.sqrt(np.mean(np.square(fidelity_matrix), axis=0))

    low, upp = np.percentile(fidelity_matrix, [10, 90], axis=0)  

    fig, ax = pl.subplots()
    if figsize is not None:
        fig.set_size_inches(figsize)

    center = ax.step(np.arange(1,len(min_fidelity)+1)-0.5, median_fidelity, color='r',
            label='median', where='post')

    ax.bar(np.arange(1,len(min_fidelity)+1)-0.5, bottom=low, height=upp-low,
            color='lightgrey', label='10th/90th percentile',
            edgecolor='lightgrey', width=1)

    lims=ax.step(np.arange(1,len(min_fidelity)+1)-0.5, min_fidelity, color='k', linestyle='--',
            label='min/max', where='post')

    ax.step(np.arange(1,len(min_fidelity)+1)-0.5, max_fidelity, color='k',
            linestyle='--', where='post')

    ax.minorticks_on()
    ax.set_xlabel('Number of PCs')
    ax.set_ylabel(ylabel)

    ax.set_xlim(1,np.shape(fidelity_matrix)[1]+1)
    ax.grid()

    leg = ax.legend(loc='lower right')


    return fig, ax

def compute_match(wave1, wave2, delta_t, low_frequency_cutoff, psd):
    """
    Helper function to put arrays wave1,2 into timeseries objects and return the
    match
    """

    ts1 = pycbc.types.TimeSeries(wave1, delta_t=delta_t)
    ts2 = pycbc.types.TimeSeries(wave2, delta_t=delta_t)

    match, idx = pycbc.filter.match(ts1, ts2,
            low_frequency_cutoff=low_frequency_cutoff, psd=psd)

    return match
    


def reconstruct_waveform(pca, betas, npcs, mtotal_ref=250.0,
        mtotal_target=250.0, fs=512, fd_interpolation=False):
    """
    Reconstruct the specified waveform (waveform_name) from the PC basis
    contatined in pca and the projection coefficients betas

    pca is the object returned by sklearn.decomposition.PCA()

    MIGHT NOT WORK ANYMORE
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

def perform_pca(data_matrix, norm=False):
    """
    Use the scikits-learn PCA tools to do PCA on a peak-aligned time-domain
    waveform catalogue
    """

    if norm:
        # normalise each row by its norm
        for row in data_matrix:
            row /= np.linalg.norm(row)

    #pca = PCA(whiten=True)
    pca = PCA(whiten=False)
    pca.fit(data_matrix)

    return pca


def reconstruct(pca, betas, npcs):
    reconstruction = np.copy(pca.mean_)

    for n in xrange(npcs):
        reconstruction+=betas[n] * pca.components_[n,:]
    return reconstruction


# *****************************************************************************
# Class Definitions 

class waveform_pca:
    """
    The full PCA result, including projection coefficients etc. Takes a
    waveform catalogue object as input.
    """

    def __init__(self, train_catalogue, test_catalogue=None):

        self.ntrain = train_catalogue.simulation_details.nsimulations
        if test_catalogue is None:
            self.ntest = None
        else:
            self.ntest = test_catalogue.simulation_details.nsimulations 

        #
        # PCA of the NR data
        #
        print "PCA of NR waveforms"


        # Should also do Fourier domain PCA here.  Remember you'll want to save
        # those PCs as time series, though.  and figure out the geometrical
        # frequencies etc.

        self.NRhplusTimeSeriesPCA = \
                perform_pca(np.real(train_catalogue.NRComplexTimeSeries),
                        norm=False)
        self.NRhcrossTimeSeriesPCA = \
                perform_pca(-1*np.imag(train_catalogue.NRComplexTimeSeries),
                        norm=True)
        self.NRAmpTimeSeriesPCA = \
                perform_pca(train_catalogue.NRAmpTimeSeries, norm=False)
        self.NRPhaseTimeSeriesPCA = \
                perform_pca(train_catalogue.NRPhaseTimeSeries)

        # Inherit useful characteristics of data
        self.NR_deltaT = train_catalogue.NR_deltaT

        #
        # PCA of the SI data (if any)
        #
        if hasattr(train_catalogue, 'ref_mass'):
            print "PCA of SI waveforms"
            
            # XXX: Could add any extra conditioning (e.g., filters) here

            self.SIhplusTimeSeriesPCA = \
                    perform_pca(np.real(train_catalogue.SIComplexTimeSeries),
                            norm=True)
            self.SIhcrossTimeSeriesPCA = \
                    perform_pca(-1*np.imag(train_catalogue.SIComplexTimeSeries),
                            norm=True)
            self.SIAmpTimeSeriesPCA = \
                    perform_pca(train_catalogue.SIAmpTimeSeries, norm=True)
            self.SIPhaseTimeSeriesPCA = \
                    perform_pca(train_catalogue.SIPhaseTimeSeries)


            # Inherit useful characteristics of data
            self.fmin = train_catalogue.simulation_details.fmin
            self.SI_deltaT = train_catalogue.SI_deltaT
            self.SI_deltaF = train_catalogue.SI_deltaF
            self.SI_flen = train_catalogue.SI_flen


        # Project the test catalogue onto the new basis formed from the
        # training catalogue
        if test_catalogue is not None:

            # XXX: First check whether we have appropriate catalogues for SI
            # decomposition.  This means we have BOTH a test and a train
            # catalogue and that they possess SI waveforms scaled to the same
            # reference mass
            if hasattr(train_catalogue, 'ref_mass') *  hasattr(test_catalogue,
                    'ref_mass') * (getattr(train_catalogue, 'ref_mass') ==
                            getattr(test_catalogue, 'ref_mass')):

                        print "Train and test catalogues have matching SI waveforms.\
Will perform SI PCA decomposition"
                
                        self.do_si_projection=True

            else:
                    print "Train and test catalogues DO NOT have matching SI waveforms.\
Will perform SI PCA decomposition (different masses)"

                    self.do_si_projection=False
                    # but we will do the NR decomposition

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
                copy.deepcopy(test_catalogue.simulation_details.simulations)

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

            if self.do_si_projection:

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

    def compute_projection_fidelity(self, psd=None):
        """
        Compute the reconstruction fidelity as a function of nPCs, for each
        waveform in the catalogue.  Optionally computes waveform match.

        match / euclidean_distance arrays are organised as
            number_of_test_waves x num_train_waves,

        where num_train_waves = number of PCs used to reconstruct the test
        waveform

        """ 

        print "Evaluating reconstruction fidelity"
 
        # Pre-allocate

        class euclidean_distances(object):
            hplus = np.zeros(shape=(self.ntest, self.ntrain))
            amplitude = np.zeros(shape=(self.ntest, self.ntrain))
            phase = np.zeros(shape=(self.ntest, self.ntrain))

        class projections(object):
            hplus = np.zeros(shape=(self.ntest, self.ntrain))
            ampphase = np.zeros(shape=(self.ntest, self.ntrain))

        if self.do_si_projection:
            class matches(object):
                hplus = np.zeros(shape=(self.ntest, self.ntrain))
                ampphase = np.zeros(shape=(self.ntest, self.ntrain))


        for w in xrange(self.ntest):
    
            # retrieve projection coefficients
            hplus_NR_betas = self.test_catalogue_data[w]['NRhplusTimeSeriesBetas']

            amp_NR_betas   = self.test_catalogue_data[w]['NRAmpTimeSeriesBetas']
            phase_NR_betas = self.test_catalogue_data[w]['NRPhaseTimeSeriesBetas']
            
            # --- The target waveforms
            # We can just use the PCs to rebuild these using the PCA() inverse
            # transform method instead of carrying all the data round
            # XXX: can only use inverse transform to go back to training data!

           # target_NR_hplus = self.NRhplusTimeSeriesPCA.inverse_transform(
           #         self.test_catalogue_data[w]['NRhplusTimeSeriesBetas'] )
           # target_NR_amp = self.NRAmpTimeSeriesPCA.inverse_transform(
           #         self.test_catalogue_data[w]['NRAmpTimeSeriesBetas'] )
           # target_NR_phase = self.NRPhaseTimeSeriesPCA.inverse_transform(
           #         self.test_catalogue_data[w]['NRPhaseTimeSeriesBetas'] )

            target_NR_hplus = reconstruct(self.NRhplusTimeSeriesPCA,
                    self.test_catalogue_data[w]['NRhplusTimeSeriesBetas'],
                    len(self.test_catalogue_data[w]['NRhplusTimeSeriesBetas']))
            target_NR_amp = reconstruct(self.NRAmpTimeSeriesPCA,
                    self.test_catalogue_data[w]['NRAmpTimeSeriesBetas'],
                    len(self.test_catalogue_data[w]['NRAmpTimeSeriesBetas']))
            target_NR_phase = reconstruct(self.NRPhaseTimeSeriesPCA,
                    self.test_catalogue_data[w]['NRPhaseTimeSeriesBetas'],
                    len(self.test_catalogue_data[w]['NRPhaseTimeSeriesBetas']))


            target_NR_ampphase = target_NR_amp*np.exp(1j*target_NR_phase)

            if self.do_si_projection:
                hplus_SI_betas = self.test_catalogue_data[w]['SIhplusTimeSeriesBetas']
                amp_SI_betas   = self.test_catalogue_data[w]['SIAmpTimeSeriesBetas']
                phase_SI_betas = self.test_catalogue_data[w]['SIPhaseTimeSeriesBetas']
                
                # The target waveform
#               target_SI_hplus = self.SIhplusTimeSeriesPCA.inverse_transform(
#                       self.test_catalogue_data[w]['SIhplusTimeSeriesBetas'] )
#               target_SI_amp = self.SIAmpTimeSeriesPCA.inverse_transform(
#                       self.test_catalogue_data[w]['SIAmpTimeSeriesBetas'] )
#               target_SI_phase = self.SIPhaseTimeSeriesPCA.inverse_transform(
#                       self.test_catalogue_data[w]['SIPhaseTimeSeriesBetas'] )

                target_SI_hplus = reconstruct(self.SIhplusTimeSeriesPCA,
                        self.test_catalogue_data[w]['SIhplusTimeSeriesBetas'],
                        len(self.test_catalogue_data[w]['SIhplusTimeSeriesBetas']))
                target_SI_amp = reconstruct(self.SIAmpTimeSeriesPCA,
                        self.test_catalogue_data[w]['SIAmpTimeSeriesBetas'],
                        len(self.test_catalogue_data[w]['SIAmpTimeSeriesBetas']))
                target_SI_phase = reconstruct(self.SIPhaseTimeSeriesPCA,
                        self.test_catalogue_data[w]['SIPhaseTimeSeriesBetas'],
                        len(self.test_catalogue_data[w]['SIPhaseTimeSeriesBetas']))



                target_SI_ampphase = target_SI_amp * np.exp(1j*target_SI_phase)

    
            for n,npcs in enumerate(xrange(1,self.ntrain+1)):
    
                reconstructed_NR_hplus = \
                        reconstruct(self.NRhplusTimeSeriesPCA, hplus_NR_betas,
                                npcs)
                reconstructed_NR_amp = \
                        reconstruct(self.NRAmpTimeSeriesPCA, amp_NR_betas, npcs)
                reconstructed_NR_phase = \
                        reconstruct(self.NRPhaseTimeSeriesPCA, phase_NR_betas,
                                npcs)

                if self.do_si_projection:
                    reconstructed_SI_hplus = \
                            reconstruct(self.SIhplusTimeSeriesPCA, hplus_SI_betas,
                                    npcs)
                    reconstructed_SI_amp = \
                            reconstruct(self.SIAmpTimeSeriesPCA, amp_SI_betas, npcs)
                    reconstructed_SI_phase = \
                            reconstruct(self.SIPhaseTimeSeriesPCA, phase_SI_betas,
                                    npcs)

                reconstructed_NR_ampphase = \
                        reconstructed_NR_amp*np.exp(1j*reconstructed_NR_phase)
 
                euclidean_distances.hplus[w,n] = euclidean_distance(
                   reconstructed_NR_hplus, target_NR_hplus)
                euclidean_distances.amplitude[w,n] =\
                   euclidean_distance(reconstructed_NR_amp, target_NR_amp)
                euclidean_distances.phase[w,n] =\
                        euclidean_distance(reconstructed_NR_amp, target_NR_amp)

                projections.hplus[w,n] = np.dot(
                        reconstructed_NR_hplus/np.linalg.norm(reconstructed_NR_hplus),
                        target_NR_hplus/np.linalg.norm(target_NR_hplus)
                        )

                projections.ampphase[w,n] = np.vdot(
                        reconstructed_NR_ampphase/np.linalg.norm(reconstructed_NR_ampphase),
                        target_NR_ampphase/np.linalg.norm(target_NR_ampphase)
                        )
    
                if self.do_si_projection:

                    reconstructed_SI_ampphase = \
                            reconstructed_SI_amp*np.exp(1j*reconstructed_SI_phase)

                    matches.hplus[w,n] = \
                            compute_match(reconstructed_SI_hplus, target_SI_hplus,
                                    delta_t=self.SI_deltaT, psd=psd,
                                    low_frequency_cutoff=self.fmin)
        
                    matches.ampphase[w,n] = \
                            compute_match(np.real(reconstructed_SI_ampphase),
                                    np.real(target_SI_ampphase), delta_t=self.SI_deltaT,
                                    psd=psd, low_frequency_cutoff=self.fmin)

        if self.do_si_projection:
            return euclidean_distances, projections, matches
        else:
            return euclidean_distances, projections



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
    SI_datalen= 4.0

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
            ref_mass=total_mass, sample_rate=sample_rate, SI_datalen=SI_datalen,
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

    



