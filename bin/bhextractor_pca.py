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
import scipy.signal as signal
import scipy.linalg
from scipy.spatial.distance import euclidean as euclidean_distance

from sklearn.decomposition import PCA 

import lal
import lalsimulation as lalsim
import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

# *****************************************************************************
# Function Definitions 

def taper_start(input_data, fs=512):
    """  
    Taper the start of the data
    """

    timeseries = lal.CreateREAL8TimeSeries('blah', 0.0, 0,
            1.0/fs, lal.StrainUnit, int(len(input_data)))
    timeseries.data.data = input_data

    lalsim.SimInspiralREAL8WaveTaper(timeseries.data,
        lalsim.SIM_INSPIRAL_TAPER_START)

    return timeseries.data.data


def window_wave(input_data):

    nonzero=np.argwhere(abs(input_data)>0)
    idx = range(nonzero[0],nonzero[-1])
    beta = 25
    win = lal.CreateTukeyREAL8Window(len(idx), beta)
    win.data.data[int(0.5*len(idx)):]=1.0


    return input_data

def highpass(timeseries, delta_t=1./512, knee=9., order=8, attn=0.1):

    tmp = pycbc.types.TimeSeries(initial_array=timeseries, delta_t=delta_t)

    return np.array(pycbc.filter.highpass(tmp, frequency=knee, filter_order=order,
            attenuation=attn).data)

def perform_pca(aligned_td_catalogue):
    """
    Use the scikits-learn PCA tools to do PCA on a peak-aligned time-domain
    waveform catalogue
    """

    print 'performing PCA'

    #pca = PCA(whiten=True)
    pca = PCA(whiten=False)
    pca.fit(aligned_td_catalogue)

    return pca

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

        if not fd_interpolation:

            #
            # TD scaling
            #

            # Procedure: increasing the mass shortens the waveform; so interpolate
            # the original to time stamps which are mtotal_target / mtotal_ref x
            # closer together.  We do this by creating a new time axis which is
            # denser and shorter and then interpolate the *original* waveform FROM
            # those new timestamps TO the *original* time stamps.
            #from matplotlib import pyplot as pl

            time = np.arange(0,len(reconstruction)/float(fs), 1./fs)
            new_time = np.arange(0, time[-1] / amp_ratio + 1./fs / amp_ratio, 1./fs
                    / amp_ratio)[:len(time)]

            oripy = pycbc.types.TimeSeries(reconstruction, delta_t=1./fs)

            reconstruction = np.interp(new_time, time, reconstruction)*amp_ratio


            # Finally, clean up by moving the waveform back to the center (this
            # is just handy for plotting)
            non_zero_idx = \
                    np.argwhere(abs(reconstruction)>1e-5*max(abs(reconstruction)))
            trunc_wav = \
                    reconstruction[non_zero_idx[0]:non_zero_idx[-1]]

            peak_idx=np.argmax(abs(trunc_wav))
            start_idx = 0.5*len(reconstruction) - peak_idx

#           reconstruction_out = np.zeros(len(reconstruction))
#           #reconstruction_out[start_idx:start_idx+len(trunc_wav)] = trunc_wav
#           # populate right side
#           reconstruction_out[0.5*len(reconstruction):\
#                   0.5*len(reconstruction)+len(trunc_wav[peak_idx:non_zero_idx[-1]])] = \
#                           trunc_wav[peak_idx:]
#
#           # populate left side
#           p = 0.5*len(reconstruction_out)-1
#           q = peak_idx-1
#           while p>=0:
#               reconstruction_out[p] = np.copy(trunc_wav[q])
#               p-=1
#               q-=1
            reconstruction_out = reconstruction


        elif fd_interpolation:

            print "FD interpolation not yet supported"
            sys.exit()

            #
            # FD scaling
            #

            time = np.arange(0,len(reconstruction)/float(fs), 1./fs)
            from matplotlib import pyplot as pl

            f,ax = pl.subplots(nrows=3)

            ax[0].plot(time,reconstruction)

            recwav = pycbc.types.TimeSeries(reconstruction, delta_t=1./fs)
            recwav_fd = recwav.to_frequencyseries()
            recwav_mag = abs(recwav_fd)
            recwav_phase = np.unwrap(np.angle(recwav_fd))

            ax[1].plot(recwav_fd.sample_frequencies, recwav_mag)

            recwav_mag = np.interp(recwav_fd.sample_frequencies.data ,
                    recwav_fd.sample_frequencies.data / amp_ratio, recwav_mag)*amp_ratio

            ax[1].plot(recwav_fd.sample_frequencies, recwav_mag)
            ax[1].set_xlim(0,128)

            ax[2].plot(recwav_fd.sample_frequencies, recwav_phase)
            ax[2].set_xlim(0,128)

            recwav_phase = np.interp(recwav_fd.sample_frequencies.data,
                    recwav_fd.sample_frequencies.data / amp_ratio,
                    recwav_phase)/amp_ratio

            ax[2].plot(recwav_fd.sample_frequencies, recwav_phase)
            ax[2].set_xlim(0,128)

            recwav_fd = pycbc.types.FrequencySeries(
                    recwav_mag*np.exp(1j*recwav_phase),
                    delta_f=recwav_fd.delta_f)

            reconstruction = recwav_fd.to_timeseries()
            
            reconstruction = pycbc.filter.highpass(reconstruction, frequency=9.,
                    filter_order=8, attenuation=0.1)

            ax[0].plot(time, reconstruction)

            rec = reconstruction.to_frequencyseries()
            recmag = abs(rec)
            recphase = np.unwrap(np.angle(rec))
            ax[1].plot(rec.sample_frequencies, recmag)
            ax[2].plot(rec.sample_frequencies, recphase)

            pl.show()
            sys.exit()

    return reconstruction_out

def compute_projection(waveform, pca_result):
    """
    Compute the projection coefficients (betas) for the given waveform onto the
    PCA result
    """

    return np.concatenate(pca_result.transform(waveform))


# *****************************************************************************
# Class Definitions 

class waveform_pca:
    """
    The full PCA result, including projection coefficients etc. Takes an aligned
    waveform catalogue object as input.
    """

    def __init__(self, training_catalogue, testing_catalogue):

        self.training_catalogue = training_catalogue
        self.testing_catalogue  = testing_catalogue

        # Do the PCA
        self.pca_plus = perform_pca(np.real(self.training_catalogue.aligned_catalogue))
        self.pca_cross = perform_pca(-1*np.imag(self.training_catalogue.aligned_catalogue))

        self.pca_amp = perform_pca(self.training_catalogue.aligned_amplitudes)
        self.pca_phase = perform_pca(self.training_catalogue.aligned_phases)

        # Project the catalogue onto the new basis
        self.project_catalogue()

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
        fcomplex_ampphase = filename+"_ampPhasePCs.dat"

        # First row contains the mean waveform
        dims = np.shape(self.pca_plus.components_)
        output_array_plus  = np.zeros(shape=(dims[0]+1,dims[1]))
        output_array_cross = np.zeros(shape=(dims[0]+1,dims[1]))

        output_array = np.zeros(shape=(dims[0]+1,dims[1]), dtype=complex)
        output_array_ampphase = np.zeros(shape=(dims[0]+1,dims[1]), dtype=complex)

        output_array_plus[0,:] = self.pca_plus.mean_
        output_array_plus[1:,:] = self.pca_plus.components_
        output_array_cross[0,:] = self.pca_cross.mean_
        output_array_cross[1:,:] = self.pca_cross.components_

        output_array[0,:] = self.pca_plus.mean_ - 1j*self.pca_cross.mean_
        output_array[1:,:] = self.pca_plus.components_ - 1j*self.pca_cross.components_

        output_array_ampphase[0,:] = \
                self.pca_amp.mean_*np.exp(1j*self.pca_phase.mean_)
        output_array_ampphase[1:,:] = \
                self.pca_amp.components_*np.exp(1j*self.pca_phase.components_)

        print "Performing data dump to %s"%fplus
        np.savetxt(fplus, output_array_plus)
        print "Performing data dump to %s"%fcross
        np.savetxt(fcross, output_array_cross)

        print "Performing data dump to %s"%fcomplex_strain
        fpcomplex_strain = open(fcomplex_strain,"wb")
        output_array.tofile(fpcomplex_strain)

        print "Performing data dump to %s"%fcomplex_ampphase
        fpcomplex_ampphase = open(fcomplex_ampphase,"wb")
        output_array_ampphase.tofile(fpcomplex_ampphase)


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

                     
class waveform_catalogue:
    """
    Object with the full waveform catalogue and all conditioning information
    """

    def __init__(self,catalogue_name='Q',fs=512,
            catalogue_len=4, mtotal_ref=250, Dist=1.):

        # ######################################################
        # Other initialisation
        self.catalogue_name = catalogue_name
        self.fs = fs
        self.catalogue_len = catalogue_len * fs # XXX
        self.mtotal_ref = mtotal_ref
        self.Dist = Dist

        #
        # Build the catalogue
        #
        self.build_catalogue()
        self.fft_catalogue()
        self.nwaves = np.shape(self.aligned_catalogue)[0]

    def fft_catalogue(self):
        """
        Add the amplitude and phase spectra of the waveforms from the TD
        catalogue.  Might ultimately be useful for F-domain PCA, but for now
        we'll just use it for plotting and diagnostics
        """

        print 'FFTing the catalogue'

        # Get dims and add some info
        example_td = pycbc.types.TimeSeries(
                np.real(self.aligned_catalogue[0,:]),
                delta_t=1.0/self.fs)

        self.sample_times = example_td.sample_times - \
                example_td.sample_times[np.argmax(example_td.data)]

        self.sample_times_ampalign = example_td.sample_times - \
                example_td.sample_times[np.argmax(self.aligned_amplitudes[0,:])]

        self.flen = len(example_td.to_frequencyseries())


        self.ampSpectraPlus = np.zeros(
                shape=(np.shape(self.aligned_catalogue)[0], self.flen))
        self.phaseSpectraPlus = np.zeros(
                shape=(np.shape(self.aligned_catalogue)[0], self.flen))
        self.ampSpectraCross = np.zeros(
                shape=(np.shape(self.aligned_catalogue)[0], self.flen))
        self.phaseSpectraCross = np.zeros(
                shape=(np.shape(self.aligned_catalogue)[0], self.flen))

        for w in xrange(len(self.waveform_names)):

            tdwavePlus = pycbc.types.TimeSeries(
                    np.real(self.aligned_catalogue[w,:]),
                    delta_t=1.0/self.fs)

            tdwaveCross = pycbc.types.TimeSeries(
                    -1*np.imag(self.aligned_catalogue[w,:]),
                    delta_t=1.0/self.fs)
            
            fdwavePlus  = tdwavePlus.to_frequencyseries()
            fdwaveCross = tdwaveCross.to_frequencyseries()

            self.ampSpectraPlus[w,:] = abs(fdwavePlus)
            self.phaseSpectraPlus[w,:] = \
                    np.unwrap(np.angle(fdwavePlus))

            self.ampSpectraCross[w,:] = abs(fdwaveCross)
            self.phaseSpectraCross[w,:] = \
                    np.unwrap(np.angle(fdwaveCross))

        self.sample_frequencies = fdwavePlus.sample_frequencies


    def build_catalogue(self):
        """
        Return a complex array with the time-domain waveforms for the requested
        catalogue
        """

        print >> sys.stderr, "WARNING: Turning off orientation and only using 2,2; orientation will be handled by LIB"

        # ######################################################
        # Some catalogue-specific stuff

        # Dictionary with the maximum waveform data lengths
        maxlengths={'Q':10959, 'HR':31238, 'RO3':16493}

        # Dictionary of NR sample rates. XXX Should probably get this from the data,
        # really
        NR_deltaT={'Q':0.155, 'HR':0.08, 'RO3':2./15}

        # Samples to discard due to NR noise
        NRD_sampls=1000 # keep this many samples after the peak
        NINSP_sampls=0 # discard this many samples from the start

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
        amplitude_catalogue=np.zeros(shape=(len(waveforms), resamp_len))
        phase_catalogue=np.zeros(shape=(len(waveforms), resamp_len))

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


                # Extract l,m from filename
                spharm_degree = int(mode.split('_')[2].replace('l',''))
                spharm_order  = int(mode.split('_')[3].replace('m',''))

                if spharm_degree==spharm_order==2:

                    # Load mode data
                    mode_data = np.loadtxt(waveform_dir+'/'+mode)

                    # Compute spin -2 weighted spherical harmonic
                    #sYlm = lal.SpinWeightedSphericalHarmonic(self.theta, self.phi, 
                    #        -2, spharm_degree, spharm_order)

                    # Orient Waveforms
                    #hplus[0:len(mode_data)]  += mode_data[:,1]*np.real(sYlm) + \
                    #        mode_data[:,2]*np.imag(sYlm)
                    #hcross[0:len(mode_data)] += mode_data[:,2]*np.real(sYlm) - \
                    #        mode_data[:,1]*np.imag(sYlm)

                    hplus[0:len(mode_data)]  += mode_data[:,1]
                    hcross[0:len(mode_data)] += mode_data[:,2]

            #
            # Clean up the waveform
            #

            # Get rid of junk
            peak_idx = np.argmax(abs(hplus))
            zero_idx = min(peak_idx + NRD_sampls, len(hplus))

            hplus[zero_idx:]  = 0.0
            hcross[zero_idx:] = 0.0

            startidx = np.argwhere(abs(hplus)>0)[0]
            hplus[startidx:startidx+NINSP_sampls] = 0.0
            hcross[startidx:startidx+NINSP_sampls] = 0.0

            #hplus[:peak_idx-NINSP_sampls] = 0.0
            #hcross[:peak_idx-NINSP_sampls] = 0.0

            # --- Windowing / tapering
            #hplus=window_wave(hplus)
            #hcross=window_wave(hcross)
            hplus = taper_start(hplus)
            hcross = taper_start(hcross)


            # --- Resampling 
            hplus  = signal.resample(hplus, resamp_len)
            hcross = signal.resample(hcross, resamp_len)

            # --- Filtering
            #hplus  = highpass(hplus)
            #hcross = highpass(hcross)

            # Store complex waveform
            catalogue[w,:] = hplus - 1j*hcross

            amplitude_catalogue[w,:] = abs(hplus - 1j*hcross)
            phase_catalogue[w,:] = np.unwrap(np.angle(hplus - 1j*hcross))
            
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Plus/Cross Catalogue Conditioning

        # Standardised catalogue (4 seconds long)
        self.aligned_catalogue = \
                np.zeros(shape=(len(waveforms),self.catalogue_len),
                        dtype=complex)

        print 'aligning peak times for +,x waves'

        # Find peak indices
        peak_idx=np.argmax(abs(catalogue),axis=0)

        # Align all waveform peaks to the 0.5 of the way through the final catalogue
        align_idx=np.floor(0.75*self.catalogue_len)

        for w in xrange(len(waveforms)):
            #print 'aligning %d of %d'%(w, len(waveforms))

            non_zero_idx = \
                    np.argwhere(abs(catalogue[w,:])>1e-3*max(abs(catalogue[w,:])))
            trunc_wav = \
                    catalogue[w,non_zero_idx[0]:non_zero_idx[-1]]

            peak_idx=np.argmax(abs(trunc_wav))
            start_idx = align_idx - peak_idx

            self.aligned_catalogue[w,start_idx:start_idx+len(trunc_wav)] = trunc_wav

            # --- Normalisation
            self.aligned_catalogue[w,:] /= \
                    np.linalg.norm(self.aligned_catalogue[w,:])


        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Amplitude / Phase Catalogue Conditioning

        # Standardised catalogue (4 seconds long)
        self.aligned_amplitudes = \
                np.zeros(shape=(len(waveforms),self.catalogue_len))

        self.aligned_phases = \
                np.zeros(shape=(len(waveforms),self.catalogue_len))

        print 'aligning peak times for +,x waves'

        # Find peak indices
        peak_idx=np.argmax(abs(amplitude_catalogue),axis=0)

        # Align all waveform peaks to the 0.5 of the way through the final catalogue
        align_idx=np.floor(0.75*self.catalogue_len)

        for w in xrange(len(waveforms)):
            #print 'aligning %d of %d'%(w, len(waveforms))

            non_zero_idx = \
                    np.argwhere(abs(amplitude_catalogue[w,:])>1e-3*max(abs(amplitude_catalogue[w,:])))
            amp_trunc_wav = \
                    amplitude_catalogue[w,non_zero_idx[0]:non_zero_idx[-1]]
            phase_trunc_wav = \
                    phase_catalogue[w,non_zero_idx[0]:non_zero_idx[-1]]

            peak_idx=np.argmax(abs(amp_trunc_wav))
            start_idx = align_idx - peak_idx

            self.aligned_amplitudes[w,start_idx:start_idx+len(amp_trunc_wav)] = amp_trunc_wav
            self.aligned_phases[w,start_idx:start_idx+len(phase_trunc_wav)] = phase_trunc_wav

            # --- Normalisation
            #self.aligned_amplitudes[w,:] /= \
            #        np.linalg.norm(self.aligned_amplitudes[w,:])
            #self.aligned_phases[w,:] /= \
            #        np.linalg.norm(self.aligned_phases[w,:])


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
            catalogue_len=4, mtotal_ref=250, Dist=1.)

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

    



