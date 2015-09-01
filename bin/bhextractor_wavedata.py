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
bhextractor_wavedata.py

Module for loading and building catalogues from the GT_BBH_BURST_CATALOG

"""

from __future__ import division

import os
import sys 
import glob

import numpy as np
import scipy.signal as signal

#import lal
#import lalsimulation as lalsim
#import pycbc.types
#import pycbc.filter
#from pycbc.psd import aLIGOZeroDetHighPower

# *****************************************************************************
global __param_names__
__param_names__ = ['D', 'q', 'a1', 'a2', 'th1L', 'th2L', 'ph1', 'ph2', 'th12',
                'thSL', 'thJL', 'Mmin30Hz', 'Mmin10Hz']

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

    nonzero=np.argwhere(abs(input_data)>1e-3*max(abs(input_data)))
    idx = range(nonzero[0],nonzero[-1])


    win = planckwin(len(idx), 0.3)
    win[0.5*len(win):] = 1.0

#   from matplotlib import pyplot as pl
#   pl.figure()
#   pl.plot(input_data[idx]/input_data[idx].max())
    input_data[idx] *= win
    #input_data *= win

#   pl.plot(input_data[idx]/input_data[idx].max())
#   pl.plot(win)
#   pl.show()
#   sys.exit()

    return input_data

def planckwin(N, epsilon):

    t1 = -0.5*N
    t2 = -0.5*N * (1.-2.*epsilon)
    t3 = 0.5*N * (1.-2.*epsilon)
    t4 = 0.5*N

    Zp = lambda t: (t2-t1)/(t-t1) + (t2-t1)/(t-t2)
    Zm = lambda t: (t3-t4)/(t-t3) + (t3-t4)/(t-t4)

    win = np.zeros(N)
    ts = np.arange(-0.5*N, 0.5*N)

    for n,t in enumerate(ts):
        if t<=t1:
            win[n] = 0.0
        elif t1<t<t2:
            win[n] = 1./(np.exp(Zp(t))+1)
        elif t2<=t<=t3:
            win[n] = 1.0
        elif t3<t<t4:
            win[n] = 1./(np.exp(Zm(t))+1)

    return win


# *****************************************************************************
# Class Definitions 

                     
class waveform_data:
    """
    Object with the full waveform catalogue for the chosen series
    """

    def __init__(self, series_names='RO3-series', Mmin30Hz=100.,
            param_bounds=None):

        # ######################################################
        # Other initialisation
        self.series_names = series_names
        self.Mmin30Hz = Mmin30Hz

        print "____"
        print "Finding matching waveforms"

        # Get all waveforms
        all_simulations = self.list_simulations(series_names)

        # Down-select to those with a desirably small minimum mass
        self.simulations = self._select_param_values(all_simulations, 'Mmin30Hz',
                [-np.inf, self.Mmin30Hz])

        # Now down-select on any other parameters
        if param_bounds is not None:
            for bound_param in param_bounds.keys():
                self.simulations = self._select_param_values(self.simulations,
                        bound_param, param_bounds[bound_param])

        print "----"
        print "Found %d waveforms matching criteria:"%(len(self.simulations))
        print "Series: ", series_names
        print "Mmin30Hz: ", Mmin30Hz
        print "Bounds: ", param_bounds
        print "____"


    def list_simulations(self,series_names,
            catdir="GT_BBH_BURST_CATALOG"):
        """
        Creates a list of simulation diectionaries which contain the locations
        of the data files and the physical parameters the requested series

        TODO: add support for selecting waveforms based on physical parameters
        """


        valid_series = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
                "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series"
                "TP-series"]

        if type(series_names)!=list:
            series_names = [series_names]

        simulations = []
        for series_name in series_names:
            if series_name not in valid_series:
                print >> sys.stderr, "ERROR: series name (%s)must be in "%(series_name), valid_series
                sys.exit()

            datadir = os.path.join(os.environ['BHEX_PREFIX'], 'data/NR_data',
                    catdir, series_name)

            readme_file = os.path.join(datadir, 'README_%s.txt'%series_name)

            simulations += self._get_series(datadir,readme_file)

        return simulations

    @staticmethod
    def _select_param_values(simulations, param, bounds):
        """
        Return the list of simulations with parameter values in [low_bound,
        upp_bound]
        """
        return [sim for sim in simulations if sim[param] >= min(bounds) and
                sim[param]<=max(bounds) ]

    @classmethod
    def _get_series(cls, datadir, readme_file):
        """
        Read the parameters from the readme file and return a list of simulation
        dictionaries
        """
        readme_data = np.loadtxt(readme_file, dtype=str)

        #param_names = ['D', 'q', 'a1', 'a2', 'th1L', 'th2L', 'ph1', 'ph2', 'th12',
        #        'thSL', 'thJL', 'Mmin30Hz', 'Mmin10Hz']

        simulations = []
        for s in xrange(len(readme_data)):

            sim = dict()

            param_vals = [float(param) for param in readme_data[s,2:]]
            # physical params
            for p,param_name in enumerate(__param_names__):
                sim[param_name] = param_vals[p]

            # book-keeping
            sim['wavename'] = readme_data[s,1]
            sim['wavefile'] = glob.glob(os.path.join(datadir,sim['wavename'],'*asc'))

            if len(sim['wavefile'])>1:
                print >> sys.stderr, "Error more than one data file in directory: %s"%(
                        os.path.join(datadir,sim['wavename']))
                sys.exit(-1)

            simulations.append(sim)

        return simulations

    def build_waveforms(self):
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
        NRD_sampls=0.2*self.fs   # keep this many samples after the peak
        NINSP_sampls=0.75*self.fs   # retain this many samples before the peak

        # XXX: we should compute this from Karan's information..

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
                    #sYlm = lal.SpinWeightedSphericalHarmonic(0, 0, 
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

            # --- Resampling 
            hplus  = signal.resample(hplus, resamp_len)
            hcross = signal.resample(hcross, resamp_len)

            # Get rid of junk at end
            peak_idx = np.argmax(abs(hplus-1j*hcross))
            zero_idx = min(peak_idx + NRD_sampls, len(hplus))

            hplus[zero_idx:]  = 0.0
            hcross[zero_idx:] = 0.0

            # Get rid of junk at start
            startidx = max(0,np.argmax(abs(hplus-1j*hcross)) - NINSP_sampls)

            hplus[:startidx] = 0.0
            hcross[:startidx] = 0.0

            # --- Windowing / tapering
            hplus=window_wave(hplus)
            hcross=window_wave(hcross)

            #hplus = taper_start(hplus)
            #hcross = taper_start(hcross)

            # --- Filtering
            #hplus  = highpass(hplus)
            #hcross = highpass(hcross)

            # Store complex waveform
            catalogue[w,:] = hplus - 1j*hcross

            amplitude_catalogue[w,:] = abs(catalogue[w,:])
            phase_catalogue[w,:] = np.unwrap(np.angle(catalogue[w,:]))
            
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Plus/Cross Catalogue Conditioning

        # Standardised catalogue 
        self.aligned_catalogue = \
                np.zeros(shape=(len(waveforms),self.catalogue_len),
                        dtype=complex)

        print 'aligning peak times for +,x waves'

        # Find peak indices
        peak_idx=np.argmax(abs(catalogue),axis=0)

        # Align all waveform peaks to the 0.5 of the way through the final catalogue
        #align_idx=np.floor(0.75*self.catalogue_len)
        align_idx=np.floor(0.5*self.catalogue_len)

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
        align_idx=np.floor(0.5*self.catalogue_len)

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
            self.aligned_amplitudes[w,:] /= \
                    np.linalg.norm(self.aligned_amplitudes[w,:])
            #self.aligned_phases[w,:] /= \
            #        np.linalg.norm(self.aligned_phases[w,:])

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

# *******************************************************************************
def main():

    print '-----'
    print 'Generating example waveform list'
    print 'waveform list is a class with a list of \
dictionaries where keys are physical parameters and file paths.  The user \
specifies single or multiple waveform series from the GT catalogue and, \
optionally, some constraints on the physical parameters.'

    # select waveforms with mass ratio >= 2
    bounds = dict()
    bounds['q'] = [2, np.inf]

    waves = waveform_data(series_names='RO3-series', Mmin30Hz=100.,
                param_bounds=bounds)

    return waves

#
# End definitions
#
if __name__ == "__main__":

    waveform_list = main()
    



