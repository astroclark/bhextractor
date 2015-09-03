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

import lal
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

def phase_of(z):
    return np.unwrap(np.angle(z))


# *****************************************************************************
# Class Definitions 

                     
class simulation_details:
    """
    The waveform catalogue for the chosen series (possibly plural) with
    user-specified parameter bounds.
    
    Example Usage:

    In [24]: bounds = dict()

    In [25]: bounds['q'] = [2, np.inf]

    In [26]: simcat = bwave.simulation_details(series_names=['Eq-series',
    'RO3-series'], param_bounds=bounds, Mmin30Hz=100)

    In [28]: for sim in simcat.simulations: print sim
    {'a1': 0.6, 'th2L': 90.0, 'D': 7.0, 'thJL': 19.8, 'th1L': 90.0, 'q': 2.5, 'th12': 180.0, 'a2': 0.6,  'wavefile': ['/home/jclark308/Projects/bhextractor/data/NR_data/GT_BBH_BURST_CATALOG/Eq-series/Eq_D7_q2.50_a0.6_ph270_m140/Strain_jinit_l2_m2_r75_Eq_D7_q2.50_a0.6_ph270_m140.asc'], 'wavename': 'Eq_D7_q2.50_a0.6_ph270_m140', 'ph2': 90.0, 'ph1': -90.0, 'Mmin30Hz': 97.3, 'Mmin10Hz': 292.0, 'thSL': 90.0}

    ... and so on ...

    """

    def __init__(self, series_names='RO3-series', Mmin30Hz=100.,
            param_bounds=None):

        # ######################################################
        # Other initialisation
        self.series_names = series_names
        self.Mmin30Hz = Mmin30Hz

        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
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

        self.nsimulations = len(self.simulations)

        print "----"
        print "Found %d waveforms matching criteria:"%(self.nsimulations)
        print "Series: ", series_names
        print "Mmin30Hz: ", Mmin30Hz
        print "Bounds: ", param_bounds
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"


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

    @staticmethod
    def _get_series(datadir, readme_file):
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

            wavename = readme_data[s,1]
            wavefile = glob.glob(os.path.join(datadir, wavename, '*asc'))

            # Check that this waveform exists
            if len(wavefile)>1:
                print >> sys.stderr, "Error, more than one data file in directory: %s"%(
                        os.path.join(datadir,wavename))
                sys.exit(-1)
            elif len(wavefile)==0:
                print >> sys.stderr, "WARNING, No file matching glob pattern: \n%s"%(
                        os.path.join(datadir,wavename, '*asc'))
                continue
            sim['wavename'] = wavename
            sim['wavefile'] = wavefile[0]

            param_vals = [float(param) for param in readme_data[s,2:]]
            # physical params
            for p,param_name in enumerate(__param_names__):
                sim[param_name] = param_vals[p]

            simulations.append(sim)

        return simulations

class waveform_catalogue:
    """
    This contains the waveform data (as +,x, Amp and phase time-series and
    F-domain waveforms) constructed from the catalogue contained in a
    simulation_details class

    """

    def __init__(self, simulation_details, ref_mass=None, distance=1,
            sample_rate=1024, datalen=4): 
        """
        Initalise with a simulation_details object and a reference mass ref_mass to
        which waveforms get scaled.  If ref_mass is None, waveforms are left in
        code units
        """

        self.simulation_details = simulation_details

        # Load in the NR waveform data
        self.load_wavedata()

        # Produce physical catalogue if reference mass specified
        if ref_mass is not None:
            # catalogue_to_SI() adds variables to self
            self.catalogue_to_SI(ref_mass, sample_rate, distance, datalen)

    def load_wavedata(self):
        """
        Load the waveform data pointed to by the simulation_details object

        NOTE: Currently set up to compute amplitude and phase timeseries from
        the resampled complex time series.  This is probably the wrong way
        round; the amp/phase representation is smoother and might be a more
        sensible thing to resample.  If weirdness is encountered, I recommend
        investigating that option as a resolution.
        """

        # Load the waveform data into lists for time, plus and cross.  When all
        # data is loaded, identify the longest waveform and insert all waveforms
        # into a numpy array of that length

        time_data  = []
        plus_data  = []
        cross_data = []

        # Build waveforms in this catalogue
        max_time=-np.inf
        self.NR_deltaT=np.inf

        for w, sim in enumerate(self.simulation_details.simulations):

            print 'Loading %s waveform'%sim['wavename']

            wavedata = np.loadtxt(sim['wavefile'])

            time_data.append(wavedata[:,0])
            plus_data.append(wavedata[:,1])
            cross_data.append(wavedata[:,2])

            # Get characteristics to standardise waveforms
            if wavedata[-1,0]>max_time:
                max_time = np.copy(wavedata[-1,0])

            if np.diff(wavedata[:,0])[0]<self.NR_deltaT:
                self.NR_deltaT = np.diff(wavedata[:,0])[0]


        # Now resample to a consistent deltaT
        time_data_resampled  = []
        plus_data_resampled  = []
        cross_data_resampled = []

        print 'Resampling to uniform rate'
        max_length = -np.inf
        for w in xrange(self.simulation_details.nsimulations): 
            deltaT = np.diff(time_data[w])[0]
            if deltaT > self.NR_deltaT:
                resamp_len = deltaT / self.NR_deltaT * len(plus_data[w])
                plus_data_resampled.append(signal.resample(plus_data[w],
                    resamp_len))
                cross_data_resampled.append(signal.resample(cross_data[w],
                    resamp_len))
            else:
                plus_data_resampled.append(np.copy(plus_data[w]))
                cross_data_resampled.append(np.copy(cross_data[w]))

            if len(plus_data_resampled[w])>max_length:
                max_length = len(plus_data_resampled[w])

        self.NRdata_len = 2*max_length

        #
        # Insert into a numpy array
        #
        self.NRComplexTimeSeries = np.zeros(shape=(len(plus_data_resampled),
            self.NRdata_len), dtype=complex)
        self.NRAmpTimeSeries = np.zeros(shape=(len(plus_data_resampled),
            self.NRdata_len), dtype=complex)
        self.NRPhaseTimeSeries = np.zeros(shape=(len(plus_data_resampled),
            self.NRdata_len), dtype=complex)

        # Alignment
        align_idx = 0.5*self.NRdata_len
        for w in xrange(self.simulation_details.nsimulations):

            wave = plus_data_resampled[w] - 1j*cross_data_resampled[w]

            peak_idx=np.argmax(abs(wave))
            start_idx = align_idx - peak_idx

            self.NRComplexTimeSeries[w,start_idx:start_idx+len(wave)] = wave
            self.NRAmpTimeSeries[w,:] = abs(self.NRComplexTimeSeries[w,:])
            self.NRPhaseTimeSeries[w,:] = phase_of(self.NRComplexTimeSeries[w,:])

        del time_data, plus_data, cross_data, plus_data_resampled, cross_data_resampled

        return 0

    # XXX: TRUNCATION HERE?

    def catalogue_to_SI(self, ref_mass, sample_rate=1024, distance=100.,
            datalen=4):
        """
        Convert waveforms in self.NRdata to physical time stamps / amplitudes
        """

        # Add physical attributes
        self.ref_mass = ref_mass
        self.sample_rate=sample_rate
        self.distance=distance

        # Time steps of NR waveforms at the reference mass: 
        SI_deltaT = ref_mass * lal.MTSUN_SI * self.NR_deltaT
        Mscale = ref_mass * lal.MRSUN_SI / (distance * 1e6 * lal.PC_SI)

        # Resample to duration in seconds (=nsamples x deltaT) x samples / second
        resamp_len = np.ceil(self.NRdata_len*SI_deltaT*self.sample_rate)

        self.SIComplexTimeSeries = np.zeros(shape=(self.simulation_details.nsimulations,
            datalen*sample_rate), dtype=complex)
        self.SIAmpTimeSeries = np.zeros(shape=(self.simulation_details.nsimulations,
            datalen*sample_rate), dtype=complex)
        self.SIPhaseTimeSeries = np.zeros(shape=(self.simulation_details.nsimulations,
            datalen*sample_rate), dtype=complex)

        for w in xrange(self.simulation_details.nsimulations):

            resampled_re = signal.resample(np.real(self.NRComplexTimeSeries[w,:]), resamp_len)
            resampled_im = signal.resample(np.imag(self.NRComplexTimeSeries[w,:]), resamp_len)

            peakidx = np.argmax(abs(resampled_re + 1j*resampled_im))
            startidx = 0.5*datalen*sample_rate - peakidx

            self.SIComplexTimeSeries[w,startidx:startidx+len(resampled_re)] = \
                    Mscale*(resampled_re + 1j*resampled_im)

            self.SIAmpTimeSeries[w,:] = abs(self.SIComplexTimeSeries[w,:])
            self.SIPhaseTimeSeries[w,:] = phase_of(self.SIComplexTimeSeries[w,:])

        return 0

#   def fft_catalogue(self):
#       """
#       Add the amplitude and phase spectra of the waveforms from the TD
#       catalogue.  Might ultimately be useful for F-domain PCA, but for now
#       we'll just use it for plotting and diagnostics
#       """
#
#       print 'FFTing the catalogue'
#
#       # Get dims and add some info
#       example_td = pycbc.types.TimeSeries(
#               np.real(self.aligned_catalogue[0,:]),
#               delta_t=1.0/self.fs)
#
#       self.sample_times = example_td.sample_times - \
#               example_td.sample_times[np.argmax(example_td.data)]
#
#       self.sample_times_ampalign = example_td.sample_times - \
#               example_td.sample_times[np.argmax(self.aligned_amplitudes[0,:])]
#
#       self.flen = len(example_td.to_frequencyseries())
#
#
#       self.ampSpectraPlus = np.zeros(
#               shape=(np.shape(self.aligned_catalogue)[0], self.flen))
#       self.phaseSpectraPlus = np.zeros(
#               shape=(np.shape(self.aligned_catalogue)[0], self.flen))
#       self.ampSpectraCross = np.zeros(
#               shape=(np.shape(self.aligned_catalogue)[0], self.flen))
#       self.phaseSpectraCross = np.zeros(
#               shape=(np.shape(self.aligned_catalogue)[0], self.flen))
#
#       for w in xrange(len(self.waveform_names)):
#
#           tdwavePlus = pycbc.types.TimeSeries(
#                   np.real(self.aligned_catalogue[w,:]),
#                   delta_t=1.0/self.fs)
#
#           tdwaveCross = pycbc.types.TimeSeries(
#                   -1*np.imag(self.aligned_catalogue[w,:]),
#                   delta_t=1.0/self.fs)
#           
#           fdwavePlus  = tdwavePlus.to_frequencyseries()
#           fdwaveCross = tdwaveCross.to_frequencyseries()
#
#           self.ampSpectraPlus[w,:] = abs(fdwavePlus)
#           self.phaseSpectraPlus[w,:] = \
#                   np.unwrap(np.angle(fdwavePlus))
#
#           self.ampSpectraCross[w,:] = abs(fdwaveCross)
#           self.phaseSpectraCross[w,:] = \
#                   np.unwrap(np.angle(fdwaveCross))
#
#       self.sample_frequencies = fdwavePlus.sample_frequencies

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

    waves = simulation_details(series_names='RO3-series', Mmin30Hz=100.,
                param_bounds=bounds)

    return waves

#
# End definitions
#
if __name__ == "__main__":

    waveform_list = main()
    



