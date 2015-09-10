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
import cPickle as pickle

import numpy as np
import scipy.signal as signal

import lal
import pycbc.types

# *****************************************************************************
global __param_names__
__param_names__ = ['D', 'q', 'a1', 'a2', 'th1L', 'th2L', 'ph1', 'ph2', 'th12',
                'thSL', 'thJL', 'Mmin30Hz', 'Mmin10Hz']

# *****************************************************************************
# Function Definitions 

def highpass(timeseries, delta_t=1./512, knee=9., order=12, attn=0.1):
    """
    Trivial interface function to make looping through catalogues neater
    """

    tmp = pycbc.types.TimeSeries(initial_array=timeseries, delta_t=delta_t)

    return np.array(pycbc.filter.highpass(tmp, frequency=knee, filter_order=order,
            attenuation=attn).data)

def window_wave(input_data):

    nonzero=np.argwhere(abs(input_data)>1e-3*max(abs(input_data)))
    idx = range(nonzero[0],nonzero[-1])
    win = planckwin(len(idx), 0.3)
    win[0.5*len(win):] = 1.0
    input_data[idx] *= win

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

    Note: will default to waveforms with min mass = 100 Msun permissable for a
    low-frequency-cutoff at 30 Hz, unless Mmin10Hz or a different minimum mass
    is defined

    """

    def __init__(self, series_names='RO3-series', Mmin30Hz=100.,
            Mmin10Hz=None, param_bounds=None):

        # ######################################################
        # Other initialisation
        self.series_names = series_names
        self.param_bounds = param_bounds

        self.Mmin10Hz = Mmin10Hz
        self.Mmin30Hz = Mmin30Hz
        if Mmin10Hz is not None:
            self.fmin = 10.0
        else:
            self.fmin = 30.0

        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print "Finding matching waveforms"

        # Get all waveforms
        self.simulations = self.list_simulations(series_names)
        self.nsimulations = len(self.simulations)

        print "----"
        print "Found %d waveforms matching criteria:"%(self.nsimulations)
        print "Series: ", series_names
        print "Bounds: ", param_bounds
        if Mmin10Hz is not None:
            print "Mmin10Hz: ", Mmin10Hz
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

            simulations += self._get_series(datadir,readme_file,series_name)

        # Down-select to those with a desirably small minimum mass
        if self.Mmin10Hz is not None:
            simulations = self._select_param_values(simulations, 'Mmin10Hz',
                    [-np.inf, self.Mmin30Hz])
        else:
            simulations = self._select_param_values(simulations, 'Mmin30Hz',
                    [-np.inf, self.Mmin30Hz])

        # Now down-select on any other parameters
        if self.param_bounds is not None:
            for bound_param in self.param_bounds.keys():
                simulations = self._select_param_values(simulations,
                        bound_param, self.param_bounds[bound_param])

        simulations = self._check_sim_unique(simulations)

        return simulations

    @staticmethod
    def _check_sim_unique(simulations):
        """
        Reduce the list of simulations to those which have unique parameter
        combinations

        XXX: Note that this (currently) takes the FIRST unique simulation; it
        doesn't care what the resolution or series was
        """
        print "Ensuring uniqueness of simulations"


        unique_simulations = list(simulations)

        physical_params = 'q', 'a1', 'a2', 'th1L', 'th2L', 'ph1', 'ph2', \
                'th12', 'thSL', 'thJL'

        param_sets = []
        for s in xrange(len(simulations)):
            param_vals = np.zeros(len(physical_params))
            for p,param_name in enumerate(physical_params):
                param_vals[p] = simulations[s][param_name]
            param_vals[np.isnan(param_vals)] = np.inf
            param_sets.append(tuple(param_vals))

        unique_param_sets = list(set(param_sets))

        for unique_param_set in unique_param_sets:

            indices = [i for i, x in enumerate(param_sets) if [x] ==
                    [unique_param_set]]

            if len(indices)>1:
                for index in indices[1:]:
                    unique_simulations.remove(simulations[index])

        return unique_simulations

    @staticmethod
    def _select_param_values(simulations, param, bounds):
        """
        Return the list of simulations with parameter values in [low_bound,
        upp_bound]
        """
        return [sim for sim in simulations if sim[param] >= min(bounds) and
                sim[param]<=max(bounds) ]

    @staticmethod
    def _get_series(datadir, readme_file, series_name):
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
            sim['series'] = series_name

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
            SI_deltaT=1./1024, SI_datalen=4, NR_deltaT=0.1, NR_datalen=10000,
            trunc_time=False): 
        """
        Initalise with a simulation_details object and a reference mass ref_mass to
        which waveforms get scaled.  If ref_mass is None, waveforms are left in
        code units

        """

        self.simulation_details = simulation_details
        self.NR_deltaT = NR_deltaT
        self.NR_datalen = NR_datalen

        self.trunc_time = trunc_time
        # Load in the NR waveform data
        self.load_wavedata()


        # Produce physical catalogue if reference mass specified
        if ref_mass is not None:
            # catalogue_to_SI() adds variables to self
            self.catalogue_to_SI(ref_mass=ref_mass, SI_deltaT=SI_deltaT,
                    distance=distance, SI_datalen=SI_datalen)


            # assign some characteristics which will be useful elsewhere (e.g.,
            # making PSDs)
            self.SI_deltaT = SI_deltaT

            example_ts = \
                    pycbc.types.TimeSeries(np.real(self.SIComplexTimeSeries[0,:]),
                            delta_t=self.SI_deltaT)
            example_fs = example_ts.to_frequencyseries()
            self.SI_deltaF = example_fs.delta_f
            self.SI_flen = len(example_fs)
            del example_ts, example_fs


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

        for w, sim in enumerate(self.simulation_details.simulations):

            print 'Loading %s waveform (%s)'%(sim['wavename'], sim['series'])

            wavedata = np.loadtxt(sim['wavefile'])

            time_data.append(wavedata[:,0])
            plus_data.append(wavedata[:,1])
            cross_data.append(wavedata[:,2])

        # Now resample to a consistent deltaT
        time_data_resampled  = []
        plus_data_resampled  = []
        cross_data_resampled = []

        print 'Resampling to uniform rate'
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

        #
        # Insert into a numpy array
        #
        NR_nsamples = self.NR_datalen / self.NR_deltaT
        self.NRComplexTimeSeries = np.zeros(shape=(len(plus_data_resampled),
            NR_nsamples), dtype=complex)
        self.NRAmpTimeSeries = np.zeros(shape=(len(plus_data_resampled),
            NR_nsamples))
        self.NRPhaseTimeSeries = np.zeros(shape=(len(plus_data_resampled),
            NR_nsamples))

        # Alignment & Normalisation
        align_idx = 0.5*NR_nsamples

        trunc_epsilon = 1e-3
        trunc_len = np.inf
        for w in xrange(self.simulation_details.nsimulations):

            wave = window_wave(plus_data_resampled[w] -
                    1j*cross_data_resampled[w])

            # Normalisation
            wave /= np.vdot(wave,wave)

            peak_idx=np.argmax(abs(wave))
            start_idx = align_idx - peak_idx

            self.NRComplexTimeSeries[w,start_idx:start_idx+len(wave)] = wave

            self.NRAmpTimeSeries[w,:] = abs(self.NRComplexTimeSeries[w,:])
            self.NRPhaseTimeSeries[w,:] = \
                    phase_of(self.NRComplexTimeSeries[w,:])

            # TRUNCATION ALGORITHM:
                # 1) Loop through waveforms, find lengths of data with
                # Amp>epsilon*max(Amp)
                # 2) construct Planck window with length = min(lengths from 1)
                # 3) Apply to all time series
            if sum(self.NRAmpTimeSeries[w,:] > \
                    trunc_epsilon*max(self.NRAmpTimeSeries[w,:])) < trunc_len:
                trunc_len = sum(self.NRAmpTimeSeries[w,:] > \
                    trunc_epsilon*max(self.NRAmpTimeSeries[w,:]))
                trunc_idx = self.NRAmpTimeSeries[w,:] > \
                    trunc_epsilon*max(self.NRAmpTimeSeries[w,:])

        if self.trunc_time:

            for w in xrange(self.simulation_details.nsimulations):

                # Truncate to the length of the shortest waveform
                self.NRComplexTimeSeries[w,trunc_idx] = \
                        window_wave(self.NRComplexTimeSeries[w,trunc_idx])
                self.NRComplexTimeSeries[w,np.invert(trunc_idx)] = 0.0

                self.NRAmpTimeSeries[w,:] = abs(self.NRComplexTimeSeries[w,:])
                self.NRPhaseTimeSeries[w,:] = \
                        phase_of(self.NRComplexTimeSeries[w,:])
                

        del time_data, plus_data, cross_data, plus_data_resampled, cross_data_resampled

        return 0

    def catalogue_to_SI(self, ref_mass, SI_deltaT=1./1024, distance=100.,
            SI_datalen=4):
        """
        Convert waveforms in self.NRdata to physical time stamps / amplitudes
        """

        print "converting NR catalogue to SI"

        # Add physical attributes
        self.ref_mass = ref_mass
        self.SI_deltaT=SI_deltaT
        self.distance=distance

        # Time steps of NR waveforms at the reference mass: 
        SI_deltaT_of_NR = ref_mass * lal.MTSUN_SI * self.NR_deltaT
        Mscale = ref_mass * lal.MRSUN_SI / (distance * 1e6 * lal.PC_SI)

        # Resample to duration in seconds (=nsamples x deltaT) x samples / second
        resamp_len = np.ceil(self.NR_datalen/self.NR_deltaT*SI_deltaT_of_NR/self.SI_deltaT)

        self.SIComplexTimeSeries = \
                np.zeros(shape=(self.simulation_details.nsimulations,
                    SI_datalen/SI_deltaT), dtype=complex)

        self.SIAmpTimeSeries = \
                np.zeros(shape=(self.simulation_details.nsimulations,
                    SI_datalen/SI_deltaT))
        self.SIPhaseTimeSeries = \
                np.zeros(shape=(self.simulation_details.nsimulations,
                    SI_datalen/SI_deltaT))

        for w in xrange(self.simulation_details.nsimulations):

            resampled_re = signal.resample(np.real(self.NRComplexTimeSeries[w,:]), resamp_len)
            resampled_im = signal.resample(np.imag(self.NRComplexTimeSeries[w,:]), resamp_len)
            wave = resampled_re + 1j*resampled_im

            idx = abs(wave)>1e-3*max(abs(wave))
            wave = window_wave(wave)
            wave = wave[idx]

            peakidx = np.argmax(abs(wave))
            startidx = 0.5*SI_datalen/SI_deltaT - peakidx

            self.SIComplexTimeSeries[w,startidx:startidx+len(wave)] = \
                    Mscale*wave

            self.SIAmpTimeSeries[w,:] = abs(self.SIComplexTimeSeries[w,:])
            self.SIPhaseTimeSeries[w,:] = phase_of(self.SIComplexTimeSeries[w,:])

        return 0

    def file_dump(self, catalogue_name=None):
        """
        Dump the catalogue data to pickle.  This will be useful for reproducing
        PCA results later with consistent catalogues.
        """

        # Need to call it something but make the user do it rather than parsing
        # all the parameters
        if catalogue_name is None:
            print >> sys.stderr, "ERROR: you must provide a name for catalogue \
file dumps"
            sys.exit(-1)


        # Output path
        try:
            data_path=os.path.join(os.environ['BHEX_PREFIX'],"data/NR_data")
        except KeyError:
            print >> sys.stderr, "BHEX_PREFIX environment variable appears to be un-set"
            sys.exit(-1)

        # Make output directory
        if not os.path.exists(data_path):
            os.makedirs(data_path)


        # Build filename
        filename = os.path.join(data_path, catalogue_name)

        print "Performing data dump to %s.pickle"%filename
        
        # Pickle me!
        pickle.dump(self, open(filename+".pickle",'wb'))


# *******************************************************************************
def main():

    print '-----'
    print 'Generating example waveform list'
    print 'waveform list is a class with a list of \
dictionaries where keys are physical parameters and file paths.  The user \
specifies single or multiple waveform series from the GT catalogue and, \
optionally, some constraints on the physical parameters.'

    # Some Example Inputs
    sample_rate = 1024
    SI_datalen= 4.0

    total_mass = 150. 
    distance=1. # Mpc

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Select and build NR waveforms from GT Burst Catalogue

    #
    # As an example, use only the non-spinning waveforms in the HR-series with
    # mass ratio < 3:
    #

    series_names = ['HR-series']

    bounds=dict()
    bounds['a1'] = [0, 0]
    bounds['a2'] = [0, 0]
    bounds['q'] = [-np.inf, 3] 

    print '~~~~~~~~~~~~~~~~~~~~~'
    print 'Selecting Simulations'
    print ''
    simulations_list = simulation_details(series_names=series_names,
            param_bounds=bounds, Mmin30Hz=total_mass)

    print '~~~~~~~~~~~~~~~~~~~~~'
    print 'Building NR catalogue'
    print ''
    NR_catalogue = waveform_catalogue(simulations_list, ref_mass=total_mass,
            sample_rate=sample_rate, SI_datalen=SI_datalen, distance=distance)


    print '~~~~~~~~~~~~~~~~~~~~~'
    print 'Plotting'
    print ''
    from matplotlib import pyplot as pl
    f, ax = pl.subplots(nrows=2, figsize=(8,6))
    ax[0].plot(NR_catalogue.NRAmpTimeSeries.T)
    ax[0].set_xlabel('Time Sample #')
    ax[0].set_ylabel(r'Amplitude')

    ax[1].plot(NR_catalogue.NRPhaseTimeSeries.T)
    ax[1].set_xlabel('Time Sample #')
    ax[1].set_ylabel('Phase')

    pl.show()


    NR_catalogue.file_dump('test')
    

    return simulations_list, NR_catalogue

#
# End definitions
#
if __name__ == "__main__":

    simulations_list, NR_catalogue = main()
    



