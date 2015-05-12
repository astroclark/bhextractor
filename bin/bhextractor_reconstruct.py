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
bhextractor_reconstruct.py

Loads a posterior samples file into a BayesPostProc results object and builds
the waveform.  If the injection is provided, matches for each sampled point are
also returned.
"""

from __future__ import division

import os
import sys 

import numpy as np
import scipy.io as sio
from matplotlib import pyplot as pl

from optparse import OptionParser

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

from pylal import bayespputils as bppu

import triangle

def parser():

    # --- Command line input
    parser = OptionParser()
    #parser.add_option("-p", "--posterior-file", default=None, type=str)
    parser.add_option("-i", "--injection-file", type=str, default=None)
    parser.add_option("-n", "--npcs", type=int, default=1)
    parser.add_option("-q", "--catalog", type=str, default="Q")
    parser.add_option("-t", "--true-vals") 

    (opts,args) = parser.parse_args()

    if len(args)==0:
        print >> sys.stderr, "Must supply a posterior samples file as a commandline argument"
        sys.exit()

    if not os.path.isfile(args[0]):
        print >> sys.stderr, "posterior samples file requested: %s does not exist"%args[0]
        sys.exit()

    return opts, args

def reconstructions(posterior, basis_functions, npcs):

    betas = np.zeros(shape=(len(posterior),npcs))
    for b in xrange(npcs):
        betas[:,b] = np.concatenate(posterior['beta'+str(b+1)].samples)

    mtotal = np.concatenate(posterior['mtotal'].samples)

    # Compute the weighted sum of basis functions
    rec_waveforms = np.zeros(shape=(len(posterior),
        np.shape(basis_functions)[0]), dtype=complex)
    for p in xrange(len(posterior)):

        # Sum the PCs
        pc_wave = np.zeros(np.shape(basis_functions)[0], dtype=complex)
        for b in xrange(npcs):
            pc_wave += basis_functions[:,b]*betas[p,b]


        # Scale for the mass
        rec_waveforms[p,:] = scale_for_mass(mtotal[p], pc_wave)

    recon_data = {'betas':betas, 'mtotal':mtotal,
            'reconstructed_fd_waves':rec_waveforms}

    return recon_data

def scale_for_mass(mtotal, pc_waveform, mtotal_ref=250., delta_t=1./2048,
        fftlen=4):

    freqs = np.arange(0, len(pc_waveform)/fftlen, 1./fftlen)

    magnitude = abs(pc_waveform)
    phase = np.unwrap(np.angle(pc_waveform))

    new_magnitude = np.interp(mtotal/mtotal_ref * freqs, freqs, magnitude)*\
            mtotal/mtotal_ref
    new_phase = np.interp(mtotal/mtotal_ref * freqs, freqs, magnitude)/\
            mtotal/mtotal_ref


#   pl.figure(figsize=(8,6))
#   pl.plot(freqs,magnitude)
#   pl.plot(freqs,new_magnitude)
#   pl.show()
#   sys.exit()

    return new_magnitude * np.exp(1j*new_phase)

def compute_match(reconstruction, injection, delta_f=0.25):

    rec_fd = pycbc.types.FrequencySeries(reconstruction, delta_f=delta_f)
    inj_fd = pycbc.types.FrequencySeries(injection, delta_f=delta_f)
    psd    = aLIGOZeroDetHighPower(len(rec_fd), delta_f=delta_f,
            low_freq_cutoff=1)

    return pycbc.filter.match(rec_fd, inj_fd, psd=psd,
            low_frequency_cutoff=10)[0]


# --- MAIN

def main():

    opts, args = parser()

    #
    # Load PC data
    #
    try:
        pcs_path=os.environ['BHEX_PREFIX']+"/data/PCA_data"
        catalog_path=os.environ['BHEX_PREFIX']+"/data/signal_data"
    except KeyError:
        print >> sys.stderr, "BHEX_PREFIX environment variable appears to be un-set"

    print >> sys.stdout, "Reading PCs file"
    pcdata = sio.loadmat("%s/%s_PCs_theta-0.mat"%(pcs_path, opts.catalog))
    basis_functions = pcdata['pcs_plus']

    catalog = sio.loadmat("%s/%s_catalogue_theta-0.mat"%(catalog_path, opts.catalog))

    
    #
    # Load and parse posterior samples file
    #
    peparser = bppu.PEOutputParser('common')
    resultsObj = peparser.parse(open(args[0], 'r'))
    posterior = bppu.Posterior(resultsObj)


    #
    # Reconstruct the waveforms from the samples
    #
    reconstruction_results = reconstructions(posterior, basis_functions,
            opts.npcs)

    #
    # Compute matches
    #
    if opts.injection_file is not None:
        print 'Loading injection file for match calculations'
        #injdata = np.loadtxt(opts.injection_file)
        #injfreq = injdata[:,0]
        #injwav  = injdata[:,1]+1j*injdata[:,2]
        injwav = \
                np.array(pycbc.types.TimeSeries(np.real(catalog['MDC_final'][:,12]),
                    delta_t=1./2048).to_frequencyseries())

        reconstruction_results['matches'] = np.zeros(len(posterior))

        print 'Computing matches: '
        for p in xrange(len(posterior)):
            print '%d of %d'%(p, len(posterior))
            reconstruction_results['matches'][p] = compute_match(\
                    reconstruction_results['reconstructed_fd_waves'][p,:],
                    injwav)

        # Get nominal reconstruction for Q1
        beta_true = [ -0.311405204456, -0.0629567276758, 0.225561207648]
        pc_wav = np.zeros(np.shape(basis_functions)[0], dtype=complex)
        for b in xrange(3):
            pc_wav += beta_true[b]*basis_functions[:,b]
 
        nomwav = scale_for_mass(400, pc_wav)
 
        print compute_match(nomwav, injwav)



    #
    # Plots
    #
    samples = np.zeros(shape=(len(posterior), opts.npcs+2))
    samples[:,0] = np.concatenate(posterior['mtotal'].samples)
    samples[:,1] = reconstruction_results['matches']
    b=1
    for n in xrange(2,opts.npcs+2):
        samples[:,n] = np.concatenate(posterior['beta'+str(b)].samples)
        b+=1

    labels=['M$_{\\rm total}$', 'Match', '$\\beta_1$', '$\\beta_2$', '$\\beta_3$',
            '$\\beta_4$', '$\\beta_5$','$\\beta_6$', '$\\beta_7$', '$\\beta_8$',
            '$\\beta_9$', '$\\beta_{10}$']


    trifig = triangle.corner(samples, labels=labels[:opts.npcs+2], 
            quantiles=[0.25, 0.5, 0.75])

    trifig.savefig('%s'%(args[0].replace('.dat','.png')))
        

    return reconstruction_results


#
# End definitions
#
if __name__ == "__main__":
        reconstruction_results = main()




