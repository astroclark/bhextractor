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
from matplotlib import pyplot as pl
import cPickle as pickle

from optparse import OptionParser

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

from pylal import bayespputils as bppu

import triangle

import bhextractor_pca as bhex

def parser():

    # --- Command line input
    parser = OptionParser()
    #parser.add_option("-p", "--posterior-file", default=None, type=str)
    parser.add_option("-i", "--injection-name", type=str, default=None)
    parser.add_option("-f", "--injection-file", type=str, default=None)
    parser.add_option("-n", "--npcs", type=int, default=1)
    parser.add_option("-q", "--catalogue", type=str, default="Q")
    parser.add_option("-m", "--target-mass", type=float, default=250.0) 

    (opts,args) = parser.parse_args()

    if opts.injection_name is None:
        print >> sys.stderr, "Require an injection name (e.g., D10_a0.0_q1.00_m103_Qs)"
        sys.exit()


    if len(args)==0:
        print >> sys.stderr, "Must supply a posterior samples file as a commandline argument"
        sys.exit()

    if not os.path.isfile(args[0]):
        print >> sys.stderr, "posterior samples file requested: %s does not exist"%args[0]
        sys.exit()

    return opts, args


def compute_match(reconstruction, injection, delta_t=1./512):

    rec_td = pycbc.types.TimeSeries(reconstruction, delta_t=delta_t)
    inj_td = pycbc.types.TimeSeries(injection, delta_t=delta_t)

    rec_fd = rec_td.to_frequencyseries()
    inj_fd = inj_td.to_frequencyseries()

    psd    = aLIGOZeroDetHighPower(len(rec_fd), delta_f=inj_fd.delta_f,
            low_freq_cutoff=1)

    return pycbc.filter.match(rec_fd, inj_fd, psd=psd,
            low_frequency_cutoff=10)[0]

def reconstructions(posterior, pcaresults, npcs, injwav):

    basis_functions = pcaresults.pca_plus.components_

    betas = np.zeros(shape=(len(posterior),npcs))

    for b in xrange(npcs):
        betas[:,b] = np.concatenate(posterior['beta'+str(b+1)].samples)

    mtotal = np.concatenate(posterior['mtotal'].samples)

    # Compute th  weighted sum of basis functions
    rec_waveforms = np.zeros(shape=(len(posterior),
        np.shape(basis_functions)[1]))
    matches = np.zeros(shape=(len(posterior)))

    for p in xrange(len(posterior)):

        print '%d of %d'%(p, len(posterior))

        ptest = bhex.reconstruct_waveform(pcaresults.pca_plus,
                        betas[p,:], npcs=npcs, mtotal_target=mtotal[p])

        rec_waveforms[p,:] = bhex.reconstruct_waveform(pcaresults.pca_plus,
                betas[p,:], npcs=npcs, mtotal_target=mtotal[p])

        matches[p] = compute_match(rec_waveforms[p,:], injwav)



    recon_data = {'betas':betas, 'mtotal':mtotal,
        'reconstructed_waves':rec_waveforms, 'matches':matches,
        'injected':injwav}

    return recon_data

# --- MAIN

def main():

    opts, args = parser()

    #
    # Load PC data
    #
    try:
        pcs_path=os.environ['BHEX_PREFIX']+"/data/PCA_data"
    except KeyError:
        print >> sys.stderr, "BHEX_PREFIX environment variable appears to be un-set"

    print >> sys.stdout, "Reading PCA file"
    pcadata = pickle.load(open("%s/%s-0.0-PCAresults.pickle"%(pcs_path,
        opts.catalogue),'r'))
    basis_functions = pcadata.pca_plus.components_
    catalogue = pcadata.catalogue.aligned_catalogue

    
    #
    # Load and parse posterior samples file
    #
    peparser = bppu.PEOutputParser('common')
    resultsObj = peparser.parse(open(args[0], 'r'))
    posterior = bppu.Posterior(resultsObj)

    #print opts.injection_file
    #sys.exit()


    #
    # Compute matches
    #
    betas = pcadata.projection_plus[opts.injection_name]

    if opts.injection_file is not None:
        print 'Loading injection file for match calculations (from LIB)'
        injwav = np.loadtxt(opts.injection_file)[:,1]


    else:
        print 'Loading injection file for match calculations (from PCs)'
        idx = pcadata.catalogue.waveform_names.index(opts.injection_name)

        injwav = bhex.reconstruct_waveform(pcadata.pca_plus, betas,
                npcs=len(pcadata.catalogue.waveform_names), mtotal_target=opts.target_mass)



    print 'Computing matches: '

    #
    # Reconstruct the waveforms from the samples
    #
    reconstruction_results = reconstructions(posterior, pcadata, opts.npcs,
            injwav)

    # Get nominal reconstruction
    nomwav = bhex.reconstruct_waveform(pcadata.pca_plus, betas, opts.npcs,
            mtotal_target=opts.target_mass)

    print compute_match(nomwav, injwav, delta_t=1.0/pcadata.catalogue.fs)

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
        

    return reconstruction_results, posterior


#
# End definitions
#
if __name__ == "__main__":
        reconstruction_results, posterior = main()




