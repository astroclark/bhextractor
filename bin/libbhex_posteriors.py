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

import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib import pyplot as pl

from optparse import OptionParser


from glue.ligolw import lsctables, table, utils, ligolw, ilwd
table.use_in(ligolw.LIGOLWContentHandler)
lsctables.use_in(ligolw.LIGOLWContentHandler)

import lal
import lalmetaio

import pycbc.types
import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

from pylal import Fr
from pylal import bayespputils as bppu
import triangle

from bhex_utils import bhex_wavedata as bwave
from bhex_utils import bhex_pca as bpca

def parser():

    # --- Command line input
    parser = OptionParser()
    parser.add_option("-f", "--sim-inspiral", type=str, default=None)
    parser.add_option("-n", "--npcs", type=int, default=1)
    parser.add_option("-t", "--event", type=int, default=0)

    (opts,args) = parser.parse_args()


    if len(args)==0:
        print >> sys.stderr, "Must supply a posterior samples file as a commandline argument"
        sys.exit()

    if not os.path.isfile(args[0]):
        print >> sys.stderr, "posterior samples file requested: %s does not exist"%args[0]
        sys.exit()


    if opts.sim_inspiral is not None:
        if not os.path.isfile(opts.sim_inspiral):
            print >> sys.stderr, "sim-inspiral file not found at: "%opts.sim_inspiral
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


def copy_sim_inspiral( row ):
    """
    Turn a lsctables.SimInspiral into a SWIG wrapped lalburst.SimInspiral
    """
    swigrow = lalmetaio.SimInspiralTable()
    for simattr in lsctables.SimInspiralTable.validcolumns.keys():
        if simattr in ["waveform", "source", "numrel_data", "taper"]:
            # unicode -> char* doesn't work
            setattr( swigrow, simattr, str(getattr(row, simattr)) )
        else:
            setattr( swigrow, simattr, getattr(row, simattr) )
    # FIXME: This doesn't get copied properly, and so is done manually here.
    swigrow.geocent_end_time = lal.LIGOTimeGPS(row.geocent_end_time, row.geocent_end_time_ns)
    return swigrow

def get_sims( sim_inspiral_file ):
    """
    Return a list of swig-wrapped siminspiral table entries
    """

    xmldoc = utils.load_filename(sim_inspiral_file,
            contenthandler=ligolw.LIGOLWContentHandler)

    sims = []
    try:
        sim_insp = table.get_table( xmldoc,
                lsctables.SimInspiralTable.tableName )
        sims.extend( map(copy_sim_inspiral, sim_insp) )
    except ValueError:
        if opts.verbose:
            print >> sys.stderr, "No SimInspiral table found, \
skipping..."


    return sims



# --- MAIN
def main():

    opts, args = parser()
    posfile = args[0]

    if opts.sim_inspiral is not None:
        #
        # Parse sim_inspiral table
        #
        sims = get_sims(opts.sim_inspiral)

        # Sim_inspiral row for this event:
        this_sim = sims[opts.event]

        # Retrieve the numrel data


    #
    # Load and parse posterior samples file
    #
    peparser = bppu.PEOutputParser('common')
    resultsObj = peparser.parse(open(posfile, 'r'))
    posterior = bppu.Posterior(resultsObj)


    return this_sim

    sys.exit()




    #
    # Load PC data
    #
    try:
        pcs_path=os.environ['BHEX_PREFIX']+"/data/PCA_data"
    except KeyError:
        print >> sys.stderr, "BHEX_PREFIX environment variable appears to be un-set"

    print >> sys.stdout, "Reading PCA file"
    pcadata = pickle.load(open("%s/%s-PCA.pickle"%(pcs_path,
        opts.catalogue),'r'))

#    basis_functions = pcadata.pca_plus.components_

    print >> sys.stdout, "Read PCA file"

    


    #
    # Compute matches
    #
    amp_betas = pcadata.projection_amp[opts.injection_name]
    phase_betas = pcadata.projection_amp[opts.injection_name]

    if opts.injection_file is not None:
        print 'Loading injection file for match calculations (from LIB)'
        injwav = np.loadtxt(opts.injection_file)[:,1]


    else:
        print 'Loading injection file for match calculations (from PCs)'

        reconstructed_amp = bhex.reconstruct_waveform(pcadata.pca_amp,
                amp_betas, len(amp_betas), mtotal_ref=opts.mtotal_ref,
                mtotal_target=opts.target_mass)
        reconstructed_phase = bhex.reconstruct_waveform(pcadata.pca_phase,
                phase_betas, len(phase_betas), mtotal_ref=opts.mtotal_ref,
                mtotal_target=opts.target_mass)

        injwav = reconstructed_amp * np.exp(1j*reconstructed_phase)


    print 'Computing matches: '

    #
    # Reconstruct the waveforms from the samples
    #
    reconstruction_results = reconstructions(posterior, pcadata, opts.npcs,
            injwav)

    # Get nominal reconstruction for this number of PCs
    reconstructed_amp = bhex.reconstruct_waveform(pcadata.pca_amp,
            amp_betas, opts.npcs, mtotal_ref=opts.mtotal_ref,
            mtotal_target=opts.target_mass)
    reconstructed_phase = bhex.reconstruct_waveform(pcadata.pca_phase,
            phase_betas, opts.npcs, mtotal_ref=opts.mtotal_ref,
            mtotal_target=opts.target_mass)
    nomwav = reconstructed_amp*np.exp(1j*reconstructed_phase)

    print compute_match(np.real(nomwav), np.real(injwav), delta_t=1.0/512)

    #
    # Plots
    #
    nother=7 # increase this for new parameters
    samples = np.zeros(shape=(len(posterior), 2*opts.npcs+nother)) #x2 for amp/phase
    samples[:,0] = np.concatenate(posterior['mtotal'].samples)
    #samples[:,0] = reconstruction_results['matches']
    samples[:,1] = reconstruction_results['matches']
    samples[:,2] = np.concatenate(posterior['time'].samples) - opts.trig_time
    samples[:,3] = np.concatenate(posterior['hrss'].samples)
    samples[:,4] = np.concatenate(posterior['theta_jn'].samples)
    samples[:,5] = np.concatenate(posterior['phi_orb'].samples)
    samples[:,6] = np.concatenate(posterior['psi'].samples)
    #samples[:,7] = np.concatenate(posterior['ra'].samples)
    #samples[:,8] = np.concatenate(posterior['dec'].samples)

    b=1
    for n in xrange(nother,opts.npcs+nother):
        samples[:,n] = np.concatenate(posterior['amp_beta'+str(b)].samples)
        b+=1

    b=1
    for n in xrange(opts.npcs+nother,2*opts.npcs+nother):
        samples[:,n] = np.concatenate(posterior['phase_beta'+str(b)].samples)
        b+=1

    labels=['M$_{\\rm total}$', 'Match', 'Time', 'hrss', 'inclination', 'phase',
            '$\\Psi$', 'A,$\\beta_1$', '$\\Phi$,$\\beta_1$', 'A,$\\beta_2$',
            '$\\Phi$,$\\beta_2$', 'A,$\\beta_3$', '$\\Phi$,$\\beta_3$']
            #'$\\Psi$', 'ra', 'dec', 'A,$\\beta_1$', '$\\Phi$,$\\beta_1$', 'A,$\\beta_2$',

    trifig = triangle.corner(samples, labels=labels[:2*opts.npcs+nother], 
            quantiles=[0.25, 0.5, 0.75])

    trifig.savefig('%s'%(args[0].replace('.dat','.png')))
        

    return reconstruction_results, posterior


#
# End definitions
#
if __name__ == "__main__":
        ans = main()




