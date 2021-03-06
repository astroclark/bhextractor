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
bhextractor_plot_catalogue.py

Construct catalogues and principal component analysis for NR BBH waveforms

This script simply builds and plots a catalogue (defined by the user in this
script)

"""

import sys
from bhex_utils import bhex_wavedata as bwave
import burst_nr_utils as bnru
from pycbc.psd import aLIGOZeroDetHighPower
import pycbc.types
import numpy as np
from matplotlib import pyplot as pl

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

SI_deltaT = 1./4096
SI_datalen= 4.0

init_total_mass=100.0
Mchirpmin30Hz = 35.0
smallest_observed_chirpmass = 25.0
distance=1. # Mpc

#
# Modify for imposing parameter bounds on the catalogue:
#
bounds = dict()
bounds['Mchirpmin30Hz'] = [-np.inf, Mchirpmin30Hz]
bounds['a1'] = [0, 0]
bounds['a2'] = [0, 0]
bounds['q'] = [1, 1]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do the work:

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
simulations = \
        bwave.simulation_details(param_bounds=bounds)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
catalogue = bwave.waveform_catalogue(simulations, ref_mass=init_total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
        trunc_time=False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting
#
# Magnitude Spectra
#
fspecSI, axspecSI = pl.subplots()

for w, wave in enumerate(catalogue.SIComplexTimeSeries):

    hplus_NR = \
            pycbc.types.TimeSeries(np.real(wave), delta_t=SI_deltaT)
    #hplus_NR.data = bwave.window_wave(hplus_NR.data)

    # Scale this to the smallest OBSERVED chirp mass
    mtot_target = bnru.mtot_from_mchirp(smallest_observed_chirpmass,
            simulations.simulations[w]['q'])
    print 'scaling to ', mtot_target
    hplus_NR = bnru.scale_wave(hplus_NR, mtot_target, init_total_mass)
    hplus_NR = pycbc.types.TimeSeries(hplus_NR, delta_t=SI_deltaT)

    Hplus_NR = hplus_NR.to_frequencyseries()

    axspecSI.loglog(Hplus_NR.sample_frequencies,
            2*abs(Hplus_NR)*np.sqrt(Hplus_NR.sample_frequencies), color='grey')

    delta_f = Hplus_NR.delta_f
    flen = len(Hplus_NR)
    psd = aLIGOZeroDetHighPower(flen, delta_f, 0.1) 
    axspecSI.loglog(psd.sample_frequencies, np.sqrt(psd), label='aLIGO',
            color='k', linestyle='--')


axspecSI.set_xlim(5, 512)
axspecSI.axvline(30, color='r')
axspecSI.minorticks_on()
axspecSI.set_xlabel('Frequency [Hz]')
axspecSI.set_ylabel('2|H$_+$($f$)|$\sqrt{f}$ & $\sqrt{S(f)}$')

pl.tight_layout()

pl.show()

sys.exit()

#
# Strain Time Series
#
fstrainNR, axstrainNR = pl.subplots(nrows=2, sharex=True)
NRtime = np.arange(0, catalogue.NR_datalen, catalogue.NR_deltaT)
NRtime -= NRtime[np.argmax(catalogue.NRAmpTimeSeries[0])]

axstrainNR[0].plot(NRtime, np.real(catalogue.NRComplexTimeSeries.T), color='grey')
axstrainNR[1].plot(NRtime, -1*np.imag(catalogue.NRComplexTimeSeries.T),
        linestyle='--', color='grey')

axstrainNR[0].set_xlim(-200, 50)
axstrainNR[1].set_xlim(-200, 50)

axstrainNR[0].minorticks_on()
axstrainNR[1].minorticks_on()

axstrainNR[0].set_ylabel(r'h$_+$(t) $\times$ D [M$_{\odot}$]')
axstrainNR[1].set_ylabel(r'h$_{\times}$(t) $\times$ D [M$_{\odot}$]')
axstrainNR[1].set_xlabel(r'Time [M$_{\odot}$]')

pl.tight_layout()
pl.subplots_adjust(hspace=0)

# ---- SI:
fstrainSI, axstrainSI = pl.subplots(nrows=2, sharex=True)

SItime = np.arange(0, SI_datalen, SI_deltaT)
SItime -= SItime[np.argmax(catalogue.SIAmpTimeSeries[0])]

axstrainSI[0].plot(SItime, np.real(catalogue.SIComplexTimeSeries.T),
        color='grey')
axstrainSI[1].plot(SItime, -1*np.imag(catalogue.SIComplexTimeSeries.T),
        linestyle='--', color='grey')

pl.suptitle(r'M$_{\rm tot}$=%d M$_{\odot}$'%total_mass)

axstrainSI[0].set_xlim(-0.5, 0.1)
axstrainSI[1].set_xlim(-0.5, 0.1)

axstrainSI[0].minorticks_on()
axstrainSI[1].minorticks_on()

axstrainSI[0].set_ylabel(r'h$_+$(t) $\times$ D [Mpc]')
axstrainSI[1].set_ylabel(r'h$_{\times}$(t) $\times$ D [Mpc]')
axstrainSI[1].set_xlabel(r'Time [s]')

pl.tight_layout()
pl.subplots_adjust(hspace=0)


#
# Amplitude / Phase Time Series
#
fampNR, axampNR = pl.subplots(nrows=2, sharex=True)

axampNR[0].plot(NRtime, catalogue.NRAmpTimeSeries.T, color='grey')
axampNR[1].plot(NRtime, catalogue.NRPhaseTimeSeries.T, linestyle='-',
        color='grey')

axampNR[0].set_xlim(-200, 50)
axampNR[1].set_xlim(-200, 50)

axampNR[0].minorticks_on()
axampNR[1].minorticks_on()

axampNR[0].set_ylabel(r'A(t) $\times$ D [M$_{\odot}$]')
axampNR[1].set_ylabel(r'$\Phi$(t)')
axampNR[1].set_xlabel(r'Time [M$_{\odot}$]')

pl.tight_layout()
pl.subplots_adjust(hspace=0)

# ---- SI:
fampSI, axampSI = pl.subplots(nrows=2, sharex=True)

axampSI[0].plot(SItime, catalogue.SIAmpTimeSeries.T, color='grey')
axampSI[1].plot(SItime, catalogue.SIPhaseTimeSeries.T, linestyle='-',
        color='grey')

axampSI[0].set_xlim(-0.5, 0.1)
axampSI[1].set_xlim(-0.5, 0.1)

axampSI[0].minorticks_on()
axampSI[1].minorticks_on()

axampSI[0].set_ylabel(r'A(t) $\times$ D [Mpc]')
axampSI[1].set_ylabel(r'$\Phi$(t)')
axampSI[1].set_xlabel(r'Time [s]')

pl.tight_layout()
pl.subplots_adjust(hspace=0)



