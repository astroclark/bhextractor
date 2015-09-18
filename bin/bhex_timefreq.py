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
from bhex_utils import bhex_pca as bpca
from pycbc.psd import aLIGOZeroDetHighPower
import pycbc.types
import numpy as np
from matplotlib import pyplot as pl
from pmns_utils import pmns_pca as ppca


pl.rcParams.update({'axes.labelsize': 18})
pl.rcParams.update({'xtick.labelsize':18})
pl.rcParams.update({'ytick.labelsize':18})
pl.rcParams.update({'legend.fontsize':18})


def tripleplot(cwt_result):

    signal = cwt_result['analysed_data']

    Z = np.copy(cwt_result['map'])

    for c in xrange(np.shape(Z)[1]):
        Z[:,c] /= np.log10(np.sqrt(cwt_result['frequencies']) / 2)

    # Open the figure
    fig, ax_cont = pl.subplots()

    #maxcol = 0.45*Z.max()
    maxcol = 0.4*Z.max()
    vmin, vmax = 0, maxcol
    collevs=np.linspace(vmin, vmax, 100)

    # Characteristic times
    tmax = signal.sample_times[np.argmax(signal.data)]

    # --- CWT
    c = ax_cont.contourf(cwt_result['times']-tmax, cwt_result['frequencies'], Z,
            levels=collevs, extend='both')

    #c.cmap.set_over('k')
    #c.set_cmap('gnuplot2_r')

    #ax_cont.set_yscale('log')
    ax_cont.set_xlabel('Time [s]')
    ax_cont.set_ylabel('Frequency [Hz]')

    #ax_cont.set_xlim(-2e-3, 1.5e-2)
    #ax_cont.set_xlim(-2e-3, 1.0e-2)

    ax_cont.minorticks_on()

    #ax_cont.invert_yaxis()
 
    fig.tight_layout()

    return fig, ax_cont


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Some useful info:
#

#valid_series = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
#        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series"
#        "TP-series"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT

SI_deltaT = 1./512
SI_datalen= 4.0

total_mass = 150.  # Generate SI waveforms at this mass
distance=1e3 # Mpc

series_names = ['RO3-series'] # (see above for valid choices)

#
# Modify for imposing parameter bounds on the catalogue:
#
#bounds=None
bounds=dict()
bounds['a1'] = [0.6, 0.6]
bounds['a2'] = [0.6, 0.6]
bounds['thJL'] =  [16.9, 16.9]
bounds['q'] = [3, np.inf] 

save_pcs = ['NRhplusTimeSeriesPCA', 'NRhcrossTimeSeriesPCA',
'NRAmpTimeSeriesPCA', 'NRPhaseTimeSeriesPCA',  'SIhplusTimeSeriesPCA',
'SIhcrossTimeSeriesPCA', 'SIAmpTimeSeriesPCA', 'SIPhaseTimeSeriesPCA']

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do the work:

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
simulations = \
        bwave.simulation_details(series_names=series_names,
                param_bounds=bounds, Mmin30Hz=total_mass)

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Building NR catalogue'
print ''
catalogue = bwave.waveform_catalogue(simulations, ref_mass=total_mass,
        SI_deltaT=SI_deltaT, SI_datalen=SI_datalen, distance=distance,
        trunc_time=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting

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

axstrainSI[0].set_ylim(-1.5e-21, 1.5e-21)
axstrainSI[1].set_ylim(-1.5e-21, 1.5e-21)

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

#
# Magnitude Spectra
#
fspecSI, axspecSI = pl.subplots()

for wave in catalogue.SIComplexTimeSeries:

    hplus_NR = \
            pycbc.types.TimeSeries(np.real(wave), delta_t=SI_deltaT)
    hplus_NR.data = bwave.window_wave(hplus_NR.data)
    Hplus_NR = hplus_NR.to_frequencyseries()

    axspecSI.loglog(Hplus_NR.sample_frequencies,
            2*abs(Hplus_NR)*np.sqrt(Hplus_NR.sample_frequencies), color='grey')

    delta_f = Hplus_NR.delta_f
    flen = len(Hplus_NR)
    psd = aLIGOZeroDetHighPower(flen, delta_f, 0.1) 
    axspecSI.loglog(psd.sample_frequencies, np.sqrt(psd), label='aLIGO',
            color='k', linestyle='--')


axspecSI.set_xlim(5, 512)
#axspecSI.axvline(30, color='r')
axspecSI.minorticks_on()
axspecSI.set_xlabel('Frequency [Hz]')
axspecSI.set_ylabel('2|H$_+$($f$)|$\sqrt{f}$ & $\sqrt{S(f)}$')

pl.tight_layout()

pl.show()

# XXX
ts = pycbc.types.TimeSeries(np.real(catalogue.SIComplexTimeSeries[0,:]),
        delta_t=1./512)
cwt_result = ppca.build_cwt(ts, mother_freq=3, fmin=4)
fig, ax_cont = tripleplot(cwt_result)
ax_cont.set_ylim(16,256)
ax_cont.set_xlim(-0.5, 0.25)

pl.show()
