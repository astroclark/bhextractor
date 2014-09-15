#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2014-2015 James Clark <james.clark@ligo.org>
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
"""

import sys
import numpy as np
import scipy.io as sio
from scipy import signal
import HHT
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as pl

#font = {'family' : 'Helvetica',
#        'weight' : 'bold',
#        'size'   : 12}
#
#matplotlib.rc('font', **font)

fig_width_pt = 170  # Get this from LaTeX using \showthe\columnwidth
#fig_height_pt = 300 
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
#fig_height_pt = fig_width_pt*golden_mean
fig_height_pt = fig_width_pt*1.4
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_height_pt*inches_per_pt
fig_size =  [fig_width,fig_height]

matplotlib.rcParams.update(
        {'axes.labelsize': 6,
        'text.fontsize':   6,  
        'legend.fontsize': 4,
        'xtick.labelsize': 4,
        'ytick.labelsize': 4,
        'text.usetex': True,
        'figure.figsize': fig_size,
        'font.family': "serif",
        'font.serif': ["Times"]
        })  

matplotlib.rcParams.update(
        {'savefig1.dpi': 200,
        'xtick.major.size':2,
        'xtick.minor.size':2,
        'ytick.major.size':2,
        'ytick.minor.size':2
        })
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

import lal

inputfile=sys.argv[1]
outname=inputfile.replace('.mat','')
outname=outname.split('/')[-1]
data=sio.loadmat(inputfile)

Mtot=250.
Dist=1.
Mscale_D = Mtot * lal.MRSUN_SI / (Dist * 1e9 * lal.PC_SI)
Mscale_T = Mtot * lal.MTSUN_SI #/ (Dist * 1e9 * lal.PC_SI)
fs=16384.

#ZDHP=np.loadtxt('/Users/jclark/Projects/BHEX/SMEE_repo/SRDs/ZERO_DET_high_P.txt')
#noise_freqs=ZDHP[:,0]
#Sf=ZDHP[:,1]**2

#time_axis=np.arange(0,1,1.0/fs)

# --- Labels for the different waveforms
Q_labels=['q=1.0','q=1.15','q=1.3','q=1.45','q=1.5','q=1.6',
        'q=1.75', 'q=1.9', 'q=2.0', 'q=2.05', 'q=2.2', 'q=2.35',
        'q=2.5']

HR_labels=['q=1.0,a=0.0','q=1.0,a=0.1','q=1.0, a=0.2','q=1.0, a=0.3', \
        'q=1.0, a=0.4', 'q=1.0, a=0.5','q=1.0, a=0.6','q=1.0, a=0.7',\
        'q=1.0, a=0.8','q=1.5, a=0.0','q=1.5, a=0.2','q=1.5, a=0.4',\
        'q=2.0, a=0.0','q=2.5, a=0.0','q=4.0, a=0.0']

RO3_labels=['q=1.5, a=0.4, $\Theta$=60$^{\circ}$', \
        'q=1.5, a=0.6, $\Theta$=45$^{\circ}$',
        'q=1.5, a=0.6, $\Theta$=60$^{\circ}$',
        'q=1.5, a=0.6, $\Theta$=90$^{\circ}$',

        'q=2.0, a=0.4, $\Theta$=60$^{\circ}$',
        'q=2.0, a=0.6, $\Theta$=45$^{\circ}$',
        'q=2.0, a=0.6, $\Theta$=60$^{\circ}$',
        'q=2.0, a=0.6, $\Theta$=90$^{\circ}$',
        'q=2.0, a=0.6, $\Theta$=135$^{\circ}$',
        'q=2.0, a=0.6, $\Theta$=180$^{\circ}$',
        'q=2.0, a=0.6, $\Theta$=270$^{\circ}$',

        'q=2.5, a=0.4, $\Theta$=45$^{\circ}$',
        'q=2.5, a=0.4, $\Theta$=60$^{\circ}$',
        'q=2.5, a=0.4, $\Theta$=90$^{\circ}$',

        'q=2.5, a=0.6, $\Theta$=45$^{\circ}$',
        'q=2.5, a=0.6, $\Theta$=60$^{\circ}$',
        'q=2.5, a=0.6, $\Theta$=90$^{\circ}$',

        'q=4.0, a=0.6, $\Theta$=45$^{\circ}$',
        'q=4.0, a=0.6, $\Theta$=60$^{\circ}$',
        'q=4.0, a=0.6, $\Theta$=90$^{\circ}$',
        ]
labels={'Q':Q_labels, 'HR':HR_labels, 'RO3':RO3_labels}

# Get catalogue label from filename XXX: careful with this (doesn't care about
# other identifiers)
catname=outname.split('/')[-1].split('_')[0]

print np.shape(data['MDC_final'])
time_axis=np.linspace(0,np.shape(data['MDC_final'])[0]/fs,np.shape(data['MDC_final'])[0])
print np.shape(time_axis)


# --- Time series
#fig,ax=pl.subplots(np.shape(data['MDC_final'])[1],figsize=(10,15),sharex='col')
fig,ax=pl.subplots(np.shape(data['MDC_final'])[1],sharex='col')
for r,row in enumerate(ax):
    #row.plot(data['MDC_final'][:,r]/max(data['MDC_final'][:,r]),)
    #row.set_ylim(-1.1,1.1)
    #row.set_yticklabels('')
    row.set_yticks(np.arange(-3e-3,3e-3+0.001,0.002))
    row.set_ylim(-3e-3-0.001,3e-3+0.001)
    #row.set_yticklabels('')
    row.plot(time_axis/Mscale_T, np.imag(data['MDC_final'][:,r])/Mscale_D,
            color='k',linewidth=0.5, label=labels[catname][r])
    #row.minorticks_on()
    #row.grid(which='major',color='grey',linestyle='-')
    #row.set_xlim(0.25/Mscale_T,0.85/Mscale_T)
    if r==np.floor(len(ax)/2): 
        #row.set_title(
        #        '%s Waveform Catalogue'%catname,
        #        weight='bold')
        row.set_ylabel('rh$_{+}$ [M$_{\odot}$]')
    #leg=row.legend(loc='upper right',frameon=True)
    #for label in leg.get_lines():
    #        label.set_linewidth(0.01)  # the legend line width

pl.xlabel('Time [M$_{\odot}$]')
pl.subplots_adjust(hspace=0.0,wspace=0.35,left=0.2,bottom=0.08,top=0.98)
#fig.tight_layout()
#fig.savefig('%s_waveforms_TD.png'%outname)
fig.savefig('%s_waveforms_TD.pdf'%outname)
#pl.show()

sys.exit()

# --- HHT
fig,ax=pl.subplots(np.shape(data['MDC_final'])[1],ncols=2,figsize=(8,20),sharex='col')
for r,row in enumerate(ax):

    # Take HHT
    time=np.arange(0,len(data['MDC_final'][:,r])/fs,1/fs)
    h=HHT.HHT(time,data['MDC_final'][:,r],N=0)

    row[0].plot(time,h[0]/max(h[0]))
    row[0].set_ylim(0,1.1)

    row[1].plot(time,h[1])
    row[1].set_ylim(0,100)

    if r==0:
        row[0].set_title('Instaneous Amplitude [arb units]')
        row[1].set_title('Instaneous Frequency [Hz]')

row[0].set_xlabel('Time [s]')
row[1].set_xlabel('Time [s]')

pl.suptitle('%s Hilbert Huang Transform'%outname)

fig.subplots_adjust(hspace=0)

pl.savefig('%s_waveforms_HHT.png'%outname)


# --- PSDs
fig,ax=pl.subplots(np.shape(data['MDC_final'])[1],figsize=(8,12),sharex='col')
for r,row in enumerate(ax):
    freq, Pxx_den = signal.periodogram(data['MDC_final'][:,r], fs)
    row.plot(freq,Pxx_den/max(Pxx_den))
    row.set_ylim(0,1.1)
    row.set_xlim(5,50)
    row.axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
pl.xlabel('Frequency [Hz]')
pl.suptitle('%s [arb y-units]'%outname)
pl.subplots_adjust(hspace=0.2)

pl.savefig('%s_waveforms_PSD.png'%outname)

# --- Whitened PSDs

# Generate PSD/ASD
psd=np.zeros(len(freq))
for i in range(len(freq)):
    psd[i]=lalsim.SimNoisePSDaLIGOZeroDetHighPower(freq[i])
zdhp_asd=np.sqrt(psd)

fig,ax=pl.subplots(np.shape(data['MDC_final'])[1],figsize=(8,12),sharex='col')
for r,row in enumerate(ax):
    freq, Pxx_den = signal.periodogram(data['MDC_final'][:,r], fs)
    row.plot(freq,Pxx_den/zdhp_asd)
    row.set_yticklabels('')
    row.set_xlim(5,100)
    row.axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
pl.xlabel('Frequency [Hz]')
pl.suptitle('Whitened %s [arb y-units]'%outname)
pl.subplots_adjust(hspace=0.2)

pl.savefig('%s_waveforms_whitenedPSD.png'%outname)

# --- HHT

# --- Complex Spectra
fig,ax=pl.subplots(np.shape(data['MDC_final'])[1],ncols=2,figsize=(8,12),sharex='col')
for r,row in enumerate(ax):

    # Take FFT
    Fspec=np.fft.fft(data['MDC_final'][:,r])
    freq=np.fft.fftfreq(len(data['MDC_final'][:,r]),1./fs)

    # Take +ve freq
    Fspec=Fspec[freq>=0]
    freq=freq[freq>=0]

    row[0].plot(freq,Fspec/max(Fspec))
    row[0].set_ylim(-1.1,1.1)
    row[0].set_xlim(0,100)
    row[0].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
    row[1].plot(freq,Fspec/max(Fspec))
    row[1].set_ylim(-1.1,1.1)
    row[1].set_xlim(0,100)
    row[0].set_yticklabels('')
    row[1].set_yticklabels('')
    row[1].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')

    if r==0:
        row[0].set_title('Re[FFT]')
        row[1].set_title('Im[FFT]')

row[0].set_xlabel('Frequency [Hz]')
row[1].set_xlabel('Frequency [Hz]')
pl.suptitle('%s [arb y-units]'%outname)

fig.subplots_adjust(hspace=0)

pl.savefig('%s_waveforms_FD.png'%outname)

# --- Whitened Complex Spectra
fig,ax=pl.subplots(np.shape(data['MDC_final'])[1],ncols=2,figsize=(8,12),sharex='col')
for r,row in enumerate(ax):

    # Take FFT
    Fspec=np.fft.fft(data['MDC_final'][:,r])
    freq=np.fft.fftfreq(len(data['MDC_final'][:,r]),1./fs)

    # Take +ve freq
    Fspec=Fspec[freq>=0]
    freq=freq[freq>=0]

    # Interpolate the noise PSD to data freqs
    zdhp_asd=np.interp(freq,noise_freqs,Sf)

    # whitened complex spectrum
    FspecWhite=Fspec/np.sqrt(zdhp_asd)

    row[0].plot(freq,np.real(FspecWhite)/max(np.real(FspecWhite)))
    row[0].set_ylim(-1.1,1.1)
    row[0].set_xlim(0,100)
    row[0].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
    row[1].plot(freq,np.imag(FspecWhite)/max(np.imag(FspecWhite)))
    row[1].set_ylim(-1.1,1.1)
    row[1].set_xlim(0,100)
    row[0].set_yticklabels('')
    row[1].set_yticklabels('')
    row[1].axvline(10,color='k',label=r'$f_{\mathrm{low}}$')

    if r==0:
        row[0].set_title('Re[FFT]')
        row[1].set_title('Im[FFT]')

row[0].set_xlabel('Frequency [Hz]')
row[1].set_xlabel('Frequency [Hz]')
pl.suptitle('%s [arb y-units]'%outname)

fig.subplots_adjust(hspace=0)

pl.savefig('%s_waveforms_whitenedFD.png'%outname)

