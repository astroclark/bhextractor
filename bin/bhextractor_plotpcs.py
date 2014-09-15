#!/usr/bin/env python

import sys
import numpy as np
import scipy.io as sio
from scipy import signal
import HHT
from matplotlib import pyplot as pl
import matplotlib

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
        'legend.fontsize': 6,
        'xtick.labelsize': 4,
        'ytick.labelsize': 4,
        'text.usetex': True,
        'figure.figsize': fig_size,
        'font.family': "serif",
        'font.serif': ["Times"]
        })  

matplotlib.rcParams.update(
        {'savefig1.dpi': 200,
        'xtick.major.size':4,
        'xtick.minor.size':4,
        'ytick.major.size':4,
        'ytick.minor.size':4
        })
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


import lal

inputfile=sys.argv[1]
outname=inputfile.replace('.mat','').split('/')[-1]
data=sio.loadmat(inputfile)

fs=16384.

#ZDHP=np.loadtxt('/Users/jclark/Projects/BHEX/SMEE_repo/SRDs/ZERO_DET_high_P.txt')
#noise_freqs=ZDHP[:,0]
#Sf=ZDHP[:,1]**2

Mtot=250.
Dist=1.
Mscale_D = Mtot * lal.MRSUN_SI / (Dist * 1e9 * lal.PC_SI)
Mscale_T = Mtot * lal.MTSUN_SI #/ (Dist * 1e9 * lal.LAL_PC_SI)

# Get catalogue label from filename XXX: careful with this (doesn't care about
# other identifiers)
catname=outname.split('/')[-1].split('_')[0]

#time_axis=np.arange(0,1,1.0/fs)
time_axis=np.arange(0, len(data['PCs_final'][:,0])/16384.0, 1.0/16384)

npcs={'Q':3, 'HR':6, 'RO3':7}

# --- Time series
#fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],figsize=(8,12),sharex='col')
fig,ax=pl.subplots(npcs[catname],sharex='col')
for r,row in enumerate(ax):
    #row.plot(data['PCs_final'][:,r]/max(data['PCs_final'][:,r]))
    row.plot(time_axis/Mscale_T,np.imag(data['PCs_final'][:,r]), color='k', linewidth=0.5)
    #row.set_ylim(-1.1,1.1)
    #row.set_yticklabels('')

    ticks=row.get_yticks()
    ticks=np.arange(-0.025,0.03,0.005)
    row.set_yticks(ticks)
    row.set_ylim(-0.03,0.03)

    #row.minorticks_on()
    #row.grid(which='major',color='grey',linestyle='-')
    ##row.set_xlim(0.25/Mscale_T,0.85/Mscale_T)
    #if r==0: 
    #    row.set_title(
    #            '%s-catalogue dominant principle components'%catname,
    #            weight='bold')

#pl.suptitle('%s [arb y-units]'%outname)
pl.xlabel('Time [M$_{\odot}$]',weight='bold')
pl.subplots_adjust(hspace=0.0,wspace=0.35,bottom=0.08,top=0.98)
#pl.savefig('%s_princcomps_TD.png'%outname)
pl.savefig('%s_princcomps_TD.pdf'%outname)

sys.exit()

# --- HHT
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],ncols=2,figsize=(8,20),sharex='col')
for r,row in enumerate(ax):

    # Take HHT
    time=np.arange(0,len(data['PCs_final'][:,r])/fs,1/fs)
    h=HHT.HHT(time,data['PCs_final'][:,r],N=0)

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

pl.savefig('%s_princcomps_HHT.png'%outname)

# --- PSDs
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],figsize=(8,12),sharex='col')
for r,row in enumerate(ax):
    freq, Pxx_den = signal.periodogram(data['PCs_final'][:,r], fs)
    row.plot(freq,Pxx_den/max(Pxx_den))
    row.set_ylim(0,1.1)
    row.set_xlim(5,50)
    row.axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
pl.xlabel('Frequency [Hz]')
pl.suptitle('%s [arb y-units]'%outname)
pl.subplots_adjust(hspace=0.2)

pl.savefig('%s_princcomps_PSD.png'%outname)

# --- Whitened PSDs
# Generate PSD/ASD
psd=np.zeros(len(freq))
for i in range(len(freq)):
    psd[i]=lalsim.SimNoisePSDaLIGOZeroDetHighPower(freq[i])
zdhp_asd=np.sqrt(psd)

fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],figsize=(8,12),sharex='col')
for r,row in enumerate(ax):
    freq, Pxx_den = signal.periodogram(data['PCs_final'][:,r], fs)
    row.plot(freq,Pxx_den/zdhp_asd)
    row.set_yticklabels('')
    row.set_xlim(5,100)
    row.axvline(10,color='k',label=r'$f_{\mathrm{low}}$')
pl.xlabel('Frequency [Hz]')
pl.suptitle('Whitened %s [arb y-units]'%outname)
pl.subplots_adjust(hspace=0.2)

pl.savefig('%s_princcomps_whitenedPSD.png'%outname)

# --- HHT

# --- Complex Spectra
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],ncols=2,figsize=(8,12),sharex='col')
for r,row in enumerate(ax):

    # Take FFT
    Fspec=np.fft.fft(data['PCs_final'][:,r])
    freq=np.fft.fftfreq(len(data['PCs_final'][:,r]),1./fs)

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

pl.savefig('%s_princcomps_FD.png'%outname)

# --- Whitened Complex Spectra
fig,ax=pl.subplots(np.shape(data['PCs_final'])[1],ncols=2,figsize=(8,12),sharex='col')
for r,row in enumerate(ax):

    # Take FFT
    Fspec=np.fft.fft(data['PCs_final'][:,r])
    freq=np.fft.fftfreq(len(data['PCs_final'][:,r]),1./fs)

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

pl.savefig('%s_princcomps_whitenedFD.png'%outname)

