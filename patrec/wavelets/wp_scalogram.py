#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.cm as cm
import pylab

import pywt

Fs=2048
x = pylab.arange(0, 1, 1. / Fs)
#data = pylab.sin((5 * 50 * pylab.pi * x ** 2))
data = pylab.sin((100 * 2*pylab.pi * x**2))

wavelet = 'db4'
level = 4
order = "freq"  # "normal"
interpolation = 'nearest'
cmap = cm.cool

# perform wavelet decomposition at levels 1->maxlevel
wp = pywt.WaveletPacket(data, wavelet, 'sym', maxlevel=level)

# Get the nodes of the decomposition at the level you specified
nodes = wp.get_level(level, order=order)
#labels = [n.path for n in nodes]
scales=len(nodes)
freqs = pywt.scal2frq(wavelet, (1+pylab.arange(0,scales))[::-1], 1./Fs)
labels = ['%.2f'%f for f in freqs]
values = pylab.array([n.data for n in nodes], 'd')
values = abs(values)**2

# plotting:

f = pylab.figure()
f.subplots_adjust(hspace=0.2, bottom=.03, left=.07, right=.97, top=.92)
pylab.subplot(2, 1, 1)
pylab.title("linchirp signal")
pylab.plot(x, data, 'b')
pylab.xlim(0, x[-1])

ax = pylab.subplot(2, 1, 2)
pylab.title("Wavelet packet coefficients at level %d" % level)
pylab.imshow(values, interpolation=interpolation, cmap=cmap, aspect="auto",
    origin="lower", extent=[0, 1, 0, max(freqs)])
    #origin="lower", extent=[0, 1, 0, len(values)])
#pylab.yticks(pylab.arange(0.5, len(labels) + 0.5), labels)
#pylab.yticks(pylab.arange(0.5, len(labels) + 0.5), labels)
#pylab.setp(ax.get_xticklabels(), visible=False)

#   pylab.figure(2)
#   pylab.specgram(data, NFFT=64, noverlap=32, cmap=cmap)
#   pylab.imshow(values, origin='upper', extent=[-1,1,-1,1],
#    interpolation='nearest')
#
pylab.show()
