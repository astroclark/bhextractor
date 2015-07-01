#!/bin/sh 
# inspinj.sh
# Makes a call to lalapps_inspinj with some standard naming conventions.
# James Clark, <james.clark@ligo.org>

# FRAMEPATH should contain the ninja frame and the ninja xml file created by
# fr_ninja.sh and xninja.sh, respectively.

gpsstart=1111178318
gpsend=1111264718

waveform="D10_a0-0_q1-00_m103_Q"
outname=Q_M-300_snr-50

lalapps_inspinj \
    --i-distr uniform --seed `lalapps_tconvert now` \
    --waveform NumRelNinja2 \
    --gps-start-time ${gpsstart} --gps-end-time ${gpsend} --time-step 128 \
    --time-interval 10 --l-distr random \
    --min-mtotal 300 --max-mtotal 300 \
    --m-distr nrwaves --f-lower 10 \
    --real8-ninja2 \
    --nr-file "${waveform}.xml"  \
    --snr-distr volume \
    --min-snr  50 --max-snr 50 \
    --ligo-psd aligopsd.txt \
    --ligo-start-freq 10 \
    --virgo-psd advirgopsd.txt \
    --virgo-start-freq 10 \
    --ifos H1,L1,V1 \
    --ninja-snr \
    --verbose \
    --output ${outname}.xml
