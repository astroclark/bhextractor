#!/bin/sh 
# inspinj.sh
# Makes a call to lalapps_inspinj with some standard naming conventions.
# James Clark, <james.clark@ligo.org>

# FRAMEPATH should contain the ninja frame and the ninja xml file created by
# fr_ninja.sh and xninja.sh, respectively.

gpsstart=1111178318
gpsend=1111178718

mass=${1}
snr=${2}

waveform="NonSpinning"
outname=${waveform}_M-${mass}_snr-${snr}

lalapps_inspinj \
    --i-distr uniform --seed `lalapps_tconvert now` \
    --waveform NumRelNinja2 \
    --gps-start-time ${gpsstart} --gps-end-time ${gpsend} --time-step 128 \
    --time-interval 10 --l-distr random \
    --min-mtotal ${mass} --max-mtotal ${mass} \
    --m-distr nrwaves --f-lower 10 \
    --real8-ninja2 \
    --nr-file "${waveform}.xml"  \
    --snr-distr volume \
    --min-snr  ${snr} --max-snr ${snr} \
    --ligo-psd aligopsd.txt \
    --ligo-start-freq 30 \
    --virgo-psd advirgopsd.txt \
    --virgo-start-freq 30 \
    --ifos H1,L1,V1 \
    --ninja-snr \
    --verbose \
    --output ${outname}.xml
