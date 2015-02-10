#!/bin/sh 
# inspinj.sh
# Makes a call to lalapps_inspinj with some standard naming conventions.
# James Clark, <james.clark@ligo.org>

# FRAMEPATH should contain the ninja frame and the ninja xml file created by
# fr_ninja.sh and xninja.sh, respectively.

gpsstart=1097965864
gpsend=1097966376

lalapps_inspinj \
    --i-distr fixed  --seed 101 \
    --fixed-inc 0.0 \
    --polarization 0.0 \
    --waveform NumRelNinja2 \
    --gps-start-time ${gpsstart} --gps-end-time ${gpsend} --time-step 128 \
    --time-interval 10 --l-distr random --d-distr volume \
    --min-distance 1000000 --max-distance 1000000 \
    --min-mtotal 250 --max-mtotal 250 \
    --m-distr nrwaves --f-lower 10 \
    --real8-ninja2 \
    --nr-file "D10_a0-0_q1-50_m103_Qs.xml"  \
    --output 1GpcTest.xml
