#!/bin/sh
# James Clark, <james.clark@ligo.org>

datadir=${PWD}
waveform_name="D10_a0-0_q1-50_m103_Qs"
outfile="${waveform_name}.xml"

echo "lalapps_ninja \ --format NINJA2 --datadir ${datadir} --outfile=${outfile} --min-mass-ratio 0 --max-mass-ratio 10 --pattern *gwf"

lalapps_ninja \
    --format NINJA2 \
    --datadir ${datadir} --outfile=${outfile} \
    --min-mass-ratio 0 --max-mass-ratio 10 \
    --pattern *gwf
