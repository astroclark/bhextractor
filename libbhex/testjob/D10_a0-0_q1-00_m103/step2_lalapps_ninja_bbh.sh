#!/bin/sh
# James Clark, <james.clark@ligo.org>

waveform_name="D10_a0-0_q1-00_m103_Q"
outfile="${waveform_name}.xml"

lalapps_ninja \
    --format NINJA2 \
    --datadir ./ --outfile=${outfile} \
    --min-mass-ratio 0 --max-mass-ratio 10 \
    --pattern *gwf
