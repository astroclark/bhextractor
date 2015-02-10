#!/bin/sh
# James Clark <james.clark@ligo.org>

# location of NINJA2-format ascii files
waveform_name="D10_a0-0_q1-50_m103_Qs"
nrpath="./"

echo "lalapps_fr_ninja --verbose --format NINJA2 --double-precision --nr-data-dir ${nrpath} --nr-meta-file ${nrpath}/${waveform_name}.bbh --output "${waveform_name}.gwf""

lalapps_fr_ninja \
    --verbose --format NINJA2 \
    --double-precision \
    --nr-data-dir ${nrpath} \
    --nr-meta-file ${nrpath}/${waveform_name}.bbh \
    --output "${waveform_name}.gwf"
