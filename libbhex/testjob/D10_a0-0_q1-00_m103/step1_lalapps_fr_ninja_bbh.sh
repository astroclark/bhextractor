#!/bin/sh
# James Clark <james.clark@ligo.org>

# location of NINJA2-format ascii files
waveform_name=D10_a0-0_q1-00_m103_Qs
nrpath="${BHEX_PREFIX}/data/NR_data/Q-series/D10_a0.0_q1.00_m103_Qs"

cp ${nrpath}/*asc .

lalapps_fr_ninja \
    --verbose --format NINJA2 \
    --double-precision \
    --nr-data-dir ./ \
    --nr-meta-file ${waveform_name}.bbh \
    --output "${waveform_name}.gwf"

rm *asc
