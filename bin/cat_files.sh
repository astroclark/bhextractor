#!/bin/bash

for sigfile in signal*
do
    awk '{print $2}' ${sigfile} > tmp
    tr '\n' ' ' < tmp >>  waveform_geo_1000.dat
    echo '\n' >>  waveform_geo_1000.dat
done
rm tmp
