#!/usr/bin/bash

mass=${1}
snr=${2}

inj_path="/home/jclark308/Projects/bhextractor/ninja_injections"
injection_file="NonSpinning_M-150_snr-15.xml"

config_file="libbhex_config_nonspinning.ini"

outdir=`echo ${injection_file} | sed 's/\.xml//g'`


lalinference_pipe \
	-I ${inj_path}/${injection_file} \
	-r ${outdir} \
	-p ${outdir}/logs \
	${config_file}
	
