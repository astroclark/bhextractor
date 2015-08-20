#!/bin/bash

ampPCs=Q-PCA_AmpPCs.dat
phasePCs=Q-PCA_PhasePCs.dat

mass=${1}
snr=${2}

original_injection_file="Q_M-${mass}_snr-${snr}.xml"
new_injection_file=`echo ${original_injection_file} | sed 's/.xml/_relpath.xml/g'`

if test -f ${new_injection_file}; then rm ${new_injection_file}; fi
cp ${original_injection_file} ${new_injection_file}

# Reset frame paths in sim_inspiral table to relative paths
for line in `ligolw_print ${original_injection_file} -t sim_inspiral -c numrel_data`
do
    replacement=`basename ${line}`
    sed -i 's|'"${line}"'|'"${replacement}"'|g' ${new_injection_file}
done

config_file="libbhex.ini"

bhex_pipe/lalinference_pipe ${config_file}\
    -r OSGtest_Q_q1-00-${mass}-${snr}\
    -p OSGtest_Q_q1-00-${mass}-${snr}/logs \
    -I ${new_injection_file}

# Copy executables to run directory
cp lalinference/* OSGtest_Q_q1-00-${mass}-${snr}
cp ${ampPCs} OSGtest_Q_q1-00-${mass}-${snr}
cp ${phasePCs} OSGtest_Q_q1-00-${mass}-${snr}

