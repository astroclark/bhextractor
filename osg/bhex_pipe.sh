#!/bin/bash

mass=${1}
snr=${2}

original_injection_file="Q_M-${mass}_snr-${snr}.xml"
new_injection_file=`echo ${original_injection_file} | sed 's/.xml/_relpath.xml/g'`

if test -f ${new_injection_file}; rm ${new_injection_file}; fi

# Reset frame paths in sim_inspiral table to relative paths
for line in `ligolw_print Q_M-300_snr-15.xml -t sim_inspiral -c numrel_data`
do
    replacement=`basename ${line}`
    sed -e 's|'"${line}"'|'"${replacement}"'|g' ${injection_file} >> ${new_injection_file}
done

config_file="libbhex.ini"

./lalinference_pipe ${config_file}\
    -r OSGtest_Q_q1-00-${mass}-${snr}\
    -p OSGtest_Q_q1-00-${mass}-${snr}/logs \
    -I ${injection_file}


