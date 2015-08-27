#!/bin/bash

ampPCs=Q-PCA_AmpPCs.dat
phasePCs=Q-PCA_PhasePCs.dat
numrelfile=D10_a0-0_q1-00_m103_Qs.gwf

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
    -p condor-log \
    -I ${new_injection_file} \
    -N ${numrelfile}

# Tarball up execution package
tar -cjf lalinference_execute.tar.bz2  \
    ${ampPCs} ${phasePCs} ${numrelfile} ${new_injection_file} \
    lalinference/* lalinference_nest.sh

# Tarball up the entire submission (since we'll probably need to build the
# analysis locally and scp)
cp lalinference_nest.sh OSGtest_Q_q1-00-${mass}-${snr}
cp lalinference/lalapps_nest2pos OSGtest_Q_q1-00-${mass}-${snr}
cp lalinference/lalapps_nest2pos.py OSGtest_Q_q1-00-${mass}-${snr}
cp lalinference/nest2pos.py OSGtest_Q_q1-00-${mass}-${snr}
mv lalinference_execute.tar.bz2 OSGtest_Q_q1-00-${mass}-${snr}
mv condor-log OSGtest_Q_q1-00-${mass}-${snr}
tar -cjf OSGtest_Q_q1-00-${mass}-${snr}.tar.bz2 \
    OSGtest_Q_q1-00-${mass}-${snr} 

# Clean up
rm ${new_injection_file}
rm -r OSGtest_Q_q1-00-${mass}-${snr}



