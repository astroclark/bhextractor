#!/bin/bash
#
# bhextractor_daggen.sh
#
# Copyright (C) 2014-2015 James Clark <clark@physics.umass.edu>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

if [ $# -ne 3 ] 
then
    echo "insufficient args"
    exit
fi

# ---------------------------------
# INPUTS

# The name of the Injection catalogue
signal=${1}

# The waveform number
waveformN=${2}

snr=${3}

# The name of the condor submission file for this run
subfile=inj-${signal}_wf-${waveformN}_snr-${snr}.sub
dagfile=inj-${signal}_wf-${waveformN}_snr-${snr}.dag
shellfile=inj-${signal}_wf-${waveformN}_snr-${snr}.sh

# --- Sanity check for existing files
if [ -f ${subfile} ] 
then
	echo "error, remove existing subfile: ${subfile}"
	exit
fi
if [ -f ${dagfile} ] 
then
	echo "error, remove existing dag file: ${dagfile}"
	exit
fi
if [ -f ${shellfile} ] 
then
	echo "error, remove existing dag file: ${shellfile}"
	exit
fi

echo '
#!/bin/bash 
####################
# BHEXTRACTOR
####################
# Run this interactively

' >> ${shellfile}


# The header for the sub file
echo "writing condor submission file: ${subfile}"
subtext="\
#########################
# BHEXTRACTOR: SUB file #
#########################

executable = `which bhextractor_nestsamp.sh`
universe   = vanilla 
arguments  = \$(macroarguments)
output     = condor_logs/bhextractor-\$(macroid)-\$(cluster)-\$(process).out
error      = condor_logs/bhextractor-\$(macroid)-\$(cluster)-\$(process).err
log        = condor_logs/bhextractor.log
getenv     = True

queue
"
echo "${subtext}" > ${subfile}

# The header for the dag file
echo "writing condor DAG file: ${dagfile}"
dagtext="\
####################
# BHEXTRACTOR: DAG #
####################
"
echo "${dagtext}"> ${dagfile}


# --- identify executables
# get full path for bhextractor_nestsamp.sh
bhexnest_exec=`which bhextractor_nestsamp.sh`
if [ -z "${bhexnest_exec}" ]
then
    echo "no bhextractor_nestsamp.sh in path, check environment"
    exit
fi
#sed -i 's|BHEXNEST_EXEC|'"${bhexnest_exec}"'|g' ${subfile}

if [ ! -d condor_logs ]
then
    echo "making condor log directory"
    mkdir condor_logs
fi

shaexec=`which shasum`
if [ -z "${shaexec}" ]
    # this works locally
then
    # this works on CIT
    shaexec=`which sha256sum`
fi

for j in `seq 1 50`
do

    # seed for this job
    seed=$((${RANDOM}+${RANDOM}))

    # arguments for this job
    jobargs="${signal} ${waveformN} ${seed} ${snr}"

	# Set up DAG lines
	jobname="bhex_${seed}"
    echo "JOB ${jobname} ${subfile}" >> ${dagfile}
    echo "VARS ${jobname} macroid=\"${jobname}\" macroarguments=\"${jobargs}\"" >> ${dagfile}
    echo "RETRY ${jobname} 3 " >> ${dagfile}
    echo "" >> ${dagfile}

    echo "# JOB ${jobname}" >> ${shellfile}
    echo """
    sh ${bhexnest_exec} ${signal} ${waveformN} ${seed} ${snr}
    """ >> ${shellfile}

done

echo "created ${dagfile}, ${subfile}, ${shellfile}"
echo "to run in condor: condor_submit_dag ${dagfile}"
echo "to run jobs sequentially: sh ${shellfile}"
echo "(Shell script is also useful for identifying and running jobs by hand)"

        
