
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
shellfile=inj-${signal}_wf-${waveformN}-${snr}.sh

echo '
#!/bin/bash 
####################
# BHEXTRACTOR
####################
# Run this interactively

' >> ${shellfile}


# The header for the sub file
echo '
####################
# BHEXTRACTOR
####################

executable = BHEXNEST_EXEC
universe   = vanilla 

output     = condor_logs/bhextractor_$(Process)_RANDOM.out
error      = condor_logs/bhextractor_$(Process)_RANDOM.err
log        = condor_logs/bhextractor_$(Process)_RANDOM.log

getenv=True
' >> ${subfile}


# get full path for bhextractor_nestsamp.sh
bhexnest_exec=`which bhextractor_nestsamp.sh`
if [ -z "${bhexnest_exec}" ]
then
    echo "no bhextractor_nestsamp.sh in path, check environment"
    exit
fi
sed -i 's|BHEXNEST_EXEC|'"${bhexnest_exec}"'|g' ${subfile}
sed -i 's|RANDOM|'"${RANDOM}"'|g' ${subfile}

if [ ! -d condor_logs ]
then
    echo "making condor log directory"
    mkdir condor_logs
fi

for j in `seq 1 50`
do

    seed=${RANDOM}

    echo """
    arguments = ${signal} ${waveformN} ${seed} ${snr}
    queue
    """ >> ${subfile}

    echo """
    sh SMEE_BBH.sh ${signal} ${waveformN} ${seed} ${snr}
    """ >> ${shellfile}

done

        
