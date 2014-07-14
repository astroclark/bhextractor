#!/bin/bash

if [ $# -ne 4 ] 
then
    echo "insufficient args"
    exit
fi

# User inputs
injection=${1}
waveformN=${2}
seed=${3}
targetSNR=${4}

# Constants
Ndet=1 # DON"T CHANGE THIS
doplots=0
dodistance=0
dotimeshift=0
ra=1.57
dec=0.81
psi=0.89
trigtime=1087588740
realnoise=0
dosky=0
Fmin=10
Fmax=200

# Make run name
run_name=Inj-${injection}_Wf-${waveformN}_Seed-${seed}

# Make directory if it doesn't already exist
echo ${PWD}/${run_name}
mkdir -p ${PWD}/${run_name}

# Your matlab path (matlab willl now know about all the codes in this directory)
export MATLABPATH=${MATLABPATH}:"${BHEX_PREFIX}/SMEE"


# ---- Loop through models
for model in Q HR RO3
do

    if [ ${model} == "Q" ]
    then
        NumPCA=2
    elif [ ${model} == "HR" ]
    then
        NumPCA=4
    elif [ ${model} == "RO3" ]
    then
        NumPCA=5
    fi  

    echo matlab -nosplash -nojvm -nodisplay -r "bhextractor_nestsamp('${run_name}','aligo','${model}','${injection}',${waveformN},${ra},${dec},${psi},${Fmin},${Fmax},${trigtime},${Ndet},${NumPCA},${seed},'SNR',${targetSNR},${realnoise},${doplot},${dosky},${dodistance},${dotimeshift})"

    matlab -nosplash -nojvm -nodisplay -r "bhextractor_nestsamp('${run_name}','aligo','${model}','${injection}',${waveformN},${ra},${dec},${psi},${Fmin},${Fmax},${trigtime},${Ndet},${NumPCA},${seed},'SNR',${targetSNR},${realnoise},${doplots},${dosky},${dodistance},${dotimeshift})"


done

