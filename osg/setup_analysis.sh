#!/usr/bin/bash

mass=${1}
snr=${2}

injection_file="Q_M-${mass}_snr-${snr}.xml"

config_file="Qs_q1-00_snr-50.ini"

./lalinference_pipe ${config_file}
