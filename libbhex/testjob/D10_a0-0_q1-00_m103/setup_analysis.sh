#!/usr/bin/bash

injection_file="Q_M-300_snr-50.xml"

config_file="Qs_q1-00_snr-50.ini"

lalinference_pipe \
	-I ${injection_file} \
	-r Q_M-300_snr-50 \
	-p Q_M-300_snr-50/logs \
	${config_file}
	
