#!/bin/bash
 
${BHEX_PREFIX}/bin/bhextractor_match.py \
    ${BHEX_PREFIX}/data/PCA_data/Q_PCs_theta-0.mat,${BHEX_PREFIX}/data/PCA_data/HR_PCs_theta-0.mat,${BHEX_PREFIX}/data/PCA_data/RO3_PCs_theta-0.mat\
    ${BHEX_PREFIX}/data/signal_data/Q_catalogue_theta-0.mat,${BHEX_PREFIX}/data/signal_data/HR_catalogue_theta-0.mat,${BHEX_PREFIX}/data/signal_data/RO3_catalogue_theta-0.mat
