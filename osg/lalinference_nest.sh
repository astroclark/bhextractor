#!/bin/bash

input_args=${@}

execute_dir="lalinference_execute"

if test ! -d ${execute_dir}; then tar xjf ${execute_dir}.tar.bz2; fi

ls -la

echo ${input_args}

if test ! -d log; then mkdir log; fi
if test ! -d posterior_samples; then mkdir posterior_samples; fi
if test ! -d engine; then mkdir engine; fi

lalinference/lalinference_nest ${input_args}

ls -la
