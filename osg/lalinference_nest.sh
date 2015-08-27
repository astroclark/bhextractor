#!/bin/bash

input_args=${@}

tar xjf lalinference_execute.tar.bz2

ls -la

echo ${input_args}

if test ! -d log; then mkdir log; fi

lalinference/lalinference_nest ${input_args}

ls -la
