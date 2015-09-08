#! /bin/sh
# Make env setup script
# James Clark, james.clark@ligo.org

# Get the location of the git repository by finding the full path of this file
# and stripping off the name
BHEX_PREFIX=`python -c "import os, sys; print os.path.realpath('${0}')" | sed 's|/setup.sh||g'`

# create an etc directory
test -d "${BHEX_PREFIX}/etc" || mkdir "${BHEX_PREFIX}/etc"


echo "# add script directory to path" > "${BHEX_PREFIX}/etc/bhex-user-env.sh"
echo "export PATH=$BHEX_PREFIX/bin:\$PATH" >> "${BHEX_PREFIX}/etc/bhex-user-env.sh"
echo "export PYTHONPATH=$BHEX_PREFIX/bin:$BHEX_PREFIX:\${PYTHONPATH}" >> "${BHEX_PREFIX}/etc/bhex-user-env.sh"
echo "# define variable for location of waveform data" >> "${BHEX_PREFIX}/etc/bhex-user-env.sh"
echo "export BHEX_PREFIX=${BHEX_PREFIX}" >> "${BHEX_PREFIX}/etc/bhex-user-env.sh"

