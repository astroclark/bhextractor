#! /bin/sh
# Make env setup script
# James Clark, james.clark@ligo.org

# Get the location of the git repository by finding the full path of this file
# and stripping off the name
BHEX_PREFIX=`python -c "import os, sys; print os.path.realpath('${0}')" | sed 's|/setup.sh||g'`
FRGETVECT_PATH="/home/jclark/opt/xpipeline/share/xpipeline/matlab/"
MATLABPATH=${MATLABPATH}:"${BHEX_PREFIX}/bin/matlab"
PYTHONPATH=${PYTHONPATH}:"${BHEX_PREFIX}/bin/hht"

echo "export BHEX_PREFIX=${BHEX_PREFIX}" > bhex_env.sh
echo "export FRGETVECT_PATH=${FRGETVECT_PATH}" >> bhex_env.sh
echo "export PATH=${PATH}:${BHEX_PREFIX}/bin" >> bhex_env.sh
echo "export MATLABPATH=${MATLABPATH}" >> bhex_env.sh
echo "export PYTHONPATH=${PYTHONPATH}" >> bhex_env.sh


