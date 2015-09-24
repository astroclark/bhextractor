#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2015-2016 James Clark <james.clark@ligo.org>
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
"""
libbhex_makeframes.py

Construct a catalogue and build the NINJA input frames for each waveform in that
catalogue
"""

import os,sys
import shutil, subprocess
from bhex_utils import bhex_wavedata as bwave
from bhex_utils import bhex_pca as bpca


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT: define the catalogue

#series_names = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
#        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series",
#        "TP-series"]

total_mass = 100

series_names = ['HRq-series'] # (see above for valid choices)

#
# Modify for imposing parameter bounds on the catalogue:
#
bounds=None

catalogue_name = 'HRq'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do the work:

print '~~~~~~~~~~~~~~~~~~~~~'
print 'Selecting Simulations'
print ''
simulation_selection = \
        bwave.simulation_details(series_names=series_names,
                param_bounds=bounds, Mmin30Hz=total_mass)


print '~~~~~~~~~~~~~~~~~~~~~~~'
print 'Constructing NINJA data'
print ''

#
# NINJA Files
#

template="""mass-ratio = {mass_ratio}
simulation-details = {wavename}
nr-group = GATech

2,2 = {wavefile}"""

"writing ninja config files for these simulations"

if os.path.exists(catalogue_name):
    print >> sys.stderr, "%s exists, remove"%catalogue_name
    sys.exit()
os.makedirs(catalogue_name)

for simulation in simulation_selection.simulations:


    wavefile = os.path.join(catalogue_name,
            os.path.basename(simulation['wavefile']))
        
    # Write the config file
    meta_file_name = wavefile.replace('asc', 'bbh')
    meta_file = open(meta_file_name, 'w')
    text = template.format(mass_ratio=simulation['q'],
            wavename=simulation['wavename'],
            wavefile=wavefile
            )
    meta_file.writelines(text)
    meta_file.close()

    # Copy the ascii file here
    shutil.copyfile(simulation['wavefile'], wavefile)

    # Make frame
    subprocess.call(['lalapps_fr_ninja', 
        "--verbose", "--format", "NINJA2", "--double-precision", 
        "--nr-data-dir", "./", "--nr-meta-file", meta_file_name,
        "--output", meta_file_name.replace('bbh','gwf')])

    # Clean up
    os.remove(wavefile)
    os.remove(meta_file_name)

# Create the NINJA xml
subprocess.call(['lalapps_ninja', "--format", "NINJA2", "--datadir",
    catalogue_name, "--outfile", catalogue_name+".xml", "--min-mass-ratio", "0",
    "--max-mass-ratio", "100", "--pattern", "*gwf"])

