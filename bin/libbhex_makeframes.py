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

# Options common to train & test data:
SI_deltaT = 1./512
SI_datalen= 4.0
total_mass = 150. 
distance=1. # Mpc

#series_names = ['HRq-series'] # (see above for valid choices)
series_names = ["Eq-series", "HRq-series", "HR-series",  "Lq-series",
        "RO3-series",  "Sq-series",  "S-series-v2",  "TP2-series",
        "TP-series"]

#
# Modify for imposing parameter bounds on the catalogue:
#
#bounds=None
bounds=dict()
bounds['a1'] = [0, 0]
bounds['a2'] = [0, 0]
#bounds['q'] = [-np.inf, 3] 

catalogue_name = 'NonSpinning'

save_pcs = ['NRhplusTimeSeriesPCA', 'NRhcrossTimeSeriesPCA',
'NRAmpTimeSeriesPCA', 'NRPhaseTimeSeriesPCA',  'SIhplusTimeSeriesPCA',
'SIhcrossTimeSeriesPCA', 'SIAmpTimeSeriesPCA', 'SIPhaseTimeSeriesPCA']

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

for simulation in simulation_selection.simulations:

    wavefile = os.path.basename(simulation['wavefile'])
        
    # Write the config file
    meta_file_name = wavefile.replace('asc', 'bbh')
    meta_file = open(meta_file_name, 'w')
    text = template.format(mass_ratio=simulation['q'],
            wavename=simulation['wavename'],
            wavefile=os.path.basename(simulation['wavefile'])
            )
    meta_file.writelines(text)
    meta_file.close()

    # Copy the ascii file here
    shutil.copyfile(simulation['wavefile'],
            os.path.basename(simulation['wavefile']))

    # Make frame
    subprocess.call(['lalapps_fr_ninja', 
        "--verbose", "--format", "NINJA2", "--double-precision", 
        "--nr-data-dir", "./", "--nr-meta-file", meta_file_name,
        "--output", meta_file_name.replace('bbh','gwf')])

    # Clean up
    os.remove(wavefile)
    os.remove(meta_file_name)

