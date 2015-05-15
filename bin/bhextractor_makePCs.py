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
bhextractor_makePCs.py

Construct catalogues and principal component analysis for NR BBH waveforms; save
the PCA results to file using the file dump method on the PCA result object
"""

import sys
import bhextractor_pca as bhex

# -------------------------------
# USER INPUT

catalogue_name=sys.argv[1]
theta=float(sys.argv[2])

# END USER INPUT
# -------------------------------

# -------------------------------
# ANALYSIS

#
# Setup and then build the catalogue
#
catalogue = bhex.waveform_catalogue(catalogue_name=catalogue_name, fs=2048,
        catalogue_len=4, mtotal_ref=250, Dist=1., theta=theta)

#
# Do the PCA
#
pca = bhex.waveform_pca(catalogue)

#
# Dump files
#
pca.file_dump()
