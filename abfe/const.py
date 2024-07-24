#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the const module that setup some constant(default) parameters.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import os

AMBERHOME = os.environ['AMBERHOME']
AMBERBIN = os.path.join(os.environ['AMBERHOME'], 'bin')
LEAPHOME = os.path.join(os.environ['AMBERHOME'], 'dat/leap')

WATER_DEFS = {
    'tip3p': {'source': 'leaprc.water.tip3p', 'box': 'TIP3PBOX'},
    'tip4pew': {'source': 'leaprc.water.tip4pew', 'box': 'TIP4PEWBOX'},
    'opc': {'source': 'leaprc.water.opc', 'box': 'OPCBOX'}
    }

PRO_DEFS = {
    'ff14SB': {'source': 'leaprc.protein.ff14SB'},
    }

PHOSAA_DEFS = {
    'phosaa14SB': {'source': 'leaprc.phosaa14SB'},
    }

LIPID_DEFS = {
    'lipid21': {'source': 'leaprc.lipid21'},
    }