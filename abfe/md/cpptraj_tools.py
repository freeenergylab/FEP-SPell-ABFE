#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the cpptraj tools module that can be used to do trajectory analysis.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import math
import numpy as np
from scipy import stats

SIX_DOFS_INPUT = """parm PRMTOP
trajin MDCRD

distance bnd_r L1 P1 out bnd_r.dat
angle bnd_alpha P1 L1 L2 out bnd_alpha.dat
angle bnd_theta P2 P1 L1 out bnd_theta.dat
dihedral bnd_gamma P1 L1 L2 L3 out bnd_gamma.dat
dihedral bnd_beta P2 P1 L1 L2 out bnd_beta.dat
dihedral bnd_phi P3 P2 P1 L1 out bnd_phi.dat
"""

def mean_std(input_file, isPeriodic=False):
    """Calculate the average and standard deviation for data
        from the cpptraj .dat file.

    Args:
        input_file (str):
            the data file
        isPeriodic (bool, optional):
            dihedral/torsion should be periodic, defaults to False.

    Returns:
        mean (float): the average value
        std (float): the standard deviation value
    """
    _values = []
    with open(input_file, 'r') as _finput:
        lines = _finput.readlines()
        for line in lines[1:]:
            _, _value = line.split()
            _values.append(float(_value))
    if isPeriodic:
        d2r = math.pi / 180.
        values = [value*d2r for value in _values] # convert degree to radian
        mean = stats.circmean(values, high=math.pi, low=-math.pi)/d2r # convert radian to degree
        std = stats.circstd(values, high=math.pi, low=-math.pi)/d2r # convert radian to degree
    else:
        values = np.array(_values)
        mean = round(values.mean(), 2)
        std = round(values.std(), 2)

    return mean, std

def autoimage(output_file, parm7, rst7):
    """Construct the cpptraj script to do autoimage.

    Args:
        output_file (str): the output cpptraj script file
        parm7 (str): the .parm7 file
        rst7 (str): the .rst7 file
    """
    with open(output_file, 'w') as _foutput:
        _foutput.write('parm %s\n' % parm7)
        _foutput.write('trajin %s\n' % rst7)
        _foutput.write('autoimage\n')
        _foutput.write('trajout %s\n' % rst7)
        _foutput.write('go\n')
        _foutput.write('quit\n')