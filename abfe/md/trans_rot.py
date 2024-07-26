#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the trans_rot module that can be used to do Boresch-style restraint
analytic correction.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
from math import pi, log, sin, sqrt

def restraint_correction(T, ktheta, theta0, kphi, kr, r0):
    """Do Boresch-style restraint analytic correction.

    Args:
        T (float): system simulation temperature
        ktheta (float): force constant for theta
        theta0 (float): equilibrium value for theta
        kphi (float): force constant for phi
        kr (float): force constant for r
        r0 (float): equilibrium value for r
 
    Returns:
        FtC0 (float): Translational components
        Fr (float): Rotational components
        FtC0 + Fr (float): Boresch analytic restraint correction value
    """
    ktheta = float(ktheta)*2.0
    theta0 = float(theta0)
    kphi   = float(kphi)*2.0
    kr     = float(kr)*2.0
    r0     = float(r0)

    kB = 0.0019872041 # Boltzmann's constant (kcal/mol/K)
    FtC0 = -(kB*T)*log(r0**2*sin(theta0/180.0*pi)\
        *sqrt((2*pi*kB*T)**3/(kr*ktheta*kphi))/1660.0)
    Fr = -(kB*T)*log(1.0/(8*pi**2)*(2*pi*kB*T/kr)**1.5)
    print('Translational: %.2f kcal/mol' % FtC0)
    print('Rotational: %.2f kcal/mol' % Fr)

    return FtC0, Fr, FtC0 + Fr

if __name__ == '__main__':
    """Just one test for restraint_correction function.
    """
    correction = restraint_correction(T=300., ktheta=100., theta0=63.0, 
                         kphi=100., kr=10., r0=9.8) 
    print('The reference value is 8.74 kcal/mol,\
          \nwhile this calculated value is %.2f kcal/mol.' % correction)