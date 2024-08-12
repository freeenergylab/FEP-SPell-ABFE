#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
===============================================================================
CopyRight (c) By freeenergylab.
@Description:
This module can be used to do Boresch-style restraint analytic correction.

@Author: Pengfei Li
@Date: May 30th, 2022
===============================================================================
"""
from math import pi, log, sin, sqrt, exp, erf

def restraint_correction(T, r0, theta0, kr, ktheta, kphi):
    """Do Boresch-style restraint analytic correction, refer to

    Calculation of Standard Binding Free Energies: Aromatic Molecules
    in the T4 Lysozyme L99A Mutant

    Yuqing Deng & Benoit Roux,
    J. Chem. Theory Comput. 2006, 2, 1255-1273

    https://pubs.acs.org/doi/epdf/10.1021/ct060037v
    
    Args:
        T (float): system simulation temperature
        r0 (float): equilibrium value for distance r
        theta0 (float): equilibrium value for angle theta
        kr (float): force constant for distance r
        ktheta (float): force constant for angle theta
        kphi (float): force constant for dihedral phi
 
    Returns:
        FtC0 + Fr (float): Boresch analytic restraint correction value
    """
    T      = float(T)
    r0     = float(r0)
    theta0 = float(theta0)
    kr     = float(kr)*2.0
    ktheta = float(ktheta)*2.0
    kphi   = float(kphi)*2.0

    V = 1660. # in A^3
    kB = 0.0019872041 # Boltzmann's constant (kcal/mol/K)
    FtC0 = -(kB*T)*log(r0**2*sin(theta0/180.0*pi)\
        *sqrt((2*pi*kB*T)**3/(kr*ktheta*kphi))/V) # Eq. 38
    Fr = -(kB*T)*log(1.0/(8*pi**2)*(2*pi*kB*T/kr)**1.5) # Eq. 40
    print('Translational: %.2f kcal/mol' % FtC0)
    print('Rotational: %.2f kcal/mol' % Fr)

    return FtC0 + Fr

def restraint_correction_schrodinger(T, r0, alpha0, theta0, gamma0, beta0, phi0,
                                     Kr, Kang, Kdihed):
    """Do Boresch-style restraint analytic correction, refer to

    Enhancing Hit Discovery in Virtual Screening through Absolute
    Protein-Ligand Binding Free-Energy Calculations

    Wei Chen & Lingle Wang et al.
    J Chem. Inf. Model. 2023, 63, 3171-3185
    
    https://pubs.acs.org/doi/10.1021/acs.jcim.3c00013
    https://chemrxiv.org/engage/chemrxiv/article-details/63a23f6116e9a872d32f81ef

    Args:
        T (float): system simulation temperature
        r0 (float): equilibrium value for r
        alpha0 (float): equilibrium value for alpha
        theta0 (float): equilibrium value for theta
        gamma0 (float): equilibrium value for gamma
        beta0 (float): equilibrium value for beta
        phi0 (float): equilibrium value for phi
        Kr (float): force constant for distance
        Kang (float): force constant for angle
        Kdihed (float): force constant for dihedral
 
    Returns:
        result (float): Boresch analytic restraint correction value
    """
    T      = float(T)
    r0     = float(r0)
    alpha0 = float(alpha0)
    theta0 = float(theta0)
    gamma0 = float(gamma0)
    beta0  = float(beta0)
    phi0   = float(phi0)
    Kr     = float(Kr)
    Kang   = float(Kang)
    Kdihed = float(Kdihed)

    V = 1660. # in A^3
    kB = 0.0019872041 # Boltzmann's constant (kcal/mol/K)
    beta = 1./(kB*T)
    Z_dist_r = r0/(2*beta*Kr)*exp(-beta*Kr*r0**2) + \
        sqrt(pi)/(4*beta*Kr*sqrt(beta*Kr))*(1+2*beta*Kr*r0**2)*(1+erf(sqrt(beta*Kr)*r0)) # Eq. 4
    Z_ang_alpha_r = sqrt(pi/(beta*Kang))*exp(-1./(4.*beta*Kang))*sin(alpha0/180.*pi) # Eq. 5
    Z_ang_theta_r = sqrt(pi/(beta*Kang))*exp(-1./(4.*beta*Kang))*sin(theta0/180.*pi) # Eq. 5
    Z_dihed_r =  sqrt(pi/(beta*Kdihed))*erf(pi*sqrt(beta*Kdihed)) # Eq. 6

    # result = -Z_ang_theta_r*Z_ang_alpha_r*kB*T*log(Z_dist_r*(Z_dihed_r**3)/(8.*(pi**2)*V)) # JCIM Eq. 3
    result = -kB*T*log(Z_dist_r*Z_ang_theta_r*Z_ang_alpha_r*(Z_dihed_r**3)/(8.*(pi**2)*V)) # ChemRxiv Eq. 3

    return result

if __name__ == '__main__':
    """Just do test for restraint_correction function.
    """
    correction = restraint_correction(
        T=298.15, ktheta=100., theta0=85.370, kphi=100., kr=10., r0=7.230,
        )
    print('The reference value is 8.99 kcal/mol, while this calculated value is %.2f kcal/mol.' % correction)