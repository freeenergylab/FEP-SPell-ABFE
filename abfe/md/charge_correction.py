#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the charge correction module that performs Poisson-Boltzmann based
corrections for charged ligand in solvated/complex using APBS software.

APBS: https://apbs.readthedocs.io/en/latest/getting/index.html

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import math
import os
import pandas as pd
import sys
import time

import griddata
import numpy as np
import parmed as pmd

import abfe.const as const
import abfe.utils.common_tools as common_tools

def align_complex(prmtop, mdcrd, solvent_mask, output_file='align.cpptraj.in'):
    """Align trajectory for charge correction calculations."""
    solute_mask = f'!{solvent_mask}'
    basename, ext = os.path.splitext(os.path.basename(mdcrd))
    new_mdcrd = f'{basename}.autoimage{ext}'
  
    with open(output_file, 'w') as _foutput:
        _foutput.write('parm %s\n' % prmtop)
        _foutput.write('trajin %s lastframe\n' % mdcrd)
        _foutput.write('autoimage anchor %s\n' % solute_mask)
        _foutput.write('center %s origin\n' % solute_mask)
        _foutput.write("principal '%s & !@H=' dorotation\n" % solute_mask)
        _foutput.write('trajout %s\n' % new_mdcrd)
        _foutput.write('go\n')
        _foutput.write('quit\n')

    cpptraj_cmd = [os.path.join(const.AMBERBIN, 'cpptraj'), '-i', output_file]
    status = os.system(' '.join(cpptraj_cmd))
    if status != 0:
        print(time.strftime("%c"))
        print("cpptraj failed to align system during charge correction!")
        sys.exit(1)

    return new_mdcrd

def compute_charge_correction(prmtop, mdcrd, lig_rname, solvent_mask, temperature,
                              epsilon_solv, wat_rname='WAT'):
    """Compute the free energy correction for ligand _decharging_ using the
    Poisson-Boltzmann (PB) method described by Rocklin et al.
    
    Args:
        prmtop (str): Amber parm7 file
        mdcrd (str): Amber trajectory file
        lig_rname (str): Ligand residue name
        solvent_mask (str): Solvent AmberMask
        temperature (float): Simulation temperature
        epsilon_solv (float): Dielectric constant of solvent
        wat_rname (str, optional): Water residue name. Defaults to 'WAT'.

    Returns:
        result (DataFrame): The charge correction free energy value from PB method.

    References:
        (1) Gabriel J. Rocklin, David L. Mobley, Ken A. Dill, and Philippe H. Hunenberger
        "Calculating the Binding Free Energies of Charged Species Based on Explicit-Solvent
        Simulations Employing Lattice-Sum Methods: An Accurate Correction Scheme for
        Electrostatic Finite-Size Effects" J. Chem. Phys. 139, 184103 (2013)

        https://aip.scitation.org/doi/pdf/10.1063/1.4826261
        https://aip.scitation.org/doi/suppl/10.1063/1.4826261

        (2) Wei Chen, Yuqing Deng, Ellery Russell, Yujie Wu, Robert Abel, and Lingle Wang
        "Accurate Calculation of Relative Binding Free Energies between Ligands with
        Different Net Charges" J. Chem. Theory Comput. 14, 6346 (2018)

        https://pubs.acs.org/doi/10.1021/acs.jctc.8b00825
    """
    mol = PBCAmberParm(prmtop, xyz=mdcrd)
    # delete copies of the ligand
    lig_indices = [r.number + 1 for r in mol.residues if r.name == lig_rname]
    if len(lig_indices) > 1:
        mol.strip(':' + ','.join([str(i) for i in lig_indices[1:]]))

    nwat = sum(1 for r in mol.residues if r.name == wat_rname)
    gamma = mol.compute_quadrupole_trace(wat_rname)
    # To round the charges since they are computed as a sum of many small numbers.
    tol = 1e-3
    lig_atoms = mol.residues[(lig_indices[0] - 1)].atoms
    lig_charge = round(sum(a.charge for a in lig_atoms) / tol)*tol
    # Strip the solvent as well as counter-ions and determine how to run APBS.
    mol.strip(solvent_mask)
    is_complex = (len(mol.residues) > 1)
    # Please refer to
    # Ref. 1 Sec. V "Procedure for appling the analytical correction scheme" part
    if not is_complex:
        rip_P = 0.0
        rip_I = compute_apbs_rip(mol, '!:%s'%lig_rname, lig_charge, temperature, epsilon_solv)
    else:
        rec_charge = round(sum(a.charge for a in mol.atoms) / tol)*tol - lig_charge
        rip_P = compute_apbs_rip(mol, ':%s'%lig_rname, rec_charge, temperature, epsilon_solv)
        rip_I = compute_apbs_rip(mol, '!:%s'%lig_rname, lig_charge, temperature, epsilon_solv)

    # Analytic correction scheme for finite-size effects.
    # According to Eq. 11 and Eq. 14 of Ref. 1, let's compute the
    # charge correction components: NET_USV, RIP, EMP(omitted), DSC.
    # As seen from Table IV of Ref. 1, EMP terms always are very small, so it can be omitted here.
    volume = mol.volume
    COULOMB = APBS.COULOMB
    xi_LS = APBS.xi_LS
    NET_USV = -(xi_LS*COULOMB*(lig_charge**2)/(2*epsilon_solv*volume**(1/3))) # Eq. 15 + Eq. 16 of Ref. 1
    RIP = ((rip_P + rip_I)*lig_charge)/volume # Eq. 17 of Ref. 1
    DSC = -2*math.pi*COULOMB*gamma*nwat*lig_charge/(3.*volume) # Eq. 35 of Ref. 1
    Total = NET_USV + RIP + DSC

    print('NET_USV: % 8.3f RIP: % 8.3f DSC: % 8.3f ChargeCorrection: % 8.2f\n'%(NET_USV, RIP, DSC, Total))

    result = pd.DataFrame({'NET_USV': [NET_USV],
                           'RIP': [RIP],
                           'DSC': [DSC],
                           'Total (kcal/mol)': [Total],
                           })
    return result

class PBCAmberParm(pmd.amber.AmberParm):
    """This class is used to handle periodic box info and others.

    Attributes: 
        volume (float): the periodic cell volume (in A^3)
    """
    def __init__(self, *args, **kwargs):
        super(PBCAmberParm, self).__init__(*args, **kwargs)

    @property
    def volume(self):
        """The volume of the periodic cell (in A^3)"""
        volume = self.box[:3].prod()
        # IFBOX: set to 1 if periodic cubic box, 2 when truncated octahedral
        if self.pointers['IFBOX'] == 1:
            return volume 
        elif self.pointers['IFBOX'] == 2:
            # Conversion factor from AmberTools/src/leap/src/leap/tools.c
            return volume*0.7698004

    def compute_orthorhombic_box(self, iso=False):
        """Compute an equivalent orthorhombic cell for the system.

        Args:
            iso (bool, optional):
                If True, return a cubic box. The default elongates the box
                to encompass the coordinates of the system.

        Returns:
            lens (ndarray): three lengths of the box edges
        """
        if self.pointers['IFBOX'] == 1: 
            return self.box[:3]
        elif self.pointers['IFBOX'] >= 2:
            if iso: return np.ones(3)*self.volume**(1/3.)
            # Rescale the cubic dimensions by considering the shape of the system.
            dims = self.coordinates.max(axis=0) - self.coordinates.min(axis=0)
            scale = dims / dims[0]
            lens = scale*(self.volume / scale.prod())**(1/3.)
            while np.any(dims > lens):
                lens *= 1.05
            return lens

    def zero_charge(self, mask):
        """Zero the charges on atoms in a given mask selection."""
        for i in pmd.amber.AmberMask(self, mask).Selected():
            self.atoms[i].charge *= 0.0

    def compute_quadrupole_trace(self, resname='WAT'):
        """Compute the quadrupole trace of a solvent molecule (in e-A^2).
        Ref. 1 Sec. II E: "...for a solvent model with a single vdW interaction site..."
        """
        mask = pmd.amber.AmberMask(self, f':{resname}')
        resid = self.atoms[next(mask.Selected())].residue.number + 1
        mask = pmd.amber.AmberMask(self, ':%d'%resid)
        centers = [i for i in mask.Selected() if self.atoms[i].epsilon > 0.0]
        if len(centers) != 1:
            error_msg = 'Can only compute quadrupole trace for water molecule with one vdW center.'
            raise RuntimeError(error_msg)
        xj = self.coordinates[centers[0]]
        gamma = 0.0
        for i in mask.Selected():
            xi = self.coordinates[i]
            gamma += self.atoms[i].charge*((xi - xj)**2).sum()

        return gamma

def compute_apbs_rip(mol, chg_mask, net_charge, temperature, epsilon_solv):
    """Compute RIP energy term using APBS.

    Args:
        mol (AmberParm object): the molecule object
        chg_mask (str): in Amber-Mask style
        net_charge (float): the net charge
        temperature (float): the simulation temperature
        epsilon_solv (float): the dielectric constant of solvent

    Returns:
        rip (float):  in kcal-A^3/mol-e (in AMBER units, not APBS units).
    """
    tmp_mol = mol[:]
    # zero out charges in chg_mask
    tmp_mol.zero_charge(chg_mask) 
    tmp_mol.save('apbs.pqr', None, True)
    del tmp_mol
    is_complex = (len(mol.residues) > 1)
    box = mol.compute_orthorhombic_box(iso=(not is_complex))

    apbs = APBS('apbs.pqr', temperature, epsilon_solv, box)
    apbs.run('apbs.in')
    rip = apbs.compute_rip('apbs.dx', net_charge)

    return rip

class APBS(object):
    """Create and launch APBS jobs

    Args: 
        temperature (float): in K
        sdie (float): the dielectric of the external medium (unitless)
        grid_spacing (float): maximum grid spacing (in A/voxel)
    """
    xi_CB = -2.380077 # for a cubic box, unitless from Sec. II D of Ref. 1
    xi_LS = -2.837297 # for a cubic box, unitless from Sec. II D of Ref. 1
    BOLTZMANN = 0.001987192 # from AMBER, in kcal/mol-K
    AMBER_ELECTROSTATIC = 18.2223 # from parmed.constants.AMBER_ELECTROSTATIC
    COULOMB = AMBER_ELECTROSTATIC**2 # in kcal-A/mol, 4*pi*epsilon0 = COULOMB in AMBER

    def __init__(self, pqr, temperature, sdie, box, grid_spacing=0.2, pdie=1.0):
        self._pqr = str(pqr) 
        self.temperature = max(0.0, float(temperature))
        self.sdie = max(1.0, float(sdie))
        self.box = np.array(box, dtype=np.float64).flatten()
        assert self.box.size == 3
        self.grid_spacing = max(0.01, float(grid_spacing))
        self.pdie = max(1.0, float(pdie))

    @property
    def dime(self):
        """Uniform grid size that gives the requested grid spacing

        For performance reasons, APBS requires grid sizes of the form 2^k + 1,
        where k is a natural number. For simplicity, the lowest value of k
        considered here is 6. The value is iteratively increased until the grid
        spacing in each direction is now greater than `grid_spacing`, 
        """
        ngrid = lambda k: 2**k + 1
        k = 6
        dime = ngrid(k)
        while not np.all((self.box / dime) <= self.grid_spacing) and k < 8:
            k += 1
            dime = ngrid(k)

        return [dime]*3

    def run(self, apbs_in, do_force=False):
        """Run APBS with the current system settings.

        Args:
            apbs_in (str): write APBS input to apbs_in file
            do_force (bool, optional): if True, compute forces
        """
        lines = [
            'read',
            '    mol pqr %s'%self._pqr,
            'end',
            '',
            'elec name potential',
            '    mg-manual',
            '    dime %s'%(' '.join(['%d'%n for n in self.dime])),
            '    glen %s'%(' '.join(['%f'%b for b in self.box])),
            '    gcent 0.0 0.0 0.0',
            '    mol 1',
            '    lpbe',
            '    bcfl mdh',
            '    pdie %f'%self.pdie,
            '    sdie %f'%self.sdie,
            '    chgm spl4',
            '    srfm smol',
            '    srad 1.4',
            '    swin 0.3',
            '    sdens 40.0',
            '    temp %f'%self.temperature,
            '    calcenergy total',
            '    calcforce %s'%('total' if do_force else 'no'),
            '    write pot dx apbs',
            'end',
            ]
        with open(apbs_in, 'w') as _apbs_in:
            _apbs_in.write('\n'.join(lines))

        cmd_list = ['apbs', apbs_in, '>', 'apbs.out']
        os.system(' '.join(cmd_list))

    def compute_rip(self, dx_file, net_charge):
        """Compute the residual integrated potential (RIP).
        Ref. 1 SI provides the code example for calculating RIP:
        https://aip.scitation.org/doi/suppl/10.1063/1.4826261
        
        The result is in kcal-A^3/mol-e (in AMBER units, not APBS units).
        """
        dxfile = griddata.Grid(dx_file)
        dV = dxfile.delta.prod()
        V = dxfile.grid.size*dV
        # Label this as B for box integral, X because it refers to the system X.
        # Ref. 1 Eq. 19
        B_X = dxfile.grid.sum()*dV*self.BOLTZMANN*self.temperature
        # Calculate the expected total potential based on the net charge QX
        # Label this as B for box integral, QX because it refers to net charge.
        # Ref. 1 Eq. 21
        B_QX = -(self.xi_CB*self.COULOMB/self.sdie*net_charge*V**(2./3.))
        # Calculate the RIP.
        # Ref. 1 Eq. 18
        RIP = B_X - B_QX

        return RIP

if __name__ == '__main__':
    """Do charge correction for ABFE prediction for a charged ligand in solvated/complex.
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--prmtop', type=str, default='')
    parser.add_argument('--mdcrd', type=str, default='')
    parser.add_argument('--solvent_mask', type=str, default=':WAT,K+,Na+,Cl-')
    parser.add_argument('--temperature', type=float, default=298.15)
    args = parser.parse_args()

    work_dir = os.path.dirname(args.mdcrd)
    with common_tools.DirManager(work_dir):
        new_mdcrd = align_complex(args.prmtop, args.mdcrd, args.solvent_mask)
        epsilon_solv = 97. # relative permittivity for TIP3P water
        dG_correction = compute_charge_correction(
            prmtop=args.prmtop, mdcrd=new_mdcrd, lig_rname='MOL',
            solvent_mask=args.solvent_mask, temperature=args.temperature,
            epsilon_solv=epsilon_solv, wat_rname='WAT'
            )
        with open('results.csv', 'w') as f:
            dG_correction.to_csv(f, float_format='%.3f')