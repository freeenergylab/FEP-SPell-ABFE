#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
===============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the alchemy md module that runs alchemical simulations.

@Author: Pengfei Li
@Date: May 30th, 2022
===============================================================================
"""
import argparse
import os
import pickle
import random
import re
import socket
import stat
import time
from collections import OrderedDict
from math import ceil

import numpy as np
import pandas as pd
import parmed as pmd

import abfe.const as const
import abfe.md.cpptraj_tools as cpptraj_tools
import abfe.md.md_tools as md_tools
import abfe.md.analytic as analytic
import abfe.utils.common_tools as ctools
from abfe.md.amber_alchemy_mdin import SOLVATED_MIN
from abfe.md.amber_alchemy_mdin import SOLVATED_HEAT, SOLVATED_PRESS
from abfe.md.amber_alchemy_mdin import SOLVATED_TI

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parse arguments for alchemy md module!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    # construction arguments
    parser.add_argument(
        '-bw', '--abfe_workdir',
        help='Workdir'
        ) 
    parser.add_argument(
        '-cn', '--complex_name',
        help='Str contains the complex name'
        )
    parser.add_argument(
        '-ln', '--ligand_name',
        help='Str contains the ligand name'
        )
    parser.add_argument(
        '-ho', '--hfe_only', action='store_true', default=False,
        help='Only perform HFE'
        )
    # input arguments
    parser.add_argument(
        '-sl', '--solvated_length', type=float, default=5.0,
        help='Length of the solvated simulations (in ns).'
        )
    parser.add_argument(
        '-cl', '--complex_length', type=float, default=5.0,
        help='Length of the complex simulations (in ns).'
        )
    parser.add_argument(
        '--restraintmask', type=str, default="!:WAT,Cl-,K+,Na+ & !@H=",
        help='specify atoms for positional restraint during min-heat-press simulations'
        )
    parser.add_argument(
        '--restraint_wt', type=float, default=5.0,
        help='force constant of positional restraint, in kcal/mol/A^2'
        )
    parser.add_argument(
        '--cpus_per_task', type=int, default=2,
        help='The number of CPU process cores'
        )
    parser.add_argument(
        '-dr', '--dry_run', action='store_true', default=False,
        help='If turned on, dry run.'
        )
    parser.add_argument(
        '--partition', type=str, default='gpuq',
        help='Slurm partition.'
        )
    parser.add_argument(
        '--exclude', type=str, default='',
        help='Slurm exclude nodes.'
        )
    parser.add_argument(
        '--gpus_per_node', type=int, default=4,
        help='Number of gpu cards per node'
        )
    # alchemical simulation protocol and control options
    parser.add_argument(
        '--abfe_solvated_vdw_lambdas', nargs='*', type=float, default=None
        )
    parser.add_argument(
        '--abfe_complex_vdw_lambdas', nargs='*', type=float, default=None
        )
    parser.add_argument(
        '--abfe_complex_restraint_lambdas', nargs='*', type=float, default=None
        )
    parser.add_argument(
        '--remd', action='store_true',
        help='Run alchemical simulations using replica exchange'
        )
    parser.add_argument(
        '--remd_numexchg', type=int, default=1000,
        help='Number of exchange attempts that will be performed between replica pairs'
        )
    parser.add_argument(
        '--temperature', type=float, default=298.15,
        help='System temperature (in K)'
        )
    parser.add_argument(
         '--heat_temps', default=[5.0, 100.0, 200.0, 298.15], type=float, nargs='*',
         help='Temperature (in K) sequence for heat protocol'
        )
    parser.add_argument(
        '--cutoff', type=float, default=9.0,
        help='Cutoff radius for long-range interactions (in A)'
        )
    parser.add_argument(
        '--timestep', type=float, default=0.002,
        help='Timestep for integral timestep in simulations (in ps)'
        )
    parser.add_argument(
        '--ntwx', type=float, default=1.,
        help='Every ntwx/timestep steps to save trajectory (in ps)'
        )
    parser.add_argument(
        '--gti_add_sc', type=int, default=1,
        help='Treatment of the non-bonded interactions of SC internal and SC-CC cross regions'
        )
    # Boresch-style restraint options
    boresch_grp = parser.add_argument_group('Boresch-style restraints')
    boresch_grp.add_argument(
        '-bc', '--bond_const', default=10.0, type=float,
        help='Force constant for bond restraint, in kcal/mol/A^2'
        )
    boresch_grp.add_argument(
        '-ac', '--angle_const', default=100.0, type=float,
        help='Force constant for angle restraint, in kcal/mol/radian^2'
        )
    boresch_grp.add_argument(
        '-dc', '--dihedral_const', default=100.0, type=float,
        help='Force constant for dihedral restraint, in kcal/mol/radian^2'
        )
    boresch_grp.add_argument(
        '-rs', '--restraint_seed', type=int, default=None,
        help='Random seed used to set atoms Boresch restraints'
        )
    boresch_grp.add_argument(
        '--lig_restraint_atoms', type=str, nargs=3,
        default=[], help='Directly select ligand atoms for restraint'
        )
    boresch_grp.add_argument(
        '--rec_restraint_atoms', type=str, nargs=3,
        default=[], help='Directly select receptor atoms for restraint'
        )

    args = parser.parse_args()
    if args.restraint_seed is None:
        args.restraint_seed = random.randint(0, 99999)
    print(f'Using restraint seed of {args.restraint_seed} for all steps.')

    args.heat_temps[-1] = args.temperature

    return args

class BoreschRestraint(object):
    """A Boresch-style restraint is a set of six harmonic restraints that orient
    a ligand (lig) relative to a receptor (rec). Each restraint corresponds to
    a set of six internal coordinates: one distance, two angles, and three dihedrals.
    """
    def __init__(self, prmtop, lig_masks, rec_masks, rst_refs, rst_wts):
        self._mol = self.prmtop_factory(prmtop)
        self._lig_masks = lig_masks
        self._rec_masks = rec_masks
        self._rst_refs = rst_refs
        self._rst_wts = rst_wts

    def prmtop_factory(self, prmtop, xyz=None):
         """Convert input to a parmed AmberParm object or else check that it
         already is an AmberParm object.
         """
         if isinstance(prmtop, str):
             mol = pmd.load_file(prmtop, xyz)
         elif isinstance(prmtop, pmd.amber.AmberParm):
             mol = prmtop
         else:
             _type = type(prmtop)
             raise ValueError(f'Cannot create AmberParm object from object of type {_type}.')
         return mol

    def _gen_rst_str(self, lig_idx):
        """Constructe multiple restraints"""
        dis_fmt = '&rst iat=%d,%d r1=0.0, r2=%.2f, r3=%.2f, r4=99.0, rk2=%.2f, rk3=%.2f/'
        ang_fmt = '&rst iat=%d,%d,%d r1=0.0, r2=%.2f, r3=%.2f, r4=180.0, rk2=%.2f, rk3=%.2f/'
        dih_fmt = '&rst iat=%d,%d,%d,%d r1=%.2f, r2=%.2f, r3=%.2f, r4=%.2f, rk2=%.2f, rk3=%.2f/'
        # atom indices
        L1, L2, L3 = self._lig_atoms(lig_idx)
        P1, P2, P3 = self._rec_atoms()
        # equilibrium values
        r0, alpha0, theta0, gamma0, beta0, phi0 = self._rst_refs
        # force constants
        KR, KALPHA, KTHETA, KGAMMA, KBETA, KPHI = self._rst_wts
        # adjusted boundaries for dihedrals (accounts for periodicity)
        dih11, dih14 = gamma0 - 180.0, gamma0 + 180.0
        dih21, dih24 = beta0 - 180.0, beta0 + 180.0
        dih31, dih34 = phi0 - 180.0, phi0 + 180.0
        # Build the list of rst strings.
        rst_list = [
            dis_fmt%(L1, P1, r0, r0, KR, KR),
            ang_fmt%(P1, L1, L2, alpha0, alpha0, KALPHA, KALPHA),
            ang_fmt%(P2, P1, L1, theta0, theta0, KTHETA, KTHETA),
            dih_fmt%(P1, L1, L2, L3, dih11, gamma0, gamma0, dih14, KGAMMA, KGAMMA),
            dih_fmt%(P2, P1, L1, L2, dih21, beta0, beta0, dih24, KBETA, KBETA),
            dih_fmt%(P3, P2, P1, L1, dih31, phi0, phi0, dih34, KPHI, KPHI)
        ]
        return '\n'.join(rst_list)

    def _lig_atoms(self, lig_idx=0):
        """Get a list of the indices of the restrained ligand atoms.

        Parameters
        ----------
        lig_idx (int, optional):
          The index of the ligand copy to extract. The first/only copy is returned by default.

        Returns
        -------
        list of int:
          The atom indices in the prmtop frame corresponding to a copy of the
          ligand atoms needed to define the restraint coordinates
        """
        indices = []
        for mask in self._lig_masks:
            mask_iter = pmd.amber.AmberMask(self._mol, mask).Selected()
            idx = [i+1 for i in mask_iter]
            indices.append(idx[lig_idx])
        return indices

    def _rec_atoms(self):
        """Get a list of the indices of the restrained receptor atoms.

        Returns
        -------
        list of int:
          The atom indices in the prmtop frame corresponding to the receptor
          atoms needed to define the restraint coordinates
        """
        indices = []
        for mask in self._rec_masks:
            mask_iter = pmd.amber.AmberMask(self._mol, mask).Selected()
            idx = [i+1 for i in mask_iter]
            indices.append(idx[0])
        return indices

def _add_boresch_restraint(args, prmtop, lig_masks, rec_masks, lig_idx, rst_refs, rst_wts):
    """Add Boresch-style restraints on protein-ligand."""
    fixed_lig = (len(args.lig_restraint_atoms) == 3)
    fixed_rec = (len(args.rec_restraint_atoms) == 3)
    if fixed_lig and fixed_rec:
        br = BoreschRestraint(
            prmtop=prmtop, lig_masks=lig_masks, rec_masks=rec_masks,
            rst_refs=rst_refs, rst_wts=rst_wts
            )
        rst = br._gen_rst_str(lig_idx=lig_idx)
    else: #TODO random ligand and random receptor
        rst = 'pass'
        print('Currently, only support manual selections on Boresch-style restraint six atoms.')

    with open('restraints.inp', 'w') as _frst:
        _frst.write(str(rst) + '\n')

if __name__ == '__main__':
    """Parse the input arguments, then run the alchemy md stage.
    """
    args = _parse_args()
    print(time.strftime("%c"))
    print(f"Alchemy MD started on host {socket.gethostname()}.")
    #==========================================================================
    # 1. SOLVATED
    #==========================================================================
    alchemy_morph = os.path.join(args.abfe_workdir, '_alchemy_morph', args.ligand_name)
    morph_pkl_file = os.path.join(alchemy_morph, 'objects.pkl')
    error_msg = "Please ensure that each ligand has one objects.pkl file!"
    assert os.path.isfile(morph_pkl_file)==True, error_msg

    lig = pickle.load(open(morph_pkl_file, "rb"))
    workdir = os.path.join(args.abfe_workdir, '_alchemy_md', args.ligand_name)
    #==========================================================================
    # 1.1 VDW
    #==========================================================================
    #==========================================================================
    # 1.1.1 input files
    #==========================================================================
    lig_md_steps = OrderedDict()
    _lambdas =[str(_lambda) for _lambda in args.abfe_solvated_vdw_lambdas]
    for _lambda in _lambdas:
        if _lambda == _lambdas[0] or _lambda == _lambdas[-1]:
            ntwx = args.ntwx
        else:
            ntwx = 0.
        lig_md_steps['min'] = SOLVATED_MIN\
            .replace('CLAMBDA', _lambda)\
            .replace('CUT', str(args.cutoff))\
            .replace('GTI_ADD_SC', str(args.gti_add_sc))\
            .replace('TIMASK1', str(lig.timask1))\
            .replace('TIMASK2', str(lig.timask2))\
            .replace('SCMASK1', str(lig.scmask1))\
            .replace('SCMASK2', str(lig.scmask2))\
            .replace('RESTRAINT_WT', str(args.restraint_wt))\
            .replace('RESTRAINTMASK', str(args.restraintmask))
        for n, (tempi, temp0) in enumerate(zip(args.heat_temps[:-1], args.heat_temps[1:]), 1):
            lig_md_steps['heat-%s'%str(n)] = SOLVATED_HEAT\
                .replace('CLAMBDA', str(_lambda))\
                .replace('DT', str(0.002))\
                .replace('CUT', str(args.cutoff))\
                .replace('TEMPI', str(tempi))\
                .replace('TEMP0', str(temp0))\
                .replace('TEMP_0', 'TEMP0')\
                .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                .replace('TIMASK1', str(lig.timask1))\
                .replace('TIMASK2', str(lig.timask2))\
                .replace('SCMASK1', str(lig.scmask1))\
                .replace('SCMASK2', str(lig.scmask2))\
                .replace('RESTRAINT_WT', str(args.restraint_wt))\
                .replace('RESTRAINTMASK', str(args.restraintmask))
            lig_md_steps['press-%s'%str(n)] = SOLVATED_PRESS\
                .replace('CLAMBDA', str(_lambda))\
                .replace('DT', str(0.002))\
                .replace('CUT', str(args.cutoff))\
                .replace('TEMP0', str(temp0))\
                .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                .replace('TIMASK1', str(lig.timask1))\
                .replace('TIMASK2', str(lig.timask2))\
                .replace('SCMASK1', str(lig.scmask1))\
                .replace('SCMASK2', str(lig.scmask2))\
                .replace('RESTRAINT_WT', str(args.restraint_wt))\
                .replace('RESTRAINTMASK', str(args.restraintmask))
        if args.remd:
            NSTLIM = str(int(args.solvated_length/args.timestep*1000./args.remd_numexchg))
        else:
            NSTLIM = str(int(args.solvated_length/args.timestep*1000))
        lig_md_steps['ti'] = SOLVATED_TI\
            .replace('CLAMBDA', str(_lambda))\
            .replace('MBAR_STATES', str(len(_lambdas)))\
            .replace('MBAR_LAMBDA', ', '.join(_lambdas))\
            .replace('NSTLIM', NSTLIM)\
            .replace('DT', str(args.timestep))\
            .replace('CUT', str(args.cutoff))\
            .replace('TEMPI', str(args.temperature))\
            .replace('TEMP0', str(args.temperature))\
            .replace('NTWX', str(int(ntwx/args.timestep)))\
            .replace('GTI_ADD_SC', str(args.gti_add_sc))\
            .replace('TIMASK1', str(lig.timask1))\
            .replace('TIMASK2', str(lig.timask2))\
            .replace('SCMASK1', str(lig.scmask1))\
            .replace('SCMASK2', str(lig.scmask2))
        if args.remd:
            lig_md_steps['ti'] = lig_md_steps['ti']\
                .replace('!gremd_acyc =', 'gremd_acyc =')\
                .replace('!numexchg =', 'numexchg =')\
                .replace('NUMEXCHG', str(args.remd_numexchg))
        _workdir = os.path.join(workdir, 'vdw', str(_lambda))
        with ctools.DirManager(_workdir):
            os.system(f"ln -s -f {lig.alchemy_steps['vdw'][0]} .")
            os.system(f"ln -s -f {lig.alchemy_steps['vdw'][1]} .")
            for step, protocol in lig_md_steps.items():
                md_tools.write_mdin_file(f'{step}.in', protocol)
    #==========================================================================
    # 1.1.2 bash files (submit jobs scripts)
    #==========================================================================
    _workdir = os.path.join(workdir, 'vdw')
    with ctools.DirManager(_workdir):
        _mdsteps = list(lig_md_steps.keys())
        _fparm7 = lig.alchemy_steps['vdw'][0]
        _frst7 = lig.alchemy_steps['vdw'][1]
        if args.remd:
            group_cmd = md_tools.write_groupfile(_mdsteps, _fparm7, _lambdas)
            with open('run.remd.sh', 'w') as _fsh:
                _fsh.write(group_cmd)
            os.chmod('run.remd.sh', stat.S_IRWXU)
            _mdsteps.remove('ti')
        _last_md_info = md_tools.write_alchemy_md_submit_sh(
            fmdsh='run.sh',
            mdsteps=_mdsteps,
            fparm7=_fparm7,
            frst7=_frst7,
            nprocess=args.cpus_per_task,
            dry_run=True,
            )
        slurm_params = {
        'job-name': args.ligand_name+'.vdw.md',
        'output': 'slurm.out',
        'error': 'slurm.err',
        'partition': args.partition,
        'nodes': 1,
        'ntasks-per-node': 1,
        'cpus-per-task': 2,
        'time': '72:00:00',
        'gpujob': True,
        'gpus': 1,
        'array': '0-%s'%(str(len(_lambdas)-1)),
        }
        exclude = args.exclude
        if exclude.strip(): 
            slurm_params.update(exclude=exclude)

        cmd_list = [
            'DIR_ARR=( %s )'%' '.join(_lambdas),
            'DIR=${DIR_ARR[$SLURM_ARRAY_TASK_ID]}',
            'cd $DIR',
            'echo INFO: ARRAY JOB ID $SLURM_ARRAY_JOB_ID',
            'echo INFO: TASK ID $SLURM_ARRAY_TASK_ID',
            'echo INFO: JOB ID $SLURM_JOBID',
            'echo INFO: WORKDIR $DIR',
            'echo INFO: HOSTNAME `hostname`',
            'echo INFO: CUDA DEVICE $CUDA_VISIBLE_DEVICES',
            '',
            'bash ../run.sh',
            '',
            'cd ..',
            ]
        cmd = '\n'.join(cmd_list)
        ctools.gen_slurm_bash('submit.slurm.sh', cmd, slurm_params)
        jobid = ''
        if not args.dry_run:
            shcmd = 'sbatch submit.slurm.sh'
            _, out, _ = ctools.command_caller(command=shcmd, shell=True)
            jobid = re.findall(r"job (.*?)\n", str(out))
            with open('slurm_jobid.log', 'w') as f:
                f.writelines(f'INFO: ARRAY JOB ID {jobid[0]}')

        if args.remd:
            slurm_params['job-name'] = args.ligand_name+'.vdw.remd.md'
            slurm_params['output'] = 'slurm_remd.out'
            slurm_params['error'] = 'slurm_remd.err'
            slurm_params['ntasks-per-node'] = args.gpus_per_node
            slurm_params['nodes'] = ceil(float(len(_lambdas))/args.gpus_per_node)
            slurm_params['gpus'] = args.gpus_per_node
            del slurm_params['array']
            if not args.dry_run:
                slurm_params['dependency'] = f'afterok:{jobid[0]}'
            ctools.gen_slurm_bash('submit.slurm.remd.sh', './run.remd.sh', slurm_params)
            if not args.dry_run:
                shcmd = 'sbatch submit.slurm.remd.sh'
                _, out, _ = ctools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                with open('slurm_jobid.log', 'w') as f:
                    f.writelines(f'INFO: ARRAY JOB ID {jobid[0]}')
    if not args.hfe_only:
        #==========================================================================
        # 2. COMPLEX
        #==========================================================================
        alchemy_morph = os.path.join(args.abfe_workdir, '_alchemy_morph', args.complex_name)
        morph_pkl_file = os.path.join(alchemy_morph, 'objects.pkl')
        error_msg = "Please ensure that each complex has one objects.pkl file!"
        assert os.path.isfile(morph_pkl_file)==True, error_msg

        comp = pickle.load(open(morph_pkl_file, "rb"))
        workdir = os.path.join(args.abfe_workdir, '_alchemy_md', args.complex_name)
        if comp.protein_type == 'soluble':
            from abfe.md.amber_alchemy_mdin import COMPLEX_MIN
            from abfe.md.amber_alchemy_mdin import COMPLEX_HEAT, COMPLEX_PRESS
            from abfe.md.amber_alchemy_mdin import COMPLEX_TI
            from abfe.md.amber_alchemy_mdin import COMPLEX_RESTRAINT_MIN
            from abfe.md.amber_alchemy_mdin import COMPLEX_RESTRAINT_HEAT, COMPLEX_RESTRAINT_PRESS
            from abfe.md.amber_alchemy_mdin import COMPLEX_RESTRAINT_TI
        elif comp.protein_type == 'membrane':
            from abfe.md.amber_alchemy_mdin_mem import COMPLEX_MIN
            from abfe.md.amber_alchemy_mdin_mem import COMPLEX_HEAT, COMPLEX_PRESS
            from abfe.md.amber_alchemy_mdin_mem import COMPLEX_TI
            from abfe.md.amber_alchemy_mdin_mem import COMPLEX_RESTRAINT_MIN
            from abfe.md.amber_alchemy_mdin_mem import COMPLEX_RESTRAINT_HEAT, COMPLEX_RESTRAINT_PRESS
            from abfe.md.amber_alchemy_mdin_mem import COMPLEX_RESTRAINT_TI
        else:
            print(f'Only support the protein type of [soluble, membrane].')
            print(f'The provided protein type is {comp.protein_type}.')
        #==========================================================================
        # 2.1 Analytic Boresch Restraint
        #==========================================================================
        _workdir = os.path.join(workdir, 'analytic')
        with ctools.DirManager(_workdir):
            last_parm7 = comp.equil_last_md_info['last_parm7']
            last_rst7 =  comp.equil_last_md_info['last_rst7']
            last_traj = comp.equil_last_md_info['last_traj']
            trajin = cpptraj_tools.SIX_DOFS_INPUT\
                .replace('PRMTOP', last_parm7)\
                .replace('MDCRD', last_traj)\
                .replace('L1', args.lig_restraint_atoms[0])\
                .replace('L2', args.lig_restraint_atoms[1])\
                .replace('L3', args.lig_restraint_atoms[2])\
                .replace('P1', args.rec_restraint_atoms[0])\
                .replace('P2', args.rec_restraint_atoms[1])\
                .replace('P3', args.rec_restraint_atoms[2])
            md_tools.write_mdin_file('traj.in', trajin)
            _cmd = ' '.join([os.path.join(const.AMBERBIN, 'cpptraj'), '-i traj.in'])
            os.system(_cmd)
            r0, r0_std = cpptraj_tools.mean_std('bnd_r.dat', isPeriodic=False)
            alpha0, alpha0_std = cpptraj_tools.mean_std('bnd_alpha.dat', isPeriodic=False)
            theta0, theta0_std = cpptraj_tools.mean_std('bnd_theta.dat', isPeriodic=False)
            gamma0, gamma0_std = cpptraj_tools.mean_std('bnd_gamma.dat', isPeriodic=True)
            beta0, beta0_std = cpptraj_tools.mean_std('bnd_beta.dat', isPeriodic=True)
            phi0, phi0_std = cpptraj_tools.mean_std('bnd_phi.dat', isPeriodic=True)
            with open('boresch_six_dofs_stats_info.dat', 'w') as _fin:
                _fin.write('%5s %9s %9s %9s %9s %9s %9s\n'%(
                           '', 'r', 'alpha', 'theta', 'gamma', 'beta', 'phi')
                           )
                _fin.write('%s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n'%(
                           'mean:', r0, alpha0, theta0, gamma0, beta0, phi0)
                           )
                _fin.write('%s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n'%(
                           ' std:', r0_std, alpha0_std, theta0_std, gamma0_std, beta0_std, phi0_std)
                           )
            temperature = args.temperature
            kr = args.bond_const
            ktheta = args.angle_const
            kphi = args.dihedral_const
            rst_correction = analytic.restraint_correction(
                T=temperature, r0=r0, theta0=theta0, kr=kr, ktheta=ktheta, kphi=kphi
                )
            results = np.array([[rst_correction]]).T
            df = pd.DataFrame(results, columns=['Total (kcal/mol)'])
            df.to_csv('results.csv', float_format='%.2f')
        #==========================================================================
        # 2.2 VDW
        #==========================================================================
        #==========================================================================
        # 2.2.1 input files
        #==========================================================================
        comp_md_steps = OrderedDict()
        _lambdas =[str(_lambda) for _lambda in args.abfe_complex_vdw_lambdas]
        rst_wts = [
            args.bond_const,
            args.angle_const, args.angle_const,
            args.dihedral_const, args.dihedral_const, args.dihedral_const
            ]
        rst_refs = [r0, alpha0, theta0, gamma0, beta0, phi0]
        for _lambda in _lambdas:
            if _lambda == _lambdas[0] or _lambda == _lambdas[-1]:
                ntwx = args.ntwx
            else:
                ntwx = 0.
            comp_md_steps['min'] = COMPLEX_MIN\
                .replace('CLAMBDA', _lambda)\
                .replace('CUT', str(args.cutoff))\
                .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                .replace('TIMASK1', str(comp.timask1))\
                .replace('TIMASK2', str(comp.timask2))\
                .replace('SCMASK1', str(comp.scmask1))\
                .replace('SCMASK2', str(comp.scmask2))\
                .replace('RESTRAINT_WT', str(args.restraint_wt))\
                .replace('RESTRAINTMASK', str(args.restraintmask))
            for n, (tempi, temp0) in enumerate(zip(args.heat_temps[:-1], args.heat_temps[1:]), 1):
                comp_md_steps['heat-%s'%str(n)] = COMPLEX_HEAT\
                    .replace('CLAMBDA', str(_lambda))\
                    .replace('DT', str(0.002))\
                    .replace('CUT', str(args.cutoff))\
                    .replace('TEMPI', str(tempi))\
                    .replace('TEMP0', str(temp0))\
                    .replace('TEMP_0', 'TEMP0')\
                    .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                    .replace('TIMASK1', str(comp.timask1))\
                    .replace('TIMASK2', str(comp.timask2))\
                    .replace('SCMASK1', str(comp.scmask1))\
                    .replace('SCMASK2', str(comp.scmask2))\
                    .replace('RESTRAINT_WT', str(args.restraint_wt))\
                    .replace('RESTRAINTMASK', str(args.restraintmask))
                comp_md_steps['press-%s'%str(n)] = COMPLEX_PRESS\
                    .replace('CLAMBDA', str(_lambda))\
                    .replace('DT', str(0.002))\
                    .replace('CUT', str(args.cutoff))\
                    .replace('TEMP0', str(temp0))\
                    .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                    .replace('TIMASK1', str(comp.timask1))\
                    .replace('TIMASK2', str(comp.timask2))\
                    .replace('SCMASK1', str(comp.scmask1))\
                    .replace('SCMASK2', str(comp.scmask2))\
                    .replace('RESTRAINT_WT', str(args.restraint_wt))\
                    .replace('RESTRAINTMASK', str(args.restraintmask))
            if args.remd:
                NSTLIM = str(int(args.complex_length/args.timestep*1000./args.remd_numexchg))
            else:
                NSTLIM = str(int(args.complex_length/args.timestep*1000))
            comp_md_steps['ti'] = COMPLEX_TI\
                .replace('CLAMBDA', str(_lambda))\
                .replace('MBAR_STATES', str(len(_lambdas)))\
                .replace('MBAR_LAMBDA', ', '.join(_lambdas))\
                .replace('NSTLIM', NSTLIM)\
                .replace('DT', str(args.timestep))\
                .replace('CUT', str(args.cutoff))\
                .replace('TEMPI', str(args.temperature))\
                .replace('TEMP0', str(args.temperature))\
                .replace('NTWX', str(int(ntwx/args.timestep)))\
                .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                .replace('TIMASK1', str(comp.timask1))\
                .replace('TIMASK2', str(comp.timask2))\
                .replace('SCMASK1', str(comp.scmask1))\
                .replace('SCMASK2', str(comp.scmask2))
            if args.remd:
                comp_md_steps['ti'] = comp_md_steps['ti']\
                    .replace('!gremd_acyc =', 'gremd_acyc =')\
                    .replace('!numexchg =', 'numexchg =')\
                    .replace('NUMEXCHG', str(args.remd_numexchg))
            _workdir = os.path.join(workdir, 'vdw', str(_lambda))
            with ctools.DirManager(_workdir):
                os.system(f"ln -s -f {comp.alchemy_steps['vdw'][0]} .")
                os.system(f"ln -s -f {comp.alchemy_steps['vdw'][1]} .")
                for step, protocol in comp_md_steps.items():
                    md_tools.write_mdin_file(f'{step}.in', protocol)
                _add_boresch_restraint(args=args, prmtop=comp.alchemy_steps['vdw'][0],
                                       lig_masks=args.lig_restraint_atoms,
                                       rec_masks=args.rec_restraint_atoms,
                                       lig_idx=0, rst_refs=rst_refs, rst_wts=rst_wts
                                       )
        if args.remd:
            _workdir = os.path.join(workdir, 'vdw')
            with ctools.DirManager(_workdir):
                _add_boresch_restraint(args=args, prmtop=comp.alchemy_steps['vdw'][0],
                                       lig_masks=args.lig_restraint_atoms,
                                       rec_masks=args.rec_restraint_atoms,
                                       lig_idx=0, rst_refs=rst_refs, rst_wts=rst_wts
                                       )
        #==========================================================================
        # 2.2.2 bash files (submit jobs scripts)
        #==========================================================================
        _workdir = os.path.join(workdir, 'vdw')
        with ctools.DirManager(_workdir):
            _mdsteps = list(comp_md_steps.keys())
            _fparm7 = comp.alchemy_steps['vdw'][0]
            _frst7 = comp.alchemy_steps['vdw'][1]
            if args.remd:
                group_cmd = md_tools.write_groupfile(_mdsteps, _fparm7, _lambdas)
                with open('run.remd.sh', 'w') as _fsh:
                    _fsh.write(group_cmd)
                os.chmod('run.remd.sh', stat.S_IRWXU)
                _mdsteps.remove('ti')
            _last_md_info = md_tools.write_alchemy_md_submit_sh(
                fmdsh='run.sh',
                mdsteps=_mdsteps,
                fparm7=_fparm7,
                frst7=_frst7,
                nprocess=args.cpus_per_task,
                dry_run=True,
                )
            slurm_params = {
            'job-name': args.complex_name+'.vdw.md',
            'output': 'slurm.out',
            'error': 'slurm.err',
            'partition': args.partition,
            'nodes': 1,
            'ntasks-per-node': 1,
            'cpus-per-task': 2,
            'time': '72:00:00',
            'gpujob': True,
            'gpus': 1,
            'array': '0-%s'%(str(len(_lambdas)-1)),
            }
            exclude = args.exclude
            if exclude.strip(): 
                slurm_params.update(exclude=exclude)
            cmd_list = [
                'DIR_ARR=( %s )'%' '.join(_lambdas),
                'DIR=${DIR_ARR[$SLURM_ARRAY_TASK_ID]}',
                'cd $DIR',
                'echo INFO: ARRAY JOB ID $SLURM_ARRAY_JOB_ID',
                'echo INFO: TASK ID $SLURM_ARRAY_TASK_ID',
                'echo INFO: JOB ID $SLURM_JOBID',
                'echo INFO: WORKDIR $DIR',
                'echo INFO: HOSTNAME `hostname`',
                'echo INFO: CUDA DEVICE $CUDA_VISIBLE_DEVICES',
                '',
                'bash ../run.sh',
                '',
                'cd ..',
                ]
            cmd = '\n'.join(cmd_list)
            ctools.gen_slurm_bash('submit.slurm.sh', cmd, slurm_params)
            jobid = ''
            if not args.dry_run:
                shcmd = 'sbatch submit.slurm.sh'
                _, out, _ = ctools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                with open('slurm_jobid.log', 'w') as f:
                    f.writelines(f'INFO: ARRAY JOB ID {jobid[0]}')

            if args.remd:
                slurm_params['job-name'] = args.complex_name+'.vdw.remd.md'
                slurm_params['output'] = 'slurm_remd.out'
                slurm_params['error'] = 'slurm_remd.err'
                slurm_params['ntasks-per-node'] = args.gpus_per_node
                slurm_params['nodes'] = ceil(float(len(_lambdas))/args.gpus_per_node)
                slurm_params['gpus'] = args.gpus_per_node
                del slurm_params['array']
                if not args.dry_run:
                    slurm_params['dependency'] = f'afterok:{jobid[0]}'
                ctools.gen_slurm_bash('submit.slurm.remd.sh', './run.remd.sh', slurm_params)
                if not args.dry_run:
                    shcmd = 'sbatch submit.slurm.remd.sh'
                    _, out, _ = ctools.command_caller(command=shcmd, shell=True)
                    jobid = re.findall(r"job (.*?)\n", str(out))
                    with open('slurm_jobid.log', 'w') as f:
                        f.writelines(f'INFO: ARRAY JOB ID {jobid[0]}')
        #==========================================================================
        # 2.3 RESTRAINT
        #==========================================================================
        #==========================================================================
        # 2.3.1 input files
        #==========================================================================
        rec_restraint_atoms = md_tools.gen_masks(org_masks=args.rec_restraint_atoms, shift=1)
        _lambdas =[str(_lambda) for _lambda in args.abfe_complex_restraint_lambdas]
        for _lambda in _lambdas:
            if _lambda == _lambdas[0] or _lambda == _lambdas[-1]:
                ntwx = args.ntwx
            else:
                ntwx = 0.
            comp_md_steps['min'] = COMPLEX_RESTRAINT_MIN\
                .replace('CLAMBDA', _lambda)\
                .replace('CUT', str(args.cutoff))\
                .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                .replace('RESTRAINT_WT', str(args.restraint_wt))\
                .replace('RESTRAINTMASK', str(args.restraintmask))
            for n, (tempi, temp0) in enumerate(zip(args.heat_temps[:-1], args.heat_temps[1:]), 1):
                comp_md_steps['heat-%s'%str(n)] = COMPLEX_RESTRAINT_HEAT\
                    .replace('CLAMBDA', str(_lambda))\
                    .replace('DT', str(0.002))\
                    .replace('CUT', str(args.cutoff))\
                    .replace('TEMPI', str(tempi))\
                    .replace('TEMP0', str(temp0))\
                    .replace('TEMP_0', 'TEMP0')\
                    .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                    .replace('RESTRAINT_WT', str(args.restraint_wt))\
                    .replace('RESTRAINTMASK', str(args.restraintmask))
                comp_md_steps['press-%s'%str(n)] = COMPLEX_RESTRAINT_PRESS\
                    .replace('CLAMBDA', str(_lambda))\
                    .replace('DT', str(0.002))\
                    .replace('CUT', str(args.cutoff))\
                    .replace('TEMP0', str(temp0))\
                    .replace('GTI_ADD_SC', str(args.gti_add_sc))\
                    .replace('RESTRAINT_WT', str(args.restraint_wt))\
                    .replace('RESTRAINTMASK', str(args.restraintmask))
            if args.remd:
                NSTLIM = str(int(args.complex_length/args.timestep*1000./args.remd_numexchg))
            else:
                NSTLIM = str(int(args.complex_length/args.timestep*1000))
            comp_md_steps['ti'] = COMPLEX_RESTRAINT_TI\
                .replace('CLAMBDA', str(_lambda))\
                .replace('MBAR_STATES', str(len(_lambdas)))\
                .replace('MBAR_LAMBDA', ', '.join(_lambdas))\
                .replace('NSTLIM', NSTLIM)\
                .replace('DT', str(args.timestep))\
                .replace('CUT', str(args.cutoff))\
                .replace('TEMPI', str(args.temperature))\
                .replace('TEMP0', str(args.temperature))\
                .replace('NTWX', str(int(ntwx/args.timestep)))
            if args.remd:
                comp_md_steps['ti'] = comp_md_steps['ti']\
                    .replace('!gremd_acyc =', 'gremd_acyc =')\
                    .replace('!numexchg =', 'numexchg =')\
                    .replace('NUMEXCHG', str(args.remd_numexchg))
            _workdir = os.path.join(workdir, 'restraint', str(_lambda))
            with ctools.DirManager(_workdir):
                os.system(f"ln -s -f {comp.alchemy_steps['restraint'][0]} .")
                os.system(f"ln -s -f {comp.alchemy_steps['restraint'][1]} .")
                for step, protocol in comp_md_steps.items():
                    md_tools.write_mdin_file(f'{step}.in', protocol)
                _add_boresch_restraint(args=args, prmtop=comp.alchemy_steps['restraint'][0],
                                       lig_masks=args.lig_restraint_atoms,
                                       rec_masks=rec_restraint_atoms,
                                       lig_idx=0, rst_refs=rst_refs, rst_wts=rst_wts
                                       )
        if args.remd:
            _workdir = os.path.join(workdir, 'restraint')
            with ctools.DirManager(_workdir):
                _add_boresch_restraint(args=args, prmtop=comp.alchemy_steps['restraint'][0],
                                       lig_masks=args.lig_restraint_atoms,
                                       rec_masks=rec_restraint_atoms,
                                       lig_idx=0, rst_refs=rst_refs, rst_wts=rst_wts
                                       )
        #==========================================================================
        # 2.3.2 bash files (submit jobs scripts)
        #==========================================================================
        _workdir = os.path.join(workdir, 'restraint')
        with ctools.DirManager(_workdir):
            _mdsteps = list(comp_md_steps.keys())
            _fparm7 = comp.alchemy_steps['restraint'][0]
            _frst7 = comp.alchemy_steps['restraint'][1]
            if args.remd:
                group_cmd = md_tools.write_groupfile(_mdsteps, _fparm7, _lambdas)
                with open('run.remd.sh', 'w') as _fsh:
                    _fsh.write(group_cmd)
                os.chmod('run.remd.sh', stat.S_IRWXU)
                _mdsteps.remove('ti')
            _last_md_info = md_tools.write_alchemy_md_submit_sh(
                fmdsh='run.sh',
                mdsteps=_mdsteps,
                fparm7=_fparm7,
                frst7=_frst7,
                nprocess=args.cpus_per_task,
                dry_run=True,
                )
            slurm_params = {
            'job-name': args.complex_name+'.restraint.md',
            'output': 'slurm.out',
            'error': 'slurm.err',
            'partition': args.partition,
            'nodes': 1,
            'ntasks-per-node': 1,
            'cpus-per-task': 2,
            'time': '72:00:00',
            'gpujob': True,
            'gpus': 1,
            'array': '0-%s'%(str(len(_lambdas)-1)),
            }
            exclude = args.exclude
            if exclude.strip(): 
                slurm_params.update(exclude=exclude)
            cmd_list = [
                'DIR_ARR=( %s )'%' '.join(_lambdas),
                'DIR=${DIR_ARR[$SLURM_ARRAY_TASK_ID]}',
                'cd $DIR',
                'echo INFO: ARRAY JOB ID $SLURM_ARRAY_JOB_ID',
                'echo INFO: TASK ID $SLURM_ARRAY_TASK_ID',
                'echo INFO: JOB ID $SLURM_JOBID',
                'echo INFO: WORKDIR $DIR',
                'echo INFO: HOSTNAME `hostname`',
                'echo INFO: CUDA DEVICE $CUDA_VISIBLE_DEVICES',
                '',
                'bash ../run.sh',
                '',
                'cd ..',
                ]
            cmd = '\n'.join(cmd_list)
            ctools.gen_slurm_bash('submit.slurm.sh', cmd, slurm_params)
            jobid = ''
            if not args.dry_run:
                shcmd = 'sbatch submit.slurm.sh'
                _, out, _ = ctools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                with open('slurm_jobid.log', 'w') as f:
                    f.writelines(f'INFO: ARRAY JOB ID {jobid[0]}')

            if args.remd:
                slurm_params['job-name'] = args.complex_name+'.restraint.remd.md'
                slurm_params['output'] = 'slurm_remd.out'
                slurm_params['error'] = 'slurm_remd.err'
                slurm_params['ntasks-per-node'] = args.gpus_per_node
                slurm_params['nodes'] = ceil(float(len(_lambdas))/args.gpus_per_node)
                slurm_params['gpus'] = args.gpus_per_node
                del slurm_params['array']
                if not args.dry_run:
                    slurm_params['dependency'] = f'afterok:{jobid[0]}'
                ctools.gen_slurm_bash('submit.slurm.remd.sh', './run.remd.sh', slurm_params)
                if not args.dry_run:
                    shcmd = 'sbatch submit.slurm.remd.sh'
                    _, out, _ = ctools.command_caller(command=shcmd, shell=True)
                    jobid = re.findall(r"job (.*?)\n", str(out))
                    with open('slurm_jobid.log', 'w') as f:
                        f.writelines(f'INFO: ARRAY JOB ID {jobid[0]}')
