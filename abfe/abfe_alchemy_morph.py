#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the alchemy morph module that prepares input files for different
alchemical transformation pathway.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""

import argparse
from collections import defaultdict
import os
import parmed as pmd
import pickle
import pytraj as pt
import shutil
import socket
import sys
import time

import abfe.utils.common_tools as ctools
import abfe.const as const

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parse arguments for alchemy morph module!',
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
        '-hr', '--hmr', action='store_true', default=False,
        help='Toggles hydrogen mass repartitioning functionality to allow for larger timestep.'
        )
    parser.add_argument(
        '-hm', '--hmass', default=3.024,
        help='New hydrogen mass'
        )

    return parser.parse_args()

if __name__ == '__main__':
    """Parse the input arguments, then run the alchemy morph stage.
    """
    args = _parse_args()
    print(time.strftime("%c"))
    print(f"Alchemy morph started on host {socket.gethostname()}.")
    #==========================================================================
    # Morph ligand
    #==========================================================================
    equil_workdir = os.path.join(args.abfe_workdir, '_equilibration', args.ligand_name)
    equil_pkl_file = os.path.join(equil_workdir, 'objects.pkl')
    error_msg = "Please ensure that each ligand has one objects.pkl file!"
    assert os.path.isfile(equil_pkl_file)==True, error_msg

    lig = pickle.load(open(equil_pkl_file, "rb"))
    workdir = os.path.join(args.abfe_workdir, '_alchemy_morph', args.ligand_name)
    alchemy_steps = defaultdict(list)
    _workdir = os.path.join(workdir, 'vdw')
    with ctools.DirManager(_workdir):
        shutil.copy(lig.equil_last_md_info['last_parm7'], 'solvated_vdw.parm7')
        shutil.copy(lig.equil_last_md_info['last_rst7'], 'solvated_vdw.rst7')
        alchemy_steps['vdw'].append(os.path.abspath('solvated_vdw.parm7'))
        alchemy_steps['vdw'].append(os.path.abspath('solvated_vdw.rst7'))
        # alchemy_mdin setup for alchemy_md
        lig.timask1='":1"'
        lig.timask2='""'
        lig.scmask1='":1"'
        lig.scmask2='""'
        # Hydrogen Mass Repartitioning
        if args.hmr:
            prmtop = alchemy_steps['vdw'][0]
            mol = pmd.load_file(prmtop, xyz=None)
            print(time.strftime("%c"))
            print(f"Use HMR for {prmtop}.")
            os.rename(prmtop, prmtop.replace('.parm7', '.nohmr.parm7'))
            pmd.tools.actions.HMassRepartition(mol, args.hmass).execute()
            mol.save(alchemy_steps['vdw'][0], overwrite=True, format='amber')
            del mol # flush the buffer

    lig.alchemy_steps = alchemy_steps
    pkl_file = os.path.join(workdir, 'objects.pkl')
    lig.pkl_file = os.path.abspath(pkl_file)
    # dump the pkl files
    with open(lig.pkl_file, 'wb') as pickf:
        pickle.dump(lig, pickf)

    if not args.hfe_only:
        #==========================================================================
        # Morph complex
        #==========================================================================
        equil_workdir = os.path.join(args.abfe_workdir, '_equilibration', args.complex_name)
        equil_pkl_file = os.path.join(equil_workdir, 'objects.pkl')
        error_msg = "Please ensure that each ligand has one objects.pkl file!"
        assert os.path.isfile(equil_pkl_file)==True, error_msg
    
        comp = pickle.load(open(equil_pkl_file, "rb"))
        workdir = os.path.join(args.abfe_workdir, '_alchemy_morph', args.complex_name)
        alchemy_steps = defaultdict(list)
        _workdir = os.path.join(workdir, 'vdw')
        with ctools.DirManager(_workdir):
            shutil.copy(comp.equil_last_md_info['last_parm7'], 'complex_vdw.parm7')
            shutil.copy(comp.equil_last_md_info['last_rst7'], 'complex_vdw.rst7')
            alchemy_steps['vdw'].append(os.path.abspath('complex_vdw.parm7'))
            alchemy_steps['vdw'].append(os.path.abspath('complex_vdw.rst7'))
            # alchemy_mdin setup for alchemy_md
            comp.timask1='":1"'
            comp.timask2='""'
            comp.scmask1='":1"'
            comp.scmask2='""'
            # Hydrogen Mass Repartitioning
            if args.hmr:
                prmtop = alchemy_steps['vdw'][0]
                mol = pmd.load_file(prmtop, xyz=None)
                print(time.strftime("%c"))
                print(f"Use HMR for {prmtop}.")
                os.rename(prmtop, prmtop.replace('.parm7', '.nohmr.parm7'))
                pmd.tools.actions.HMassRepartition(mol, args.hmass).execute()
                mol.save(alchemy_steps['vdw'][0], overwrite=True, format='amber')
                del mol # flush the buffer

        _workdir = os.path.join(workdir, 'restraint')
        with ctools.DirManager(_workdir):
            prmtop = comp.equil_last_md_info['last_parm7']
            inpcrd = comp.equil_last_md_info['last_rst7']
            traj = pt.load(inpcrd, prmtop)
            _traj = pt.strip(traj, ':MOL')
            pt.write_traj('protein.pdb', _traj, overwrite=True)
            _traj = pt.strip(traj, '!:MOL')
            pt.write_traj('vacuum.mol2', _traj, overwrite=True)
            shutil.copy(comp.tleap_complex_in, 'tleap.complex.in')
            box_info = ctools.parse_pdb('protein.pdb')
            # box_info: X/Y/Z dimension box length, box alpha/beta/gamma angles
            ctools.parse_leapin(
                leapin='tleap.complex.in',
                box_info=box_info,
                leapout='tleap.in',
                duplicate=True
                )
            tleap_cmd = [f"{os.path.join(const.AMBERBIN, 'tleap')}", '-f', \
                'tleap.in', '> /dev/null 2>&1']
            tleap_status = os.system(' '.join(tleap_cmd))
            if tleap_status != 0:
                print(time.strftime("%c"))
                print("tleap failed for alchemy morph 'restraint' step!")
                sys.exit(1)
            # When using solvateOct (a truncated octahedron box),
            # 'set Unit box {x, y, z}' in tleap only affect the box dimension info
            # in parm7 and rst7 files. The box alpha/beta/gamms angles will be
            # accordingly set to 90 degrees (when x=y=z in truncated octahedron box),
            # which is wrong and will cause failed simualtions at next step.
            # Thus, we need to use ChBox to specify the box angles and also
            # re-specify the box dimensions from 'box_info' inheriting from
            # '_equilibration' step.
            #
            # Although ChBox only modifies the rst7 files, it should be enough.
            # Because the box info in parm7 ('%FLAG BOX_DIMENSIONS') are redundant.
            # Ths box info in rst7 should be enough. Please see the documentation:
            # http://ambermd.org/FileFormats.php#topo.cntrl
            #
            # BTW: The IFBOX in parm7 does not matter here. 
            ChBox_cmd = [f"{os.path.join(const.AMBERBIN, 'ChBox')}",
                        '-c complex.rst7 ', '-o complex.rst7', 
                        f'-al {box_info[4]}', f'-bt {box_info[5]}', f'-gm {box_info[6]}',
                        f'-X {box_info[1]}', f'-Y {box_info[2]}', f'-Z {box_info[3]}',
                        '> /dev/null 2>&1']
            print(f"ChBox command is {' '.join(ChBox_cmd)}.")
            ChBox_status = os.system(' '.join(ChBox_cmd))
            if ChBox_status != 0:
                print(time.strftime("%c"))
                print("ChBox failed for alchemy morph 'restraint' step!")
                sys.exit(1)
            os.system('mv complex.parm7 complex_restraint.parm7')
            os.system('mv complex.rst7 complex_restraint.rst7')
            alchemy_steps['restraint'].append(os.path.abspath('complex_restraint.parm7'))
            alchemy_steps['restraint'].append(os.path.abspath('complex_restraint.rst7'))
            # Hydrogen Mass Repartitioning
            if args.hmr:
                prmtop = alchemy_steps['restraint'][0]
                mol = pmd.load_file(prmtop, xyz=None)
                print(time.strftime("%c"))
                print(f"Use HMR for {prmtop}.")
                os.rename(prmtop, prmtop.replace('.parm7', '.nohmr.parm7'))
                pmd.tools.actions.HMassRepartition(mol, args.hmass).execute()
                mol.save(alchemy_steps['restraint'][0], overwrite=True, format='amber')
                del mol # flush the buffer
    
        comp.alchemy_steps = alchemy_steps
        pkl_file = os.path.join(workdir, 'objects.pkl')
        comp.pkl_file = os.path.abspath(pkl_file)
        # dump the pkl files
        with open(comp.pkl_file, 'wb') as pickf:
            pickle.dump(comp, pickf)