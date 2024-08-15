#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the equilibration module that runs equilibrations of complex and
ligand in solvated environment.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import argparse
import pickle
import os
import socket
import time
from collections import OrderedDict

import abfe.md.md_tools as md_tools
import abfe.utils.common_tools as ctools
from abfe.md.amber_mdin import *

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parse arguments for equilibration module!',
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
    # simulation arguments
    parser.add_argument(
        '-dr', '--dry_run', action='store_true', default=False,
        help='If turned on, dry run.'
        )
    parser.add_argument(
        '-cl', '--complex_length', default=0.5, type=float,
        help='Length (in ns) of the solvated complex equilibration'
        )
    parser.add_argument(
        '-sl', '--solvated_length', default=0.5, type=float,
        help='Length (in ns) of the solvated ligand equilibration'
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

    args = parser.parse_args()

    args.heat_temps[-1] = args.temperature

    return args

if __name__ == '__main__':
    """Parse the input arguments, then run the equilibration stage.
    """
    args = _parse_args()
    print(time.strftime("%c"))
    print(f"Equilibration started on host {socket.gethostname()}.")
    #==========================================================================
    # Equilibrate the solvated ligand
    #==========================================================================
    topo_workdir = os.path.join(args.abfe_workdir, '_topology', args.ligand_name)
    topo_pkl_file = os.path.join(topo_workdir, 'objects.pkl')
    error_msg = "Please ensure that each ligand has one objects.pkl file!"
    assert os.path.isfile(topo_pkl_file)==True, error_msg

    lig = pickle.load(open(topo_pkl_file, "rb"))
    workdir = os.path.join(args.abfe_workdir, '_equilibration', args.ligand_name)
    sh_filename = 'submit.sh'
    lig_equil_steps = OrderedDict()
    lig_equil_steps['min-1'] = SOLVATED_MIN_1\
        .replace('CUT', str(args.cutoff))\
        .replace('RESTRAINT_WT', str(args.restraint_wt))\
        .replace('RESTRAINTMASK', str(args.restraintmask))
    lig_equil_steps['min-2'] = SOLVATED_MIN_2\
        .replace('CUT', str(args.cutoff))
    for n, (tempi, temp0) in enumerate(zip(args.heat_temps[:-1], args.heat_temps[1:]), 1):
        lig_equil_steps['heat-%s'%str(n)] = SOLVATED_HEAT\
            .replace('DT', str(0.002))\
            .replace('CUT', str(args.cutoff))\
            .replace('TEMPI', str(tempi))\
            .replace('TEMP0', str(temp0))\
            .replace('TEMP_0', 'TEMP0')\
            .replace('RESTRAINT_WT', str(args.restraint_wt))\
            .replace('RESTRAINTMASK', str(args.restraintmask))
        lig_equil_steps['press-%s'%str(n)] = SOLVATED_PRESS\
            .replace('DT', str(0.002))\
            .replace('CUT', str(args.cutoff))\
            .replace('TEMP0', str(temp0))\
            .replace('RESTRAINT_WT', str(args.restraint_wt))\
            .replace('RESTRAINTMASK', str(args.restraintmask))
    lig_equil_steps['relax'] = SOLVATED_RELAX\
        .replace('NSTLIM', str(int(args.solvated_length/args.timestep*1000)))\
        .replace('DT', str(args.timestep))\
        .replace('CUT', str(args.cutoff))\
        .replace('TEMPI', str(args.temperature))\
        .replace('TEMP0', str(args.temperature))
    with ctools.DirManager(workdir):
        for step, protocol in lig_equil_steps.items():
            md_tools.write_mdin_file(f'{step}.in', protocol)
        os.system(f"ln -s -f {lig.parm7_file} .")
        os.system(f"ln -s -f {lig.rst7_file} .")
        last_md_info = md_tools.write_md_submit_sh(
            fmdsh=sh_filename,
            mdsteps=list(lig_equil_steps.keys()),
            fparm7=lig.parm7_file,
            frst7=lig.rst7_file,
            nprocess=args.cpus_per_task,
            dry_run=args.dry_run
            )
        if not last_md_info:
            print(time.strftime("%c"))
            error_msg = f"Equilibration for {lig} failed."
            raise RuntimeError(error_msg)
        else:
            print(time.strftime("%c"))
            print(f"Equilibration for {lig} finished.")

    lig.equil_last_md_info = last_md_info
    lig.equil_dir = os.path.abspath(workdir)
    pkl_file = os.path.join(workdir, 'objects.pkl')
    lig.pkl_file = os.path.abspath(pkl_file)
    # dump the pkl files
    with open(lig.pkl_file, 'wb') as pickf:
        pickle.dump(lig, pickf)

    if not args.hfe_only:
        #==========================================================================
        # Equilibrate the solvated complex
        #==========================================================================
        topo_workdir = os.path.join(args.abfe_workdir, '_topology', args.complex_name)
        topo_pkl_file = os.path.join(topo_workdir, 'objects.pkl')
        error_msg = "Please ensure that each complex has one objects.pkl file!"
        assert os.path.isfile(topo_pkl_file)==True, error_msg
    
        comp = pickle.load(open(topo_pkl_file, "rb"))
        workdir = os.path.join(args.abfe_workdir, '_equilibration', args.complex_name)
        comp.equil_dir = os.path.abspath(workdir)
        sh_filename = 'submit.sh'
        comp_equil_steps = OrderedDict()
        comp_equil_steps['min-1'] = COMPLEX_MIN_1\
            .replace('CUT', str(args.cutoff))\
            .replace('RESTRAINT_WT', str(args.restraint_wt))\
            .replace('RESTRAINTMASK', str(args.restraintmask))
        comp_equil_steps['min-2'] = COMPLEX_MIN_2\
            .replace('CUT', str(args.cutoff))
        for n, (tempi, temp0) in enumerate(zip(args.heat_temps[:-1], args.heat_temps[1:]), 1):
            comp_equil_steps['heat-%s'%str(n)] = COMPLEX_HEAT\
                .replace('DT', str(0.002))\
                .replace('CUT', str(args.cutoff))\
                .replace('TEMPI', str(tempi))\
                .replace('TEMP0', str(temp0))\
                .replace('TEMP_0', 'TEMP0')\
                .replace('RESTRAINT_WT', str(args.restraint_wt))\
                .replace('RESTRAINTMASK', str(args.restraintmask))
            comp_equil_steps['press-%s'%str(n)] = COMPLEX_PRESS\
                .replace('DT', str(0.002))\
                .replace('CUT', str(args.cutoff))\
                .replace('TEMP0', str(temp0))\
                .replace('RESTRAINT_WT', str(args.restraint_wt))\
                .replace('RESTRAINTMASK', str(args.restraintmask))
        comp_equil_steps['relax'] = COMPLEX_RELAX\
            .replace('NSTLIM', str(int(args.complex_length/args.timestep*1000)))\
            .replace('DT', str(args.timestep))\
            .replace('CUT', str(args.cutoff))\
            .replace('TEMPI', str(args.temperature))\
            .replace('TEMP0', str(args.temperature))
        with ctools.DirManager(workdir):
            for step, protocol in comp_equil_steps.items():
                md_tools.write_mdin_file(f'{step}.in', protocol)
            os.system(f"ln -s -f {comp.parm7_file} .")
            os.system(f"ln -s -f {comp.rst7_file} .")
            last_md_info = md_tools.write_md_submit_sh(
                fmdsh=sh_filename,
                mdsteps=list(comp_equil_steps.keys()),
                fparm7=comp.parm7_file,
                frst7=comp.rst7_file,
                nprocess=args.cpus_per_task,
                dry_run=args.dry_run
                )
            if not last_md_info:
                print(time.strftime("%c"))
                error_msg = f"Equilibration for {comp} failed."
                raise RuntimeError(error_msg)
            else:
                print(time.strftime("%c"))
                print(f"Equilibration for {comp} finished.")
    
            comp.equil_last_md_info = last_md_info
            comp.equil_dir = os.path.abspath(workdir)
            pkl_file = os.path.join(workdir, 'objects.pkl')
            comp.pkl_file = os.path.abspath(pkl_file)
            with open(comp.pkl_file, 'wb') as pickf:
                pickle.dump(comp, pickf)
            lig.complex = comp
            with open(lig.pkl_file, 'wb') as pickf2:
                pickle.dump(lig, pickf2)