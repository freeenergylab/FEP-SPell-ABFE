#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
===============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the alchemy analysis module that does alchemical analysis.

@Author: Pengfei Li
@Date: Dec. 18, 2023
===============================================================================
"""
import argparse
import os
import socket
import time
import numpy as np
import pandas as pd
import pickle
from collections import defaultdict
from glob import glob

from abfe.md.charge_correction import align_complex
from abfe.md.charge_correction import compute_charge_correction
from abfe.utils.check_md_jobs import JobManager
from abfe.utils.common_tools import DirManager
from abfe.utils.alchemlyb_dUoverlap import alchemlyb_analysis

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parse arguments for alchemy analysis module!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    # construction arguments
    parser.add_argument(
        '-bw', '--abfe_workdir',
        help='Workdir'
        ) 
    parser.add_argument(
        '-ln', '--ligand_name',
        help='Str contains the ligand name'
        )
    parser.add_argument(
        '-cn', '--complex_name',
        help='Str contains the complex name'
        )
    parser.add_argument(
        '-ho', '--hfe_only', action='store_true', default=False,
        help='Only perform HFE'
        )
    # input arguments
    parser.add_argument(
        '-m', '--methods', default=['MBAR', 'TI'], type=str, nargs='*',
        help='Alchemical analysis methods.'
        )
    parser.add_argument(
        '--temperature', type=float, default=298.15,
        help='System temperature (in K)'
        )
    parser.add_argument(
        '-sm', '--solvent_mask', type=str, default=':WAT,K+,Na+,Cl-',
        help='Str contains the non-solute molecule mask'
        )
    parser.add_argument(
        '-es', '--epsilon_solv', type=float, default=97.,
        help='The relative permittivity for TIP3P water'
        )

    return parser.parse_args()

if __name__ == '__main__':
    """Parse the input arguments, then run the alchemy analysis stage.
    """
    args = _parse_args()
    print(time.strftime("%c"))
    print(f"Alchemy analysis started on host {socket.gethostname()}.")

    abfe_workdir = args.abfe_workdir
    lig_name = args.ligand_name
    comp_name = args.complex_name
    methods = args.methods
    T = args.temperature
    solvent_mask = args.solvent_mask
    epsilon_solv = args.epsilon_solv
    #=========================================================================
    # Analyze solvated/vdw
    #=========================================================================
    solv_check_obj = JobManager(workdir=abfe_workdir, sys_name=lig_name, stage='.', step='vdw')
    solv_check_obj.seek_jobid()
    if solv_check_obj.flag == True:
        data_dir = os.path.join(abfe_workdir, '_alchemy_md', lig_name, 'vdw')
        save_dir = os.path.join(abfe_workdir, '_alchemy_analysis', lig_name, 'vdw')
        ti_out_files = glob(os.path.join(data_dir, '[01]*/ti*.out'))
        ti_out_files.sort()
        with DirManager(save_dir):
            df_solvated_vdw = alchemlyb_analysis(ti_out_files, methods, T, units='kcal/mol')
            with open('results.csv', 'w') as f:
                df_solvated_vdw.to_csv(f)
    else:
        print(f'Failed in doing alchemlyb analysis for solvated/vdw.')

    # Doing charge correction for charged ligand.
    charge_correction_dir=defaultdict(list)
    equil_workdir = os.path.join(abfe_workdir, '_equilibration', lig_name)
    equil_pkl_file = os.path.join(equil_workdir, 'objects.pkl')
    lig = pickle.load(open(equil_pkl_file, "rb"))
    if lig.net_charge and solv_check_obj.flag == True:
        solv_parm=os.path.join(abfe_workdir, '_alchemy_morph', lig_name, 'vdw', 'solvated_vdw.parm7')
        solv_crd=os.path.join(abfe_workdir, '_alchemy_md', lig_name, 'vdw', '0.0', 'ti.nc')
        ## solvated leg
        _dir = os.path.dirname(solv_crd)
        charge_correction_dir[lig_name].append(_dir)
        with DirManager(_dir):
            new_mdcrd = align_complex(solv_parm, solv_crd, solvent_mask)
            dG_correction = compute_charge_correction(
                prmtop=solv_parm, mdcrd=new_mdcrd, lig_rname='MOL',
                solvent_mask=solvent_mask, temperature=T,
                epsilon_solv=epsilon_solv, wat_rname='WAT'
                )
            with open('results.csv', 'w') as f:
                dG_correction.to_csv(f, float_format='%.3f')
    elif not lig.net_charge:
        print(f'Ligand {lig_name} is neutral: do not need to do charge correction for solvated leg.')
    else:
        print(f'Charge correction for solvated leg failed, due to failure in solvated/vdw md.')

    if not args.hfe_only:
        #=========================================================================
        # Analyze complex/vdw
        #=========================================================================
        comp_check_obj = JobManager(workdir=abfe_workdir, sys_name=comp_name, stage='.', step='vdw')
        comp_check_obj.seek_jobid()
        if comp_check_obj.flag == True:
            data_dir = os.path.join(abfe_workdir, '_alchemy_md', comp_name, 'vdw')
            save_dir = os.path.join(abfe_workdir, '_alchemy_analysis', comp_name, 'vdw')
            ti_out_files = glob(os.path.join(data_dir, '[01]*/ti*.out'))
            ti_out_files.sort()
            with DirManager(save_dir):
                df_complex_vdw = alchemlyb_analysis(ti_out_files, methods, T, units='kcal/mol')
                with open('results.csv', 'w') as f:
                    df_complex_vdw.to_csv(f)
        else:
            print(f'Failed in doing alchemlyb analysis for complex/vdw.')

        # Doing charge correction for charged ligand.
        if lig.net_charge and comp_check_obj.flag == True:
            ## complex leg
            comp_parm=os.path.join(abfe_workdir, '_alchemy_morph', comp_name, 'vdw', 'complex_vdw.parm7')
            comp_crd=os.path.join(abfe_workdir, '_alchemy_md', comp_name, 'vdw', '0.0', 'ti.nc')
            _dir = os.path.dirname(comp_crd)
            charge_correction_dir[lig_name].append(_dir)
            with DirManager(_dir):
                new_mdcrd = align_complex(comp_parm, comp_crd, solvent_mask)
                dG_correction = compute_charge_correction(
                    prmtop=comp_parm, mdcrd=new_mdcrd, lig_rname='MOL',
                    solvent_mask=solvent_mask, temperature=T,
                    epsilon_solv=epsilon_solv, wat_rname='WAT'
                    )
                with open('results.csv', 'w') as f:
                    dG_correction.to_csv(f, float_format='%.3f')
        elif not lig.net_charge:
            print(f'Ligand {lig_name} is neutral: do not need to do charge correction for complex leg.')
        else:
            print(f'Charge correction for complex leg failed, due to failure in complex/vdw md.')

        #======================================================================
        # Analyze complex/restraint
        #======================================================================
        comp_rst_check_obj = JobManager(workdir=abfe_workdir, sys_name=comp_name, stage='.', step='restraint')
        comp_rst_check_obj.seek_jobid()
        if comp_rst_check_obj.flag == True:
            data_dir = os.path.join(abfe_workdir, '_alchemy_md', comp_name, 'restraint')
            save_dir = os.path.join(abfe_workdir, '_alchemy_analysis', comp_name, 'restraint')
            ti_out_files = glob(os.path.join(data_dir, '[01]*/ti*.out'))
            ti_out_files.sort()
            with DirManager(save_dir):
                df_complex_restraint = alchemlyb_analysis(ti_out_files, methods, T, units='kcal/mol')
                with open('results.csv', 'w') as f:
                    df_complex_restraint.to_csv(f)
        else:
            print(f'Failed in doing alchemlyb analysis for complex/restraint.')

    #=========================================================================
    # Summarize the final results.
    #=========================================================================
    dG_components_result=pd.DataFrame()
    dG_result=pd.DataFrame()

    _dir = os.path.join(abfe_workdir, '_alchemy_analysis')
    f_solvated_vdw = os.path.join(_dir, lig_name, 'vdw', 'results.csv')
    df_solvated_vdw = pd.read_csv(f_solvated_vdw, index_col=[0])
    if lig_name in charge_correction_dir.keys():
        f_solvated_cc = os.path.join(charge_correction_dir[lig_name][0], 'results.csv')
        df_solvated_cc = pd.read_csv(f_solvated_cc, index_col=[0])
    else:
        df_solvated_cc = pd.DataFrame()
        df_solvated_cc.loc[0, 'Total (kcal/mol)'] = 0.0
    for method in methods:
        dG=-1.0*(df_solvated_vdw.loc[method, 'dF']-df_solvated_cc.loc[0, 'Total (kcal/mol)'])
        dG_std=df_solvated_vdw.loc[method, 'ddF']
        dG_components_result.loc[lig_name, f'{method}_dG_solvated_vdw']=df_solvated_vdw.loc[method, 'dF']
        dG_components_result.loc[lig_name, f'{method}_dG_solvated_cc']=df_solvated_cc.loc[0, 'Total (kcal/mol)']
        dG_result.loc[lig_name, method+'_dG']=dG
        dG_result.loc[lig_name, method+'_dG_std']=dG_std

    if not args.hfe_only:
        f_complex_vdw = os.path.join(_dir, comp_name, 'vdw', 'results.csv')
        f_complex_restraint = os.path.join(_dir, comp_name, 'restraint', 'results.csv')
        from_dir = os.path.join(abfe_workdir, '_alchemy_md', comp_name, 'analytic')
        with DirManager(os.path.join(_dir, comp_name)):
            os.system(f'cp -r {from_dir} .')
        f_analytic = os.path.join(_dir, comp_name, 'analytic', 'results.csv')

        df_complex_vdw = pd.read_csv(f_complex_vdw, index_col=[0])
        df_complex_restraint = pd.read_csv(f_complex_restraint, index_col=[0])
        df_analytic = pd.read_csv(f_analytic, index_col=[0])
        if lig_name in charge_correction_dir.keys():
            f_solvated_cc = os.path.join(charge_correction_dir[lig_name][0], 'results.csv')
            df_solvated_cc = pd.read_csv(f_solvated_cc, index_col=[0])
            f_complex_cc = os.path.join(charge_correction_dir[lig_name][1], 'results.csv')
            df_complex_cc = pd.read_csv(f_complex_cc, index_col=[0])
        else:
            df_solvated_cc = pd.DataFrame()
            df_complex_cc = pd.DataFrame()
            df_solvated_cc.loc[0, 'Total (kcal/mol)'] = 0.0
            df_complex_cc.loc[0, 'Total (kcal/mol)'] = 0.0
        for method in methods:
            ## dG(final) = (dG(ligand)-dG_charge_correction(ligand)]-[dG(complex)-dG_charge_correction(complex)]+dG(restraint)+dG(Translational)+dG(Rotational)
            dG = (df_solvated_vdw.loc[method, 'dF'] - df_solvated_cc.loc[0, 'Total (kcal/mol)']) \
                - (df_complex_vdw.loc[method, 'dF'] - df_complex_cc.loc[0, 'Total (kcal/mol)']) \
                + df_complex_restraint.loc[method, 'dF'] \
                + df_analytic.loc[0, 'Total (kcal/mol)']
            dG_std = np.sqrt(df_solvated_vdw.loc[method, 'ddF']**2 \
                + df_complex_vdw.loc[method, 'ddF']**2 \
                + df_complex_restraint.loc[method, 'ddF']**2
                )
            dG_components_result.loc[lig_name, f'{method}_dG_solvated_vdw']=df_solvated_vdw.loc[method, 'dF']
            dG_components_result.loc[lig_name, f'{method}_dG_solvated_cc']=df_solvated_cc.loc[0, 'Total (kcal/mol)']
            dG_components_result.loc[lig_name, f'{method}_dG_complex_vdw']=df_complex_vdw.loc[method, 'dF']
            dG_components_result.loc[lig_name, f'{method}_dG_complex_cc']=df_complex_cc.loc[0, 'Total (kcal/mol)']
            dG_components_result.loc[lig_name, f'{method}_dG_complex_restraint']=df_complex_restraint.loc[method, 'dF']
            dG_components_result.loc[lig_name, f'{method}_dG_analytic']=df_analytic.loc[0, 'Total (kcal/mol)']
            dG_result.loc[lig_name, method+'_dG'] = dG
            dG_result.loc[lig_name, method+'_dG_std'] = dG_std

    with DirManager(abfe_workdir):
        dG_components_result.to_csv('dG_components_results.csv', float_format='%.2f', mode='a', header=True)
        dG_result.to_csv('dG_results.csv', float_format='%.2f', mode='a', header=True)