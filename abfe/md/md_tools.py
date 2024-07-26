#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
===============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the md tools module that can be used to help construct AMBER simulation
input protocol.

@Author: Pengfei Li
@Date: May 30th, 2022
===============================================================================
"""
import os
import re
import stat

import abfe.const as const

def write_mdin_file(fmdin, protocol):
    """Write md cotrol parameter protocol file.

    Args:
        fmdin (str): A mdin file name.
        protocol (str): One line simulation protocol.
    """
    with open(fmdin, 'w') as _fmdin:
        _fmdin.write(protocol)

def write_md_submit_sh(fmdsh, mdsteps, fparm7, frst7, dry_run=False):
    """Write md submit bash file.

    Args:
        fmdsh (str): A md bash file name.
        mdsteps (list): A list contains several md steps.
        fparm7 (str): initial parm7 filename
        frst7 (str): initial rst7 filename
        dry_run (bool): False: submit; True: only slurm bash script
    """
    header = ['#!/usr/bin/env bash']
    with open(fmdsh, 'w') as _fmdsh:
        _fmdsh.write('\n'.join(header) + '\n')
        for j, step in enumerate(mdsteps):
            _fmdsh.write('\n')
            if j==0 and step.startswith('min'):
                min_cmd = [
                    os.path.join(const.AMBERBIN, 'pmemd'),
                    '-O', f'-i {step}.in', f'-p {fparm7}',
                    f'-c {frst7}', f'-ref {frst7}',
                    f'-o {step}.out', f'-r {step}.rst7',
                    f'-inf {step}.info'
                    ]
                _fmdsh.write(' '.join(min_cmd) + '\n')
            elif j==1 and step.startswith('min'):
                min_cmd = [
                    os.path.join(const.AMBERBIN, 'pmemd'),
                    '-O', f'-i {step}.in', f'-p {fparm7}',
                    f'-c {mdsteps[j-1]}.rst7', f'-ref {mdsteps[j-1]}.rst7',
                    f'-o {step}.out', f'-r {step}.rst7',
                    f'-inf {step}.info'
                    ]
                _fmdsh.write(' '.join(min_cmd) + '\n')
            else:
                if step.startswith('press-1'):
                    md_engine = os.path.join(const.AMBERBIN, 'pmemd')
                else:
                    md_engine = os.path.join(const.AMBERBIN, 'pmemd.cuda')
                md_cmd = [
                    md_engine,
                    '-O', f'-i {step}.in', f'-p {fparm7}',
                    f'-c {mdsteps[j-1]}.rst7', f'-ref {mdsteps[j-1]}.rst7',
                    f'-o {step}.out', f'-r {step}.rst7', f'-x {step}.nc',
                    f'-inf {step}.info'
                    ]
                _fmdsh.write(' '.join(md_cmd) + '\n')

    last_parm7 = os.path.abspath(f'{fparm7}')
    last_rst7 = os.path.abspath(f'{mdsteps[-1]}.rst7')
    last_traj = os.path.abspath(f'{mdsteps[-1]}.nc')
    last_md_info = {'last_parm7': last_parm7, 'last_rst7': last_rst7, 'last_traj': last_traj}
    try:
        os.chmod(fmdsh, stat.S_IRWXU)
        if dry_run:
            status = 'dry_run'
        else: 
            os.system(f'bash {fmdsh}')
            status = 'run'
        last_md_info['status'] = status
        return last_md_info
    except Exception:
        return False

def write_alchemy_md_submit_sh(fmdsh, mdsteps, fparm7, frst7, dry_run=False):
    """Write alchemy md submit bash file.

    Args:
        fmdsh (str): A md bash file name.
        mdsteps (list): A list contains several md steps.
        fparm7 (str): initial parm7 filename
        frst7 (str): initial rst7 filename
        dry_run (bool): False: submit; True: only slurm bash script
    """
    header = ['#!/usr/bin/env bash']
    with open(fmdsh, 'w') as _fmdsh:
        _fmdsh.write('\n'.join(header) + '\n')
        for j, step in enumerate(mdsteps):
            _fmdsh.write('\n')
            if j==0 and step.startswith('min'):
                md_cmd = [
                    os.path.join(const.AMBERBIN, 'pmemd'),
                    '-O', f'-i {step}.in', f'-p {fparm7}',
                    f'-c {frst7}', f'-ref {frst7}',
                    f'-o {step}.out', f'-r {step}.rst7', f'-x {step}.nc',
                    f'-inf {step}.info'
                    ]
                _fmdsh.write(' '.join(md_cmd) + '\n')
            else:
                md_cmd = [
                    os.path.join(const.AMBERBIN, 'pmemd.cuda'),
                    '-O', f'-i {step}.in', f'-p {fparm7}',
                    f'-c {mdsteps[j-1]}.rst7', f'-ref {mdsteps[j-1]}.rst7',
                    f'-o {step}.out', f'-r {step}.rst7', f'-x {step}.nc',
                    f'-inf {step}.info'
                    ]
                _fmdsh.write(' '.join(md_cmd) + '\n')

    last_parm7 = os.path.abspath(f'{fparm7}')
    last_rst7 = os.path.abspath(f'{mdsteps[-1]}.rst7')
    last_traj = os.path.abspath(f'{mdsteps[-1]}.nc')
    last_md_info = {'last_parm7': last_parm7, 'last_rst7': last_rst7, 'last_traj': last_traj}
    try:
        os.chmod(fmdsh, stat.S_IRWXU)
        if dry_run:
            status = 'dry_run'
        else: 
            os.system(f'bash {fmdsh}')
            status = 'run'
        last_md_info['status'] = status
        return last_md_info
    except Exception:
        return False

def gen_masks(org_masks, shift):
    """Shift Amber masks.

    Args:
        org_masks (list): a list of Amber masks
        shift (int): shift the residue number index

    Returns:
        new_masks (list): a list of index-shifted Amber masks
    """    
    new_masks = []
    for org_mask in org_masks:
        org_tmp=re.search(r':(.*?)@', org_mask).group(1)
        new_tmp=str(int(org_tmp)+shift)
        new_mask = org_mask.replace(org_tmp, new_tmp)
        new_masks.append(new_mask)
    return new_masks

def write_groupfile(mdsteps, fparm7, lambdas):
    """Generate groupfile for REMD simulations.

    Args:
        mdsteps (list): contains simulation steps
        fparm7 (str): initial parm7 filename
        lambdas (list): contains the replicas (lambdas here)

    Returns:
        ' '.join(cmd) (str): one command line
    """    
    step = mdsteps[-1]
    with open(f'{step}.groupfile', 'w') as _file:
        for _lambda in lambdas:
            md_cmd = ['-O', f'-i {_lambda}/{step}.in', f'-p {fparm7}',
                      f'-c {_lambda}/{mdsteps[-2]}.rst7',
                      f'-ref {_lambda}/{mdsteps[-2]}.rst7',
                      f'-o {_lambda}/{step}.out', f'-r {_lambda}/{step}.rst7',
                      f'-x {_lambda}/{step}.nc', f'-inf {_lambda}/{step}.info'
                      ]
            _file.write(' '.join(md_cmd) + '\n')
    cmd = ['mpirun', '--oversubscribe', '-np', str(len(lambdas)),
           os.path.join(const.AMBERBIN, 'pmemd.cuda.MPI'),
           '-rem 3', f'-remlog {step}.remd.log',
           '-ng', str(len(lambdas)), '-groupfile', f'{step}.groupfile'
           ]

    return ' '.join(cmd)