#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the main module that controls the whole abfe jobs workflow.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import argparse
import os
import re
import sys
import time
import warnings
from glob import glob

import abfe.const as const
import abfe.utils.common_tools as common_tools

def _parse_args():
    """Parse arguments from the config file.

    Returns:
        args (dict): contains the arguments from the config file
    """    
    parser = argparse.ArgumentParser(
        description='Parse one input configuration file!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        '-i', '--inconfig', 
        help='One configuration file contains all input parameters.'
        )
    pargs = parser.parse_args()
    args = common_tools.read_config(pargs.inconfig)
    args['main']['inconfig'] = os.path.abspath(pargs.inconfig)

    for stage, params in args.items():
        args[stage] = common_tools.DotDict(params)
    args = common_tools.DotDict(args)

    assert args.main.workflow_type in const.WORKFLOW_TYPES, \
        f"Cannot support the user-defined workflow type: {args.main.workflow_type}!"

    return args

class abfeWorkflow(object):
    """This class is designed to control the whole abfe jobs workflow.
    """
    def __init__(self, args):
        self.args=args
        self.hfe_only = args.main.workflow_type == 'ahfe'
        self._get_workdir()
        if not self.hfe_only:
            self._find_protein_dir()
            self._find_cofactor_dir()
        else:
            self.protein_dir = ''
            self.cofactor_dirs = []
            self.cofactor_names = []
            self.cofactor_charges = []
        self._find_ligand_dir()

    def _get_workdir(self):
        """Get working directory.
        """
        self.abfe_workdir=os.getcwd()

    def _find_protein_dir(self):
        """Find the protein directory.
        """
        protein_dir = os.path.join(self.abfe_workdir, 'proteins', self.args.topology.protein_name)
        if not os.path.isdir(protein_dir):
            print(f"The protein '{self.args.topology.protein_name}' dir is not found in 'proteins' dir")
            sys.exit(1)
        else:
            self.protein_dir = os.path.abspath(protein_dir)

    def _find_cofactor_dir(self):
        """Find the cofactor directory.
        """
        self.cofactor_dirs = []
        self.cofactor_names = []
        for name in self.args.topology.cofactor_names:
            _dir = os.path.join(self.abfe_workdir, 'cofactors', f'{name}')
            self.cofactor_dirs.append(_dir)
            self.cofactor_names.append(name)

    def _find_ligand_dir(self):
        """Find the ligand directory.
        """
        self.ligand_names = []
        _ligandsin = os.path.join(self.abfe_workdir, self.args.main.ligands_in)
        with open(_ligandsin, 'r') as ligandsin:
            for n, line in enumerate(ligandsin.readlines()):
                try:
                    line = re.sub(r'\s+', '', line)
                    lig_name = line
                    self.ligand_names.append(lig_name)
                except ValueError:
                    warn_msg = f"Cannot parse line {n+1} in {self.args.main.ligands_in}:\n"
                    warn_msg += line.rstrip()
                    warnings.warn(warn_msg)
                    continue
        self.ligand_names = list(set(self.ligand_names))

        self.skipped_ligand_names = []
        self.ligand_dirs = {}
        for ln in self.ligand_names:
            _dir = os.path.join(self.abfe_workdir, 'ligands', f'{ln}')
            lig_dirs = [d for d in glob(_dir) if os.path.isdir(d)]
            if not lig_dirs:
                self.skipped_ligand_names.append(ln)
                print(f"No dir for {ln} is found in 'ligands'.")
                print(f"Skipping {ln}.")
                continue
            self.ligand_dirs[ln] = lig_dirs[0]

        self.skipped_ligand_names = set(self.skipped_ligand_names)
        with common_tools.DirManager(os.path.join(self.abfe_workdir, '_logs')):
            _outfl = os.path.join(self.abfe_workdir, '_logs', 'skipped_ligands.out')
            with open(_outfl, 'w') as outfl:
                for lig_name in self.skipped_ligand_names:
                    outfl.write(f"Ligand {lig_name} is skipped.")

        if not self.ligand_names:
            print(time.strftime("%c"))
            print("All ligands have been skipped. The abfe workflow normally terminated.")
            sys.exit(1)

    def run_topology(self):
        """Perform the topology stage run.
        """
        _common_cmd_list = [
            'python',
            os.path.join(os.environ['ABFE_PKGPATH'], 'abfe', 'abfe_topology.py')
            ]
        self.topology_jobids = {}
        for lig_name, lig_dir in self.ligand_dirs.items():
            topology_args = dict(self.args.topology)
            self.ligand_dir = lig_dir
            topology_args.update(
                abfe_workdir=self.abfe_workdir,
                ligand_dir=self.ligand_dir,
                protein_dir=self.protein_dir,
                cofactor_dirs=self.cofactor_dirs,
                hfe_only=self.hfe_only,
                )
            sh_dir = os.path.join(self.abfe_workdir, '_topology', lig_name)
            slurm_params = {'job-name': lig_name+'.topology',
                            'output': 'slurm.out',
                            'error': 'slurm.err',
                            'partition': self.args.main.slurm['partition'],
                            'nodes': 1,
                            'ntasks-per-node': 1,
                            'cpus-per-task': 2,
                            'time': self.args.main.slurm['time'][0],
                            'gpujob': False,
                            }
            exclude =  self.args.main.slurm['exclude']
            if exclude.strip(): 
                slurm_params.update(exclude=exclude)

            with common_tools.DirManager(sh_dir):
                cmd = common_tools.dict_to_cmd(topology_args)
                cmd = ' '.join(_common_cmd_list) + ' ' + cmd
                common_tools.gen_slurm_bash(os.path.join(sh_dir, 'slurm.sh'), cmd, slurm_params)
                shcmd = 'sbatch slurm.sh'
                _, out, _ = common_tools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                self.topology_jobids[lig_name] = str(jobid[0])

    def run_equilibration(self):
        """Perform the equilibration stage run.
        """
        _common_cmd_list = [
            'python',
            os.path.join(os.environ['ABFE_PKGPATH'], 'abfe', 'abfe_equilibration.py')
            ]
        protein_name = self.args.topology.protein_name
        # submit equilibration setup jobs
        self.equilibration_jobids = {}
        for lig_name, _lig_dir in self.ligand_dirs.items():
            comp_name = f'{protein_name}_{lig_name}'
            equil_args = dict(self.args.equilibration)
            equil_args.update(
                ligand_name=lig_name,
                complex_name=comp_name,
                abfe_workdir=self.abfe_workdir,
                dry_run=self.args.equilibration.dry_run,
                restraintmask=self.args.equilibration.restraintmask,
                restraint_wt=self.args.equilibration.restraint_wt,
                temperature=self.args.main.mdcntrl['temperature'],
                cutoff=self.args.main.mdcntrl['cutoff'],
                timestep=self.args.equilibration.timestep,
                cpus_per_task=self.args.main.slurm['cpus-per-task'],
                hfe_only=self.hfe_only,
                )
            slurm_params = {
                'job-name': lig_name+'.equilibration',
                'output': 'slurm.out',
                'error': 'slurm.err',
                'partition': self.args.main.slurm['partition'],
                'nodes': 1,
                'ntasks-per-node': 1,
                'cpus-per-task': self.args.main.slurm['cpus-per-task'],
                'time': self.args.main.slurm['time'][1],
                'gpujob': True,
                'gpus': 1,
                'dependency': f'afterok:{self.topology_jobids[lig_name]}'
                }
            exclude =  self.args.main.slurm['exclude']
            if exclude.strip(): 
                slurm_params.update(exclude=exclude)

            sh_dir = os.path.join(self.abfe_workdir, '_equilibration', lig_name)
            with common_tools.DirManager(sh_dir):
                cmd = common_tools.dict_to_cmd(equil_args)
                cmd = ' '.join(_common_cmd_list) + ' ' + cmd
                common_tools.gen_slurm_bash(os.path.join(sh_dir, 'slurm.sh'), cmd, slurm_params)
                shcmd = 'sbatch slurm.sh'
                _, out, _ = common_tools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                self.equilibration_jobids[lig_name] = str(jobid[0])

    def run_alchemy_morph(self):
        """Perform alchemy morph stage run.
        """
        _common_cmd_list = [
            'python',
            os.path.join(os.environ['ABFE_PKGPATH'], 'abfe', 'abfe_alchemy_morph.py')
            ]
        protein_name = self.args.topology.protein_name
        self.alchemy_morph_jobids = {}
        for lig_name, lig_dir in self.ligand_dirs.items():
            comp_name = f'{protein_name}_{lig_name}'
            alchemy_morph_args = dict(self.args.alchemy_morph)
            alchemy_morph_args.update(
                ligand_name=lig_name,
                complex_name=comp_name,
                abfe_workdir=self.abfe_workdir,
                hfe_only=self.hfe_only,
                )
            sh_dir = os.path.join(self.abfe_workdir, '_alchemy_morph', lig_name)
            slurm_params = {'job-name': lig_name+'.alchemy_morph',
                            'output': 'slurm.out',
                            'error': 'slurm.err',
                            'partition': self.args.main.slurm['partition'],
                            'nodes': 1,
                            'ntasks-per-node': 1,
                            'cpus-per-task': 2,
                            'time': self.args.main.slurm['time'][2],
                            'gpujob': False,
                            'dependency': f'afterok:{self.equilibration_jobids[lig_name]}'
                            }
            exclude =  self.args.main.slurm['exclude']
            if exclude.strip(): 
                slurm_params.update(exclude=exclude)

            with common_tools.DirManager(sh_dir):
                cmd = common_tools.dict_to_cmd(alchemy_morph_args)
                cmd = ' '.join(_common_cmd_list) + ' ' + cmd
                common_tools.gen_slurm_bash(os.path.join(sh_dir, 'slurm.sh'), cmd, slurm_params)
                shcmd = 'sbatch slurm.sh'
                _, out, _ = common_tools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                self.alchemy_morph_jobids[lig_name] = str(jobid[0])

    def run_alchemy_md(self):
        """Perform alchemy md stage run.
        """
        _common_cmd_list = [
            'python',
            os.path.join(os.environ['ABFE_PKGPATH'], 'abfe', 'abfe_alchemy_md.py')
            ]
        protein_name = self.args.topology.protein_name
        self.alchemy_md_jobids = {}
        for lig_name, _lig_dir in self.ligand_dirs.items():
            comp_name = f'{protein_name}_{lig_name}'
            alchemy_md_args= dict(self.args.alchemy_md)
            alchemy_md_args.update(
                hfe_only=self.hfe_only,
                ligand_name=lig_name,
                complex_name=comp_name,
                abfe_workdir=self.abfe_workdir,
                lig_restraint_atoms=self.args.alchemy_md.lig_restraint_atoms[lig_name],
                rec_restraint_atoms=self.args.alchemy_md.rec_restraint_atoms[lig_name],
                dry_run=self.args.alchemy_md.dry_run,
                temperature=self.args.main.mdcntrl['temperature'],
                cutoff=self.args.main.mdcntrl['cutoff'],
                timestep=self.args.alchemy_md.timestep,
                ntwx=self.args.main.mdcntrl['ntwx'],
                gti_add_sc=self.args.main.mdcntrl['gti_add_sc'],
                partition=self.args.main.slurm['partition'],
                exclude=self.args.main.slurm['exclude'],
                gpus_per_node=self.args.main.slurm['gpus-per-node'],
                )
            slurm_params = {
                'job-name': lig_name+'.alchemy_md',
                'output': 'slurm.out',
                'error': 'slurm.err',
                'partition': self.args.main.slurm['partition'],
                'nodes': 1,
                'ntasks-per-node': 1,
                'cpus-per-task': self.args.main.slurm['cpus-per-task'],
                'time': self.args.main.slurm['time'][3],
                'gpujob': True,
                'gpus': 1,
                'dependency': f'afterok:{self.alchemy_morph_jobids[lig_name]}'
                }
            exclude =  self.args.main.slurm['exclude']
            if exclude.strip(): 
                slurm_params.update(exclude=exclude)

            sh_dir = os.path.join(self.abfe_workdir, '_alchemy_md', lig_name)
            with common_tools.DirManager(sh_dir):
                cmd = common_tools.dict_to_cmd(alchemy_md_args)
                cmd = ' '.join(_common_cmd_list) + ' ' + cmd
                common_tools.gen_slurm_bash(os.path.join(sh_dir, 'slurm.sh'), cmd, slurm_params)
                shcmd = 'sbatch slurm.sh'
                _, out, _ = common_tools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                self.alchemy_md_jobids[lig_name] = str(jobid[0])

    def run_alchemy_analysis(self):
        """Perform alchemy analysis stage run.
        """
        _common_cmd_list = [
            'python',
            os.path.join(os.environ['ABFE_PKGPATH'], 'abfe', 'abfe_alchemy_analysis.py')
            ]
        protein_name = self.args.topology.protein_name
        self.alchemy_analysis_jobids = {}
        for lig_name, lig_dir in self.ligand_dirs.items():
            comp_name = f'{protein_name}_{lig_name}'
            alchemy_analysis_args = dict(self.args.alchemy_analysis)
            alchemy_analysis_args.update(
                ligand_name=lig_name,
                complex_name=comp_name,
                abfe_workdir=self.abfe_workdir,
                hfe_only=self.hfe_only,
                temperature=self.args.main.mdcntrl['temperature'],
                )
            sh_dir = os.path.join(self.abfe_workdir, '_alchemy_analysis', lig_name)
            slurm_params = {'job-name': lig_name+'.alchemy_analysis',
                            'output': 'slurm.out',
                            'error': 'slurm.err',
                            'partition': self.args.main.slurm['partition'],
                            'nodes': 1,
                            'ntasks-per-node': 1,
                            'cpus-per-task': 2,
                            'time': self.args.main.slurm['time'][4],
                            'gpujob': False,
                            'dependency': f'afterok:{self.alchemy_md_jobids[lig_name]}'
                            }
            exclude =  self.args.main.slurm['exclude']
            if exclude.strip(): 
                slurm_params.update(exclude=exclude)

            with common_tools.DirManager(sh_dir):
                cmd = common_tools.dict_to_cmd(alchemy_analysis_args)
                cmd = ' '.join(_common_cmd_list) + ' ' + cmd
                common_tools.gen_slurm_bash(os.path.join(sh_dir, 'slurm.sh'), cmd, slurm_params)
                shcmd = 'sbatch slurm.sh'
                _, out, _ = common_tools.command_caller(command=shcmd, shell=True)
                jobid = re.findall(r"job (.*?)\n", str(out))
                self.alchemy_analysis_jobids[lig_name] = str(jobid[0])

if __name__ == '__main__':
    """Parse the arguments of config.yaml, then run the whole abfe workflow.
    """
    args = _parse_args()
    workflow = abfeWorkflow(args)
    if 'topology' in args.main.stages or 'all' in args.main.stages:
        workflow.run_topology()
    if 'equilibration' in args.main.stages or 'all' in args.main.stages:
        workflow.run_equilibration()
    if 'alchemy_morph' in args.main.stages or 'all' in args.main.stages:
        workflow.run_alchemy_morph()
    if 'alchemy_md' in args.main.stages or 'all' in args.main.stages:
        workflow.run_alchemy_md()
    if 'alchemy_analysis' in args.main.stages or 'all' in args.main.stages:
        workflow.run_alchemy_analysis()