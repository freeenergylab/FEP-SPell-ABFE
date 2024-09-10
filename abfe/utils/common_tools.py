#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the common_tools module that provides some common tools.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import os
import pickle
import stat
import subprocess
import time
import warnings
import yaml

import abfe.const as const

def read_config(config_file):
    """Read the configuration file for the workflow and parse it into arg=value pairs.

    Args:
        config_file (str): Path of the configuration file.

    Returns:
        A dictionary contains the parameter arg=value pairs.
    """
    try: 
        with open(config_file, 'r') as infile:
            params_dict = yaml.load(infile, Loader=yaml.FullLoader)
        return params_dict
    except:
        error_msg = f"Cannot read options from the config file {config_file}."
        raise RuntimeError(error_msg)

class DotDict(dict):
    """Transform a dictionary such that the values of the dict can be accessed using dot representation.
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

class DirManager(object):
    """Manage enter and exit of one directory.

    Args:
        dst (str): Path of the dir to manage.
        verbose (bool): If True, print messages of the current actions.
    """
    def __init__(self, dst, verbose=False):
        self.topdir = os.getcwd()
        self.dst = os.path.join(self.topdir, dst)
        self.verbose = verbose

    def __enter__(self):
        if not os.access(self.dst, os.F_OK):
            if self.verbose:
                print(time.strftime("%c"))
                print(f"Creating directory {self.dst}")
            try:
                os.makedirs(self.dst)
            except OSError as ose:
                print(f"Failed to create the folder {self.dst}: {str(ose)}")
                time.sleep(5)
        if self.verbose:
            print(time.strftime("%c"))
            print(f"Entering directory {self.dst}")
        os.chdir(self.dst)
        return

    def __exit__(self, type, value, trace):
        if self.verbose:
            print(time.strftime("%c"))
            print("Entering directory {self.topdir}")
        os.chdir(self.topdir)
        return

def dict_to_cmd(args):
    """Transform one arg=value dictionary to a command-line option string.

    Args:
        args (dict): Dictionary whose key-value pairs represent command-line argument names and values.
    
    Returns:
        A string reproduces the string one will type when running the script.
    """
    cmd = []
    for key, value in args.items():
        arg = f'--{str(key)}'
        if type(value) == bool:
            if value == True:
                cmd.append(arg)
            else:
                continue
        elif type(value) in [list, set, tuple]:
            cmd.append(arg)
            cmd.extend(str(v) for v in value)
        elif type(value) in [dict]:
            cmd.append(arg)
            for k, v in value.items():
                cmd.extend([str(k), str(v)])
        else:
            cmd.extend((arg, '"' + str(value) + '"'))

    return ' '.join(cmd)

def gen_bash_file(fsh, cmd):
    """Generate and chmod the bash file.

    Args:
        fsh (str): A shell file name.
        cmd (str): A command line.
    """
    cmd_list = ['#!/usr/bin/env bash']
    with open(fsh, 'w') as _fsh:
        _fsh.write('\n'.join(cmd_list) + '\n')
        _fsh.write('\n')
        _fsh.write(cmd + '\n')
    os.chmod(fsh, stat.S_IRWXU)

def gen_slurm_bash(fsh, cmd, slurm_params):
    """Generate and chmod the slurm bash file.

    Args:
        fsh (str): A shell file name.
        cmd (str): A command line.
    """
    cmd_list = ['#!/usr/bin/env bash']
    cmdline_list = []

    if 'job-name' in list(slurm_params.keys()):
        jobname = slurm_params['job-name']
        cmd_list.append(f'#SBATCH --job-name="{jobname}"')
        cmdline_list.append(f'--job-name "{jobname}"')
    if 'output' in list(slurm_params.keys()):
        output = slurm_params['output']
        cmd_list.append(f'#SBATCH --output="{output}"')
        cmdline_list.append(f'--output "{output}"')
    if 'error' in list(slurm_params.keys()):
        error = slurm_params['error']
        cmd_list.append(f'#SBATCH --error="{error}"')
        cmdline_list.append(f'--error "{error}"')
    if 'partition' in list(slurm_params.keys()):
        partition = slurm_params['partition']
        cmd_list.append(f'#SBATCH --partition={partition}')
        cmdline_list.append(f'--partition {partition}')
    if 'nodes' in list(slurm_params.keys()):
        nodes = slurm_params['nodes']
        cmd_list.append(f'#SBATCH --nodes={nodes}')
        cmdline_list.append(f'--nodes {nodes}')
    if 'ntasks-per-node' in list(slurm_params.keys()):
        ntaskspernode = slurm_params['ntasks-per-node']
        cmd_list.append(f'#SBATCH --ntasks-per-node={ntaskspernode}')
        cmdline_list.append(f'--ntasks-per-node {ntaskspernode}')
    if 'cpus-per-task' in list(slurm_params.keys()):
        cpuspertask   = slurm_params['cpus-per-task']
        cmd_list.append(f'#SBATCH --cpus-per-task={cpuspertask}')
        cmdline_list.append(f'--cpus-per-task {cpuspertask}')
    if 'time' in list(slurm_params.keys()):
        time = slurm_params['time']
        cmd_list.append(f'#SBATCH --time={time}')
        cmdline_list.append(f'--time {time}')
    if 'dependency' in list(slurm_params.keys()):
        dependency = slurm_params['dependency']
        cmd_list.append(f'#SBATCH --dependency={dependency}')
        cmdline_list.append(f'--dependency {dependency}')
    if 'array' in list(slurm_params.keys()):
        array = slurm_params['array']
        cmd_list.append(f'#SBATCH --array={array}')
        cmdline_list.append(f'--array {array}')
    if 'exclude' in list(slurm_params.keys()):
        exclude = slurm_params['exclude']
        cmd_list.append(f'#SBATCH --exclude={exclude}')
        cmdline_list.append(f'--exclude {exclude}')
    if slurm_params['gpujob']:
        gpus = slurm_params['gpus']
        cmd_list.append(f'#SBATCH --gres=gpu:{gpus}')
        cmdline_list.append(f'--gres gpu:{gpus}')
    with open(fsh, 'w') as _fsh:
        _fsh.write('\n'.join(cmd_list) + '\n')
        _fsh.write('\n')
        _fsh.write(cmd + '\n')
    os.chmod(fsh, stat.S_IRWXU)

    return cmdline_list

def guess_file_info(filepath):
    """Guess the name of the molecule and the format based on filepath.

    Args:
        filepath (str): Path of the file.

    Returns:
        molname (str): Guessed name of the molecule.
        ext (str): Guessed format of the file.
    """
    abs_path = os.path.abspath(filepath)
    molname = os.path.basename(os.path.dirname(abs_path))
    filename = os.path.basename(abs_path)
    _, ext = os.path.splitext(filename)
    ext = ext.replace('.', '')
    return molname, ext

def get_tleap_input(component_list, output_prefix, water_model, boxtype,
                    boxbuffer, neutralize=True, cation='Na+', anion='Cl-',
                    num_cations=0, num_anions=0, phosaa_ff=None, lipid_ff=None,
                    protein_type='soluble', box_info=[]):
    """Generate an input cmd list for tleap.

    Args:
        component_list (list):
            List of Protein/Ligand/Cofactors objects which have already loaded force
            field information. The minimal requirement is an attribute `mol2_file`
            _or_ `coord_file` (these may be None) and an attribute dict `ff_info`
            with possible keys 'source', 'lib', and 'dat'.
        output_prefix (str):
            Prefix for the output parm7 and rst7 files.
        water_model (str):
            {`tip3p`, `opc`, `tip4pew`, None}
            Name of the water model to use. If None, this is skipped.
        boxtype (str):
            {`rect`, `oct`, `nosolv`}
            Add a rectangular or truncated octhedtral waterbox, unless `nosolv` is
            given, in which case this can used for no periodic boundary.
        boxbuffer (float):
            Buffer for the solvation box (in A).
        neutralize (bool, optional):
            If True, neutralize the system with ions and/or add extra ions (see
            `num_cations/`num_anions`)
        cation (str, optional):
            {`K+`, `Na+`}
            Name of the cation to use for neutralization (if neutralize is True)
        anion (str, optional):
            {'Cl-'}
            Name of the anion to use for neutralization (if neutralize is True)
        num_cations (int, optional):
            Number of cations to add. If zero, only enough to neutralize the system
            will be added (can also be zero)
        num_anions (int, optional):
            Number of anions to add. If zero, only enough to neutralize the system
            will be added (can also be zero)
        protein_type (str):
            {`soluble`, `membrane`}
        box_info (list):
            The box information,
            e.g. CRYST1   81.971   82.209   88.609  90.00  90.00  90.00

    Returns:
        leap_cmds (list):
            Commands to be sent to tleap. These should be written to separate lines
            in a text file in order to be correctly interpreted.
    """
    output_prefix = str(output_prefix)
    pdb_nowat = f'{output_prefix}_nowat.pdb'
    pdb = f'{output_prefix}.pdb'
    prmtop = f'{output_prefix}.parm7'
    restrt = f'{output_prefix}.rst7'

    leap_cmds = []
    if water_model is not None:
        error_msg = "Only support TIP3P, OPC, TI4PEW watermodel."
        assert water_model in const.WATER_DEFS, error_msg
        leap_cmds.append('source %s' % os.path.join(const.LEAPHOME, 'cmd',
            const.WATER_DEFS[water_model]['source'])
            )
    water_box = const.WATER_DEFS[water_model]['box']

    # Lists for tracking which files have already been loaded.
    used_files = {'source': [], 'lib': [], 'dat': []}
    # Commands for loading different file types.
    load_cmd = {'source':'source', 'lib':'loadOff', 'dat':'loadAmberParams'}
    comp_names = []
    for n, comp in enumerate(component_list):
        cname = f'Comp{n}'
        for k, v in used_files.items():
            if k in comp.ff_info and comp.ff_info[k] not in v:
                leap_cmds.append(f'{load_cmd[k]} {comp.ff_info[k]}')
                used_files[k].append(comp.ff_info[k])
        if hasattr(comp, 'mol2_file') and comp.mol2_file is not None:
            leap_cmds.append(f'{cname} = loadMol2 {comp.mol2_file}')
            comp_names.append(cname)
        elif comp.coord_file is not None:
            if protein_type == 'membrane':
                if lipid_ff is not None:
                    error_msg = "Only support lipid21 model."
                    assert lipid_ff in const.LIPID_DEFS, error_msg
                    leap_cmds.append('source %s' % os.path.join(const.LEAPHOME, 'cmd',
                        const.LIPID_DEFS[lipid_ff]['source']))
                else:
                    print(f'Please specify a force field for lipid membrane.')
            if phosaa_ff is not None:
                error_msg = "Only support phosaa14SB and phosaa19SB model."
                assert phosaa_ff in const.PHOSAA_DEFS, error_msg
                leap_cmds.append('source %s' % os.path.join(const.LEAPHOME, 'cmd',
                    const.PHOSAA_DEFS[phosaa_ff]['source'])
                    )
            leap_cmds.append(f'{cname} = loadPdb {comp.coord_file}')
            comp_names.append(cname)

    leap_cmds.append('system = combine {%s}' % ' '.join(comp_names))
    leap_cmds.append(f'savePdb system {pdb_nowat}')

    boxbuffer = max(0.0, float(boxbuffer))
    num_cations = max(0, int(num_cations))
    num_anions = max(0, int(num_anions))

    if protein_type == 'soluble' or output_prefix == 'solvated':
        if boxtype == 'rect':
            leap_cmds.append(f'solvateBox system {water_box} {boxbuffer}')
        elif boxtype == 'oct':
            leap_cmds.append(f'solvateOct system {water_box} {boxbuffer}')
        else:
            print('Only support a rectangular or truncated octhedtral waterbox.')

        if neutralize:
            leap_cmds.append(f'addions system {cation} 0')
            leap_cmds.append(f'addions system {anion} 0')

        if not num_cations == 0 or not num_anions == 0:
            leap_cmds.append(f'addions system {cation} {num_cations} {anion} {num_anions}')
    elif protein_type == 'membrane' and output_prefix == 'complex':
        leap_cmds.append('set system box {%s %s %s}'%(box_info[1], box_info[2], box_info[3]))

    leap_cmds.append(f'saveAmberParm system {prmtop} {restrt}')
    leap_cmds.append(f'savePdb system {pdb}')
    leap_cmds.append('quit')

    return leap_cmds

def command_caller(command, shell=False):
    """Execute one Linux command line.

    Args:
        command (string): Linux command line
        shell (bool, optional): Defaults to False. If true, the command will be executed through the shell.

    Returns:
        sp.returncode (int)
        out.decode() (string)
        err (bytes)
    """
    sp = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=shell)
    out, err = sp.communicate()
    if sp.returncode:
        print(f"Return code: {sp.returncode} Error message: {err}")
    return sp.returncode, out.decode(), err

def parse_pdb(fpdb):
    """Parse the pdb file.

    Args:
        fpdb (str): pdb filename

    Returns:
        box_info (list): contains the box information
    """    
    with open(fpdb, 'r') as _fpdb:
        line = _fpdb.readline()
        box_info = line.split()

    return box_info

def parse_leapin(leapin, box_info, leapout, duplicate=True,
                 stage=None, resWATpdb=None, leap_cmd=None):
    """Parse the tleap input file and dump one tleap output file.

    Args:
        leapin (string): tleap input file name
        box_info (list): contains the box information
        leapout (string): tleap output file name
        others: other functions for some specific purposes
    """    
    new_cmds = []
    with open(leapin, 'r') as _fleapin:
        lines = _fleapin.readlines()
        for line in lines:
            if "loadMol2" in line:
                mol2_cmd = line.split()
                mol2_cmd[-1] = "vacuum.mol2" + "\n"
                line = ' '.join(mol2_cmd)
                #in solvated stage, the is no Comp1 in tleap.in
                if stage=='solvated':
                    line += 'Comp1 = loadPdb protein.pdb\n'
                if leap_cmd is not None:
                    for cmd in leap_cmd:
                        line += cmd+'\n'
                if resWATpdb is not None:
                    CompWats = []
                    for i, pdb in enumerate(resWATpdb):
                        line += 'CompWat%s = loadPdb %s\n'%(str(i), pdb)
                        CompWats.append('CompWat%s'%str(i))
            elif "loadPdb" in line:
                pdb_cmd = line.split()
                pdb_cmd[-1] = "protein.pdb" + "\n"
                line = ' '.join(pdb_cmd)
            elif "combine" in line:
                #in complex_restratin stage, it needs two copies of ligand
                if duplicate:
                    line = line.replace("Comp0", "Comp0 Comp0")
                #in solvated stage, the is no Comp1 in tleap.in
                if stage=='solvated':
                    line = line.replace("Comp0", "Comp0 Comp1")
                if resWATpdb is not None:
                    line = line.replace("Comp0 Comp1", "Comp0 Comp1 %s"
                                        %(' '.join(CompWats)))
            elif "solvateBox" in line or "solvateOct" in line:
                line = 'set system box {%s %s %s} \n' % \
                    (box_info[1], box_info[2], box_info[3])
            elif 'addions' in line:
                continue
            elif 'savePdb' in line:
                continue
            new_cmds.append(line)
    with open(leapout, 'w') as _leapout:
        _leapout.writelines(new_cmds)

def pkl_check(obj, pklfile, check_keys):
    """Check whether the object stored in the pklfile can be considered as
    equivalent to the obj based on criteria provided by check_keys.

    Args:
        obj (object): The reference object.
        pklfile (str): Path of the pklfile.
        check_keys (list): List of keys representing the attributes of the objects to check

    Returns:
        A bool value representing whether the pkl file passes the check.
    """
    if not file_check(pklfile):
        return False

    try:
        with open(pklfile, 'rb') as pickl:
            stored_obj = pickle.load(pickl)
    except Exception:
        warnings.warn(f"Failed to read pkl file {os.path.abspath(pklfile)}.")
        return False

    ref_dict, work_dict = obj.__dict__, stored_obj.__dict__
    for ckey in check_keys:
        if ckey not in ref_dict.keys():
            warnings.warn("'%s' is not an attribute of %s." \
                          % (ckey, obj))
            continue
        elif ckey not in ref_dict.keys():
            warnings.warn("'%s' is not an attribute of %s loaded from %s." \
                          % (ckey, stored_obj, os.path.abspath(pklfile)))
            continue
        elif ref_dict[ckey] == work_dict[ckey]:
            continue
        else:
            return False
    return True

def file_check(filepath, size_threshold=1e3):
    """Check the validity of a file based on its size.
    """
    if os.access(filepath, os.R_OK):
        return (os.path.getsize(filepath) >= size_threshold)
    return False