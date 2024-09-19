#!/usr/bin/python3
# coding: utf-8
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This script is used to check whether the MD simulations terminated normally.

@Author: Writted by Tingting Pu
         Modified by Pengfei Li
@Date: Dec. 18, 2023
==============================================================================
"""
import os
import time
from collections import defaultdict
from glob import glob

import abfe.utils.common_tools as ctools

class JobManager(object):
    def __init__(self, workdir, sys_name, stage, step):
        self.workdir = workdir
        self.sys_name = sys_name
        self.stage = stage
        self.step = step

    def seek_jobid(self):
        '''Seek jobid for each sys_name/stage/step in _alchemy_md
        '''
        print(f'Seeking jobid for {self.sys_name}/{self.stage}/{self.step}...')
        _log_file = os.path.join(self.workdir, '_alchemy_md', self.sys_name,
                                 self.stage, self.step, 'slurm_jobid.log')
        while not os.path.exists(_log_file):
            print(f'The {_log_file} has not been generated yet.')
            time.sleep(5)
        else:
            with open(_log_file, 'r') as file:
                first_line = file.readline()
                if 'INFO: ARRAY JOB ID' in first_line:
                    job_id = (first_line.split()[-1].strip())
                    print(f'INFO: ARRAY JOB ID is {job_id}')
                    self._check_squeue(job_id)
                else:
                    print(f'The {_log_file} file is an empty file.')

    def _get_state_with_sacct(self, job_id):
        """Get job state with sacct command.

        Args:
            job_id (str): the slurm job id

        Returns:
            job_states (dict): the job states
        """
        job_states = defaultdict()
        cmd = 'sacct -j '
        cmd += job_id
        cmd += ' -n --format=jobid%30,node,jobname%60,start,end,state%30'
        with os.popen(cmd, 'r') as f:
            for line in f:
                line = line.strip()
                if 'batch' not in line:
                    _line = line.split()
                    jobid = str(_line[0])
                    job_states[jobid] = str(_line[5])

        return job_states

    def _check_squeue(self, job_id):
        '''Read job information for each jobid.
        '''
        while True:
            job_states = self._get_state_with_sacct(job_id)
            job_states_list = list(job_states.values())
            status = job_states_list.count('COMPLETED') == len(job_states_list)
            if status or 'FAILED' in job_states_list:
                break
            time.sleep(60)

        _base_dir = os.path.join(self.workdir, '_alchemy_md', self.sys_name,
                                 self.stage, self.step)
        print(f'Checking all ti.out information for {_base_dir} ...')
        ti_out_files = glob(os.path.join(_base_dir, '[01]*/ti*.out'))
        ti_out_files.sort()
        flags = []
        for file_path in ti_out_files:
            _cmd = f'tail -n 5 {file_path}'
            _, _out, _ = ctools.command_caller(command=_cmd, shell=True)
            _line = ''.join(_out).strip().split("\n")
            if _line[-1].startswith('|  Total wall time:'):
                flags.append(True)
                print(f'Total wall time found in {file_path}.')
            else:
                flags.append(False)
                print(f'The {file_path} may fail during simulation.')

        if False in flags:
            self.flag = False
        else:
            self.flag = True