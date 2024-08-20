#!/usr/bin/env bash
# -*- coding:utf-8 -*-
#==============================================================================
#CopyRight (c) By freeenergylab.
#@Description:
#This is one shell script that can setup the necessary environments, then
#submit ABFE jobs.
#
#@Author: Pengfei Li
#@Date: April 9th, 2024
#==============================================================================

if [ -f "/etc/profile" ] ; then
    source /etc/profile
fi
#==============================================================================
# Setup User-defined Variables
#==============================================================================
module purge
module load slurm/slurm/20.02.7
module load cuda12.0/toolkit/12.0.1
module load openmpi/5.0.3
module load anaconda3/FEP-SPell-ABFE
module load amber/amber24_ambertools24
module load apbs/3.4.1

export CONFIG_FILE='config.yaml'
export ABFE_PKGPATH=$HOME/software/freeenergylab/FEP-SPell-ABFE/dev
export PYTHONPATH=$ABFE_PKGPATH:$PYTHONPATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
#==============================================================================
# Setup Python Requirement
#==============================================================================
MIN_REQUIRED_PYTHON_VERSION="3.9"
parsed_MIN_VERSION=$(echo "${MIN_REQUIRED_PYTHON_VERSION//./}")
# check if python is installed
if ! command -v python >/dev/null 2>&1; then
    echo "Error: Python is not installed."
fi
python_executable=$(command -v python)
# get the actual Python version
python_version=$(python --version 2>&1 | awk '{print $2}')
parsed_python_version=$(echo "${python_version//./}")
# check if the version is greater than the minimum
if [[ $parsed_python_version -lt $parsed_MIN_VERSION ]]; then
    echo "Error: Python version $python_version is not compatible."
    echo "Python $MIN_REQUIRED_PYTHON_VERSION or higher is required."
else
    echo "Python $python_version detected (compatible)"
fi
#==============================================================================
# Submit ABFE Jobs
#==============================================================================
if [ -d "_logs" ]; then
    echo "Delete the previous _logs directory..."
    rm -rf _logs
fi
if [ -f $CONFIG_FILE ] ; then
    export PYTHON_EXE=$python_executable
    $PYTHON_EXE $ABFE_PKGPATH/abfe/abfe_main.py -i $CONFIG_FILE
    echo "Used python is <$PYTHON_EXE>."
    echo "Job submission is finished sucessfully now."
else
    echo "Job submission quits here."
    echo "Please provide one configure file named <$CONFIG_FILE>."
fi
