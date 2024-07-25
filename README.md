# FEP-SPell-ABFE
FEP-SPell-ABFE: An Open Source Alchemical Absolute Binding Free Energy Calculation Workflow for Drug Discovery

# Notable features:
- Support for performing not only ABFE but also AHFE (absolute hydration free energy) prediction when turn on "hfe_only".
- Support AlChemical Enhanced Sampling (so-called ACES: one REST2-like method) by setting "gti_add_sc=5/6".
- Support Replica-exchange method for FEP production simulations, but this requires one lambda window simulation to occupy one GPU card.
- Support for specifying ion concentration condition.
- Support for doing charge correction for charged ligands based on Poisson-Boltzmann algorithm using APBS software.
- Support for running one or more stages ("topology", "equilibration", "alchemy_morph", "alchemy_md", "alchemy_analysis") in one submit.
- Support local parsimonious: the simulations of stages ("equilibration", "alchemy_md") can be skipped when turn on "dry_run".
- Support for generating alchemical analysis results in .svg plot.
- Support for including 'cofactors' in the system.

# Basic usage:
Before launching ABFE jobs, your working directory should include the following files:
- A `proteins` directory with structure `proteins/[protein_name]/protein.pdb`, when "hfe_only=False".
- A `ligands` directory with structure `ligands/[ligand_name]/ligand.sdf(or ligand.mol2)`.
- A `cofactors` directory with structure `cofactors/[cofactor_name]/cofactor.sdf (or cofactor.mol2)`, if cofactors are present. `[cofactor_name]` should be three-capital-letter style, e.g. CFA, CFB et al.
- A `ligands.in` text file determining which compounds will be submitted to calculate ABFE.
- A `config.yaml` yaml file. Check out the example file: `testing/abfe_testing/config.yaml`.
- A `sbmitBFE.sh` bash file. Check out the example file: `testing/abfe_testing/submitBFE.sh`.

After the above files are prepared well, please submit your jobs using the following command:
```
cd Your_Working_Directory
submitBFE.sh -i config.yaml # submitBFE.sh should be modified according to your environment
```

# Advanced usage:
Users can fine-tune the behavior of this workflow by editing the `config.yaml` file. Please note the comments in the `config.yaml` file.

# Dependencies:
This workflow mainly depends on the packages including Slurm, CUDA, OpenMPI, Amber, Anaconda3 and APBS. With the environment-modules's help, you can use `module load` command to load the dependent softwares, like:
```
module load slurm/slurm/20.02.7
module load cuda12.0/toolkit/12.0.1
module load openmpi/5.0.3
module load amber24_ambertools24
module load anaconda3/FEP-SPell-ABFE
module load apbs/3.4.1
```
## Anaconda3 environment depolyment refers to:
```
 1. cd Your_Anaconda3_Installation_PATH # e.g. cd $HOME/software/anaconda3/2024.06
 2. wget https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh
 3. chmod +x Anaconda3-2024.06-1-Linux-x86_64.sh
 4. ./Anaconda3-2024.06-1-Linux-x86_64.sh -u
 5. source $HOME/software/anaconda3/2024.06/bin/activate base
 6. conda create --name FEP-SPell-ABFE --clone base
 7. conda activate FEP-SPell-ABFE
 8. conda env create --name FEP-SPell-ABFE --file=environment_abfe.yml
 9. pip install rdkit
10. pip install parmed

###############################################################################

module load anaconda3/FEP-SPell-ABFE

# FEP-SPell-ABFE modulefile, like:
"""
conflict anaconda3

setenv ANACONDA3HOME "$env(HOME)/software/anaconda3/2024.06"
prepend-path PATH "$env(ANACONDA3HOME)/envs/FEP-SPell-ABFE/bin"
prepend-path PYTHONPATH "$env(ANACONDA3HOME)/envs/FEP-SPell-ABFE/lib/python3.12/site-packages"
"""
```
## Amber and AmberTools software installation refers to:
```
chmod +x amber.sh
###############################################################################

module load amber24_ambertools24

# amber24_ambertools24 modulefile, like:
"""
conflict amber24

setenv AMBERHOME "$env(HOME)/software/amber/amber24_ambertools24/amber24"

system "$env(AMBERHOME)/amber.sh"
"""
```
## APBS software installation refers to:
```
 1. wget https://github.com/Electrostatics/apbs/releases/download/v3.4.1/APBS-3.4.1.Linux.zip
 2. unzip APBS-3.4.1.Linux.zip
 3. module load apbs/3.4.1
    (configure modulefile like:
     setenv APBSHOME "[Path to install APBS]/APBS-3.4.1.Linux"
     prepend-path PATH "$env(APBSHOME)/bin"
     prepend-path LD_LIBRARY_PATH "$env(APBSHOME)/lib"
     )
```