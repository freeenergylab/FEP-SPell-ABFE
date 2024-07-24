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
- A 'proteins' directory with structure `proteins/[protein_name]/protein.pdb`, when "hfe_only=False".
- A 'ligands' directory with structure `ligands/[ligand_name]/ligand.sdf(or ligand.mol2)`.
- A 'cofactors' directory with structure `cofactors/[cofactor_name]/cofactor.sdf (or cofactor.mol2)`, if cofactors are present.
- A configuration yaml format file. Check out the example file: `testing/abfe_testing/config.yaml`.
- A sbmitBFE.sh bash shell file. Check out the example file: `testing/abfe_testing/submitBFE.sh`.

After the above files are prepared well, please submit your jobs using the following command:
```
cd Your_Working_Directory
submitBFE.sh -i config.yaml # submitBFE.sh should be modified according to your environment
```

# Advanced usage:
Users can fine-tune the behavior of this workflow by editing the `config.yaml` file. Please note the comments in the `config.yaml` file.

# Dependencies:
This workflow mainly depends on the following packages:
```
module load slurm/slurm/20.02.7
module load cuda12.0/toolkit/12.0.1
module load openmpi/5.0.3
module load amber24_ambertools24
module load anaconda3/FEP-SPell-ABFE
module load apbs/3.4.1
```
Anaconda3 conda environment depolyment refers to:
```
conda env create --name FEP-SPell-ABFE --file=environment_abfe.yml
```
APBS software installation refers to:
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