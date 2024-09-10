#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
===============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the biochemsystem module that prepares the biochemical systems.

@Author: Pengfei Li
@Date: May 30th, 2022
===============================================================================
"""
import os
import sys
import time

import parmed as pmd
import rdkit.Chem as Chem

import abfe.const as const
import abfe.utils.common_tools as ctools

class _BCS_Base(object):
    """Base class for constructing BioChemSystem.
    """
    def __init__(self, coord_file):
        self.coord_file = os.path.abspath(coord_file)
        self.molname, self.fmt = ctools.guess_file_info(self.coord_file)
        self.resname = None # needs to be initilized in Child class
        self.net_charge = None # needs to be initilized in Child class
        self.ff_info = {}
        self.use_gaff2 = False

    def _initialize_rmol(self):
        """Initialzing rdkit rmol object for small molecules.

        Raises:
            ValueError: only support mol2 or sdf formats
        """
        if self.fmt == 'mol2':
            try:
                self.rmol = Chem.MolFromMol2File(self.coord_file, removeHs=False)
                Chem.Kekulize(self.rmol, clearAromaticFlags=False)
            except Exception:
                self.rmol = Chem.MolFromMol2File(self.coord_file, removeHs=False, sanitize=False)
                Chem.Kekulize(self.rmol, clearAromaticFlags=False)
        elif self.fmt == 'sdf':
            try:
                self.rmol = Chem.SDMolSupplier(self.coord_file, removeHs=False)[0]
                Chem.Kekulize(self.rmol, clearAromaticFlags=False)
            except Exception:
                self.rmol = Chem.SDMolSupplier(self.coord_file, removeHs=False, sanitize=False)[0]
                Chem.Kekulize(self.rmol, clearAromaticFlags=False)
        else:
            error_msg = "Currently only support mol2 or sdf formats."
            raise ValueError(error_msg)

    def _add_ff_info(self, **kwargs):
        """Add force field parameter files for tleap into ff_info dict.

        Raises:
            KeyError: the ff files should be recognized by tleap
        """
        for key, path in kwargs.items():
            if key not in ['source', 'lib', 'dat']:
                error_msg = "Only source/lib/dat keys are supported in tleap."
                raise KeyError(error_msg)
            self.ff_info[key] = path

    def _load_ff(self, ff='gaff2', workdir=None):
        """Load the force field parameters for small molecules, e.g. Amber GAFF2 as default.
        """
        if self.use_gaff2:
            print(f"{self} has been parameterized by GAFF2.")
            return

        if workdir is None:
            workdir = os.path.dirname(self.coord_file)

        with ctools.DirManager(workdir):
            print(time.strftime("%c"))
            print(f"Started generating {ff} parameters for {self}.")
            self.mol2_file = os.path.abspath('vacuum.mol2')
            self.frcmod_file = os.path.abspath('vacuum.frcmod')
            if not os.access(self.mol2_file, os.R_OK):
                # 0). generate the ac mol2 file with GAFF2 atom-type from the input sdf/mol2 file
                ac_cmd = [
                    f"{os.path.join(const.AMBERBIN, 'antechamber')}",
                    '-i', f'{self.coord_file}',
                    '-fi', f'{self.fmt}',
                    '-o', f'{self.mol2_file}',
                    '-fo', 'mol2',
                    '-at', f'{ff}',
                    '-c', 'bcc',
                    '-an', 'n',
                    '-nc', f'{self.net_charge}',
                    '-rn', self.resname,
                    '-dr', 'no',
                    '-pf', 'y', '>', '/dev/null', '2>&1']
                ac_status = os.system(' '.join(ac_cmd))
                ctools.gen_bash_file('antechamber_0.sh', ' '.join(ac_cmd))
                if ac_status != 0:
                    print(time.strftime("%c"))
                    print(f"Antechamber failed for {self}!")
                    sys.exit(1)

                # 1). generate the ac mol2 file with sybyl atom-type from the above mol2 file
                ac_sybyl_cmd = [
                    f"{os.path.join(const.AMBERBIN, 'antechamber')}",
                    '-i', f'{self.mol2_file}',
                    '-fi', 'mol2',
                    '-o', f'{os.path.splitext(self.mol2_file)[0]}_sybyl.mol2',
                    '-fo', 'mol2',
                    '-at', 'sybyl',
                    '-j', '1',
                    '-an', 'n',
                    '-nc', f'{self.net_charge}',
                    '-rn', self.resname,
                    '-dr', 'no',
                    '-pf', 'y', '>', '/dev/null', '2>&1']
                ac_sybyl_status = os.system(' '.join(ac_sybyl_cmd))
                ctools.gen_bash_file('antechamber_sybyl_1.sh', ' '.join(ac_sybyl_cmd))
                if ac_sybyl_status != 0:
                    # 2). if 1). failed, with adding -c bcc, then retry to generate the ac mol2 file with sybyl atom-type again
                    ac_sybyl_cmd = [
                        f"{os.path.join(const.AMBERBIN, 'antechamber')}",
                        '-i', f'{self.coord_file}',
                        '-fi', f'{self.fmt}',
                        '-o', f'{os.path.splitext(self.mol2_file)[0]}_sybyl.mol2',
                        '-fo', 'mol2',
                        '-at', 'sybyl',
                        '-c', 'bcc',
                        '-j', '1',
                        '-an', 'n',
                        '-nc', f'{self.net_charge}',
                        '-rn', self.resname,
                        '-dr', 'no',
                        '-pf', 'y', '>', '/dev/null', '2>&1']
                    ac_sybyl_status = os.system(' '.join(ac_sybyl_cmd))
                    ctools.gen_bash_file('antechamber_sybyl_2.sh', ' '.join(ac_cmd))
                    if ac_sybyl_status != 0:
                        print(time.strftime("%c"))
                        print(f"Antechamber failed for {self}!")
                        sys.exit(1)
            print(time.strftime("%c"))
            print(f"Finished loading {ff} parameters for {self}.")

            # Run parmchk2 to generate the frcmod file
            if not os.access(self.frcmod_file, os.R_OK):
                parmchk2_cmd = [
                    f"{os.path.join(const.AMBERBIN, 'parmchk2')}",
                    '-i', f'{self.mol2_file}',
                    '-f', 'mol2',
                    '-p', f"{os.path.join(const.LEAPHOME, 'parm', f'{ff}.dat')}",
                    '-o', f'{self.frcmod_file}']
                parmchk2_status = os.system(' '.join(parmchk2_cmd))
                if parmchk2_status != 0:
                    print(time.strftime("%c"))
                    print(f"Parmchk2 failed for {self}!")
                    sys.exit(1)
        self.leaprc = os.path.join(const.LEAPHOME, 'cmd', f'leaprc.{ff}')
        self._add_ff_info(source=self.leaprc, dat=self.frcmod_file)
        if ff == 'gaff2': 
            self.use_gaff2 = True

    def _create_lib(self, workdir):
        """Create .lib force field file for Cofactor.
        """
        with ctools.DirManager(workdir):
            self.lib_file = os.path.abspath('vacuum.lib')
            if not os.access(self.lib_file, os.R_OK):
                leapin =  f"source {self.leaprc}\n"
                leapin += f"loadAmberParams {self.frcmod_file}\n"
                leapin += f"{self.resname} = loadMol2 {self.mol2_file}\n"
                leapin += f"saveOff {self.resname} {self.lib_file}\n"
                leapin += f"quit"

                with open('tleap_lib.in', 'w') as infile:
                    infile.write(leapin)
                tleap_cmd = [f"{os.path.join(const.AMBERBIN, 'tleap')}", '-f', 'tleap_lib.in', '> /dev/null 2>&1']
                tleap_status = os.system(' '.join(tleap_cmd))
                if tleap_status != 0:
                    print(time.strftime("%c"))
                    print(f"tleap failed for {self}!")
                    sys.exit(1)
        self._add_ff_info(lib=self.lib_file)

class Cofactor(_BCS_Base):
    """Child class for constructing cofactor in complex system.
    """
    def __init__(self, coord_file):
        super(Cofactor, self).__init__(coord_file)
        self._initialize_rmol()
        self.net_charge = int(Chem.GetFormalCharge(self.rmol))
        self.resname = self.molname
        self.prepared = False

    def __repr__(self):
        return f"<Cofactor {self.molname}>"

    def prep_cofactor_ff(self, ff, workdir):
        """Prepare the force field for cofactor.
        """
        self._load_ff(ff=ff, workdir=workdir)
        self._create_lib(workdir=workdir)
        self.prepared = True

class Ligand(_BCS_Base):
    """Child class for constructing Ligand in complex system.
    """
    def __init__(self, coord_file):
        super(Ligand, self).__init__(coord_file)
        self._initialize_rmol()
        self.net_charge = int(Chem.GetFormalCharge(self.rmol))
        self.resname = 'MOL'
        self.prepared = False

    def __repr__(self):
        return f"<Ligand {self.molname}>"

    def prep_ligand_ff(self, ff, workdir, water_model, boxtype, boxbuffer, 
                       neutralize=True, hmr=False, hmass=3.024, ionconc=0.0):
        """Prepare the force field for ligand.
        """
        self.parm7_file = os.path.abspath(os.path.join(workdir, 'solvated.parm7'))
        self.rst7_file = os.path.abspath(os.path.join(workdir, 'solvated.rst7'))
        self.topo_dir = os.path.abspath(workdir)
        self.water_model = water_model
        self.boxtype = boxtype
        self.boxbuffer = boxbuffer

        self._load_ff(ff=ff, workdir=workdir)

        leapin = ctools.get_tleap_input(
            component_list=[self],
            output_prefix='solvated',
            water_model=water_model,
            boxtype=boxtype,
            boxbuffer=boxbuffer,
            neutralize=neutralize,
            )
        with ctools.DirManager(workdir):
            with open('tleap.solvated.in', 'w') as outfile:
                outfile.write('\n'.join(leapin))
            self.tleap_solvated_in = os.path.abspath('tleap.solvated.in')
            tleap_cmd = [f"{os.path.join(const.AMBERBIN, 'tleap')}", '-f', 'tleap.solvated.in', '> /dev/null 2>&1']
            tleap_status = os.system(' '.join(tleap_cmd))
            if tleap_status != 0:
                print(time.strftime("%c"))
                print(f"tleap failed for {self}!")
                sys.exit(1)

        # Add certain salt concentration
        if not ionconc == 0.0:
            prmtop = self.parm7_file
            print(time.strftime("%c"))
            print(f"Add ionconc for {prmtop}.")
            struct = pmd.load_file(prmtop, xyz=None)
            nwaters = len(struct[':WAT'].residues)
            density = 55.5
            nions = round(int(nwaters)/density*ionconc)

            leapin = ctools.get_tleap_input(
                component_list=[self],
                output_prefix='solvated',
                water_model=water_model,
                boxtype=boxtype,
                boxbuffer=boxbuffer,
                neutralize=neutralize,
                num_cations=nions,
                num_anions=nions,
                )
            with ctools.DirManager(workdir):
                with open('tleap.solvated.ionconc.in', 'w') as outfile:
                    outfile.write('\n'.join(leapin))
                self.tleap_solvated_in = os.path.abspath('tleap.solvated.ionconc.in')
                tleap_cmd = [f"{os.path.join(const.AMBERBIN, 'tleap')}", '-f', 'tleap.solvated.ionconc.in', '> /dev/null 2>&1']
                tleap_status = os.system(' '.join(tleap_cmd))
                if tleap_status != 0:
                    print(time.strftime("%c"))
                    print(f"tleap failed for {self}!")
                    sys.exit(1)        

        if not ctools.file_check(self.parm7_file, size_threshold=1e3):
            print(time.strftime("%c"))
            print(f"tleap failed to prepare AmberParm system for {self}.")
            sys.exit(1)

        # Hydrogen Mass Repartitioning
        # https://pubs.acs.org/doi/abs/10.1021/ct5010406
        if hmr:
            prmtop = self.parm7_file
            # mol = pmd.amber.AmberParm(prmtop)
            mol = pmd.load_file(prmtop, xyz=None)
            print(time.strftime("%c"))
            print(f"Use HMR for {prmtop}.")
            os.rename(prmtop, prmtop.replace('.parm7', '.nohmr.parm7'))
            pmd.tools.actions.HMassRepartition(mol, hmass).execute()
            # TODO Bugfix
            # mol.save will add the following section into parm7 file:
            # """
            # %FLAG FORCE_FIELD_TYPE
            # %FORMAT(i2,a78)
            # 1        CHARMM force field: No FF information parsed...
            # """
            # Then, this kind of parm7 after HMR operation will turn on
            # CHARMM CMAP info, which will be not incompatible with TI:
            # """
            # ERROR: TI is incompatible with CHARMM for now.
            # """
            # Currently, If turn on HMR, only ff14SB can be used here.
            mol.save(self.parm7_file, overwrite=True, format='amber')
            del mol # flush the buffer

        print(time.strftime("%c"))
        print(f"Finished preparing AmberParm system for {self}.")
        self.prepared = True

class Protein(_BCS_Base):
    """Child class for constructing Protein in complex system.
    """
    def __init__(self, coord_file):
        super(Protein, self).__init__(coord_file)

    def __repr__(self):
        return f"<Protein {self.molname}>"

    def add_ff(self, source, lib=None, dat=None):
        self._add_ff_info(source=source)
        if lib:
            self._add_ff_info(lib=lib)
        if dat:
            self._add_ff_info(dat=dat)

class Complex(object):
    """Complex class for preparing the protein-ligand-cofactors complex system.
    """
    def __init__(self, protein, ligand, cofactors):
        self.protein = protein
        self.ligand = {}
        self.ligand[ligand.molname] = ligand
        ligand.complex = self
        self.cofactors = {c.molname: c for c in cofactors}
        self.components = [ligand] + [protein] + cofactors
        self.compname = f'{protein.molname}_{ligand.molname}'
        self.prepared = False

    def __repr__(self):
        return f'<Complex {self.compname}>'

    def prep_complex_ff(self, workdir, protein_ff, water_model, boxtype,
                        boxbuffer, neutralize=True, hmr=False, hmass=3.024,
                        phosaa_ff=None, lipid_ff=None, ionconc=0.0,
                        protein_type='soluble', box_info=[]):
        """Prepare the force field for complex system.
        """
        self.parm7_file = os.path.abspath(os.path.join(workdir, 'complex.parm7'))
        self.rst7_file = os.path.abspath(os.path.join(workdir, 'complex.rst7'))
        self.topo_dir = os.path.abspath(workdir)

        self.protein_ff = protein_ff
        self.phosaa_ff = phosaa_ff
        self.lipid_ff = lipid_ff
        self.water_model = water_model
        self.boxtype = boxtype
        self.boxbuffer = boxbuffer
        self.protein_type = protein_type
        self.protein.add_ff(source=os.path.join(const.LEAPHOME, 'cmd', const.PRO_DEFS[protein_ff]['source']))

        for _, cof in self.cofactors.items():
            if not (cof.use_gaff2):
                error_msg = f"Please parameterize {cof} first."
                raise RuntimeError(error_msg)
        for _, lig in self.ligand.items():
            if not (lig.use_gaff2):
                error_msg = f"Please parameterize {lig} first."
                raise RuntimeError(error_msg)

        leapin = ctools.get_tleap_input(
            self.components, 'complex',
            water_model, boxtype, boxbuffer,
            neutralize=neutralize, phosaa_ff=phosaa_ff,
            lipid_ff=lipid_ff, protein_type=protein_type,
            box_info=box_info,
            )

        with ctools.DirManager(workdir):
            with open('tleap.complex.in', 'w') as outfile:
                outfile.write('\n'.join(leapin))
            self.tleap_complex_in = os.path.abspath('tleap.complex.in')
            tleap_cmd = [f"{os.path.join(const.AMBERBIN, 'tleap')}", '-f', \
                'tleap.complex.in', '> /dev/null 2>&1']
            tleap_status = os.system(' '.join(tleap_cmd))
            if tleap_status != 0:
                print(time.strftime("%c"))
                print(f"tleap failed for {self}!")
                sys.exit(1)

        # Add certain salt concentration
        if not ionconc == 0.0 and protein_type=='soluble':
            prmtop = self.parm7_file
            print(time.strftime("%c"))
            print(f"Add ionconc for {prmtop}.")
            struct = pmd.load_file(prmtop, xyz=None)
            nwaters = len(struct[':WAT'].residues)
            density = 55.5
            nions = round(int(nwaters)/density*ionconc)

            leapin = ctools.get_tleap_input(
                self.components, 'complex',
                water_model, boxtype, boxbuffer,
                neutralize=neutralize,  phosaa_ff=phosaa_ff,
                lipid_ff=lipid_ff, num_cations=nions, num_anions=nions
                )
            with ctools.DirManager(workdir):
                with open('tleap.complex.ionconc.in', 'w') as outfile:
                    outfile.write('\n'.join(leapin))
                self.tleap_complex_in = os.path.abspath('tleap.complex.ionconc.in')
                tleap_cmd = [f"{os.path.join(const.AMBERBIN, 'tleap')}", '-f', \
                    'tleap.complex.ionconc.in', '> /dev/null 2>&1']
                tleap_status = os.system(' '.join(tleap_cmd))
                if tleap_status != 0:
                    print(time.strftime("%c"))
                    print(f"tleap failed for {self}!")
                    sys.exit(1)  
        
        if not ctools.file_check(self.parm7_file, size_threshold=1e3):
            print(time.strftime("%c"))
            print(f"tleap failed to prepare AmberParm system for {self}.")
            sys.exit(1)

        # Hydrogen Mass Repartitioning
        if hmr:
            prmtop = self.parm7_file
            mol = pmd.load_file(prmtop, xyz=None)
            print(time.strftime("%c"))
            print(f"Use HMR for {prmtop}.")
            os.rename(prmtop, prmtop.replace('.parm7', '.nohmr.parm7'))
            pmd.tools.actions.HMassRepartition(mol, hmass).execute()
            mol.save(self.parm7_file, overwrite=True, format='amber')
            del mol # flush the buffer

        print(time.strftime("%c"))
        print(f"Finished preparing AmberParm system for {self}.")
        self.prepared = True