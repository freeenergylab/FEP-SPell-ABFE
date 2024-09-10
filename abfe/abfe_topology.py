#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
==============================================================================
CopyRight (c) By freeenergylab.
@Description:
This is the topology module that constructs initial parm7 and rst7 files.

@Author: Pengfei Li
@Date: May 30th, 2022
==============================================================================
"""
import argparse
import os
import pickle
import socket
import time

import abfe.const as const
import abfe.biochemsystems as bcs
import abfe.utils.common_tools as ctools

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parse arguments for topology module!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    # construction arguments
    parser.add_argument(
        '-bw', '--abfe_workdir',
        help='Workdir'
        )
    parser.add_argument(
        '-cd', '--cofactor_dirs', nargs='*', default=[],
        help='List contains the cofactors, e.g. cofactors/cofactor_name'
        )
    parser.add_argument(
        '-ld', '--ligand_dir',
        help='Str contains the ligand name, e.g. ligands/ligand_name'
        )
    parser.add_argument(
        '-pd', '--protein_dir',
        help='Str contains the protein name, e.g. proteins/protein_name'
        )
    parser.add_argument(
        '-ho', '--hfe_only', action='store_true', default=False,
        help='Only perform HFE'
        )

    # input arguments
    parser.add_argument(
        '-pn', '--protein_name', nargs='*', default='',
        help='Str contains the protein name'
        )
    parser.add_argument(
        '-cn', '--cofactor_names', nargs='*', default=[],
        help='List contains the cofactor names'
        )
    parser.add_argument(
        '-pt', '--protein_type', default='soluble', choices=['soluble', 'membrane'],
        help='Str contains the protein type'
        )
    parser.add_argument(
        '-bt', '--box_type', default='rect', choices=['rect', 'oct'],
        help='Water box model'
        )
    parser.add_argument(
        '-lbb', '--lig_box_buffer', default=18.0,
        help='Ligand box buffer'
        )
    parser.add_argument(
        '-cbb', '--com_box_buffer', default=8.0,
        help='Complex box buffer'
        )
    parser.add_argument(
        '-nl', '--neutralize', action='store_true'
        )
    parser.add_argument(
        '-ic', '--ionconc', default=0.0, type=float,
        help='Ion concentration'
        )
    parser.add_argument(
        '-hr', '--hmr', action='store_true', default=False,
        help='Toggles hydrogen mass repartitioning functionality \
            to allow for larger timestep'
        )
    parser.add_argument(
        '-hm', '--hmass', default=3.024,
        help='New hydrogen mass'
        )

    # force field arguments
    parser.add_argument(
        '-cf', '--cofactor_ff', default='gaff2', choices=['gaff2'],
        help='Force field for the cofactor'
        )
    parser.add_argument(
        '-lf', '--ligand_ff', default='gaff2', choices=['gaff2'],
        help='Force field for the ligand'
        )
    parser.add_argument(
        '-pf', '--protein_ff', default='ff14SB', choices=const.PRO_DEFS.keys(),
        help='Force field for the protein'
        )
    parser.add_argument(
        '-sf', '--phosaa_ff', default='phosaa14SB', choices=const.PHOSAA_DEFS.keys(),
        help='Force field for the phosphorylated amino acid'
        )
    parser.add_argument(
        '-df', '--lipid_ff', default='lipid21',
        help='Force field for the lipid'
        )
    parser.add_argument(
        '-wf', '--water_ff', default='tip3p', choices=const.WATER_DEFS.keys(),
        help='Force field for the water'
        )

    return parser.parse_args()

if __name__ == "__main__":
    """Parse the input arguments, then run the topology stage.
    """
    args = _parse_args()
    print(time.strftime("%c"))
    print(f"Topology started on host {socket.gethostname()}.")

    # Construct the ligand object. One ligand for each job.
    _lig = os.path.join(args.ligand_dir, 'ligand.sdf')
    if not os.path.exists(_lig):
        _lig = os.path.join(args.ligand_dir, 'ligand.mol2')
    try:
        ligand = bcs.Ligand(_lig)
    except:
        error_msg = f"Error reading {_lig}."
        raise RuntimeError(error_msg)

    # Prepare the solvated ligand.
    lig_workdir = os.path.join(args.abfe_workdir, '_topology', ligand.molname)
    if args.ligand_ff == 'gaff2':
        ligand.prep_ligand_ff(
            ff=args.ligand_ff,
            workdir=lig_workdir,
            water_model=args.water_ff,
            boxtype=args.box_type,
            boxbuffer=args.lig_box_buffer,
            neutralize=args.neutralize,
            hmr=args.hmr,
            hmass=args.hmass,
            ionconc=args.ionconc,
            )
    else:
        error_msg = "Ligand can only be parameterized with gaff2 force field."
        raise ValueError(error_msg)

    if not args.hfe_only:
        # Construct the protein object. One protein for each job.
        _pro = os.path.join(args.protein_dir, 'protein.pdb')
        print(f'The protein type is {args.protein_type}.')
        if args.protein_type == "membrane":
            print(f"Read the box information from the {_pro} file.")
            box_info = ctools.parse_pdb(_pro)
            error_msg = f"Please provide the box info for membrane protein, "
            error_msg += f"the first line of {_pro} file starts with CRYST. e.g. "
            error_msg += f"CRYST1   81.971   82.209   88.609  90.00  90.00  90.00"
            assert box_info[0].startswith('CRYST'), error_msg
        else:
            box_info = []
        try:
            protein = bcs.Protein(_pro)
        except:
            error_msg = f"Error reading {_pro}."
            raise RuntimeError(error_msg)

        # Construct the cofactor objects. Each protein may contain multiple different cofactors.
        cofactors = []
        for cof_dir in args.cofactor_dirs:
            _cof = os.path.join(cof_dir, 'cofactor.sdf')
            if not os.path.exists(_cof):
                _cof = os.path.join(cof_dir, 'cofactor.mol2')
            try:
                cof = bcs.Cofactor(_cof)
                cofactors.append(cof)
            except:
                error_msg = f"Error reading {_cof}."
                raise RuntimeError(error_msg)

        # Construct the complex object.
        try:
            comp = bcs.Complex(protein, ligand, cofactors)
        except:
            error_msg = f"Error constructing complex: {protein}-{ligand}-{cofactors}."
            raise RuntimeError(error_msg)
        # Load cofactors with gaff2.
        for cof_name, cof in comp.cofactors.items():
            cof_workdir = os.path.join(args.abfe_workdir, '_topology', comp.compname, cof_name)
            cof.prep_cofactor_ff(ff=args.cofactor_ff, workdir=cof_workdir)

        # Prepare the solvated complex.
        comp_workdir = os.path.join(args.abfe_workdir, '_topology', comp.compname)
        lipid_ff = args.lipid_ff if len(args.lipid_ff.strip()) > 0 else None
        phosaa_ff = args.phosaa_ff if len(args.phosaa_ff.strip()) > 0 else None
        comp.prep_complex_ff(
            workdir=comp_workdir,
            protein_ff=args.protein_ff,
            water_model=args.water_ff,
            boxtype=args.box_type,
            boxbuffer=args.com_box_buffer,
            neutralize=args.neutralize,
            phosaa_ff=phosaa_ff,
            lipid_ff=lipid_ff,
            hmr=args.hmr,
            ionconc=args.ionconc,
            protein_type=args.protein_type,
            box_info=box_info,
            )

        # Dump the pkl files.
        comp.pkl_file = os.path.abspath(os.path.join(comp_workdir, 'objects.pkl'))
        with open(comp.pkl_file, 'wb') as pickf:
            pickle.dump(comp, pickf)

    # Dump the pkl files.
    ligand.pkl_file = os.path.abspath(os.path.join(lig_workdir, 'objects.pkl'))
    with open(ligand.pkl_file, 'wb') as pickf:
        pickle.dump(ligand, pickf)