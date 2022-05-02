#!/usr/bin/python3

# igfold >= 0.0.8
# run:
# python ig_fold.py out_file heavy_chain light_chain --rosetta=0
# python ig_fold.py igfold.pdb QVQLVESGGGVVQPGRSLRLS QLVLTQSPSASASLGASVKL

import argparse
from igfold import IgFoldRunner

use_rosetta = False

def get_args():
    """Gets command line arguments"""
    desc = ('''
        Folding with IgFold.
        ''')
    ### positional args
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("out", type=str, help="""
        Output Antibody PDB.
    """)
    parser.add_argument("heavy", type=str, help="""
        heavy chain.
    """)
    parser.add_argument("light", type=str, help="""
        light chain.
    """)
    ### optional args
    parser.add_argument("--rosetta", type=int, default=0, help="""
        Use rosetta for refinement (1) or OpenMM (0).
    """)
    parser.add_argument("--renum", type=int, default=0, help="""
        Send predicted structure to AbNum server for Chothia renumbering (1)
    """)
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()

    use_rosetta = args.rosetta
    if use_rosetta == 1:
        from igfold import IgFoldRunner, init_pyrosetta
        init_pyrosetta()
        use_openmm = False
    else:
        use_openmm = True

    if args.renum == 1:
        renumber = True
    else:
        renumber = False

    sequences = {
        "H": args.heavy,
        "L": args.light,
    }

    igfold = IgFoldRunner()
    out = igfold.fold(
        args.out,  # Output PDB file
        sequences=sequences,  # Antibody sequences
        do_refine=True,  # Refine the antibody structure with PyRosetta or OpenMM
        do_renum=renumber,  # Send predicted structure to AbNum server for Chothia renumbering
	    use_openmm=use_openmm,
    )
