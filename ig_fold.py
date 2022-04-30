#!/usr/bin/python3

# run:
# python ig_fold.py out_file heavy_chain light_chain
# python ig_fold.py igfold.pdb "QVQLVESGGGVVQPGRSLRLS" "QLVLTQSPSASASLGASVKL"

import argparse


openmm = True

if openmm is True:
    from igfold import IgFoldRunner
else:
    from igfold import IgFoldRunner, init_pyrosetta


def get_args():
    """Gets command line arguments"""
    desc = ('''
        Dock and Score.
        ''')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("out", type=str, help="""
        Output Antibody PDB.
    """)
    parser.add_argument("h", type=str, help="""
        heavy chain.
    """)
    parser.add_argument("l", type=str, help="""
        light chain.
    """)
    return parser.parse_args()

if __name__ == '__main__':
    if openmm is False:
        init_pyrosetta()

    args = get_args()

    sequences = {
        "H": args.h,
        "L": args.l,
    }

    igfold = IgFoldRunner()
    out = igfold.fold(
        args.out,  # Output PDB file
        sequences=sequences,  # Antibody sequences
        do_refine=True,  # Refine the antibody structure with PyRosetta
        do_renum=True,  # Send predicted structure to AbNum server for Chothia renumbering
	    use_openmm=openmm,
    )
