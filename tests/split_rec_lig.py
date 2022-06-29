import sys

inp_file = sys.argv[1]
out_dir  = sys.argv[2]
out_file_rec = f"{out_dir}/rec.pdb"
out_file_lig = f"{out_dir}/lig.pdb"

# #pdbl.retrieve_pdb_file(infile[0], pdir=infile[1], file_format="pdb")
# structure = parser.get_structure('kir', inp_file)
#
# for model in structure:
#
#     for chain in structure.get_chains():
#         io.set_structure(chain)
#         io.save(chain.get_id() + ".pdb")


def split():
    """ Split a pdb file in different pdb files by chains: [H,L] and [A]
    """

    f = open(inp_file, "r")
    pdblines = f.readlines()
    f.close()

    # Chains selection :
    rec = []
    lig = []

    # Write the different files
    for line in pdblines:
        if line[0:4] == 'ATOM' or line[0:3] == 'TER':
            if line[21] == 'H' or line[21] == 'L':
                rec.append(line)
            else:
                # line[17:21] = 'LIG'
                lig.append(line)

    # fix TER if missing
    if rec[-1][0:3] != 'TER':
        rec.append('TER')
    rec.append('END')

    if lig[-1][0:3] != 'TER':
        lig.append('TER')
    lig.append('END')

    frec = open(out_file_rec, 'w')
    flig = open(out_file_lig, 'w')
    frec.writelines(rec)
    flig.writelines(lig)
    frec.close()
    flig.close()


if __name__ == '__main__':
    split()