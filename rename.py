import sys
from Bio.PDB import PDBList, PDBIO, PDBParser

"""
Renames chains in PDB file
1. Change the dict renames below
2. python3 rename.py data/7mem_spike_tmp.pdb data/7mem_spike_aligned.pdb
"""

pdbl = PDBList()
io = PDBIO()
parser = PDBParser(PERMISSIVE=1)

inp_file = sys.argv[1]
#infile = .path.split(inp_file)
out_file = sys.argv[2]

#pdbl.retrieve_pdb_file(infile[0], pdir=infile[1], file_format="pdb")
structure = parser.get_structure('kir', inp_file)

# renames = {
#     "B": "H",
#     "C": "L",
# }
renames = {
    "A": "H",
    "g": "A",
    "f": "A"
}

for model in structure:
    for chain in model:
        old_name = chain.get_id()
        new_name = renames.get(old_name)
        if new_name:
            # print(f"Renaming {inp_file} chain {old_name} to {new_name}")
            chain.id = new_name
        # else:
            # print(f"Keeping {inp_file} chain name {old_name}")

io.set_structure(structure)
io.save(out_file)
