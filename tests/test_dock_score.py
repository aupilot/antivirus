import subprocess

from Bio.PDB import *
import numpy as np

block_distance = 38.0 / 3


def t(p, q, r):
    x = p-q
    return np.dot(r-q, x)/np.dot(x, x)


def d(p, q, r):
    return np.linalg.norm(t(p, q, r)*(p-q)+q-r)


def save_blocking_positions_pdb(pdb_file):
    ff_h = open("block-H.txt", "wt")
    ff_l = open("block-L.txt", "wt")

    last_coords = {}
    p = PDBParser()
    structure = p.get_structure('Ab', pdb_file)
    # bld = open('test.bld', 'wt')    # arrows file to plot in ChimeraX to test the selection

    # we're expecting one model with 2 chains H and L
    if len(structure.child_list) != 1:
        raise Exception('The PDB must have only one model with Ab')
    for model in structure:
        if len(model.child_list) != 2:
            raise Exception('The Ab Model must have two chains')
        for chain in model:
            # find the last element.
            *_, last = chain
            last_coords[chain.id] = last['CA'].get_coord()
        midpoint = (last_coords['H'] + last_coords['L']) / 2.

        for chain in model:
            for residue in chain:
                atom_coord = residue['CA'].get_coord()
                # distance between the line connecting end residues and the current atom
                distance = d(last_coords['H'], last_coords['L'], atom_coord)
                if distance < block_distance:
                    # print(f'block {chain.id} {residue.id[1]} {distance}')
                    # we create a bld file to indicate the blocked atoms in ChimeraX image
                    # bld.write(f'.arrow {midpoint[0]} {midpoint[1]} {midpoint[2]} {atom_coord[0]} {atom_coord[1]} {atom_coord[2]}\n')
                    if chain.id == 'H':
                        ff_h.write(f"{residue.id[1]},")
                    elif chain.id == 'L':
                        ff_l.write(f"{residue.id[1]},")
    # bld.close()
    ff_h.close()
    ff_l.close()


if __name__ == '__main__':

    dla_threshold = 0.05
    mega_type = 0

    # -0.52415 (Vox size 1.0A) / -0.62479
    # RECEPTOR = "7e3c_Fv.pdb"
    # LIGAND = "7e3c_spike.pdb"

    # crashes on vina
    # RECEPTOR = "7lm9_Fv.pdb"
    # LIGAND = "7lm9_spike.pdb"

    # no solution (Vox size 1.0A) / no solution
    # RECEPTOR = "6zer_Fv.pdb"
    # LIGAND = "6zer_spike.pdb"

    # -0.36079 (Vox size 1.0A) / -0.24376 (Vox size 0.8A)
    # RECEPTOR = "7ps7_Fv.pdb"
    # LIGAND = "7ps7_spike.pdb"

    # -0.44919 (Vox size 1.0A) / -0.73982 (Vox size 0.8A)
    RECEPTOR = "7urs_Fv.pdb"
    LIGAND = "7urs_spike.pdb"

    # pdb errors on naccess
    # RECEPTOR = "2dd8_Fv.pdb"
    # LIGAND = "2dd8_spike.pdb"

    save_blocking_positions_pdb(f"{RECEPTOR}")

    output = subprocess.run(
        ["./test_dock_score.sh", f"/workdir/{RECEPTOR}", f"/workdir/{LIGAND}", f"{dla_threshold}", f"{mega_type}"],
        capture_output=False, check=True)  # these paths are inside the container!
