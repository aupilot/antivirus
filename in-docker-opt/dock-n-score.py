#!/usr/bin/python3

import os
import sys
import subprocess
import argparse
from os import listdir

import numpy as np
import random


##### this script sits inside a container

n_scored_decoys = 8
records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
receptor_chains = {'H', 'L', 'M', 'N'}

# all decoys with the score lower than this will not be scored
dla_ranker_threshold = 0.05

def megadock(receptor=None, legand=None):
    # read the parameters for blocking H
    blk_file = open("/workdir/data/block-H.txt")
    blk = blk_file.readline().rstrip('\n')
    tmp_name_H = f"/opt/var/{random.randint(10000, 99999)}.pdb"

    output = subprocess.run(["/opt/block3.py", receptor, "H", blk],  check=True,
                            universal_newlines = True,
                            stdout = subprocess.PIPE)

    tmp_file = open(tmp_name_H, 'w')
    buffer = output.stdout.splitlines()
    tmp_file.writelines("%s\n" % line for line in buffer)
    tmp_file.close()

    blk_file.close()

    # read the parameters for blocking L
    blk_file = open("/workdir/data/block-L.txt")
    blk = blk_file.readline().rstrip('\n')
    tmp_name_L = f"/opt/var/{random.randint(10000, 99999)}.pdb"

    output = subprocess.run(["/opt/block3.py", tmp_name_H, "L", blk],  check=True,
                            universal_newlines = True,
                            stdout = subprocess.PIPE)

    tmp_file = open(tmp_name_L, 'w')
    buffer = output.stdout.splitlines()
    tmp_file.writelines("%s\n" % line for line in buffer)
    tmp_file.close()

    blk_file.close()


    output = subprocess.run(["/opt/megadock-4.1.1/megadock-gpu","-R", tmp_name_L, "-L", legand, "-r 3600", "-v 1.0", "-o", "/opt/var/dock.out" ], capture_output=False, check=True)
    output = subprocess.run(["/opt/megadock-4.1.1/mega-post.sh", receptor, legand], capture_output=False, check=True)

    os.remove(tmp_name_L)
    os.remove(tmp_name_H)


def vina_score():
    tmp_receptor = "/opt/var/receptor.pdb"
    tmp_ligand = "/opt/var/ligand.pdb"

    decoy_list = listdir("/opt/var/good_decoys")

    scores = []
    decoy_no = 0
    for dec in decoy_list:
        if not dec.endswith("pdb"):
            continue
        if not dec.startswith("decoy"):
            continue

        decoy_no += 1
    # for d in range(n_scored_decoys):
    #     pdbfh = open(f"/opt/var/decoys/decoy.{d+1}.pdb", 'r')
        pdbfh = open(dec, 'r')
        chain_data = {}  # {chain_id: lines}
        prev_chain = None
        for line in pdbfh:
            if line.startswith(records):
                if line.startswith('TER'):
                    chain_data[line_chain].append(line)
                else:
                    line_chain = line[21]
                    if line_chain != prev_chain:
                        if line_chain not in chain_data:
                            chain_data[line_chain] = []
                        prev_chain = line_chain
                    chain_data[line_chain].append(line)
        pdbfh.close()

        lines = list()
        for chain_id in receptor_chains:
            try:
                lines = lines + chain_data[chain_id]
            except:
                ValueError(f"No chain {chain_id} in PDB file")

        with open(tmp_receptor, 'w') as fh:
            fh.write(''.join(lines))

        lines = list()
        for chain_id in sorted(chain_data.keys()):
            if chain_id in receptor_chains:
                continue
            lines = lines + chain_data[chain_id]

        if len(lines) == 0:
            Exception(ValueError(f"No ligand chains found in PDB file"))

        with open(tmp_ligand, 'w') as fh:
            fh.write(''.join(lines))

        #output = subprocess.run(["./run_adtools.sh"], capture_output=True, check=True)
        output = subprocess.run(["/opt/run_obabel.sh"], capture_output=True, check=True)

        out_text = output.stdout.decode().split("\n")
        affinity_lines  = [match for match in out_text if "Affinity" in match]

        affinity_txt = affinity_lines[0]
        print(f"Dec {decoy_no} {affinity_txt}")

        scores.append(float(affinity_txt.split()[1]))

    if len(scores) == 0:
        return 10.0
    else:
        return np.min(np.array(scores))


def get_args():
    """Gets command line arguments"""
    # project_path = os.path.abspath(os.path.join(deepab.__file__, "../.."))

    desc = ('''
        Dock and Score.
        ''')
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("receptor",
                        type=str,
                        help="""
        Antibody PDB.
    """)

    parser.add_argument("legand",
                        type=str,
                        help="""
        Spike PDB.
    """)


    # now = str(datetime.now().strftime('%y-%m-%d_%H:%M:%S'))

    # default_pred_dir = os.path.join(project_path, "pred_{}".format(now))
    # parser.add_argument("--pred_dir",
    #                     type=str,
    #                     default=default_pred_dir,
    #                     help="Directory where results should be saved.")

    # parser.add_argument("--use_gpu",
    #                     default=False,
    #                     action="store_true",
    #                     help="Run model prediction on GPU.")

    return parser.parse_args()


if __name__ == '__main__':

    args = get_args()

    # check if files exist
    if not os.path.exists(args.receptor):
        print("Receptor PDB does not exist!")
        exit(1)
    if not os.path.exists(args.legand):
        print("Legand PDB does not exist!")
        exit(1)

    os.makedirs("/opt/var",exist_ok=True)

    print('Docking:')
    print(args.receptor)
    print(args.legand)
    megadock(args.receptor, args.legand)

    print("Re-score with DLA")
    output = subprocess.run(["python3","/opt/DLA-Ranker/dla_ranker.py"], capture_output=False, check=True)

    print('Scoring with Vina:')
    min_score = vina_score()

    print(f"The best score is {min_score}")
