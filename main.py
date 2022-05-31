import argparse
import os
import random
import subprocess
import time
from multiprocessing import Pool
import cma
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
import matplotlib
from Bio.PDB import *
from numpy.linalg import norm

from embed import Embed
from helper import highlight_differences, insert_spacers, align

matplotlib.use('TKAgg')
from ttictoc import tic, toc
# from igfold import IgFoldRunner, init_pyrosetta
# from igfold import IgFoldRunner
from prody import parsePDB, writePDB



###################################
# These are set thru arguments now
dla_threshold = 0.06
mega_type = 1           # 0 - original, 1 - kir's, ...
use_rosetta = 0
renumber = 0
#############################3#####

# !!! we need to manually roughly align the spike over a typical Ab Fv
spike = "7k9i_spike_aligned.pdb"
spike2 = "7mem_spike_aligned.pdb"

aligned_over = "alignment.pdb"

# the max distance between the atom and the line from H-end and L-end is about 38Å
# so, we will block the atoms that are closer than about 1/3 of that
block_distance = 38.0 / 3   # Å


# start point Ab and the max length of the Ab chains including spacers and paddings
# TODO: use ANARCI to insert spacers to fasta
starting_point = "7lm9-Fv.fasta"

# we calc this as max len for now in main()
chain_max_length = 0

movie_cnt = 0

# folded ab name.
ig_fold_pdb = "ig_fold.pdb"
ig_fold_aligned_pdb = "ig_fold_aligned.pdb"

# we split the sequence to cdr/framework regions. we won't  optimise on constant framework
# TODO: automate cdr detection ??
# all give different CDR breakdown!
# http://dunbrack2.fccc.edu/PyIgClassify/User/UserPdb.aspx
# http://cao.labshare.cn/AbRSA/download.php
# https://github.com/mit-ll/Insilico_Ab_Variant_Generator/blob/main/scripts/parse_region.py
# inserts:
# http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/


# embedder class. keep global
emb = Embed("esm1b_t33_650M_UR50S")

# residue_letters = [
#     "-",  # spacer
#     "A",
#     "C",
#     "D",
#     "E",
#     "F",
#     "G",
#     "H",
#     "I",
#     "K",
#     "L",
#     "M",
#     "N",
#     "P",
#     "Q",
#     "R",
#     "S",
#     "T",
#     "V",
#     "W",
#     "Y",
# ]


# TODO: learn better embedding with lower dimensions and smooth space. Perhaps 2-3 layer net?
residue_embedding = np.eye(21, 21, dtype=int)

global_best_score = 999.0


# input arg: dict of Seq 'H', 'L'
# len must be less than 1022 (ESM embedding limit)
def seq2np(seq_dict):
    return emb.embed(seq_dict)


# we do not separate CDRs here. We rely on the embedding as it includes the language model.
# TODO: add auto CDR detection and candidates pruning inside the optimiser to increase the prob changes in CRDs rather than elsewhere??
def np2full_seq(eee):
    seqs = emb.de_embed(eee)

    H = seqs['H']
    L = seqs['L']

    # remove spacers if any
    H = H.replace("-", "")
    L = L.replace("-", "")

    return H, L


# do not remove "-" for demonstration only
def np2seq_show(eee):
    seqs = emb.de_embed(eee)
    H = seqs['H']
    L = seqs['L']
    return H, L


# https://stackoverflow.com/questions/39840030/distance-between-point-and-a-line-from-two-points
def t(p, q, r):
    x = p-q
    return np.dot(r-q, x)/np.dot(x, x)

def d(p, q, r):
    return np.linalg.norm(t(p, q, r)*(p-q)+q-r)

# Automatically mark blocked residues based on distance from the line between end of H to the end of L.
# This requires pdb file as input.
# It saves results only to one location. Thus, if we want split between threads, it must be updated accordingly
def save_blocking_positions_pdb(pdb_file):
    ff_h = open("data/block-H.txt", "wt")
    ff_l = open("data/block-L.txt", "wt")

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


def fake_fitness(arg):
    return random.random()


def dock_score(thread_no: int):
    average_score = 0.0
    best_score = 999.0
    for i in range(5):
        output = subprocess.run(["./run_score.sh", f"/workdir/th.{thread_no}/renamed_{i}.pdb", "/workdir/" + spike],
                                capture_output=True, check=True)  # these paths are inside the container!
        score = float(output.stdout.split()[-1])
        average_score = average_score + score
        if score < best_score:
            best_score = score
    average_score = average_score / 5.
    return (best_score, average_score)


# with containerised Rosetta or OpenMM refinement running in docker container
def ig_fold_docker(thread_no: int, xx):
    HL = np2full_seq(xx)

    try:
        output = subprocess.run(
            ["./run_igfold.sh", f"./data/th.{thread_no}/{ig_fold_pdb}", HL[0], HL[1], f"--rosetta={use_rosetta}",
             f"--renum={renumber}"],
            capture_output=True, check=True)
    except:
        # we run it again in the case of failure. Specifically designed to workaround OpenMM NaN exception for fuck knows reason
        print("Exception occured in run_igfold.sh. Trying to re-run again")
        output = subprocess.run(
            ["./run_igfold.sh", f"./data/th.{thread_no}/{ig_fold_pdb}", HL[0], HL[1], f"--rosetta={use_rosetta}",
             f"--renum={renumber}"],
            capture_output=True, check=True)

    # now we need to rotate the model in PDB to have it aligned always the same way
    align(reference=aligned_over, sample=f"./data/th.{thread_no}/{ig_fold_pdb}", output=f"./data/th.{thread_no}/{ig_fold_aligned_pdb}")

    # once folded, we need to select the residues to block. They will be the same most the time for all threads,
    # so we don't bother to run it multiple times. Only in thread 0
    if thread_no == 0:
        save_blocking_positions_pdb(f"./data/th.{thread_no}/{ig_fold_aligned_pdb}")


# input - list of 2 samples
def double_fun_igfold(X):
    global global_best_score
    global movie_cnt
    tic()

    ##### these 2 calls can be run in parallel??
    thread_numbers = (0, 1)
    with Pool(2) as pool:
        pool.starmap(ig_fold_docker, zip(thread_numbers, X))
    # thread_no = 0
    # ig_fold_openmm(thread_no, X[thread_no])
    # thread_no = 1
    # ig_fold_openmm(thread_no, X[thread_no])

    ##### the following must run in sequence!
    best_score = [10., 10.]
    for thread_no in [0,1]:
        # run docking/score for AF thread 0
        output = subprocess.run(
            ["./run_score.sh", f"/workdir/th.{thread_no}/{ig_fold_aligned_pdb}", "/workdir/" + spike, f"{dla_threshold}",
             f"{mega_type}"],
            capture_output=True, check=True)  # these paths are inside the container!
        best_score[thread_no] = float(output.stdout.split()[-1])
        # optionally copy the best pdb to save it
        if best_score[thread_no] < global_best_score:
            shutil.copy(f"./data/th.{thread_no}/{ig_fold_aligned_pdb}", f"./data/best_{movie_cnt:02d}.pdb")
            global_best_score = best_score[thread_no]
            movie_cnt += 1

    # # run docking/score for AF thread 1
    # thread_no = 1
    # output = subprocess.run(
    #     ["./run_score.sh", f"/workdir/th.{thread_no}/{ig_fold_aligned_pdb}", "/workdir/" + spike, f"{dla_threshold}",
    #      f"{mega_type}"],
    #     capture_output=True, check=True)  # these paths are inside the container!
    # best_score_1 = float(output.stdout.split()[-1])
    # # optionally copy the best pdb to save it
    # if best_score_1 < global_best_score:
    #     shutil.copy(f"./data/th.{thread_no}/{ig_fold_aligned_pdb}", f"./data/best_{movie_cnt:02d}.pdb")
    #     global_best_score = best_score_1
    #     movie_cnt += 1

    print(f"Scores: {best_score[0]:.4f} {best_score[1]:.4f}, The best: {global_best_score:.4f}, Time: {toc():.1f}")

    return (best_score[0], best_score[1])


def get_args():
    """Gets command line arguments"""
    desc = ('''
        Ab optimiser.
        ''')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("dla", type=float, default=0.06, help="""
        DLA-Ranker thereshold.
    """)
    parser.add_argument("--mega", type=int, default=0, help="""
        Use kir optimised Megadock (1) or not (0).
    """)
    parser.add_argument("--rosetta", type=int, default=0, help="""
        Use Rosetta for IgFold refinement (1) or OpenMM (0).
    """)
    parser.add_argument("--renum", type=int, default=0, help="""
        Send predicted structure to AbNum server for Chothia renumbering (1)
    """)
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    dla_threshold = args.dla
    mega_type = args.mega
    use_rosetta = args.rosetta
    renumber = args.renum

    print(time.asctime())

    start_seq = {}
    for record in SeqIO.parse(starting_point, "fasta"):
        start_seq[record.id] = record.seq
        if len(record) > chain_max_length:
            chain_max_length = len(record)
        print(f"Sequence length of {record.id}: {len(record)}")

    print("Adding spacers with ANARCI")
    start_seq_spacers = insert_spacers(start_seq)

    if not os.path.exists("data"):
        os.mkdir("data")

    x0 = seq2np(start_seq_spacers).flatten()

    fun = double_fun_igfold

    plot_avg = []
    plot_min = []
    sigma0 = 0.05  # initial standard deviation to sample new solutions - should be ~ 1/4 of range???

    # # cfun = cma.ConstrainedFitnessAL(fun, constraints)  # unconstrained function with adaptive Lagrange multipliers
    es = cma.CMAEvolutionStrategy(x0, sigma0,
                                  inopts={
                                      'ftarget': -5.0,
                                      'popsize': 18,        # must be even as we return pairs of target solutions
                                      'maxiter': 28,
                                      'bounds': [-80., 80.],
                                      'verb_time': 0,
                                      'verb_disp': 500,
                                      'seed': 3}, )

    while not es.stop():
        X = es.ask()  # sample len(X) candidate solutions

        V = []
        for i in range(len(X) // 2):
            x2 = (X[i * 2], X[i * 2 + 1])
            v2 = fun(x2)
            V.append(v2[0])
            V.append(v2[1])

        plot_avg.append(np.array(V).mean())
        plot_min.append(np.array(V).min())

        es.tell(X, V)
        # es.tell(X, [fun(x) for x in X])

        # compare orig and current
        sq1 = start_seq_spacers['H'] + ' ' + start_seq_spacers['L']
        sq2 = np2seq_show(es.result.xfavorite)[0] + ' ' + np2seq_show(es.result.xfavorite)[1]
        print(sq1)
        print(sq2)
        print(highlight_differences(sq1, sq2))

        es.logger.add()  # for later plotting
        es.disp()

    es.result_pretty()
    # === es.result explanations: https://github.com/CMA-ES/pycma/blob/9b4eb5450c020ac99d637780a42c39788f1e1045/cma/evolution_strategy.py#L977

    # res = np2full_seq(es.result.xbest)
    # the result as the mean of the metamodel gaussian rather than the best candidate:
    res = es.result.xfavorite
    seq = np2full_seq(res)

    sequences = [SeqRecord(Seq(seq[0]), id='H', description="Optimised Ab H"),
                 SeqRecord(Seq(seq[1]), id="L", description="Optimised Ab L")]
    SeqIO.write(sequences, "best.fasta", "fasta")

    print(seq)
    print(f"The best score: {es.result.fbest}, iterations: {es.result.iterations}")
    print(time.asctime())

    plot = np.array([plot_avg, plot_min]).T
    matplotlib.pyplot.plot(plot)
    matplotlib.pyplot.show(block=True)

    es.plot()
    matplotlib.pyplot.show(block=True)
