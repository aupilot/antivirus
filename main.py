import argparse
import os
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
from embed_t5 import EmbedT5
from helper import highlight_differences, insert_spacers, align
from ttictoc import tic, toc
# from prody import parsePDB, writePDB

try:
    matplotlib.use('TKAgg')
except:
    print("Headless mode")



###################################
# These are overwritten thru arguments now
dla_threshold = 0.045
mega_type = 1           # 0 - original, 1 - kir's mega, ...
use_rosetta = 0         # 0 - use OpenMM, 1 - use rosetta. OpenMM seems better
renumber = 0
#############################3#####

# !!! we need to manually roughly align the spike over a typical Ab Fv. The chain must be renamed to "A" with rename.py
# spike_list = ["7e3c_spike_aligned.pdb", "7mem_spike_aligned.pdb"]
# spike_list = ["1sy6_epitope_aligned.pdb", "6jxr_epitope_aligned.pdb"]
spike_list = ["1sy6_epitope_aligned.pdb"]
aligned_over = "alignment.pdb"

# the max distance between the atom and the line from H-end and L-end is about 38Å
# so, we will block the atoms that are closer than about 1/3 of that
block_distance = 38.0 * 0.4   # Å , 0.4 from bottom does not add to megadock score. This prevents docking to the bottom


# start point Ab. The chains must be named H and L (no extra shit in names!)
# use http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/ to cut Fv region from a longer fasta sequence
# TODO: automate cutting the Ab Fv
# pdb_tofasta.py -multi '/home/kir/Documents/Antivirus/tests/7ps7_Fv.pdb' > /home/kir/Documents/Antivirus/7ps7.fasta
# starting_point = "2dd8_Fv.fasta"
# starting_point = "7e3c-Fv.fasta"
# starting_point = "1sy6-Fv.fasta"
starting_point = "7ps7-Fv.fasta"

# we calc this as max len for now in main()
# chain_max_length = 0

movie_cnt = 0

# folded ab name.
ig_fold_pdb = "ig_fold.pdb"
ig_fold_aligned_pdb = "ig_fold_aligned.pdb"

# we split the sequence to cdr/framework regions. we won't  optimise on constant framework
# automate cdr detection - all give different CDR breakdown!
# http://dunbrack2.fccc.edu/PyIgClassify/User/UserPdb.aspx
# http://cao.labshare.cn/AbRSA/download.php
# https://github.com/mit-ll/Insilico_Ab_Variant_Generator/blob/main/scripts/parse_region.py
# inserts:
# http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/

# embedder class. keep global
global emb

global_best_score = 999.0


# input arg: dict of Seq 'H', 'L'
# len must be less than 1022 (ESM embedding limit)
# TODO: Add some noise to embedding to make it less certain?
def seq2np(seq_dict):
    return emb.embed(seq_dict)


# we do not separate CDRs here. We rely on the embedding as it includes the language model.
# TODO: add auto CDR detection and candidates pruning inside the optimiser to increase the prob changes in CRDs rather than elsewhere??
def np2full_seq(eee):
    batch_size = 2
    eee = np.reshape(eee, (batch_size, -1, emb.latent_size))

    seqs = emb.de_embed(eee)

    H = seqs['H']
    L = seqs['L']

    # remove spacers if any
    H = H.replace("-", "")
    L = L.replace("-", "")
    H = H.replace("_", "")
    L = L.replace("_", "")

    return H, L


# do not remove "-" for demonstration only
def np2seq_show(eee):
    batch_size = 2
    eee = np.reshape(eee, (batch_size, -1, emb.latent_size))
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

    ##### these 2 calls can be run in parallel
    thread_numbers = (0, 1)
    # TODO: ocasionally we can get the sequence cut short. We need to avoid calling if_fold_docker in this cases as it will crash.
    with Pool(2) as pool:
        pool.starmap(ig_fold_docker, zip(thread_numbers, X))

    ##### the following must run in sequence! Perhaps unless a better computer
    lspl = len(spike_list)
    best_score_sep = np.zeros((lspl, len(thread_numbers)))
    best_score = [0.,0.]
    for thread_no in thread_numbers:
        # run docking/score for all spikes in the list
        for s, spike in enumerate(spike_list):
            # run dock & score to (multiple) paratopes
            output = subprocess.run(
                ["./run_score.sh", f"/workdir/th.{thread_no}/{ig_fold_aligned_pdb}", "/workdir/" + spike, f"{dla_threshold}",
                 f"{mega_type}", f"{score_system}"],
                capture_output=True, check=True)  # these paths are inside the container!
            best_score_sep[s,thread_no] = float(output.stdout.split()[-1])                          # TODO: perhaps there is a better way of combining scores

        best_score[thread_no] = np.sum(best_score_sep[:,thread_no]) / lspl                          # ----

        # optionally copy the best pdb to save it
        if best_score[thread_no] < global_best_score:
            shutil.copy(f"./data/th.{thread_no}/{ig_fold_aligned_pdb}", f"./data/best_{movie_cnt:02d}.pdb")
            global_best_score = best_score[thread_no]
            movie_cnt += 1

    for i in range(lspl):
        print(f"Target {i}: {best_score_sep[i,0]:.3f} {best_score_sep[i,1]:.3f}")
    print(f"The best: {global_best_score:.3f}, Time: {toc():.1f}")

    return best_score[0], best_score[1]


def get_args():
    """Gets command line arguments"""
    desc = ('''
        Ab optimiser.
        ''')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--sigma", type=float, default=0.05, help="""
        CMA-ES sigma.
    """)
    parser.add_argument("--popsize", type=int, default=18, help="""
        CMA-ES population size. Recommended 4 + int(3 * np.log(N)). Must be even as we run 2 threads
    """)
    parser.add_argument("--emb", type=int, default=0, help="""
        Embedding type: 0-> Facebook ESM, 1-> ProtTrans prot_t5_xl_half_uniref50-enc, 2-> ProtTrans prot_t5_xxl_uniref50, 3-> prot_t5_xl_bfd
    """)
    parser.add_argument("--maxiter", type=int, default=21, help="""
        DLA-Ranker thereshold.
    """)
    parser.add_argument("--dla", type=float, default=0.06, help="""
        DLA-Ranker thereshold. Set to 0.0 to disable DLA 
    """)
    parser.add_argument("--mega", type=int, default=0, help="""
        Use Kir's optimised Megadock (1) or the original one (0).
    """)
    parser.add_argument("--rosetta", type=int, default=0, help="""
        Use Rosetta for IgFold refinement (1) or OpenMM (0).
    """)
    parser.add_argument("--renum", type=int, default=0, help="""
        Send predicted structure to AbNum server for Chothia renumbering (1)
    """)
    parser.add_argument("--score", type=int, default=0, help="""
        Type of scoring. 0 - Vina, 1 - OnionNet. Onion could be used for screening, but does not look good for optimising. Also very slow
    """)

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    dla_threshold = args.dla
    mega_type = args.mega
    use_rosetta = args.rosetta
    renumber = args.renum
    score_system = args.score

    print(time.asctime())

    if args.emb == 0:
        print('Using ESM embedding esm1v_t33_650M_UR90S_5')
        emb = Embed("esm1v_t33_650M_UR90S_5")
    elif args.emb == 1:
        print('Using ProtT5 prot_t5_xl_half_uniref50-enc embedding ')
        emb = EmbedT5('prot_t5_xl_half_uniref50-enc')
    elif args.emb == 2:
        print('Using ProtT5 prot_t5_xxl_uniref50 embedding')
        emb = EmbedT5('prot_t5_xxl_uniref50')
    elif args.emb == 3:
        print('Using ProtT5 prot_t5_xl_bfd embedding')
        emb = EmbedT5('prot_t5_xl_bfd')
    else:
        raise Exception('Unknown Embedding type')

    start_seq = {}
    for record in SeqIO.parse(starting_point, "fasta"):
        start_seq[record.id] = record.seq
        # if len(record) > chain_max_length:
        #     chain_max_length = len(record)
        print(f"Sequence length of {record.id}: {len(record)}")

    print("Adding spacers with ANARCI")
    # it seems that adding spacers confuses much the encoding. So, we should either use moderate schema (e.g. Chothia) or no spacers at all??
    start_seq_spacers = insert_spacers(start_seq)

    if not os.path.exists("data"):
        os.mkdir("data")

    x0 = seq2np(start_seq_spacers).flatten()

    fun = double_fun_igfold

    print("Embedding/de-embedding sample:")
    highligted, changes = highlight_differences(str(start_seq['H']), np2seq_show(x0)[0])
    print(f"H ({changes:02d}): {highligted}")
    highligted, changes = highlight_differences(str(start_seq['L']), np2seq_show(x0)[1])
    print(f"L ({changes:02d}): {highligted}")

    plot_avg = []
    plot_min = []

    # with T5 embedding we must use sigma ~0.1, while with ECM its ~0.5
    sigma0 = args.sigma  # initial standard deviation to sample new solutions - should be ~ 1/4 of range???

    if args.emb == 0:
        limit = 80.0
    else:
        limit = 2.0
    # # cfun = cma.ConstrainedFitnessAL(fun, constraints)  # unconstrained function with adaptive Lagrange multipliers
    es = cma.CMAEvolutionStrategy(x0, sigma0,
                                  inopts={
                                      'CMA_diagonal': True,
                                      'ftarget': -10.0,                 # unreacheable
                                      'popsize': args.popsize,          # must be even as we return pairs of target solutions
                                      'maxiter': args.maxiter,
                                      'bounds': [-limit, limit],
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

        # compare orig and current. Make sure that the original FASTA has chains named H and L
        sq1h = start_seq_spacers['H'] # + ' ' + start_seq_spacers['L']
        sq2h = np2seq_show(es.result.xfavorite)[0] # + ' ' + np2seq_show(es.result.xfavorite)[1]
        sq1l = start_seq_spacers['L']
        sq2l = np2seq_show(es.result.xfavorite)[1]
        print(sq1h + " " + sq1l)
        # print(sq2h + " " + sq2l)
        highligted_h, changes_h = highlight_differences(sq1h, sq2h)
        highligted_l, changes_l = highlight_differences(sq1l, sq2l)
        print(f"{highligted_h} {highligted_l} Diff: {changes_h}+{changes_l}" )

        es.logger.add()  # for later plotting
        es.disp()
        # es.plot()

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
    np.savetxt("plot.csv", plot, delimiter=",")

    matplotlib.pyplot.plot(plot)
    matplotlib.pyplot.show(block=True)

    # es.plot()
    # matplotlib.pyplot.show(block=True)
