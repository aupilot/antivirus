import os
import random
import subprocess

import cma
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.SeqUtils import seq3
from Bio.SeqRecord import SeqRecord
from scipy.special import softmax
import shutil

import matplotlib
matplotlib.use('TKAgg')


spike = "7cr5_SPIKE.pdb"

# use Fv only 7cr5
# we can align with ANARCY - add spacers. Do we want it?
# http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/
initial_H = "QVQLVESGGGVVQPGRSLRLSC AASGFTFSSYIMH WVRQAPGKGLEWVA VISYDGSNEA YADSVKGRFTISRDNSKNTLYLQMSSLRAEDTGVYYC ARETGDYSSSWYDS WGRGTLVTVSS"
initial_L = "QLVLTQSPSASASLGASVKLTC TLSSGHSNYAIA WHQQQPEKGPRYLM KVNSDGSHTKGD GIPDRFSGSSSGAERYLTISSLQSEDEADYYC QTWGTGIQV FGGGTKLTVL"

# we split the sequence to cdr/framework regions. we won't  optimise on constant framework
# TODO: automate cdr detection
# https://github.com/mit-ll/Insilico_Ab_Variant_Generator/blob/main/scripts/parse_region.py
framework_H1 = "QVQLVESGGGVVQPGRSLRLSC"
cdr_H1 = "AASGFTFSSYIMH"
# cdr_H1 = "AASGFTF----SSYIMH"
framework_H2 = "WVRQAPGKGLEWVA"
cdr_H2 = "VISYDGSNEA"
# cdr_H2 = "VISYD--GSNEA"
framework_H3 = "YADSVKGRFTISRDNSKNTLYLQMSSLRAEDTGVYYC"
cdr_H3 = "ARETGDYSSSWYDS"
framework_H4 = "WGRGTLVTVSS"

framework_L1 = "QLVLTQSPSASASLGASVKLTC"
cdr_L1 = "TLSSGHSNYAIA"
# cdr_L1 = "TLSSGHS-----NYAIA"
framework_L2 = "WHQQQPEKGPRYLM"
cdr_L2 = "KVNSDGSHTKGD"
# cdr_L2 = "KVNSD---GSHTKGD"
framework_L3 = "GIPDRFSGSSSGAERYLTISSLQSEDEADYYC"
cdr_L3 = "QTWGTGIQV"
# cdr_L3 = "QTWGT----GIQV"
framework_L4 = "FGGGTKLTVL"


residue_letters = [
    "-",    # spacer
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


# h_len_opt = len(cdr_H1) + len(cdr_H2) + len(cdr_H3)
# l_len_opt = len(cdr_L1) + len(cdr_L2) + len(cdr_L3)

# TODO: learn better embedding with lower dimensions and smooth space. Perhaps 2-3 layer net?
residue_embedding = np.eye(21,21,dtype=int)

def residue2vector(res):
    return residue_embedding[residue_letters.index(res)]

def np2seq(emb):
    seq = ""
    for a in emb:
        new_letter = residue_letters[softmax(a).argmax()]
        if new_letter != '-':
            seq = seq + new_letter
    return seq

def seq2np(seq):
    out = []
    for res in seq:
        out.append(residue2vector(res))
    return  np.array(out)

def hl2np(h, l):
    h_np = seq2np(h)
    l_np = seq2np(l)
    return np.vstack((h_np, l_np)).flatten()

def np2full_seq(emb):
    hl_np = np.reshape(emb, (-1,21))
    hl_seq = np2seq(hl_np)

    ptr1 = len(cdr_H1)
    ptr2 = ptr1 + len(cdr_H2)
    ptr3 = ptr2 + len(cdr_H3)
    H = framework_H1 + hl_seq[0:ptr1] + \
        framework_H2 + hl_seq[ptr1:ptr2] + \
        framework_H3 + hl_seq[ptr2:ptr3] + framework_H4

    ptr0 = ptr3
    ptr1 = ptr0 + len(cdr_L1)
    ptr2 = ptr1 + len(cdr_L2)
    ptr3 = ptr2 + len(cdr_L3)
    L = framework_L1 + hl_seq[ptr0:ptr1] + \
        framework_L2 + hl_seq[ptr1:ptr2] + \
        framework_L3 + hl_seq[ptr2:ptr3] + framework_L4

    return H,L


def test_seq():
    my_seq = "GATCG"
    embedded = seq2np(my_seq)
    seq = np2seq(embedded)
    print(seq)

def test_full_seq():
    print(np2full_seq(seq2np(cdr_H1 + cdr_H2 + cdr_H3 + cdr_L1 + cdr_L2 + cdr_L3)))


# def test_initial():
#     x0 = hl2np(initial_H, initial_L)
#     H,L = np2hl(x0)
#     print(H)
#     print(L)

def check_stop(a):
    # how do we stop the optimiser correctly??
    if a.countiter > 10000:
        print(a.countiter)

def get_fitness(x):
    HL = np2full_seq(x)
    sequences = [SeqRecord(Seq(HL[0]), id='H', description="Optimiser Sample H"), SeqRecord(Seq(HL[1]), id="L", description="Optimiser Sample L")]
    SeqIO.write(sequences, "./data/fitness.fasta", "fasta")

    # run AlphaFold
    output = subprocess.run(["./run_alpha.sh", "./data/fitness.fasta"], capture_output=True, check=True)

    # copy alphafold results to data dir
    os.system("cp /tmp/alphafold/fitness/renamed_* ./data/")

    # run docking/score
    best_score = 999
    for i in range(5):
        output = subprocess.run(["./run_score.sh", f"/workdir/renamed_{i}.pdb", "/workdir/" + spike], capture_output=True, check=True)
        score = float(output.stdout.split()[-1])
        if score < best_score:
            best_score = score

    print(f"Best score: {best_score}")

    return best_score

if __name__ == '__main__':
    # test_full_seq()
    # exit()

    os.mkdir("data")

    x0 = seq2np(cdr_H1 + cdr_H2 + cdr_H3 + cdr_L1 + cdr_L2 + cdr_L3).flatten()

    fun = get_fitness
    # fun = cma.ff.rosen  # we could use `functools.partial(cma.ff.elli, cond=1e4)` to change the condition number to 1e4

    sigma0 = 0.3  # initial standard deviation to sample new solutions - should be ~ 1/4 of range

    # # cfun = cma.ConstrainedFitnessAL(fun, constraints)  # unconstrained function with adaptive Lagrange multipliers
    # es = cma.CMAEvolutionStrategy(x0, sigma0)
    #
    # while not es.stop():
    #     X = es.ask()  # sample len(X) candidate solutions
    #     es.tell(X, [fun(x) for x in X])
    #     # cfun.update(es)
    #     es.logger.add()  # for later plotting
    #     es.disp()
    # es.result_pretty()

    res, es = cma.fmin2(fun, x0, sigma0, callback=check_stop,
                   options={
                            'ftarget': -3.0,
                            'popsize': 5,
                            'maxiter': 12,
                            'bounds': [-0.1, 1.1],
                            'verb_time':0,
                            'verb_disp': 500,
                            'seed': 3},
                   restarts=1)


    print(np2full_seq(res))

    es.plot()

    es.result_pretty()