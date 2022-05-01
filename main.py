import os
import random
import subprocess
import time
from multiprocessing import Pool
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
from ttictoc import tic,toc
# from igfold import IgFoldRunner, init_pyrosetta
from igfold import IgFoldRunner
from prody import parsePDB, writePDB

spike = "7cr5_SPIKE.pdb"
ig_fold_pdb = "ig_fold.pdb"

###############################
dla_threshold = 0.06
#############################3#

# use Fv only 7cr5
# we can align with ANARCY - add spacers. Do we want it?
# http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/
# initial_H = "QVQLVESGGGVVQPGRSLRLSC AASGFTFSSYIMH WVRQAPGKGLEWVA VISYDGSNEA YADSVKGRFTISRDNSKNTLYLQMSSLRAEDTGVYYC ARETGDYSSSWYDS WGRGTLVTVSS"
# initial_L = "QLVLTQSPSASASLGASVKLTC TLSSGHSNYAIA WHQQQPEKGPRYLM KVNSDGSHTKGD GIPDRFSGSSSGAERYLTISSLQSEDEADYYC QTWGTGIQV FGGGTKLTVL"

# we split the sequence to cdr/framework regions. we won't  optimise on constant framework
# TODO: automate cdr detection ??
# https://github.com/mit-ll/Insilico_Ab_Variant_Generator/blob/main/scripts/parse_region.py
framework_H1 = "QVQLVESGGGVVQPGRSLRLSC"   # 1-22
# cdr_H1 = "AASGFTFSSYIMH"
cdr_H1 = "AASGFTF----SSYIMH"
framework_H2 = "WVRQAPGKGLEWVA"
# cdr_H2 = "VISYDGSNEA"
cdr_H2 = "VISYD--GSNEA"
framework_H3 = "YADSVKGRFTISRDNSKNTLYLQMSSLRAEDTGVYYC"
cdr_H3 = "ARETGDYSSSWYDS"
framework_H4 = "WGRGTLVTVSS"

framework_L1 = "QLVLTQSPSASASLGASVKLTC"
cdr_L1 = "TLSSGHSNYAIA"
# cdr_L1 = "TLSSGHS-----NYAIA"
framework_L2 = "WHQQQPEKGPRYLM"
# cdr_L2 = "KVNSDGSHTKGD"
cdr_L2 = "KVNSD---GSHTKGD"
framework_L3 = "GIPDRFSGSSSGAERYLTISSLQSEDEADYYC"
# cdr_L3 = "QTWGTGIQV"
cdr_L3 = "QTWGT----GIQV"
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

# list of pieces to block when docking - four from H and four from L
# as long as the length of the Ab varies, we need to compare the new sequence every time with the framework pieces to block
block= ["LVESGGGVVQPGRSLR", "RQAPGKGLEW", "LQMSSLRAEDTGVYYC", "GTLVTV",
        "LTQSPSASASLG", "QQPEKGPR", "SSSGAERYLT", "FGGGTK"]
# block= ["LVESGGGVVQPGRSLRL", "WVRQAPGKGLEWV", "YADSVKGR", "GTLVTVSS",
#         framework_L1, framework_L2, framework_L3, framework_L4]

# TODO: learn better embedding with lower dimensions and smooth space. Perhaps 2-3 layer net?
residue_embedding = np.eye(21,21,dtype=int)


global_best_score = 999.0

def residue2vector(residue):
    return residue_embedding[residue_letters.index(residue)]


def np2seq(emb):
    seq = ""
    for a in emb:
        new_letter = residue_letters[softmax(a).argmax()]
        # if new_letter != '-':
        seq = seq + new_letter
    return seq


def seq2np(seq):
    out = []
    for res in seq:
        out.append(residue2vector(res))
    return np.array(out)


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

    # remove spacers
    H = H.replace("-", "")
    L = L.replace("-", "")

    return H, L


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
    SeqIO.write(sequences, f"./data/{get_fitness.n}_fitness.fasta", "fasta")
    get_fitness.n = get_fitness.n+1

    # run AlphaFold. To run multiple AFs set the thread_no to different ints!
    thread_no = 0
    output = subprocess.run(["./run_alpha.sh", "./data/fitness.fasta", f"{thread_no}"], capture_output=True, check=True)

    # copy alphafold results to data dir
    os.makedirs(f"./data/th.{thread_no}/", exist_ok=True)
    os.system(f"cp -f /tmp/alphafold/th.{thread_no}/fitness.{thread_no}/renamed_* ./data/th.{thread_no}/")

    # run docking/score
    average_score = 0
    best_score = 999
    for i in range(5):
        output = subprocess.run(["./run_score.sh", f"/workdir/th.{thread_no}/renamed_{i}.pdb", "/workdir/" + spike], capture_output=True, check=True)  # these paths are inside the container!
        score = float(output.stdout.split()[-1])
        average_score = average_score + score
        if score < best_score:
            best_score = score
    average_score = average_score / 5

    print(f"Best score: {best_score}, Average score {average_score}")

    # return average_score #best_score
    return best_score


# create a file with receptor residues to block (for MEGADOCK)
# we assume that the block seqs always exist!
def save_blocking_positions(sequence_H, sequence_L):
    ff = open("data/block-H.txt", "wt")
    # ff.write("H ")
    for iii in range(4):
        ptr_from = sequence_H.find(block[iii])
        if ptr_from < 0: raise Exception('Blocking sequence problem! (H)')
        ptr_to = ptr_from + len(block[iii])
        ff.write(f"{ptr_from+1}-{ptr_to},")
    ff.close()

    ff = open("data/block-L.txt", "wt")
    # ff.write("L ")
    for iii in range(4):
        ptr_from = sequence_L.find(block[iii+4])
        if ptr_from < 0: raise Exception('Blocking sequence problem! (L)')
        ptr_to = ptr_from + len(block[iii+4])
        ff.write(f"{ptr_from+1}-{ptr_to},")
    ff.close()


def fake_fitness(arg):
    return random.random()


def run_alpha(thread_no: int):
    # prepare thread
    HL = np2full_seq(X[thread_no])
    sequences = [SeqRecord(Seq(HL[0]), id='H', description="Optimiser Sample H"), SeqRecord(Seq(HL[1]), id="L", description="Optimiser Sample L")]
    SeqIO.write(sequences, f"./data/fitness.{thread_no}.fasta", "fasta")
    # SeqIO.write(sequences, f"./data/{get_fitness.n}_fitness.{thread_no}.fasta", "fasta")        # save temporary results

    # run AlphaFold.
    output = subprocess.run(["./run_alpha.sh", f"./data/fitness.{thread_no}.fasta", f"{thread_no}"], capture_output=True, check=True)

    # copy alphafold results to data dir
    os.makedirs(f"./data/th.{thread_no}/", exist_ok=True)
    os.system(f"cp -f /tmp/alphafold/th.{thread_no}/fitness.{thread_no}/renamed_* ./data/th.{thread_no}/")


def dock_score(thread_no: int):
    average_score = 0.0
    best_score = 999.0
    for i in range(5):
        output = subprocess.run(["./run_score.sh", f"/workdir/th.{thread_no}/renamed_{i}.pdb", "/workdir/" + spike], capture_output=True, check=True)  # these paths are inside the container!
        score = float(output.stdout.split()[-1])
        average_score = average_score + score
        if score < best_score:
            best_score = score
    average_score = average_score / 5.
    return (best_score, average_score)


# input - list of 2 samples
# folding with alphafold
def double_fun(X):
    # return (random.random(), random.random())     # test
    global global_best_score
    tic()

    ##### these 2 calls must be run in parallel
    # thread_no = 0
    # run_alpha(thread_no)
    # thread_no = 1
    # run_alpha(thread_no)

    # TODO: как X попадает в run_alpha() ??? какая-то хуйня
    thread_numbers = (0, 1)
    with Pool(2) as pool:
        pool.map(run_alpha, thread_numbers)
        # pool.close()  # do we need close/join when using context?
        # pool.join()

    ##### the following must run in sequence!

    # with Pool(2) as pool:
    #     result = pool.map(dock_score, thread_numbers)
    #
    # best_score_0, average_score_0 = result[0]
    # best_score_1, average_score_1 = result[1]

    # run docking/score for AF thread 0
    thread_no = 0
    average_score = 0.
    best_score = 999.
    best_score_idx = 999
    for i in range(5):
        output = subprocess.run(["./run_score.sh", f"/workdir/th.{thread_no}/renamed_{i}.pdb", "/workdir/" + spike], capture_output=True, check=True)  # these paths are inside the container!
        score = float(output.stdout.split()[-1])
        average_score = average_score + score
        if score < best_score:
            best_score = score
            best_score_idx = i
    average_score = average_score / 5.
    best_score_0 = best_score
    average_score_0 = average_score

    # optionally copy the best pdb to save it
    if best_score_idx != 999 and best_score < global_best_score:
        shutil.copy(f"./data/th.{thread_no}/renamed_{best_score_idx}.pdb", "./data/best.pdb")
        global_best_score = best_score

    # run docking/score for AF thread 1
    thread_no = 1
    average_score = 0.
    best_score = 999.
    best_score_idx = 999
    for i in range(5):
        output = subprocess.run(["./run_score.sh", f"/workdir/th.{thread_no}/renamed_{i}.pdb", "/workdir/" + spike], capture_output=True, check=True)  # these paths are inside the container!
        score = float(output.stdout.split()[-1])
        average_score = average_score + score
        if score < best_score:
            best_score = score
            best_score_idx = i
    average_score = average_score / 5.
    best_score_1 = best_score
    average_score_1 = average_score

    # optionally copy the best pdb to save it
    if best_score_idx != 999 and best_score < global_best_score:
        shutil.copy(f"./data/th.{thread_no}/renamed_{best_score_idx}.pdb", "./data/best.pdb")
        global_best_score = best_score

    print(f"Best 0,1: {best_score_0:.4f}, {best_score_1:.4f}, Average 0,1: {average_score_0:.4f} {average_score_1:.4f}, {toc():.2f}")

    # TODO: Find out what is better average or best
    return (average_score_0, average_score_1)
    # return (best_score_0, best_score_1)


# with locally installed Rosetta refinements
def ig_fold_rosetta(thread_no: int, xx):
    HL = np2full_seq(xx)
    sequences = {
        "H": HL[0],
        "L": HL[1],
    }

    igfold = IgFoldRunner()
    igfold.fold(
        f"./data/th.{thread_no}/{ig_fold_pdb}",            # Output PDB file
        sequences=sequences,    # Antibody sequences
        do_refine=True,         # Refine the antibody structure with PyRosetta
        do_renum=True,          # Send predicted structure to AbNum server for Chothia renumbering
    )


# with conteinerised OpenMM refinement
def ig_fold_openmm(thread_no: int, xx):
    HL = np2full_seq(xx)
    output = subprocess.run(["./run_igfold.sh", f"./data/th.{thread_no}/{ig_fold_pdb}", HL[0], HL[1]], capture_output=True, check=True)


# input - list of 2 samples
def double_fun_igfold(X):
    global global_best_score
    tic()

    ##### these 2 calls can be run in parallel (см ниже)
    # thread_no = 0
    # ig_fold(thread_no, X[thread_no])
    # thread_no = 1
    # ig_fold(thread_no, X[thread_no])
    thread_numbers = (0, 1)
    # with Pool(2) as pool:
    #     pool.starmap(ig_fold_rosetta, zip(thread_numbers, X))
    with Pool(2) as pool:
        pool.starmap(ig_fold_openmm, zip(thread_numbers, X))

    ##### the following must run in sequence!
    # run docking/score for AF thread 0
    thread_no = 0
    output = subprocess.run(["./run_score.sh", f"/workdir/th.{thread_no}/{ig_fold_pdb}", "/workdir/" + spike, f"{dla_threshold}"],
                                capture_output=True, check=True)  # these paths are inside the container!
    best_score_0 = float(output.stdout.split()[-1])

    # optionally copy the best pdb to save it
    if best_score_0 < global_best_score:
        shutil.copy(f"./data/th.{thread_no}/{ig_fold_pdb}", "./data/best.pdb")
        global_best_score = best_score_0

    # run docking/score for AF thread 1
    thread_no = 1
    output = subprocess.run(["./run_score.sh", f"/workdir/th.{thread_no}/{ig_fold_pdb}", "/workdir/" + spike, f"{dla_threshold}"],
                                capture_output=True, check=True)  # these paths are inside the container!
    best_score_1 = float(output.stdout.split()[-1])

    # optionally copy the best pdb to save it
    if best_score_1 < global_best_score:
        shutil.copy(f"./data/th.{thread_no}/{ig_fold_pdb}", "./data/best.pdb")
        global_best_score = best_score_1

    print(f"Scores: {best_score_0:.4f} {best_score_1:.4f}, The best: {global_best_score:.4f}, {toc()}")

    return (best_score_0, best_score_1)


if __name__ == '__main__':
    # test_full_seq()
    # exit()
    print(time.asctime())
    # init_pyrosetta()

    if not os.path.exists("data"):
        os.mkdir("data")

    x0 = seq2np(cdr_H1 + cdr_H2 + cdr_H3 + cdr_L1 + cdr_L2 + cdr_L3).flatten()

    # before docking we create 2 files with residues to block on receptor side. These should be unchangeable residues
    HL = np2full_seq(x0)
    save_blocking_positions(HL[0], HL[1])

    # fun = get_fitness
    # get_fitness.n = 0
    # fun = fake_fitness
    # fun = double_fun
    fun = double_fun_igfold
    plot_avg = []
    plot_min = []
    sigma0 = 0.15  # initial standard deviation to sample new solutions - should be ~ 1/4 of range

    # # cfun = cma.ConstrainedFitnessAL(fun, constraints)  # unconstrained function with adaptive Lagrange multipliers
    es = cma.CMAEvolutionStrategy(x0, sigma0,
                        inopts={
                            'ftarget': -3.0,
                            'popsize': 18,
                            'maxiter': 10,
                            'bounds': [-0.1, 1.1],
                            'verb_time': 0,
                            'verb_disp': 500,
                            'seed': 3},)

    while not es.stop():
        X = es.ask()  # sample len(X) candidate solutions

        V = []
        for i in range(len(X)//2):
            x2 = (X[i*2], X[i*2+1])
            v2 = fun(x2)
            V.append(v2[0])
            V.append(v2[1])

        plot_avg.append(np.array(V).mean())
        plot_min.append(np.array(V).min())

        es.tell(X, V)
        # es.tell(X, [fun(x) for x in X])

        es.logger.add()  # for later plotting
        es.disp()

    es.result_pretty()
    # === es.result explanations: https://github.com/CMA-ES/pycma/blob/9b4eb5450c020ac99d637780a42c39788f1e1045/cma/evolution_strategy.py#L977

    # res = np2full_seq(es.result.xbest)
    # the result as the mean of the metamodel gaussian rather than the best candidate:
    res = es.result.xfavorite
    seq = np2full_seq(res)

    sequences = [SeqRecord(Seq(seq[0]), id='H', description="Optimised Ab H"), SeqRecord(Seq(seq[1]), id="L", description="Optimised Ab L")]
    SeqIO.write(sequences, "best.fasta", "fasta")

    print(seq)
    print(f"The best score: {es.result.fbest}, iterations: {es.result.iterations}")
    print(time.asctime())

    plot = np.array([plot_avg, plot_min]).T
    matplotlib.pyplot.plot(plot)
    matplotlib.pyplot.show(block=True)

    es.plot()
    matplotlib.pyplot.show(block=True)

