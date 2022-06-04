import subprocess
import time
from multiprocessing import Pool
from prody import parsePDB, writePDB
from ttictoc import tic, toc
from Bio.Seq import Seq
import Bio.PDB
import difflib

ig_fold_pdb = "ig_fold.pdb"


# with conteinerised Rosetta or OpenMM refinement running in docker-anarci container
def ig_fold_sequence(thread_no: int, sequence, use_rosetta=0, do_renum=0):
    output = subprocess.run(
        ["./run_igfold.sh", f"./data/th.{thread_no}/{ig_fold_pdb}", sequence["H"], sequence["L"], f"--rosetta={use_rosetta}",
         f"--renum={do_renum}"],
        capture_output=True, check=True)


# input: list of 2 sequences, each seq is a dict with "H" and "L"
def fold_n_score2(sequences, spike, mega_type=0, dla_threshold=0.06, rosetta=0, renum=0):
    tic()

    ##### these 2 calls can be run in parallel
    thread_numbers = (0, 1)
    rosetta_n = [rosetta]*2
    renum_n = [renum]*2
    with Pool(2) as pool:
        pool.starmap(ig_fold_sequence, zip(thread_numbers, sequences, rosetta_n, renum_n))
    # thread_no = 0
    # ig_fold_sequence(thread_no, sequences[thread_no])
    # thread_no = 1
    # ig_fold_sequence(thread_no, sequences[thread_no])

    ##### the following must run sequentially!
    # run docking/score
    scores = []
    for thread_no in range(2):
        output = subprocess.run(
            ["./run_score.sh", f"/workdir/th.{thread_no}/{ig_fold_pdb}", "/workdir/" + spike, f"{dla_threshold}",
             f"{mega_type}"],
            capture_output=True, check=True)  # these paths are inside the container!
        scores.append(float(output.stdout.split()[-1]))

    # run docking/score for thread 1
    # thread_no = 1
    # output = subprocess.run(
    #     ["./run_score.sh", f"/workdir/th.{thread_no}/{ig_fold_pdb}", "/workdir/" + spike, f"{dla_threshold}",
    #      f"{mega_type}"],
    #     capture_output=True, check=True)  # these paths are inside the container!
    # best_score_1 = float(output.stdout.split()[-1])

    print(f"Scores:"+ str(scores) + f" Time: {toc():.1f}")

    return scores


# def read_fasta(fasta_file):
#     buf = []
#     with open(fasta_file, "r") as infile:
#         for line_idx, line in enumerate(infile):
#             if line.startswith(">"):  # label line
#                 line = line[1:].strip()
#                 if len(line) > 0:
#                     cur_seq_label = line
#                 else:
#                     cur_seq_label = f"seqnum{line_idx:09d}"
#             else:  # sequence line
#                 buf.append(line.strip())
#     return buf


# highligh differences in color

# red = lambda text: f"\033[38;2;255;0;0m{text}\033[38;2;255;255;255m"
# green = lambda text: f"\033[38;2;0;255;0m{text}\033[38;2;255;255;255m"
# blue = lambda text: f"\033[38;2;0;0;255m{text}\033[38;2;255;255;255m"
# white = lambda text: f"\033[38;2;255;255;255m{text}\033[38;2;255;255;255m"

red = lambda text: f"\033[38;2;255;0;0m{text}\033[m"
green = lambda text: f"\033[38;2;0;255;0m{text}\033[m"
blue = lambda text: f"\033[38;2;0;0;255m{text}\033[m"
white = lambda text: f"\033[38;2;255;255;255m{text}\033[m"
native = lambda text: f"\033[m{text}"

def highlight_differences(old, new):
    result = ""
    codes = difflib.SequenceMatcher(a=old, b=new).get_opcodes()
    for code in codes:
        if code[0] == "equal":
            result += native(old[code[1]:code[2]])
        elif code[0] == "delete":
            result += red(old[code[1]:code[2]])
        elif code[0] == "insert":
            result += green(new[code[3]:code[4]])
        elif code[0] == "replace":
            result += blue(new[code[3]:code[4]])

    return result



def insert_spacers(sequences):
    """
    use ANARCI MSA to insert spacers
    https://github.com/oxpig/ANARCI
    http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/

    Aha:
DIQLTQSPSFLSASVGDRVTITCRAS--QGIS------NYLAWYQQKPGKAPKLLIYA--------ASTLQTGVPSRFSGSGSG--TEFTLTISSLQPEDFATYYCQQLNSY----------------------PPLTFGGGSRVEIK
EVLLVES-GGGLVQPGGSLRLSCAASG-LTVSS-----NYMTWVRQAPGKGLEWVSVIYS----GGSTFYADSVKGRCTISRHNSKNTLYLQMNSLRAEDTAVYYCARDLDY---------------------YGMDVWGQGTTVTVSS
    Chothia:
DIQLTQSPSFLSASVGDRVTITCRASQGISNYLAWYQQKPGKAPKLLIYAASTLQTGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQQLNSYPPLTFGGGSRVEIK
EVLLVESGGGLVQPGGSLRLSCAASGLTVSSNYMTWVRQAPGKGLEWVSVIYSGGSTFYADSVKGRCTISRHNSKNTLYLQMNSLRAEDTAVYYCARDLDYYGMDVWGQGTTVTVSS
    IGMT:
DIQLTQSPSFLSASVGDRVTITCRASQGI------SNYLAWYQQKPGKAPKLLIYAA-------STLQTGVP-SRFSGSG--SGTEFTLTISSLQPEDFATYYCQQLNS---YPPLTFGGGSRVEIK
EVLLVESGG-GLVQPGGSLRLSCAASGLTV----SSNYMTWVRQAPGKGLEWVSVIYSG---GSTFYADSVK-GRCTISRHNSKNTLYLQMNSLRAEDTAVYYCARDLDY--YGMDVWGQGTTVTVSS
    Kabat:
DIQLTQSPSFLSASVGDRVTITCRASQGISNYLAWYQQKPGKAPKLLIYAASTLQTGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQQLNSYPPLTFGGGSRVEIK
EVLLVESGGGLVQPGGSLRLSCAASGLTVSSNYMTWVRQAPGKGLEWVSVIYSGGSTFYADSVKGRCTISRHNSKNTLYLQMNSLRAEDTAVYYCARDLDYYGMDVWGQGTTVTVSS
    Martin:
DIQLTQSPSFLSASVGDRVTITCRASQGISNYLAWYQQKPGKAPKLLIYAASTLQTGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQQLNSYPPLTFGGGSRVEIK
EVLLVESGGGLVQPGGSLRLSCAASGLTVSSNYMTWVRQAPGKGLEWVSVIYSGGSTFYADSVKGRCTISRHNSKNTLYLQMNSLRAEDTAVYYCARDLDYYGMDVWGQGTTVTVSS
    Wolfguy:
DIQLTQSPSFLSASVGDRVTITCRASQGISNYLAWYQQKPGKAPKLLIYAASTLQTGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQQLNSYPPLTFGGGSRVEIK
EVLLVESGGGLVQPGGSLRLSCAASGLTVSSNYMTWVRQAPGKGLEWVSVIYSGGSTFYADSVKGRCTISRHNSKNTLYLQMNSLRAEDTAVYYCARDLDYYGMDVWGQGTTVTVSS

    # check
print(highlight_differences(np2seq_show(x0)[0], str(start_seq['H'])))
    """

    # schema = 'a'  # Aho
    # schema = 'i'  # IMGT
    schema = 'c' # Chothia
    # schema = 'k' # kabat
    # schema = 'm' # martin
    # schema = 'w' # wolfguy

    new_sequences = {}
    for key, seq in sequences.items():
        str_seq = str(seq)
        output = subprocess.run(["docker", "run", "-t", "--rm", r"8kir8/anarci:220524.1", "ANARCI", "-s", f"{schema}", "-i", f"{str_seq}"], capture_output=True, check=True)
        # "docker run -it --rm  8kir8/anarci:220524.1 ANARCI -s a -i EVQLVQSGAEVKKPGESLKISCQGSGYSFTSYWIGWVRQMPGKGLEWMGIIYPGESDTRYSSSFQGHVTISADKSISTAYLQWSSLKASDTAMYYCARIRGVYSSGWIGGDYWGQGTLVTVSS"
        new_str_seq = ""
        out_text = output.stdout.decode().split("\n")
        for line in out_text:
            if line[0] == '#':
                continue
            if line[0] == 'H' or line[0] == 'L':
                new_str_seq += line[10]
            if line[0] == '/':
                break
        new_sequences[key] = Seq(new_str_seq)
        print(new_str_seq)
    return new_sequences


# read sample pdb file, align agaist the reference and save to output
def align(reference, sample, output):
    # Select what residues we wish to align. Ideally, we need to renumber?
    start_id = 1
    end_id = 100
    atoms_to_be_aligned = range(start_id, end_id + 1)

    pdb_parser = Bio.PDB.PDBParser(QUIET=False)

    ref_structure = pdb_parser.get_structure("reference", reference)
    sample_structure = pdb_parser.get_structure("sample", sample)

    # Use the first model in the pdb-files for alignment (we only have one)
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]

    # Make a list of the CA atoms we wish to align
    ref_atoms = []
    sample_atoms = []

    # Iterate of all chains in the model in order to find all residues
    for ref_chain in ref_model:
        # Iterate of all residues in each model in order to find proper atoms
        for ref_res in ref_chain:
            # Check if residue number is in the list
            if ref_res.get_id()[1] in atoms_to_be_aligned:
                # Append CA atom to list
                ref_atoms.append(ref_res['CA'])

    # Do the same for the sample structure
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
                sample_atoms.append(sample_res['CA'])

    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    # print(super_imposer.rms)
    if super_imposer.rms > 20.0:     # could be as much as 8.6 for some Abs that don't fit CoV2
        raise Exception(f"Something is wrong with PDB alignment - RMS is too high {super_imposer.rms}")

    # Save the aligned version of a pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure)
    io.save(output)

