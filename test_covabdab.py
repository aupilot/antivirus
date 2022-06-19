import random
import time
import pandas as pd
from helper import fold_n_score2


"""
get bind/no_bind from covabdab and compare to our pipeline results (fold/bind/score)
"""


CovAbDab_db_file = "/media/kir/Work/Datasets/CoV-AbDab/CoV-AbDab_200422.csv"


binds_to = ["SARS-CoV2_WT"]
binds_to_exclude = "weak"

# epitopes = ["S; RBD", "S; S2"]
epitopes = ["S; RBD"]

samples = None                      # if None use the min of pos/neg for both pos and neg. Otherwise use this number for both

# files are inside the "data" dir
# spike = "7k9i_spike.pdb"            # wild type SARS-CoV2 spike protein S1. Don't use Amber numbering!
spike = "7e3c_spike_aligned.pdb"            # wild type SARS-CoV2 spike protein S1. Don't use Amber numbering!
mega_type = 1                       # 0 - orig, 1 - kir
dla_threshold = 0.07                # drop the dockings with lower score as non-natural (from DLA-Ranker perspective)
rosetta = 0                         # don't use rosetta. OpenMM is way better
renum = 0

result_file = "results.csv"

random.seed(666)

if __name__ == '__main__':
    print(time.asctime())
    df = pd.read_csv(CovAbDab_db_file)
    # exclude nanobodies
    df = df.loc[df["Ab or Nb"] == "Ab"]
    # remove lines without sequence or with full Fab (long) sequences
    df = df.loc[(df["VH or VHH"].str.len() < 150) & (df["VH or VHH"].str.len() > 100)]
    df = df.loc[(df["VL"].str.len() < 150) & (df["VL"].str.len() > 100)]

    # == positive samples: bings to what we needed
    # we only need entries that binds to specific virus and not weekly
    df_pos = df.loc[df["Binds to"].isin(binds_to) + df["Binds to"].str.contains("weak") == False]
    # choose specific epitopes we're interested in
    df_pos = df_pos.loc[df_pos["Protein + Epitope"].isin(epitopes)]
    orig_pos_n = len(df_pos)

    # == negative samples
    df_neg = df.loc[df["Doesn't Bind to"].isin(binds_to)]

    result = []

    # print(df.to_string())

    #pick same number pos as neg randomly
    df_pos = df_pos.sample(n=len(df_neg))
    print(f"OrigPos: {orig_pos_n}, NewPos: {len(df_pos)}, Neg: {len(df_neg)}")

    #### now start folding-docking-scoring
    if samples is None:
        samples = len(df_neg)
    else:
        df_pos = df_pos.sample(n=samples)
        df_neg = df_neg.sample(n=samples)

    df_pos = df_pos.reset_index()
    df_neg = df_neg.reset_index()

    result_df = pd.DataFrame({
        "Name": pd.Series(dtype='str'),
        # "Epitope": pd.Series(dtype='str'),
        "Binds": pd.Series(dtype='bool'),
        "Score": pd.Series(dtype='float'),
        "H": pd.Series(dtype='str'),
        "L": pd.Series(dtype='str'),

    })

    score_pos = 0
    score_neg = 0

    for (idx, s_pos), (_, s_neg) in zip(df_pos.iterrows(), df_neg.iterrows()):
        pos_seq = {}
        pos_seq['H'] = s_pos['VH or VHH']
        pos_seq['L'] = s_pos['VL']
        neg_seq = {}
        neg_seq['H'] = s_neg['VH or VHH']
        neg_seq['L'] = s_neg['VL']
        sequence = (pos_seq, neg_seq)
        scores = fold_n_score2(sequence, spike, mega_type, dla_threshold=dla_threshold, rosetta=rosetta, renum=renum)
        # scores = [0.2, 0.3]

        # res_p = pd.Series({'Name': s_pos['Name'], "Binds": True,  'Score': scores[0], 'H': s_pos['VH or VHH'], 'L': s_pos['VL']}) #fking pandas depricated append!
        # res_n = pd.Series({'Name': s_neg['Name'], "Binds": False, 'Score': scores[1], 'H': s_neg['VH or VHH'], 'L': s_neg['VL']})
        res_p = pd.DataFrame(
            {'Name': s_pos['Name'], "Binds": True, 'Score': scores[0], 'H': s_pos['VH or VHH'], 'L': s_pos['VL']},
            index=[0])
        res_n = pd.DataFrame(
            {'Name': s_neg['Name'], "Binds": False, 'Score': scores[1], 'H': s_neg['VH or VHH'], 'L': s_neg['VL']},
            index=[0])

        # result_df = result_df.append(res_p, ignore_index=True)
        # result_df = result_df.append(res_n, ignore_index=True)
        result_df = pd.concat([result_df, res_p, res_n])

    result_df.to_csv(result_file)

    sp = result_df[result_df['Binds'] == True]['Score'].to_numpy()
    sn = result_df[result_df['Binds'] == False]['Score'].to_numpy()

    print(f"Positive score: {sp[sp < 10.].mean()}   Negative: {sn[sn < 10.].mean()} ")

    # https://dzone.com/articles/correlation-between-categorical-and-continuous-var-1
    df = pd.read_csv("results.csv")
    df_filtered = df[df['Score'] < 10]
    print(df_filtered[['Binds','Score']].corr())

    print(time.asctime())
