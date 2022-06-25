000
7e3c-Fv.fasta
esm1v_t33_650M_UR90S_5
'popsize': 12,
'maxiter': 31,
--sigma=0.5 --dla=0.045 --mega=1
-0.68662

000
7e3c-Fv.fasta
esm1v_t33_650M_UR90S_5
'popsize': 18,
'maxiter': 21,
-0.70718

001
7e3c-Fv.fasta
esm1v_t33_650M_UR90S_5
popsize': 18,
'maxiter': 41,
-0.69678

002
7e3c-Fv.fasta
esm1v_t33_650M_UR90S_5
popsize': 12,
'maxiter': 61,
-0.68264,

003
7e3c-Fv.fasta
esm1v_t33_650M_UR90S_5
popsize': 24,
'maxiter': 15,
-0.57133

004
7e3c-Fv.fasta
esm1v_t33_650M_UR90S_5
'popsize': 18,
'maxiter': 21,
Testing different renumbering (with spacers) - IGMT
-0.6507

005
7e3c-Fv.fasta
esm1v_t33_650M_UR90S_5
'popsize': 18,
'maxiter': 21,
sigma 1.0
-0.59776

006
--sigma=0.5 --dla=0.06
-0.64767


=====================================
starting from 2dd8 -- 41 iteration
2dd8_Fv.fasta
esm1v_t33_650M_UR90S_5
'popsize': 18,
'maxiter': 41,
sigma 0.5
-0.6932

======================================
1. Testng different encodings:

esm1v_t33_650M_UR90S_5
7urs-Fv.fasta
7urs_spike.pdb
alignment.pdb (from 7e3c?)
popsize: 18
maxiter: 21
--sigma=0.5 --dla=0.045 --mega=1
-0.46956/-0.473
---------------------------
Хреновый
esm1b_t33_650M_UR50S
-0.40559
---------------------------
esm1v_t33_650M_UR90S_4
-0.48246 -0.43762
---------------------------
esm1v_t33_650M_UR90S_1
crashed because it's truncated generated sequence WTF?
---------------------------
esm1v_t33_650M_UR90S_2
crashed because it's truncated generated sequence WTF?

esm1v_t33_650M_UR90S_3
crashed because it's truncated generated sequence WTF?

-----------------------------
81 iteration - does not help

---------------------------
dla 0.03  sigma 0.5     -0.48465
dla 0.045 sigma 0.7     -0.46893
dla 0.045 sigma 0.5     -0.4519
dla 0.045 sigma 0.2     -0.43734
dla 0.045 sigma 0.5 IMGT -0.46189




=====================================





текущие:

mine: new embedding test

Dmy kir1: compare with old embedding





================================
идеи
re-implementation of trRossetta for folding.
another embedding - berta
