#!/bin/bash

# -0.91540
#obabel 7e3c_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 7e3c_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

# -0.46298
#obabel 7ps7_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 7ps7_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

# -0.96857
#obabel 7urs_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 7urs_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

# -0.61630
#obabel 6zer_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 6zer_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

#-2.54717
#obabel 2dd8_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 2dd8_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

#-0.93759 (kcal/mol)
#obabel 6yla_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 6yla_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

# -0.91550 (kcal/mol)
#obabel 4xak_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 4xak_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

# -0.69870
#obabel 6c6z_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 6c6z_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

# broken ligand
#obabel 1sy6_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 1sy6_spike.pdb -O /tmp/ligand.pdbqt -l 0 -p 7.4

# -0.98265 (kcal/mol) (after removing the useless γ chain)
#obabel 1sy6_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 1sy6_epitope.pdb -O /tmp/ligand.pdbqt -l 0 -p 7.4


# 2.69287 (kcal/mol). The other epitope of the same epsilon does not match well
#obabel 1sy6_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
#obabel 6jxr_epitope_aligned_to_1sy6.pdb -O /tmp/ligand.pdbqt -l 0 -p 7.4

python3 split_rec_lig.py i_dec_1.03.pdb ./
obabel rec.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
obabel lig.pdb -O /tmp/ligand.pdbqt -l 0 -p 7.4

vina --receptor /tmp/receptor.pdbqt --ligand /tmp/ligand.pdbqt  --score_only

