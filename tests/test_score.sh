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
obabel 2dd8_Fv.pdb -xr -O /tmp/receptor.pdbqt -p 7.4
obabel 2dd8_spike.pdb -O /tmp/ligand.pdbqt -p 7.4

vina --receptor /tmp/receptor.pdbqt --ligand /tmp/ligand.pdbqt  --score_only