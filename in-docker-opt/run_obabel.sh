#!/bin/bash


# https://www.cheminformania.com/ligand-docking-with-smina/
# https://open-babel.readthedocs.io/en/latest/FileFormats/AutoDock_PDBQT_format.html
# https://buildmedia.readthedocs.org/media/pdf/open-babel/latest/open-babel.pdf

# original openbabel DOES NOT WORK! Fix with BRANCH atom number spacing 5


cd /tmp

#obabel /tmp/receptor.pdb -xr -O /tmp/receptor.pdbqt 
#obabel /tmp/ligand.pdb -O /tmp/ligand.pdbqt
#obabel /tmp/receptor.pdb -xr -O /tmp/receptor.pdbqt --addpolarh 
#obabel /tmp/ligand.pdb -O /tmp/ligand.pdbqt --addpolarh
obabel /opt/var/receptor.pdb -xr -O /tmp/receptor.pdbqt -p 7.4 
obabel /opt/var/ligand.pdb -O /tmp/ligand.pdbqt -p 7.4

/usr/local/bin/vina --receptor /tmp/receptor.pdbqt --ligand /tmp/ligand.pdbqt  --score_only
#~/Apps/smina-code/build/smina --receptor /tmp/receptor.pdbqt --ligand /tmp/ligand.pdbqt  --score_only
