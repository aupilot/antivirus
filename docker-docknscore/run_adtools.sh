#!/bin/bash

source ~/Apps/ADFRsuite-1.0/bin/adsenv.sh

cd /tmp
/home/kir/Apps/ADFRsuite-1.0/CCSBpckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l /tmp/ligand.pdb -o /tmp/ligand.pdbqt
/home/kir/Apps/ADFRsuite-1.0/CCSBpckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r /tmp/receptor.pdb -o /tmp/receptor.pdbqt

vina --receptor /tmp/receptor.pdbqt --ligand /tmp/ligand.pdbqt  --score_only
~/Apps/smina-code/build/smina --receptor /tmp/receptor.pdbqt --ligand /tmp/ligand.pdbqt  --score_only

