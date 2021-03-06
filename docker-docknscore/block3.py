#!/usr/bin/python3
#
# Copyright (C) 2020 Tokyo Institute of Technology
#
# =====================================================================
#
#   Software Name : MEGADOCK (blocking residues)
#
#   Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
#
#   Last update: December 3, 2014
#
# =====================================================================
import os
import sys
import csv

if (len(sys.argv) != 4):
    print('Usage: $ %s [pdbfile] [chain] [target residue list]' % sys.argv[0])
    print()
    print(' e.g.) $ %s 1gcq_r.pdb B 182-186,189,195-198,204 > blocked.pdb' % sys.argv[0])
    print()
    print('Note: Target residue list is separated by commas and no spaces.')
    print('      You can also use hyphen \"-\": \"182-186\" means blocking residues of 182, 183, ..., 186.')
    print('      Blocked residues are substituted for \'BLK\'.')
    print('      Updated PDB coordinates are written to \"standard output\".')
    quit(1)

ch = sys.argv[2]

reslist = []
resarg = sys.argv[3].split(",")
for r in resarg:
    if "-" in r:
        rseq = r.split("-")
        for i in range(int(rseq[0]), int(rseq[1])+1):
            reslist.append(str(i))
    else:
        reslist.append(r)

fp = open(sys.argv[1], "r")
try:
    for l in fp.readlines():
        l = l.strip()

        if l[0:4] != "ATOM" and l[0:6] != "HETATM":
            print(l)
            continue

        if l[21] != ch:
            print(l)
            continue

        ll = l
        if ll[22:26].strip() not in reslist:
            print(l)
            continue

        print(l[0:16] + " BLK" + l[20:])

finally:
    fp.close()
