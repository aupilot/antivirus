#!/bin/bash

# ./mega-post.sh fab.pdb epitope.pdb

if [ "$#" -ne 2 ]; then
    echo -e "\033[1mUsage:\033[0m mega-post.sh fab.pdb epitope.pdb"
else

mkdir -p /opt/var/decoys

for i in `seq 1 6`;
	do /opt/megadock-kir/decoygen /opt/var/decoys/lig.${i}.pdb $2 /opt/var/dock.out $i;
		cat $1 /opt/var/decoys/lig.${i}.pdb > /opt/var/decoys/decoy.${i}.pdb;
	done
fi
