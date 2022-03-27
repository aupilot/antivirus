#!/bin/bash

# run_score.sh /workdir/alphafold/fitness/ranked_0.pdb /workdir/6yla_SPIKE.pdb

cd data

#source /home/kir/Apps/alphafold/venv/bin/activate
docker run -t --gpus all -v $(pwd):/workdir --workdir /workdir molecule:0322.1 /usr/bin/python3 /opt/dock-n-score.py $1 $2
#docker ps

