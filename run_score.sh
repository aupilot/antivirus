#!/bin/bash

# run_score.sh /workdir/alphafold/fitness/ranked_0.pdb /workdir/6yla_SPIKE.pdb

# test from "data" folder:
# docker run -it --gpus all --rm -v $(pwd):/workdir --workdir /workdir 8kir8/molecule:0426.2 /usr/bin/python3 /opt/dock-n-score.py /workdir/best.pdb /workdir/7cr5_SPIKE.pdb 0.06
# docker run -it --gpus all --rm -v $(pwd):/workdir --workdir /workdir 8kir8/molecule:0501.1 /usr/bin/python3 /opt/dock-n-score.py /workdir/best.pdb /workdir/7cr5_SPIKE.pdb 0.06 1

cd data/

# docker login -u 8kir8
# pwd: Supplier1
#source /home/kir/Apps/alphafold/venv/bin/activate
docker run -t --rm --gpus all -v $(pwd):/workdir --workdir /workdir 8kir8/molecule:0501.1 /usr/bin/python3 /opt/dock-n-score.py $1 $2 $3 $4
#docker ps
