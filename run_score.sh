#!/bin/bash

# run_score.sh /workdir/alphafold/fitness/ranked_0.pdb /workdir/6yla_SPIKE.pdb 0.06 1

# test from "data" folder:
# docker run -it --gpus all --rm -v $(pwd):/workdir --workdir /workdir 8kir8/molecule:0501.1 /usr/bin/python3 /opt/dock-n-score.py /workdir/best_01.pdb /workdir/7e3c_spike_alighned.pdb 0.06 1 1

cd data/

# docker login -u 8kir8
# pwd: Supplier1
docker run -t --rm --gpus all -v $(pwd):/workdir --workdir /workdir 8kir8/score:latest /usr/bin/python3 /opt/dock-n-score.py $1 $2 $3 $4 $5
# Dmytro's
#docker run -t --rm --gpus '"device=0"' -v $(pwd):/workdir --workdir /workdir 8kir8/score:latest /usr/bin/python3 /opt/dock-n-score.py $1 $2 $3 $4 $5
#docker run -t --rm --gpus '"device=1"' -v $(pwd):/workdir --workdir /workdir 8kir8/score:latest /usr/bin/python3 /opt/dock-n-score.py $1 $2 $3 $4 $5

# inside the container
#/usr/bin/python3 /opt/dock-n-score.py /workdir/best_07.pdb /workdir/7cr5_SPIKE.pdb 0.06 1