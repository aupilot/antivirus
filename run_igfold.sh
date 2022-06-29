#!/bin/bash

# docker login -u 8kir8
docker run -t --rm --gpus all -v $(pwd):/workdir --workdir /workdir 8kir8/igfold2:latest /usr/bin/python3 /opt/IgFold/ig_fold.py $1 $2 $3 $4
# Dmytro's
#docker run -t --rm --gpus '"device=0"' -v $(pwd):/workdir --workdir /workdir 8kir8/igfold2:latest /usr/bin/python3 /opt/IgFold/ig_fold.py $1 $2 $3 $4
#docker run -t --rm --gpus '"device=1"' -v $(pwd):/workdir --workdir /workdir 8kir8/igfold2:latest /usr/bin/python3 /opt/IgFold/ig_fold.py $1 $2 $3 $4
