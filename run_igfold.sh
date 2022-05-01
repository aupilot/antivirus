#!/bin/bash

# docker login -u 8kir8
docker run -t --rm --gpus all -v $(pwd):/workdir --workdir /workdir 8kir8/igfold:0502.1 /usr/bin/python3 /opt/IgFold/ig_fold.py $1 $2 $3
