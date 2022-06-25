#!/bin/bash

docker run -t --rm --gpus all -v $(pwd):/workdir --workdir /workdir 8kir8/score:latest /usr/bin/python3 /opt/dock-n-score.py $1 $2 $3 $4
