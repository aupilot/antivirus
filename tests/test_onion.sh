#!/bin/bash

docker run -t --rm --gpus all -v $(pwd):/workdir --workdir /workdir 8kir8/score:latest /usr/bin/python3 /opt/onionnet/kir01/onion_score.sh $1