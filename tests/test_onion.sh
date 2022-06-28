#!/bin/bash

docker run -t --rm --gpus all -v $(pwd):/workdir --workdir /workdir 8kir8/score:latest /bin/bash /opt/onionnet/kir01/onion_score.sh /workdir/e_dec_0.94.pdb
