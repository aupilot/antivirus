#!/bin/bash

# usage:
# run_alpha fasta_file thread_no

source /home/kir/Apps/alphafold/venv/bin/activate
#cd /home/kir/Apps/alphafold

mkdir -p /tmp/alphafold/th.$2

#python3 /home/kir/Apps/alphafold/docker-anarci/run_docker.py --fasta_paths=$1 --max_template_date=2020-05-14 --model_preset=multimer --db_preset=reduced_dbs --data_dir=/media/kir/HDDData/Datasets/Alphafold

# in production - set num_multimer_predictions_per_model=3
python3 /home/kir/Apps/alphafold/docker/run_docker.py --fasta_paths=$1 --num_multimer_predictions_per_model=1 --max_template_date=2022-05-01 --model_preset=multimer --db_preset=reduced_dbs --output_dir=/tmp/alphafold/th.$2 --data_dir=/media/kir/HDDData/Datasets/Alphafold
#--use_precomputed_msas=true

python3 ./rename.py /tmp/alphafold/th.$2/fitness.$2/ranked_0.pdb /tmp/alphafold/th.$2/fitness.$2/renamed_0.pdb
python3 ./rename.py /tmp/alphafold/th.$2/fitness.$2/ranked_1.pdb /tmp/alphafold/th.$2/fitness.$2/renamed_1.pdb
python3 ./rename.py /tmp/alphafold/th.$2/fitness.$2/ranked_2.pdb /tmp/alphafold/th.$2/fitness.$2/renamed_2.pdb
python3 ./rename.py /tmp/alphafold/th.$2/fitness.$2/ranked_3.pdb /tmp/alphafold/th.$2/fitness.$2/renamed_3.pdb
python3 ./rename.py /tmp/alphafold/th.$2/fitness.$2/ranked_4.pdb /tmp/alphafold/th.$2/fitness.$2/renamed_4.pdb
