#!/bin/bash

source /home/kir/Apps/alphafold/venv/bin/activate
#cd /home/kir/Apps/alphafold

#python3 /home/kir/Apps/alphafold/docker/run_docker.py --fasta_paths=$1 --max_template_date=2020-05-14 --model_preset=multimer --db_preset=reduced_dbs --data_dir=/media/kir/HDDData/Datasets/Alphafold
python3 /home/kir/Apps/alphafold/docker/run_docker.py --fasta_paths=$1 --num_multimer_predictions_per_model=1 --max_template_date=2022-03-01 --use_precomputed_msas=true --model_preset=multimer --db_preset=reduced_dbs --data_dir=/media/kir/HDDData/Datasets/Alphafold
#--use_precomputed_msas=true
#--num_multimer_predictions_per_model=1
python3 ./rename.py /tmp/alphafold/fitness/ranked_0.pdb /tmp/alphafold/fitness/renamed_0.pdb
python3 ./rename.py /tmp/alphafold/fitness/ranked_1.pdb /tmp/alphafold/fitness/renamed_1.pdb
python3 ./rename.py /tmp/alphafold/fitness/ranked_2.pdb /tmp/alphafold/fitness/renamed_2.pdb
python3 ./rename.py /tmp/alphafold/fitness/ranked_3.pdb /tmp/alphafold/fitness/renamed_3.pdb
python3 ./rename.py /tmp/alphafold/fitness/ranked_4.pdb /tmp/alphafold/fitness/renamed_4.pdb
