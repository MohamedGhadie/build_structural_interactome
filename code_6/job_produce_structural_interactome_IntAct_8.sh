#!/bin/bash
#SBATCH --account=ctb-yxia
#SBATCH --time=7-00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/project/ctb-yxia/ghadie84/edgotype_fitness_effect/code_6/%x-%j.out
#SBATCH --job-name=ghadie84_produce_structural_interactome_IntAct_8

source ~/venv/bin/activate
python produce_structural_interactome_IntAct_8.py