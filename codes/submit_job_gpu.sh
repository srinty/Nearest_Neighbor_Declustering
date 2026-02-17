#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00
#SBATCH --gres=gpu:4
#SBATCH --partition=agpuq

source ~/.bashrc
conda activate pygmt

python run_NND_all.py