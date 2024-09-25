#!/bin/bash
#SBATCH -t 7:59:00
#SBATCH --mem=60GB
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -e %j.err
#SBATCH -o %j.out
#SBATCH --priority=TOP

#burst
#SBATCH --account=chunlidai
#SBATCH --qos=chunlidai-b

module unload matlab
module load matlab/2023b  # projcrs, only works for matlab Since R2020b
exec matlab < Tilemain_nov.m


