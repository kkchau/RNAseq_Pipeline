#!/usr/bin/env bash

#SBATCH --mem=40G
#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH -J testSTARindex

conda activate rnaseq
snakemake
