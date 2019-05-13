#!/usr/bin/env bash

#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --mem=15G
#SBATCH -o log/MainSnakemake.out

WORK_DIR="$(pwd)"

for workflow in workflows/*/; do
    cd $workflow
    snakemake --unlock
    cd $WORK_DIR
done

snakemake -j 999 \
    --cluster-config /home/kkchau/projects/RNAseq_Pipeline/cluster.yaml \
    --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem}"
