#!/bin/bash
#
#SBATCH --job-name=wasp
#SBATCH --output=out.cwas.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=medium
#SBATCH --mem=8G
#SBATCH -t 2-0

srun snakemake -s cwas.snakefile -j 100 -pr --latency-wait 60 --use-conda --rerun-incomplete --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {threads}" 

