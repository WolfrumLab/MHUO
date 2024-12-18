#!/bin/bash

#SBATCH --job-name="sc_snp_scaling"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

set +e

# Replace SNAME with sample name, and SNUM with the number of individuals in the pool
SNAME=""
SNUM=""

snakemake -p \
	--cores 20 \
	--use-conda \
	-s ~/path/to/sc_snp.smk \
	"${SNAME}"_"${SNUM}N_genome1K.phase3.SNP_AF5e2.chr1toX.hg38.demux" \
	--configfile config.yml
