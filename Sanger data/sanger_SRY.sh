#!/bin/bash

#SBATCH --job-name=SRY_phylogeny
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

set -eo pipefail

eval "$(/nfs/home/yunovikova/.local/bin/micromamba shell hook --shell=bash)"
micromamba activate project

set -u

BASE="/mnt/tank/scratch/yunovikova/mitogenomes/sanger_2026"
SRY_DIR="${BASE}/SRY"
ORIGINAL_FASTA="${SRY_DIR}/original/SRY.fasta"
RESULTS_DIR="${SRY_DIR}/complete_pipeline"
MAFFT_DIR="${RESULTS_DIR}/01_mafft"
GBLOCKS_DIR="${RESULTS_DIR}/02_gblocks"
IQTREE_DIR="${RESULTS_DIR}/03_iqtree"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$MAFFT_DIR" "$GBLOCKS_DIR" "$IQTREE_DIR"

# 1. Alignment with MAFFT L-INS-i
mafft \
    --localpair \
    --maxiterate 1000 \
    --thread "$THREADS" \
    "$ORIGINAL_FASTA" \
    > "${MAFFT_DIR}/SRY_aligned.fasta"

# 2. Removal of ambiguously aligned regions with Gblocks (optional)
GBLOCKS_INPUT="${GBLOCKS_DIR}/SRY_aligned.fasta"
cp "${MAFFT_DIR}/SRY_aligned.fasta" "$GBLOCKS_INPUT"

N_SEQUENCES="$(grep -c '^>' "$GBLOCKS_INPUT")"
MIN_SEQUENCES=$((N_SEQUENCES / 2 + 1))

Gblocks "$GBLOCKS_INPUT" \
    -t=d \
    -b1="$MIN_SEQUENCES" \
    -b2="$MIN_SEQUENCES" \
    -b3=10 \
    -b4=5 \
    -b5=h

mv "${GBLOCKS_INPUT}-gb.fa" \
    "${GBLOCKS_DIR}/SRY_gblocks.fasta"

# 3. Maximum-likelihood trees with model selection and UFBoot2
iqtree3 \
    -s "${MAFFT_DIR}/SRY_aligned.fasta" \
    -m MFP \
    -B 1000 \
    -nm 5000 \
    -T AUTO \
    -ntmax "$THREADS" \
    -seed 20260802 \
    --prefix "${IQTREE_DIR}/SRY_mafft"

iqtree3 \
    -s "${GBLOCKS_DIR}/SRY_gblocks.fasta" \
    -m MFP \
    -B 1000 \
    -nm 5000 \
    -T AUTO \
    -ntmax "$THREADS" \
    -seed 20260802 \
    --prefix "${IQTREE_DIR}/SRY_gblocks"
