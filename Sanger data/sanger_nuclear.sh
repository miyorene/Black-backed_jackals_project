#!/bin/bash

#SBATCH --job-name=nuclear_phylogeny
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

set -eo pipefail

eval "$(/nfs/home/yunovikova/.local/bin/micromamba shell hook --shell=bash)"
micromamba activate project

set -u

BASE="/mnt/tank/scratch/yunovikova/sanger"
NUCLEAR_DIR="${BASE}/nuclear"
ORIGINAL_DIR="${NUCLEAR_DIR}/original"
RESULTS_DIR="${NUCLEAR_DIR}/complete_pipeline"
MAFFT_DIR="${RESULTS_DIR}/01_mafft"
GBLOCKS_DIR="${RESULTS_DIR}/02_gblocks"
CONCAT_DIR="${RESULTS_DIR}/03_concatenated"
IQTREE_DIR="${RESULTS_DIR}/04_iqtree"
CONCATENATE_SCRIPT="${BASE}/concatenate_alignments.py"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

LOCI=(ACT FBN1 GHR MDH1 PRKDC TJP1 TRSP TTR1)

mkdir -p "$MAFFT_DIR" "$GBLOCKS_DIR" "$CONCAT_DIR" "$IQTREE_DIR"

# 1. Alignment of each locus with MAFFT L-INS-i
for LOCUS in "${LOCI[@]}"; do
    mafft \
        --localpair \
        --maxiterate 1000 \
        --thread "$THREADS" \
        "${ORIGINAL_DIR}/${LOCUS}.fasta" \
        > "${MAFFT_DIR}/${LOCUS}_aligned.fasta"
done

# 2. Removal of ambiguously aligned regions with Gblocks (optional)
for LOCUS in "${LOCI[@]}"; do
    GBLOCKS_INPUT="${GBLOCKS_DIR}/${LOCUS}_aligned.fasta"
    cp "${MAFFT_DIR}/${LOCUS}_aligned.fasta" "$GBLOCKS_INPUT"

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
        "${GBLOCKS_DIR}/${LOCUS}_gblocks.fasta"
done

# 3. Concatenation and creation of IQ-TREE partition files
MAFFT_ALIGNMENTS=()
GBLOCKS_ALIGNMENTS=()

for LOCUS in "${LOCI[@]}"; do
    MAFFT_ALIGNMENTS+=("${LOCUS}=${MAFFT_DIR}/${LOCUS}_aligned.fasta")
    GBLOCKS_ALIGNMENTS+=("${LOCUS}=${GBLOCKS_DIR}/${LOCUS}_gblocks.fasta")
done

python3 "$CONCATENATE_SCRIPT" \
    --output-fasta "${CONCAT_DIR}/nuclear_mafft_concatenated.fasta" \
    --output-partitions "${CONCAT_DIR}/nuclear_mafft.partitions" \
    --sample-prefix-order LMS,LMM,LAD \
    "${MAFFT_ALIGNMENTS[@]}"

python3 "$CONCATENATE_SCRIPT" \
    --output-fasta "${CONCAT_DIR}/nuclear_gblocks_concatenated.fasta" \
    --output-partitions "${CONCAT_DIR}/nuclear_gblocks.partitions" \
    --sample-prefix-order LMS,LMM,LAD \
    "${GBLOCKS_ALIGNMENTS[@]}"

# 4. Maximum-likelihood trees with per-locus model selection and UFBoot2
iqtree3 \
    -s "${CONCAT_DIR}/nuclear_mafft_concatenated.fasta" \
    -p "${CONCAT_DIR}/nuclear_mafft.partitions" \
    -m MFP \
    -B 1000 \
    -T AUTO \
    -ntmax "$THREADS" \
    -seed 20260802 \
    --prefix "${IQTREE_DIR}/nuclear_mafft"

iqtree3 \
    -s "${CONCAT_DIR}/nuclear_gblocks_concatenated.fasta" \
    -p "${CONCAT_DIR}/nuclear_gblocks.partitions" \
    -m MFP \
    -B 1000 \
    -T AUTO \
    -ntmax "$THREADS" \
    -seed 20260802 \
    --prefix "${IQTREE_DIR}/nuclear_gblocks"
