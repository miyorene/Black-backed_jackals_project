#!/bin/bash

#SBATCH --job-name=mitogenome_15_loci_gblocks
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# Analysis configuration
PROJECT_DIR="/mnt/tank/scratch/yunovikova/mitogenomes/mitogenomes"
PIPELINE_DIR="${PROJECT_DIR}/mitogenome_15_loci_pipeline"
FASTA_DIR="${PROJECT_DIR}/fasta"
GFF_DIR="${PROJECT_DIR}/gff"
RESULTS_DIR="${PROJECT_DIR}/phylogeny_15_loci_gblocks"
LOCI=(ATP6 ATP8 COX1 COX2 COX3 CYTB ND1 ND2 ND3 ND4 ND4L ND5 ND6 12S 16S)

export PATH="/nfs/home/yunovikova/micromamba/envs/project/bin:${PATH}"

PREP_DIR="${RESULTS_DIR}/01_prepared"
EXTRACTED_DIR="${RESULTS_DIR}/02_extracted_loci"
ALIGNED_DIR="${RESULTS_DIR}/03_mafft_loci"
GBLOCKS_DIR="${RESULTS_DIR}/04_gblocks_loci"
DATASETS_DIR="${RESULTS_DIR}/05_concatenated"
TREES_DIR="${RESULTS_DIR}/06_tree"

mkdir -p "$PREP_DIR" "$EXTRACTED_DIR" "$ALIGNED_DIR" "$GBLOCKS_DIR" \
    "$DATASETS_DIR" "$TREES_DIR"

echo "[1/7] Prepare FASTA and BED files from the individual GFF annotations"
python3 "${PIPELINE_DIR}/prepare_loci.py" prepare \
    --fasta-dir "$FASTA_DIR" \
    --gff-dir "$GFF_DIR" \
    --output-dir "$PREP_DIR"

echo "[2/7] Extract the 13 PCGs, 12S and 16S"
bedtools getfasta \
    -fi "${PREP_DIR}/combined_mitogenomes.fasta" \
    -bed "${PREP_DIR}/selected_15_loci.bed" \
    -s -nameOnly \
    -fo "${PREP_DIR}/extracted_15_loci.fasta"

python3 "${PIPELINE_DIR}/prepare_loci.py" split \
    --input-fasta "${PREP_DIR}/extracted_15_loci.fasta" \
    --samples-file "${PREP_DIR}/expected_samples.txt" \
    --output-dir "$EXTRACTED_DIR"

echo "[3/7] Align each locus with MAFFT L-INS-i"
for locus in "${LOCI[@]}"; do
    echo "  MAFFT: $locus"
    mafft --localpair --maxiterate 1000 --thread 16 \
        "${EXTRACTED_DIR}/${locus}.fasta" \
        > "${ALIGNED_DIR}/${locus}.fasta"
done

echo "[4/7] Filter each locus alignment with Gblocks"
N_SAMPLES="$(wc -l < "${PREP_DIR}/expected_samples.txt")"
MAJORITY=$((N_SAMPLES / 2 + 1))

for locus in "${LOCI[@]}"; do
    echo "  Gblocks: $locus"
    input_name="${locus}.mafft.fasta"
    cp "${ALIGNED_DIR}/${locus}.fasta" "${GBLOCKS_DIR}/${input_name}"

    (
        cd "$GBLOCKS_DIR"
        Gblocks "$input_name" \
            -t=d -b1="$MAJORITY" -b2="$MAJORITY" -b3=8 -b4=5 -b5=h
        mv "${input_name}-gb.fa" "${locus}.fasta"
        rm "$input_name" "${input_name}-gb.html"
    )
done

echo "[5/7] Concatenate the 15 Gblocks alignments and create partitions"
python3 "${PIPELINE_DIR}/alignment_tools.py" concat \
    --input-dir "$GBLOCKS_DIR" \
    --samples-file "${PREP_DIR}/expected_samples.txt" \
    --output-fasta "${DATASETS_DIR}/mitogenome_gblocks.fasta" \
    --output-nexus "${DATASETS_DIR}/mitogenome_gblocks.nex" \
    --output-partitions "${DATASETS_DIR}/mitogenome_gblocks.partitions" \
    --output-audit "${DATASETS_DIR}/mitogenome_gblocks_partition_audit.tsv"

echo "[6/7] Infer the partitioned maximum-likelihood tree"
iqtree3 \
    -s "${DATASETS_DIR}/mitogenome_gblocks.fasta" \
    -p "${DATASETS_DIR}/mitogenome_gblocks.partitions" \
    -m MFP \
    -B 1000 \
    -T 16 \
    -seed 20260802 \
    --prefix "${TREES_DIR}/mitogenome_gblocks" \
    --redo

echo "[7/7] Root the tree with UCI_1"
python3 "${PIPELINE_DIR}/alignment_tools.py" root \
    --input-tree "${TREES_DIR}/mitogenome_gblocks.treefile" \
    --output-tree "${TREES_DIR}/mitogenome_gblocks.rooted.treefile" \
    --outgroup UCI_1

echo "Pipeline completed successfully."
echo "Rooted tree: ${TREES_DIR}/mitogenome_gblocks.rooted.treefile"
