#!/usr/bin/env python3

"""Concatenate aligned FASTA files and create an IQ-TREE partition file."""

import argparse
import re
from pathlib import Path


def read_fasta(path: Path):
    records = {}
    order = []
    sample = None
    sequence = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if sample is not None:
                if sample in records:
                    raise ValueError(f"Duplicate sample ID '{sample}' in {path}")
                records[sample] = "".join(sequence).upper()
                order.append(sample)

            sample = line[1:].split()[0]
            sequence = []
        else:
            if sample is None:
                raise ValueError(f"Sequence found before the first FASTA header in {path}")
            sequence.append("".join(line.split()))

    if sample is not None:
        if sample in records:
            raise ValueError(f"Duplicate sample ID '{sample}' in {path}")
        records[sample] = "".join(sequence).upper()
        order.append(sample)

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Sequences in {path} do not have equal lengths")

    return records, order, lengths.pop()


def parse_alignment(specification: str):
    if "=" not in specification:
        raise argparse.ArgumentTypeError(
            f"Alignment must be specified as NAME=FASTA: {specification}"
        )

    name, fasta_path = specification.split("=", 1)
    name = name.strip()
    fasta_path = fasta_path.strip()

    if not name or not fasta_path:
        raise argparse.ArgumentTypeError(
            f"Alignment must be specified as NAME=FASTA: {specification}"
        )

    return name, Path(fasta_path)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Concatenate aligned FASTA files by sample ID and create an "
            "IQ-TREE partition file. Missing loci are filled with '?'."
        )
    )
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-partitions", required=True, type=Path)
    parser.add_argument("--datatype", default="DNA")
    parser.add_argument(
        "--sample-prefix-order",
        help=(
            "Optional comma-separated order of sample ID prefixes. "
            "Samples are otherwise kept in their first-seen FASTA order."
        ),
    )
    parser.add_argument(
        "alignments",
        nargs="+",
        type=parse_alignment,
        metavar="NAME=FASTA",
        help="Partition name and aligned FASTA path, in concatenation order",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    alignments = []
    sample_order = []
    known_samples = set()
    partition_names = set()

    for partition_name, fasta_path in args.alignments:
        if partition_name in partition_names:
            raise ValueError(f"Duplicate partition name: {partition_name}")
        partition_names.add(partition_name)

        records, order, length = read_fasta(fasta_path)
        alignments.append((partition_name, records, length))

        for sample in order:
            if sample not in known_samples:
                known_samples.add(sample)
                sample_order.append(sample)

    if args.sample_prefix_order:
        prefixes = [
            prefix.strip()
            for prefix in args.sample_prefix_order.split(",")
            if prefix.strip()
        ]
        prefix_rank = {prefix: rank for rank, prefix in enumerate(prefixes)}

        def natural_parts(value):
            return tuple(
                int(part) if part.isdigit() else part.lower()
                for part in re.split(r"(\d+)", value)
            )

        def sample_sort_key(sample):
            prefix = sample.split("_", 1)[0]
            return (
                prefix_rank.get(prefix, len(prefix_rank)),
                natural_parts(sample),
            )

        sample_order.sort(key=sample_sort_key)

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_partitions.parent.mkdir(parents=True, exist_ok=True)

    with args.output_fasta.open("w", encoding="utf-8") as output_fasta:
        for sample in sample_order:
            concatenated_sequence = "".join(
                records.get(sample, "?" * length)
                for _, records, length in alignments
            )
            output_fasta.write(f">{sample}\n{concatenated_sequence}\n")

    start = 1
    with args.output_partitions.open("w", encoding="utf-8") as output_partitions:
        for partition_name, _, length in alignments:
            end = start + length - 1
            output_partitions.write(
                f"{args.datatype}, {partition_name} = {start}-{end}\n"
            )
            start = end + 1


if __name__ == "__main__":
    main()
