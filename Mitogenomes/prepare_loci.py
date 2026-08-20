#!/usr/bin/env python3
"""Prepare and split 13-PCG + 12S + 16S mitochondrial loci.

The script deliberately uses each sample's own GFF coordinates.  It creates a
combined FASTA and BED6 file for ``bedtools getfasta`` and, after extraction,
splits the result into one FASTA file per locus.
"""

from __future__ import annotations

import argparse
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


LOCI = (
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6",
    "12S",
    "16S",
)
PCG_LOCI = frozenset(LOCI[:13])
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas")
GFF_SUFFIXES = (".gff", ".gff3")


def die(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def natural_key(value: str):
    return tuple(int(piece) if piece.isdigit() else piece.lower()
                 for piece in re.split(r"(\d+)", value))


def discover_files(directory: Path, suffixes: tuple[str, ...]) -> dict[str, Path]:
    if not directory.is_dir():
        die(f"input directory does not exist: {directory}")
    found: dict[str, Path] = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.suffix.lower() not in suffixes:
            continue
        if path.stem in found:
            die(f"multiple files have sample stem {path.stem!r} in {directory}")
        found[path.stem] = path
    if not found:
        die(f"no input files found in {directory}")
    return found


def read_single_fasta(path: Path) -> tuple[str, str]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    pieces: list[str] = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(pieces).upper()))
                name = line[1:].split()[0]
                pieces = []
            elif name is None:
                die(f"sequence occurs before the first header in {path}")
            else:
                pieces.append(re.sub(r"\s+", "", line))
    if name is not None:
        records.append((name, "".join(pieces).upper()))
    if len(records) != 1:
        die(f"{path} must contain exactly one FASTA record; found {len(records)}")
    name, sequence = records[0]
    if not sequence:
        die(f"empty sequence in {path}")
    invalid = sorted(set(sequence) - DNA_ALPHABET)
    if invalid:
        die(f"invalid DNA symbols in {path}: {''.join(invalid)}")
    return name, sequence


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key.strip()] = unquote(value.strip())
    return attributes


def canonical_locus(feature_type: str, attributes: dict[str, str]) -> str | None:
    name = attributes.get("Name", "").strip()
    if feature_type == "CDS":
        cleaned = re.sub(r"\s+CDS$", "", name, flags=re.IGNORECASE).upper()
        aliases = {
            "COI": "COX1",
            "COII": "COX2",
            "COIII": "COX3",
            "CO1": "COX1",
            "CO2": "COX2",
            "CO3": "COX3",
            "CYTB": "CYTB",
            "COB": "CYTB",
        }
        cleaned = aliases.get(cleaned, cleaned)
        return cleaned if cleaned in PCG_LOCI else None
    if feature_type == "rRNA":
        compact = re.sub(r"[^A-Za-z0-9]", "", name).upper()
        if compact in {"12S", "12SRRNA", "RRN12", "RNR1"}:
            return "12S"
        if compact in {"16S", "16SRRNA", "RRN16", "RNR2"}:
            return "16S"
    return None


def parse_target_features(path: Path, sequence_length: int):
    features: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open(encoding="utf-8-sig") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                die(f"{path}:{line_number}: expected 9 GFF3 columns")
            seqid, source, feature_type, start_text, end_text, score, strand, phase, attr_text = fields
            locus = canonical_locus(feature_type, parse_attributes(attr_text))
            if locus is None:
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                die(f"{path}:{line_number}: non-integer coordinates")
            if start < 1 or end < start or end > sequence_length:
                die(
                    f"{path}:{line_number}: invalid {locus} coordinates "
                    f"{start}-{end} for a {sequence_length}-bp FASTA"
                )
            if strand not in {"+", "-"}:
                die(f"{path}:{line_number}: {locus} has invalid strand {strand!r}")
            features[locus].append(
                {
                    "seqid": seqid,
                    "source": source,
                    "type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "line": line_number,
                }
            )

    selected: dict[str, dict[str, object]] = {}
    duplicate_notes: dict[str, str] = {}
    for locus in LOCI:
        candidates = features.get(locus, [])
        unique: dict[tuple[object, ...], dict[str, object]] = {}
        for item in candidates:
            key = (item["seqid"], item["start"], item["end"], item["strand"])
            unique.setdefault(key, item)
        if not unique:
            die(f"{path}: required feature {locus} is missing")
        if len(unique) > 1:
            details = ", ".join(
                f"{item['start']}-{item['end']}({item['strand']})"
                for item in unique.values()
            )
            die(f"{path}: multiple different annotations found for {locus}: {details}")
        selected[locus] = next(iter(unique.values()))
        if len(candidates) > 1:
            duplicate_notes[locus] = f"collapsed_{len(candidates)}_identical_GFF_records"
    return selected, duplicate_notes


def write_wrapped_fasta(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for position in range(0, len(sequence), width):
        handle.write(sequence[position:position + width] + "\n")


def command_prepare(args: argparse.Namespace) -> None:
    fasta_dir = Path(args.fasta_dir).resolve()
    gff_dir = Path(args.gff_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = discover_files(fasta_dir, FASTA_SUFFIXES)
    gff_files = discover_files(gff_dir, GFF_SUFFIXES)
    fasta_only = sorted(set(fasta_files) - set(gff_files), key=natural_key)
    gff_only = sorted(set(gff_files) - set(fasta_files), key=natural_key)
    if fasta_only or gff_only:
        die(
            "FASTA/GFF sample stems do not match; "
            f"FASTA only={fasta_only or 'none'}, GFF only={gff_only or 'none'}"
        )

    samples = sorted(fasta_files, key=natural_key)
    sequence_data: dict[str, tuple[str, str]] = {}
    feature_data: dict[str, dict[str, dict[str, object]]] = {}
    duplicate_data: dict[str, dict[str, str]] = {}

    for sample in samples:
        original_id, sequence = read_single_fasta(fasta_files[sample])
        sequence_data[sample] = (original_id, sequence)
        selected, duplicates = parse_target_features(gff_files[sample], len(sequence))
        feature_data[sample] = selected
        duplicate_data[sample] = duplicates

    medians: dict[str, float] = {}
    for locus in LOCI:
        medians[locus] = statistics.median(
            int(feature_data[sample][locus]["end"])
            - int(feature_data[sample][locus]["start"]) + 1
            for sample in samples
        )

    combined_fasta = output_dir / "combined_mitogenomes.fasta"
    bed_path = output_dir / "selected_15_loci.bed"
    audit_path = output_dir / "feature_audit.tsv"
    samples_path = output_dir / "expected_samples.txt"

    with combined_fasta.open("w") as out:
        for sample in samples:
            write_wrapped_fasta(out, sample, sequence_data[sample][1])

    with samples_path.open("w") as out:
        out.write("\n".join(samples) + "\n")

    with bed_path.open("w") as bed, audit_path.open("w") as audit:
        audit.write(
            "sample\tlocus\tfeature_type\tgff_seqid\tgff_start\tgff_end\t"
            "bed_start\tbed_end\tstrand\tlength\tmedian_length\tstatus\tnote\t"
            "fasta_header\tfasta_length\tgff_file\n"
        )
        for sample in samples:
            original_id, sequence = sequence_data[sample]
            for locus in LOCI:
                item = feature_data[sample][locus]
                start = int(item["start"])
                end = int(item["end"])
                length = end - start + 1
                median = medians[locus]
                threshold = max(15.0, median * 0.03)
                notes: list[str] = []
                status = "OK"
                if abs(length - median) > threshold:
                    status = "WARNING"
                    notes.append("length_outlier_gt_max(15bp,3pct)")
                if original_id != sample:
                    notes.append(f"FASTA_header_normalized_from_{original_id}")
                if str(item["seqid"]) != sample:
                    notes.append(f"GFF_seqid_normalized_from_{item['seqid']}")
                if locus in duplicate_data[sample]:
                    notes.append(duplicate_data[sample][locus])
                bed.write(
                    f"{sample}\t{start - 1}\t{end}\t{sample}|{locus}\t0\t{item['strand']}\n"
                )
                audit.write(
                    f"{sample}\t{locus}\t{item['type']}\t{item['seqid']}\t{start}\t{end}\t"
                    f"{start - 1}\t{end}\t{item['strand']}\t{length}\t{median:g}\t"
                    f"{status}\t{';'.join(notes) if notes else '.'}\t{original_id}\t"
                    f"{len(sequence)}\t{gff_files[sample]}\n"
                )

    print(f"Prepared {len(samples)} samples x {len(LOCI)} loci")
    print(f"Combined FASTA: {combined_fasta}")
    print(f"BED6: {bed_path}")
    print(f"Audit: {audit_path}")


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    pieces: list[str] = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(pieces).upper()))
                name = line[1:].strip()
                pieces = []
            elif name is None:
                die(f"sequence occurs before header in {path}")
            else:
                pieces.append(re.sub(r"\s+", "", line))
    if name is not None:
        records.append((name, "".join(pieces).upper()))
    if not records:
        die(f"no FASTA records in {path}")
    return records


def decode_bedtools_name(header: str) -> tuple[str, str]:
    token = header.split()[0]
    token = re.sub(r"\([+-]\)$", "", token)
    token = token.split("::", 1)[0]
    if "|" not in token:
        die(f"cannot parse bedtools FASTA header: {header!r}")
    sample, locus = token.rsplit("|", 1)
    if locus not in LOCI:
        die(f"unexpected locus {locus!r} in bedtools FASTA header {header!r}")
    return sample, locus


def read_expected_samples(path: Path) -> list[str]:
    samples = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not samples or len(samples) != len(set(samples)):
        die(f"invalid expected-sample file: {path}")
    return samples


def command_split(args: argparse.Namespace) -> None:
    input_fasta = Path(args.input_fasta).resolve()
    samples = read_expected_samples(Path(args.samples_file).resolve())
    expected = set(samples)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, dict[str, str]] = {locus: {} for locus in LOCI}

    for header, sequence in read_fasta_records(input_fasta):
        sample, locus = decode_bedtools_name(header)
        if sample not in expected:
            die(f"unexpected sample {sample!r} in {input_fasta}")
        if sample in grouped[locus]:
            die(f"duplicate extracted sequence for {sample}|{locus}")
        invalid = sorted(set(sequence) - DNA_ALPHABET)
        if invalid:
            die(f"invalid DNA symbols for {sample}|{locus}: {''.join(invalid)}")
        grouped[locus][sample] = sequence

    for locus in LOCI:
        missing = sorted(expected - set(grouped[locus]), key=natural_key)
        extra = sorted(set(grouped[locus]) - expected, key=natural_key)
        if missing or extra:
            die(f"{locus}: missing={missing or 'none'}, extra={extra or 'none'}")
        with (output_dir / f"{locus}.fasta").open("w") as out:
            for sample in samples:
                write_wrapped_fasta(out, sample, grouped[locus][sample])

    print(f"Split {len(samples) * len(LOCI)} sequences into {len(LOCI)} locus FASTAs")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="validate FASTA/GFF and create combined FASTA + BED6")
    prepare.add_argument("--fasta-dir", required=True)
    prepare.add_argument("--gff-dir", required=True)
    prepare.add_argument("--output-dir", required=True)
    prepare.set_defaults(func=command_prepare)

    split = subparsers.add_parser("split", help="split bedtools getfasta output into one FASTA per locus")
    split.add_argument("--input-fasta", required=True)
    split.add_argument("--samples-file", required=True)
    split.add_argument("--output-dir", required=True)
    split.set_defaults(func=command_split)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
