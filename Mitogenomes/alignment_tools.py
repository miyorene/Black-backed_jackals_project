#!/usr/bin/env python3
"""Concatenate, partition, and root mitochondrial alignments."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


LOCI = (
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "12S", "16S",
)
ALIGNMENT_ALPHABET = frozenset("ACGTRYSWKMBDHVN-?")


def die(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def read_samples(path: Path) -> list[str]:
    samples = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not samples or len(samples) != len(set(samples)):
        die(f"invalid expected-sample file: {path}")
    return samples


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        die(f"duplicate FASTA ID {name!r} in {path}")
                    records[name] = "".join(pieces).upper()
                name = line[1:].split()[0]
                pieces = []
            elif name is None:
                die(f"sequence occurs before first header in {path}")
            else:
                pieces.append(re.sub(r"\s+", "", line))
    if name is not None:
        if name in records:
            die(f"duplicate FASTA ID {name!r} in {path}")
        records[name] = "".join(pieces).upper()
    if not records:
        die(f"no sequences in {path}")
    for sample, sequence in records.items():
        invalid = sorted(set(sequence) - ALIGNMENT_ALPHABET)
        if invalid:
            die(f"invalid alignment symbols for {sample!r} in {path}: {''.join(invalid)}")
    return records


def validate_one(path: Path, expected_samples: list[str]) -> tuple[dict[str, str], int]:
    records = read_fasta(path)
    expected = set(expected_samples)
    observed = set(records)
    if observed != expected:
        die(
            f"taxa mismatch in {path}: missing={sorted(expected - observed) or 'none'}, "
            f"extra={sorted(observed - expected) or 'none'}"
        )
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        die(f"unaligned FASTA {path}; sequence lengths={sorted(lengths)}")
    length = lengths.pop()
    if length == 0:
        die(f"zero-length alignment: {path}")
    return records, length


def write_wrapped_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    with path.open("w") as out:
        for name, sequence in records:
            out.write(f">{name}\n")
            for position in range(0, len(sequence), width):
                out.write(sequence[position:position + width] + "\n")


def command_concat(args: argparse.Namespace) -> None:
    input_dir = Path(args.input_dir).resolve()
    samples = read_samples(Path(args.samples_file).resolve())
    output_fasta = Path(args.output_fasta).resolve()
    output_partitions = Path(args.output_partitions).resolve()
    for path in (output_fasta, output_partitions):
        path.parent.mkdir(parents=True, exist_ok=True)

    alignments: list[tuple[str, dict[str, str], int]] = []
    for locus in LOCI:
        records, length = validate_one(input_dir / f"{locus}.fasta", samples)
        alignments.append((locus, records, length))

    combined_records = [
        (sample, "".join(records[sample] for _, records, _ in alignments))
        for sample in samples
    ]
    write_wrapped_fasta(output_fasta, combined_records)

    total_length = len(combined_records[0][1])
    start = 1
    with output_partitions.open("w") as partitions:
        for locus, _, length in alignments:
            end = start + length - 1
            partitions.write(f"DNA, {locus} = {start}-{end}\n")
            start = end + 1

    if start - 1 != total_length:
        die("internal partition-length error")
    print(f"Concatenated {len(LOCI)} loci: {len(samples)} taxa x {total_length} sites")


@dataclass(eq=False)
class ParsedNode:
    label: str = ""
    length: str | None = None
    children: list["ParsedNode"] = field(default_factory=list)


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.position = 0

    def peek(self) -> str:
        return self.text[self.position] if self.position < len(self.text) else ""

    def skip_space_and_comments(self) -> None:
        while True:
            while self.peek().isspace():
                self.position += 1
            if self.peek() != "[":
                return
            depth = 0
            while self.position < len(self.text):
                char = self.text[self.position]
                self.position += 1
                if char == "[":
                    depth += 1
                elif char == "]":
                    depth -= 1
                    if depth == 0:
                        break
            if depth != 0:
                die("unterminated Newick comment")

    def parse_label(self) -> str:
        self.skip_space_and_comments()
        if self.peek() in {"'", '"'}:
            quote = self.peek()
            self.position += 1
            pieces: list[str] = []
            while self.position < len(self.text):
                char = self.text[self.position]
                self.position += 1
                if char == quote:
                    if self.peek() == quote:
                        pieces.append(quote)
                        self.position += 1
                        continue
                    return "".join(pieces)
                pieces.append(char)
            die("unterminated quoted Newick label")
        start = self.position
        while self.peek() and self.peek() not in ",():;[]" and not self.peek().isspace():
            self.position += 1
        return self.text[start:self.position]

    def parse_length(self) -> str | None:
        self.skip_space_and_comments()
        if self.peek() != ":":
            return None
        self.position += 1
        self.skip_space_and_comments()
        start = self.position
        while self.peek() and self.peek() not in ",();[]" and not self.peek().isspace():
            self.position += 1
        value = self.text[start:self.position]
        if not value:
            die("empty branch length in Newick tree")
        try:
            float(value)
        except ValueError:
            die(f"invalid Newick branch length {value!r}")
        return value

    def parse_subtree(self) -> ParsedNode:
        self.skip_space_and_comments()
        if self.peek() == "(":
            self.position += 1
            children = [self.parse_subtree()]
            while True:
                self.skip_space_and_comments()
                if self.peek() == ",":
                    self.position += 1
                    children.append(self.parse_subtree())
                elif self.peek() == ")":
                    self.position += 1
                    break
                else:
                    die(f"unexpected Newick token at position {self.position}")
            label = self.parse_label()
            length = self.parse_length()
            return ParsedNode(label=label, length=length, children=children)
        label = self.parse_label()
        if not label:
            die(f"missing leaf label at Newick position {self.position}")
        return ParsedNode(label=label, length=self.parse_length())

    def parse(self) -> ParsedNode:
        root = self.parse_subtree()
        self.skip_space_and_comments()
        if self.peek() == ";":
            self.position += 1
        self.skip_space_and_comments()
        if self.position != len(self.text):
            die(f"unexpected content after Newick tree at position {self.position}")
        return root


@dataclass
class Edge:
    length: str | None
    support: str = ""


def quote_label(label: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_.|/-]+", label):
        return label
    return "'" + label.replace("'", "''") + "'"


def half_length(length: str | None) -> str | None:
    if length is None:
        return None
    value = float(length) / 2.0
    return f"{value:.12g}"


def command_root(args: argparse.Namespace) -> None:
    input_tree = Path(args.input_tree).resolve()
    output_tree = Path(args.output_tree).resolve()
    outgroup = args.outgroup
    parsed_root = NewickParser(input_tree.read_text()).parse()

    adjacency: dict[int, dict[int, Edge]] = defaultdict(dict)
    leaf_names: dict[int, str] = {}
    next_id = 0

    def add_node(node: ParsedNode, parent_id: int | None = None) -> int:
        nonlocal next_id
        node_id = next_id
        next_id += 1
        adjacency.setdefault(node_id, {})
        if not node.children:
            leaf_names[node_id] = node.label
        if parent_id is not None:
            support = node.label if node.children else ""
            edge = Edge(node.length, support)
            adjacency[parent_id][node_id] = edge
            adjacency[node_id][parent_id] = edge
        for child in node.children:
            add_node(child, node_id)
        return node_id

    add_node(parsed_root)
    matches = [node_id for node_id, name in leaf_names.items() if name == outgroup]
    if len(matches) != 1:
        die(f"outgroup {outgroup!r} occurs {len(matches)} times in {input_tree}")
    out_id = matches[0]
    if len(adjacency[out_id]) != 1:
        die(f"outgroup {outgroup!r} is not a terminal taxon")
    neighbor = next(iter(adjacency[out_id]))
    original_edge = adjacency[out_id].pop(neighbor)
    adjacency[neighbor].pop(out_id)

    new_root = next_id
    adjacency[new_root] = {}
    split = half_length(original_edge.length)
    edge_to_outgroup = Edge(split, "")
    edge_to_ingroup = Edge(split, original_edge.support)
    adjacency[new_root][out_id] = edge_to_outgroup
    adjacency[out_id][new_root] = edge_to_outgroup
    adjacency[new_root][neighbor] = edge_to_ingroup
    adjacency[neighbor][new_root] = edge_to_ingroup

    def render(node_id: int, parent_id: int | None) -> str:
        children = [candidate for candidate in adjacency[node_id] if candidate != parent_id]
        if children:
            body = "(" + ",".join(render(child, node_id) for child in children) + ")"
        else:
            body = quote_label(leaf_names[node_id])
        if parent_id is not None:
            edge = adjacency[node_id][parent_id]
            if children and edge.support:
                body += quote_label(edge.support)
            if edge.length is not None:
                body += ":" + edge.length
        return body

    output_tree.parent.mkdir(parents=True, exist_ok=True)
    output_tree.write_text(render(new_root, None) + ";\n")
    print(f"Rooted {input_tree.name} with outgroup {outgroup}: {output_tree}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    concat = subparsers.add_parser("concat")
    concat.add_argument("--input-dir", required=True)
    concat.add_argument("--samples-file", required=True)
    concat.add_argument("--output-fasta", required=True)
    concat.add_argument("--output-partitions", required=True)
    concat.set_defaults(func=command_concat)

    root = subparsers.add_parser("root")
    root.add_argument("--input-tree", required=True)
    root.add_argument("--output-tree", required=True)
    root.add_argument("--outgroup", required=True)
    root.set_defaults(func=command_root)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
