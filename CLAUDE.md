# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

```bash
# Build all binaries
cargo build --release

# Run a specific binary
cargo run --release --bin breakpoints -- compute <gff_file> [--no-dup] [--seqid2genome <file>]
cargo run --release --bin breakpoints -- synteny <original_gff> <blocks_gff> [--seqid2genome <file>] [--breakpoints <file>] [--output <folder>] [--extend <bp>] [--supporting-elements]
cargo run --release --bin gff -- info <gff_file> [-a]
cargo run --release --bin gff -- seq <gff_file> [--seqid2genome <file>]
cargo run --release --bin gff -- block <original_gff> <new_gff> [--seqid2genome <file>] [--extend <bp>]
cargo run --release --bin gff -- partition <original_gff> <new_gff> [--seqid2genome <file>] [--extend <bp>]
cargo run --release --bin seqid2genome -- <inputs...>

# Run tests
cargo test

# Check without building
cargo check
```

## Architecture

This is a Rust crate (`allgasnobreakpoints`) for analyzing genomic breakpoints from GFF3 files containing syntenic markers (type `SO:0000856`).

### Core Data Structures (src/gff.rs)

- `Seq`: Stores markers for a single sequence (chromosome/scaffold) with parallel vectors for marker IDs, strands, start/end positions
- `Genomes`: `HashMap<String, HashMap<String, Seq>>` - nested map of genome -> seqid -> sequence data
- `Record`: Parsed GFF3 record with seqid, coordinates, ID, strand, and genome attribution

### Binaries

- **breakpoints**: Has two subcommands:
  - `compute`: Computes breakpoints from a single GFF across all sequences (parallel pairwise comparison). Outputs pairs of marker IDs that form breakpoint adjacencies to stdout.
  - `synteny`: Evaluates synteny block construction by comparing original GFF vs blocks GFF. Computes a "cover" relationship where block markers contain original markers based on coordinate overlap. Creates an output folder (default: "synteny_check") containing three files: `false_positive.txt`, `false_negative.txt`, and `true_positive.txt`. Accepts optional `--breakpoints` file to provide precomputed breakpoints instead of computing from the original GFF. The `--extend` option (default: 0) extends the end coordinate of each marker by the specified number of base pairs, allowing more original markers to be included in blocks. The `--supporting-elements` flag adds breakpoints to false negatives if either marker doesn't appear in any new marker (no supporting elements). True positives are all original breakpoints (no duplicates) that are not false negatives.
- **gff**: Utility with subcommands:
  - `info`: GFF statistics (genomes, seqids, types, sorting status, contiguity). Use `-a` or `--all` for detailed output.
  - `seq`: Converts GFF to a tab-separated format: `genome\tseqid\tmarker_list`
  - `block`: Takes two GFFs (original and new) and computes which markers from the original are contained in each marker of the new GFF based on coordinate overlap. Outputs blocks as: `genome\tseqid\t[marker1,marker2,...][marker3,...]...`. The `--extend` option (default: 0) extends the end coordinate of each new marker by the specified number of base pairs.
  - `partition`: Computes the partition of original markers into new markers. The `--extend` option (default: 0) extends the end coordinate of each new marker by the specified number of base pairs.
- **seqid2genome**: Generates seqid-to-genome mapping from GFF, FASTA, or list-of-FASTA files.

### Key Algorithms

The breakpoint computation (`breakpoints.rs:compute_breakpoints`) uses canonical signed adjacencies: each adjacency between consecutive markers is normalized so `(a,b)` and its reverse-complement form `(-b,-a)` map to the same canonical form. Breakpoints are identified as symmetric differences between adjacency sets of sequence pairs projected to their common marker IDs.

### File Formats

- **seqid2genome file**: `genome: seqid1, seqid2, ...` (one genome per line, colon-separated)
- **GFF requirements**: Expects `ID=<int>` and optionally `genome=<name>` in attributes column; filters to type `SO:0000856`
