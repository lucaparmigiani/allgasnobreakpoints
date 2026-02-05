# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

```bash
# Build all binaries
cargo build --release

# Run a specific binary
cargo run --release --bin breakpoints -- <gff_file> [--no-dup] [--seqid2genome <file>]
cargo run --release --bin breakpoints -- <original_gff> <new_gff> [--no-dup] [--seqid2genome <file>] [--breakpoints <file>]
cargo run --release --bin gff -- info <gff_file>
cargo run --release --bin gff -- seq <gff_file> [--seqid2genome <file>]
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

- **breakpoints**: With one GFF, computes breakpoints across all sequences (parallel pairwise comparison). Outputs pairs of marker IDs that form breakpoint adjacencies. With two GFFs, compares original vs new, computing a "cover" relationship where new markers contain original markers based on coordinate overlap. Accepts optional `--breakpoints` file to provide precomputed breakpoints instead of computing from the original GFF.
- **gff**: Utility with subcommands:
  - `info`: GFF statistics (genomes, seqids, types, sorting status, contiguity)
  - `seq`: Converts GFF to a tab-separated format: `genome\tseqid\tmarker_list`
- **seqid2genome**: Generates seqid-to-genome mapping from GFF, FASTA, or list-of-FASTA files.

### Key Algorithms

The breakpoint computation (`breakpoints.rs:compute_breakpoints`) uses canonical signed adjacencies: each adjacency between consecutive markers is normalized so `(a,b)` and its reverse-complement form `(-b,-a)` map to the same canonical form. Breakpoints are identified as symmetric differences between adjacency sets of sequence pairs projected to their common marker IDs.

### File Formats

- **seqid2genome file**: `genome: seqid1, seqid2, ...` (one genome per line, colon-separated)
- **GFF requirements**: Expects `ID=<int>` and optionally `genome=<name>` in attributes column; filters to type `SO:0000856`
