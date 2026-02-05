# All Gas No Breakpoints!

Utiliy tools to extract statics from GFF files regarding rearrangements/breakpoints. 

## Install

```bash
cargo install --release .
```

## File formats

Your GFF3 needs:
- Feature type `SO:0000856` (syntenic markers)
- `ID=123` attribute with numeric IDs
- Optional `genome=species_name` attribute

Example:
```
chr1  .  SO:0000856  1000  2000  .  +  .  ID=42;genome=species_A
```

Mapping files look like:
```
genome_A: chr1, chr2, chr3
genome_B: scaffold_1, scaffold_2
```

## Tools

### breakpoints

The main event. Finds breakpoints between sequences.

**Single GFF mode** - compare all sequences against each other:
```bash
breakpoints markers.gff3 > breakpoints.txt
```

**Two GFF mode** - see how refined markers map to original ones:
```bash
breakpoints original.gff3 refined.gff3 --breakpoints original_bp.txt
```

Options:
- `--seqid2genome <file>` - provide genome/chromosome mapping
- `--no-dup` - keep duplicate marker IDs instead of removing them
- `--breakpoints <file>` - use precomputed breakpoints (two GFF mode only)

### gff

Quick stats about your GFF file:
```bash
gff info markers.gff3
```

Shows genomes, chromosomes, feature types, whether it's sorted, etc.

### gff2seq

Convert GFF to simple tab-separated format:
```bash
gff2seq markers.gff3 > sequences.tsv
```

Output: `genome  seqid  id1+,id2-,id3-,id4+...`

Negative IDs mean reverse strand.

### seqid2genome

Generate a mapping file from your inputs:
```bash
seqid2genome file.gff > mapping.txt
#or
seqid2genome f1.fasta f2.fasta ... > mapping.txt
```

Works with GFF3, FASTA, or text files listing FASTA paths.

## Quick start

```bash
# Check what you're working with
gff info mydata.gff

# (optional) Get genome mappings
seqid2genome mydata.gff > mapping.txt

# Find breakpoints
breakpoints mydata.gff --seqid2genome mapping.txt > breaks.txt
```
