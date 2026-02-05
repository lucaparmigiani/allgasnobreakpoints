use anyhow::{Context, Result};
use clap::Parser;
use allgasnobreakpoints::gff::{load_genomes, load_seqid2genome, Genomes};
use std::io::{self, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    file_gff: String,
    #[arg(long)]
    seqid2genome: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let seqid_to_genome = match args.seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    let (mut genomes, _seqid_to_genome_out) = load_genomes(&args.file_gff, seqid_to_genome.as_ref(), true)
        .with_context(|| "Failed to parse the GFF")?;

    write_genomes(&mut genomes)?;

    Ok(())
}

fn write_genomes(
    genomes: &mut Genomes,
) -> Result<()> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    let mut genome_names: Vec<String> = genomes.keys().cloned().collect();
    genome_names.sort();
    for genome in genome_names {
        if let Some(seqs) = genomes.get(&genome) {
            for (seqid, seq) in seqs.iter() {
                write!(out, "{genome}\t{seqid}\t")?;
                let mut first = true;
                for &(id, strand) in &seq.markers {
                    if !first {
                        out.write_all(b",")?;
                    }
                    first = false;
                    write!(out, "{id}{strand}")?;
                }
                out.write_all(b"\n")?;
            }
        }
    }
    Ok(())
}
