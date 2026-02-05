use anyhow::{anyhow, Context, Result};
use clap::Parser;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input files (GFF, FASTA, or list-of-FASTA files). Type inferred by extension.
    #[arg(value_name = "INPUT", num_args = 1..)]
    inputs: Vec<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let mut seqid_to_genome: HashMap<String, String> = HashMap::new();
    let mut genome_to_seqids: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();

    for input in &args.inputs {
        match classify_input(input)? {
            InputKind::Gff => {
                load_gff_genomes(input, &mut seqid_to_genome, &mut genome_to_seqids)?;
            }
            InputKind::Fasta => {
                load_fasta_genome(input, &mut seqid_to_genome, &mut genome_to_seqids)?;
            }
            InputKind::FastaList => {
                let fasta_paths = read_fasta_list(input)?;
                for fasta in fasta_paths {
                    load_fasta_genome(&fasta, &mut seqid_to_genome, &mut genome_to_seqids)?;
                }
            }
        }
    }

    for (genome, seqids) in genome_to_seqids {
        let joined = seqids.into_iter().collect::<Vec<_>>().join(",");
        println!("{genome}: {joined}");
    }

    Ok(())
}

fn read_fasta_list(path: &Path) -> Result<Vec<PathBuf>> {
    let fh = File::open(path).with_context(|| format!("Opening {path:?}"))?;
    let br = BufReader::new(fh);
    let mut out = Vec::new();
    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {path:?}, line {i}"))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        for token in line.split_whitespace() {
            out.push(PathBuf::from(token));
        }
    }
    Ok(out)
}

#[derive(Debug, Copy, Clone)]
enum InputKind {
    Gff,
    Fasta,
    FastaList,
}

fn classify_input(path: &Path) -> Result<InputKind> {
    let (ext, ext2) = file_extensions_lower(path);

    if is_gff_ext(&ext, &ext2) {
        return Ok(InputKind::Gff);
    }
    if is_fasta_ext(&ext, &ext2) {
        return Ok(InputKind::Fasta);
    }

    Ok(InputKind::FastaList)
}

fn file_extensions_lower(path: &Path) -> (Option<String>, Option<String>) {
    let ext = path.extension().and_then(|s| s.to_str()).map(|s| s.to_lowercase());
    if ext.as_deref() == Some("gz") {
        let ext2 = path
            .file_stem()
            .and_then(|s| s.to_str())
            .and_then(|s| Path::new(s).extension())
            .and_then(|s| s.to_str())
            .map(|s| s.to_lowercase());
        return (ext, ext2);
    }
    (ext, None)
}

fn is_gff_ext(ext: &Option<String>, ext2: &Option<String>) -> bool {
    matches!(ext.as_deref(), Some("gff") | Some("gff3"))
        || (ext.as_deref() == Some("gz")
            && matches!(ext2.as_deref(), Some("gff") | Some("gff3")))
}

fn is_fasta_ext(ext: &Option<String>, ext2: &Option<String>) -> bool {
    let fasta = matches!(
        ext.as_deref(),
        Some("fa") | Some("fna") | Some("fasta") | Some("faa") | Some("fas")
    );
    if fasta {
        return true;
    }
    ext.as_deref() == Some("gz")
        && matches!(
            ext2.as_deref(),
            Some("fa") | Some("fna") | Some("fasta") | Some("faa") | Some("fas")
        )
}

fn load_gff_genomes(
    path: &Path,
    seqid_to_genome: &mut HashMap<String, String>,
    genome_to_seqids: &mut BTreeMap<String, BTreeSet<String>>,
) -> Result<()> {
    let fh = File::open(path).with_context(|| format!("Opening {path:?}"))?;
    let br = BufReader::new(fh);
    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {path:?}, line {i}"))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            return Err(anyhow!("Malformed GFF line (need 9 columns): {line}"));
        }
        let seqid = cols[0].to_string();
        if let Some(genome) = parse_genome_attribute(cols[8]) {
            insert_mapping(&seqid, &genome, seqid_to_genome, genome_to_seqids)?;
        }
    }
    Ok(())
}

fn parse_genome_attribute(attributes: &str) -> Option<String> {
    for field in attributes.split(';') {
        if let Some(rest) = field.strip_prefix("genome=") {
            if !rest.is_empty() {
                return Some(rest.to_string());
            }
        }
    }
    None
}

fn load_fasta_genome(
    path: &Path,
    seqid_to_genome: &mut HashMap<String, String>,
    genome_to_seqids: &mut BTreeMap<String, BTreeSet<String>>,
) -> Result<()> {
    let genome = genome_name_from_path(path)?;
    let fh = File::open(path).with_context(|| format!("Opening {path:?}"))?;
    let br = BufReader::new(fh);
    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {path:?}, line {i}"))?;
        if let Some(seqid) = parse_fasta_header(&line) {
            insert_mapping(&seqid, &genome, seqid_to_genome, genome_to_seqids)?;
        }
    }
    Ok(())
}

fn parse_fasta_header(line: &str) -> Option<String> {
    if !line.starts_with('>') {
        return None;
    }
    let header = line[1..].trim();
    if header.is_empty() {
        return None;
    }
    Some(
        header
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string(),
    )
}

fn genome_name_from_path(path: &Path) -> Result<String> {
    let file_name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow!("Cannot derive genome name from path {path:?}"))?;
    if file_name.is_empty() {
        return Err(anyhow!("Cannot derive genome name from path {path:?}"));
    }
    Ok(file_name.to_string())
}

fn insert_mapping(
    seqid: &str,
    genome: &str,
    seqid_to_genome: &mut HashMap<String, String>,
    genome_to_seqids: &mut BTreeMap<String, BTreeSet<String>>,
) -> Result<()> {
    if let Some(existing) = seqid_to_genome.get(seqid) {
        if existing != genome {
            return Err(anyhow!(
                "Seqid {seqid} is assigned to multiple genomes: {existing} and {genome}"
            ));
        }
        return Ok(());
    }

    seqid_to_genome.insert(seqid.to_string(), genome.to_string());
    genome_to_seqids
        .entry(genome.to_string())
        .or_default()
        .insert(seqid.to_string());
    Ok(())
}
