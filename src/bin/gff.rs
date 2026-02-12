use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use colored::Colorize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use allgasnobreakpoints::gff::*;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Info {
        file_gff: String,
        #[arg(long, short = 'a')]
        all: bool,
    },
    Seq {
        file_gff: String,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
        #[arg(long)]
        no_dup: bool,
    },
    Dup {
        file_gff: String,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
    },
    Find {
        file_gff: String,
        id: usize,
        #[arg(long, conflicts_with = "genome")]
        seqid: Option<String>,
        #[arg(long, conflicts_with = "seqid")]
        genome: Option<String>,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
    },
    IsDup {
        file_gff: String,
        id: usize,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
    },
    Sample {
        file_gff: String,
        n: usize,
    },
    Block {
        file_gff_original: String,
        file_gff_new: String,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
        /// Remove duplicated elements from the original GFF
        #[arg(long)]
        no_dup_element: bool,
        /// Remove duplicated markers from the new GFF
        #[arg(long)]
        no_dup_marker: bool,
        /// Extend marker end coordinates by this many base pairs (default: 0)
        #[arg(long, default_value = "0")]
        extend: u64,
    },
    Partition {
        file_gff_original: String,
        file_gff_new: String,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
        #[arg(long)]
        test: bool,
        /// Extend marker end coordinates by this many base pairs (default: 0)
        #[arg(long, default_value = "0")]
        extend: u64,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Info { file_gff, all} => info(&file_gff, all),
        Command::Seq { file_gff, seqid2genome, no_dup } => seq(&file_gff, seqid2genome, no_dup),
        Command::Dup { file_gff, seqid2genome } => dup(&file_gff, seqid2genome),
        Command::Find { file_gff, id, seqid, genome, seqid2genome } => {
            find(&file_gff, id, seqid, genome, seqid2genome)
        }
        Command::IsDup { file_gff, id, seqid2genome } => {
            is_dup(&file_gff, id, seqid2genome)
        }
        Command::Sample { file_gff, n } => {
            sample(&file_gff, n)
        }
        Command::Block { file_gff_original, file_gff_new, seqid2genome, no_dup_element, no_dup_marker, extend } => {
            block(&file_gff_original, &file_gff_new, seqid2genome, no_dup_element, no_dup_marker, extend)
        }
        Command::Partition { file_gff_original, file_gff_new, seqid2genome, test, extend } => {
            partition(&file_gff_original, &file_gff_new, seqid2genome, test, extend)
        }
    }
}

fn info(path: &str, all: bool) -> Result<()> {
    let stats = scan_gff(path).with_context(|| "Failed to parse the GFF")?;

    println!("{:<25} {}", "data_lines:", stats.data_lines);
    println!("{:<25} {}", "genomes:", stats.genomes.len());
    println!("{:<25} {}", "seqids:", stats.seqids.len());
    println!("{:<25} {}", "types_count:", stats.types.len());
    println!("{:<25} {}", "types:", join_sorted(&stats.types));

    let sorted_status = if stats.sorted_by_start {
        "true".green().to_string()
    } else {
        "false".red().to_string()
    };
    println!("{:<25} {}", "sorted_by_start:", sorted_status);

    // genomes_contiguous is not applicable if any line is missing a genome
    let contiguous_status = if stats.missing_genome_count > 0 {
        "-".to_string()
    } else if stats.genomes_contiguous {
        "true".green().to_string()
    } else {
        "false".red().to_string()
    };
    println!("{:<25} {}", "genomes_contiguous:", contiguous_status);

    // seqid_multi_genome is not applicable if all lines are missing genome
    let seqid_multi = if stats.missing_genome_count == stats.data_lines {
        "-".to_string()
    } else if stats.seqid_multi_genome_count == 0 {
        format!("{}", stats.seqid_multi_genome_count).green().to_string()
    } else {
        format!("{}", stats.seqid_multi_genome_count).red().to_string()
    };
    println!("{:<25} {}", "seqid_multi_genome:", seqid_multi);

    let missing_genome = if stats.missing_genome_count == 0 {
        format!("{}", stats.missing_genome_count).green().to_string()
    } else {
        format!("{}", stats.missing_genome_count).red().to_string()
    };
    println!("{:<25} {}", "missing_genome:", missing_genome);

    let missing_id = if stats.missing_id_count == 0 {
        format!("{}", stats.missing_id_count).green().to_string()
    } else {
        format!("{}", stats.missing_id_count).red().to_string()
    };
    println!("{:<25} {}", "missing_id:", missing_id);

    let malformed = if stats.malformed_lines == 0 {
        format!("{}", stats.malformed_lines).green().to_string()
    } else {
        format!("{}", stats.malformed_lines).red().to_string()
    };
    println!("{:<25} {}", "malformed_lines:", malformed);

    if all {
        println!("{:<25} {}", "unique_genomes:", join_sorted(&stats.genomes));
        println!("{:<25} {}", "unique_seqids:", join_sorted(&stats.seqids));
        println!("seqids_per_genome:");
        let mut genomes_sorted: Vec<String> = stats.seqids_per_genome.keys().cloned().collect();
        genomes_sorted.sort();
        for genome in genomes_sorted {
            let count = stats
                .seqids_per_genome
                .get(&genome)
                .map(|s| s.len())
                .unwrap_or(0);
            println!("  {:<23} {}", genome, count);
        }
    }

    Ok(())
}

struct InfoStats {
    genomes: HashSet<String>,
    seqids: HashSet<String>,
    types: HashSet<String>,
    seqids_per_genome: HashMap<String, HashSet<String>>,
    sorted_by_start: bool,
    genomes_contiguous: bool,
    missing_genome_count: u64,
    missing_id_count: u64,
    malformed_lines: u64,
    seqid_multi_genome_count: u64,
    data_lines: u64,
}

fn scan_gff(path: &str) -> Result<InfoStats> {
    let fh = File::open(path).with_context(|| format!("Opening {path}"))?;
    let br = BufReader::new(fh);

    let mut genomes: HashSet<String> = HashSet::new();
    let mut seqids: HashSet<String> = HashSet::new();
    let mut types: HashSet<String> = HashSet::new();
    let mut seqids_per_genome: HashMap<String, HashSet<String>> = HashMap::new();
    let mut last_start: HashMap<String, u64> = HashMap::new();
    let mut sorted_by_start = true;
    let mut genomes_contiguous = true;
    let mut finished_genomes: HashSet<String> = HashSet::new();
    let mut current_genome: Option<String> = None;

    let mut missing_genome_count = 0u64;
    let mut missing_id_count = 0u64;
    let mut malformed_lines = 0u64;
    let mut seqid_genomes: HashMap<String, HashSet<String>> = HashMap::new();
    let mut data_lines = 0u64;

    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {path}, line {i}"))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        data_lines += 1;
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            malformed_lines += 1;
            continue;
        }

        let seqid = cols[0].to_string();
        seqids.insert(seqid.clone());
        types.insert(cols[2].to_string());

        let start: u64 = match cols[3].parse() {
            Ok(v) => v,
            Err(_) => {
                malformed_lines += 1;
                continue;
            }
        };

        if let Some(prev) = last_start.get(&seqid) {
            if start < *prev {
                sorted_by_start = false;
            }
        }
        last_start.insert(seqid.clone(), start);

        let attributes = cols[8];
        let mut genome: Option<String> = None;
        let mut has_id = false;
        for field in attributes.split(';') {
            if let Some(rest) = field.strip_prefix("genome=") {
                if !rest.is_empty() {
                    genome = Some(rest.to_string());
                }
            } else if field.starts_with("ID=") {
                has_id = true;
            }
        }

        if !has_id {
            missing_id_count += 1;
        }

        if let Some(genome_val) = genome {
            seqid_genomes
                .entry(seqid.clone())
                .or_default()
                .insert(genome_val.clone());
            genomes.insert(genome_val.clone());
            seqids_per_genome
                .entry(genome_val.clone())
                .or_default()
                .insert(seqid.clone());

            match current_genome.as_deref() {
                None => current_genome = Some(genome_val),
                Some(curr) if curr == genome_val => {}
                Some(curr) => {
                    finished_genomes.insert(curr.to_string());
                    if finished_genomes.contains(&genome_val) {
                        genomes_contiguous = false;
                    }
                    current_genome = Some(genome_val);
                }
            }
        } else {
            missing_genome_count += 1;
        }
    }

    Ok(InfoStats {
        genomes,
        seqids,
        types,
        seqids_per_genome,
        sorted_by_start,
        genomes_contiguous,
        missing_genome_count,
        missing_id_count,
        malformed_lines,
        seqid_multi_genome_count: seqid_multi_genome_count(&seqid_genomes),
        data_lines,
    })
}

fn seqid_multi_genome_count(map: &HashMap<String, HashSet<String>>) -> u64 {
    map.values().filter(|s| s.len() > 1).count() as u64
}

fn join_sorted(set: &HashSet<String>) -> String {
    let mut v: Vec<String> = set.iter().cloned().collect();
    v.sort();
    v.join(",")
}

fn seq(file_gff: &str, seqid2genome: Option<PathBuf>, no_dup: bool) -> Result<()> {
    let seqid_to_genome = match seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    let (mut genomes, _seqid_to_genome_out) = load_genomes(file_gff, seqid_to_genome.as_ref(), true)
        .with_context(|| "Failed to parse the GFF")?;

    if no_dup {
        let dup_ids = duplicated_ids_by_genome(&genomes);
        remove_duplicates(&mut genomes, &dup_ids);
    }

    write_genomes(&mut genomes)?;

    Ok(())
}

fn write_genomes(genomes: &mut Genomes) -> Result<()> {
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
                        out.write(b",")?;
                    }
                    first = false;
                    write!(out, "{id}{strand}")?;
                }
                out.write(b"\n")?;
            }
        }
    }
    Ok(())
}

fn find(
    file_gff: &str,
    target_id: usize,
    seqid_filter: Option<String>,
    genome_filter: Option<String>,
    seqid2genome: Option<PathBuf>,
) -> Result<()> {
    let seqid_to_genome = match seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    let fh = File::open(file_gff).with_context(|| format!("Opening {file_gff}"))?;
    let br = BufReader::new(fh);
    let mut out = io::BufWriter::new(io::stdout().lock());

    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {file_gff}, line {i}"))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let seqid = cols[0];
        let attributes = cols[8];

        // Parse attributes to get ID and genome
        let mut id: Option<usize> = None;
        let mut genome: Option<String> = None;

        for field in attributes.split(';') {
            if let Some(rest) = field.strip_prefix("ID=") {
                if let Ok(parsed_id) = rest.parse::<usize>() {
                    id = Some(parsed_id);
                }
            } else if let Some(rest) = field.strip_prefix("genome=") {
                genome = Some(rest.to_string());
            }
        }

        // If no genome in GFF, try to get it from seqid2genome mapping
        if genome.is_none() {
            if let Some(ref map) = seqid_to_genome {
                genome = map.get(seqid).cloned();
            }
        }

        // Check if this line matches our criteria
        let id_matches = id == Some(target_id);
        if !id_matches {
            continue;
        }

        let seqid_matches = seqid_filter.as_ref().map_or(true, |filter| seqid == filter);
        let genome_matches = genome_filter.as_ref().map_or(true, |filter| {
            genome.as_ref().map_or(false, |g| g == filter)
        });

        if seqid_matches && genome_matches {
            writeln!(out, "{}", line)?;
        }
    }

    Ok(())
}

fn sample(file_gff: &str, n: usize) -> Result<()> {
    use rand::seq::index::sample as random_sample;

    // First pass: count data lines (fast - no allocations)
    let fh = File::open(file_gff).with_context(|| format!("Opening {file_gff}"))?;
    let br = BufReader::new(fh);

    let mut total_data_lines = 0usize;
    for line in br.split(b'\n') {
        let line = line?;
        if !line.starts_with(b"#") && !line.is_empty() {
            total_data_lines += 1;
        }
    }

    let sample_size = n.min(total_data_lines);
    let mut rng = rand::thread_rng();

    let mut selected_indices: Vec<usize> = random_sample(&mut rng, total_data_lines, sample_size)
        .into_iter()
        .collect();
    selected_indices.sort_unstable();

    // Second pass: output headers and selected data lines
    let fh = File::open(file_gff)?;
    let br = BufReader::new(fh);
    let mut out = io::BufWriter::new(io::stdout().lock());

    let mut data_line_idx = 0usize;
    let mut next_idx_pos = 0usize;

    for line in br.split(b'\n') {
        let line = line?;

        if line.starts_with(b"#") {
            out.write_all(&line)?;
            out.write_all(b"\n")?;
        } else if !line.is_empty() {
            // Check if this is the next line we want
            if next_idx_pos < selected_indices.len() && data_line_idx == selected_indices[next_idx_pos] {
                out.write_all(&line)?;
                out.write_all(b"\n")?;
                next_idx_pos += 1;
            }
            data_line_idx += 1;
        }
    }

    Ok(())
}

fn is_dup(file_gff: &str, target_id: usize, seqid2genome: Option<PathBuf>) -> Result<()> {
    let seqid_to_genome = match seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    let (genomes, _seqid_to_genome_out) = load_genomes(file_gff, seqid_to_genome.as_ref(), false)
        .with_context(|| "Failed to parse the GFF")?;

    let mut genomes_with_dup: Vec<String> = Vec::new();

    for (genome_name, seqs) in &genomes {
        let mut count = 0;

        for seq in seqs.values() {
            for &(id, _) in &seq.markers {
                if id == target_id {
                    count += 1;
                    if count > 1 {
                        genomes_with_dup.push(genome_name.clone());
                        break;
                    }
                }
            }
            if count > 1 {
                break;
            }
        }
    }

    if !genomes_with_dup.is_empty() {
        genomes_with_dup.sort();
        let mut out = io::BufWriter::new(io::stdout().lock());
        for genome in genomes_with_dup {
            write!(out, " {}", genome)?;
        }
        writeln!(out)?;
    }

    Ok(())
}

fn dup(file_gff: &str, seqid2genome: Option<PathBuf>) -> Result<()> {
    let seqid_to_genome = match seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    let (genomes, _seqid_to_genome_out) = load_genomes(file_gff, seqid_to_genome.as_ref(), false)
        .with_context(|| "Failed to parse the GFF")?;

    let mut duplicates: HashMap<usize, Vec<String>> = HashMap::new();

    for (genome_name, seqs) in &genomes {
        let mut seen: HashSet<usize> = HashSet::new();
        let mut genome_dups: HashSet<usize> = HashSet::new();

        for seq in seqs.values() {
            for &(id, _) in &seq.markers {
                if !seen.insert(id) {
                    genome_dups.insert(id);
                }
            }
        }

        for &dup_id in &genome_dups {
            duplicates
                .entry(dup_id)
                .or_default()
                .push(genome_name.clone());
        }
    }

    let mut out = io::BufWriter::new(io::stdout().lock());
    let mut ids: Vec<usize> = duplicates.keys().copied().collect();
    ids.sort();

    for id in ids {
        if let Some(genome_list) = duplicates.get(&id) {
            let mut sorted_genomes = genome_list.clone();
            sorted_genomes.sort();
            write!(out, "{}", id)?;
            for genome in sorted_genomes {
                write!(out, " {}", genome)?;
            }
            writeln!(out)?;
        }
    }

    Ok(())
}

fn block(
    file_gff_original: &str,
    file_gff_new: &str,
    seqid2genome: Option<PathBuf>,
    no_dup_element: bool,
    no_dup_marker: bool,
    extend: u64,
) -> Result<()> {
    let seqid_to_genome = match seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    let (mut genomes_original, seqid_to_genome) =
        load_genomes(file_gff_original, seqid_to_genome.as_ref(), true)
            .with_context(|| "Failed to parse the original GFF")?;

    if no_dup_element {
        let dup_orig = duplicated_ids_by_genome(&genomes_original);
        remove_duplicates(&mut genomes_original, &dup_orig);
    }

    let (mut genomes_new, _seqid_to_genome) =
        load_genomes(file_gff_new, Some(&seqid_to_genome), true)
            .with_context(|| "Failed to parse the new GFF")?;

    if no_dup_marker {
        let dup_new = duplicated_ids_by_genome(&genomes_new);
        remove_duplicates(&mut genomes_new, &dup_new);
    }

    let genomes_new_blocks = compute_genome_blocks(&genomes_original, &genomes_new, extend);

    print_genomes_new_blocks(&genomes_new_blocks)?;

    Ok(())
}

fn partition(
    file_gff_original: &str,
    file_gff_new: &str,
    seqid2genome: Option<PathBuf>,
    test: bool,
    extend: u64,
) -> Result<()> {
    let seqid_to_genome = match seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    let (genomes_original, seqid_to_genome) =
        load_genomes(file_gff_original, seqid_to_genome.as_ref(), true)
            .with_context(|| "Failed to parse the original GFF")?;

    let (genomes_new, _seqid_to_genome) =
        load_genomes(file_gff_new, Some(&seqid_to_genome), true)
            .with_context(|| "Failed to parse the new GFF")?;

    let genomes_new_blocks = compute_genome_blocks(&genomes_original, &genomes_new, extend);

    let partition = compute_partition(&genomes_new, &genomes_new_blocks);
    
    let mut out = io::BufWriter::new(io::stdout().lock());

    if test {
        let mut reverse_map: HashMap<usize, HashSet<usize>> = HashMap::new();
        for (new_marker_id, original_marker_ids) in partition.iter().enumerate() {
            for &original_id in original_marker_ids {
                reverse_map
                    .entry(original_id)
                    .or_insert_with(HashSet::new)
                    .insert(new_marker_id);
            }
        }

        let mut multi_partition_elements: Vec<usize> = reverse_map
            .iter()
            .filter(|(_, parts)| parts.len() > 1)
            .map(|(&id, _)| id)
            .collect();
        multi_partition_elements.sort();

        for original_id in multi_partition_elements {
            if let Some(parts) = reverse_map.get(&original_id) {
                write!(out, "{}->", original_id)?;
                let mut sorted_parts: Vec<usize> = parts.iter().copied().collect();
                sorted_parts.sort();
                for part_id in sorted_parts {
                    write!(out, " {}", part_id)?;
                }
                writeln!(out)?;
            }
        }
    } else {
        for (new_marker_id, original_marker_ids) in partition.iter().enumerate() {
            if !original_marker_ids.is_empty() {
                write!(out, "{}:", new_marker_id)?;
                let mut sorted_ids: Vec<usize> = original_marker_ids.iter().copied().collect();
                sorted_ids.sort();
                for original_id in sorted_ids {
                    write!(out, " {}", original_id)?;
                }
                writeln!(out)?;
            }
        }
    }

    Ok(())
}

fn print_genomes_new_blocks(genomes_new_blocks: &GenomesBlocks) -> Result<()> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    let mut genome_names: Vec<String> = genomes_new_blocks.keys().cloned().collect();
    genome_names.sort();
    for genome in genome_names {
        if let Some(seqs) = genomes_new_blocks.get(&genome) {
            for (seqid, blocks) in seqs.iter() {
                write!(out, "{genome}\t{seqid}\t")?;
                for block in blocks {
                    if !block.is_empty() {
                        let mut first = true;
                        write!(out, "[")?;
                        for (id, strand) in block.iter() {
                            if !first {
                                out.write(b",")?;
                            }
                            first = false;
                            write!(out, "{id}{strand}")?;
                        }
                        write!(out, "]")?;
                    }
                }
                out.write(b"\n")?;
            }
        }
    }
    Ok(())
}
