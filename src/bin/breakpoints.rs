use anyhow::{anyhow, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use allgasnobreakpoints::gff::*;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
struct FalsePositiveDebugInfo {
    genome_a: String,
    seqid_a: String,
    genome_b: String,
    seqid_b: String,
    marker_x: (usize, char),
    marker_y: (usize, char),
    block_x: Vec<(usize, char)>,
    block_y: Vec<(usize, char)>,
}

struct DebugContext<'a> {
    genome_a: &'a str,
    seqid_a: &'a str,
    genome_b: &'a str,
    seqid_b: &'a str,
}

#[derive(ValueEnum, Clone, Debug)]
enum FnMode {
    /// Compute false negatives from partition
    Partition,
    /// Compute false negatives from blocks
    Blocks,
}

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Compute breakpoints from a GFF
    Compute {
        file_gff: String,
        #[arg(long)]
        no_dup: bool,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
    },
    /// Evaluate synteny block construction
    Synteny {
        file_gff_original: String,
        file_gff_blocks: String,
        #[arg(long)]
        seqid2genome: Option<PathBuf>,
        #[arg(long)]
        breakpoints: Option<PathBuf>,
        /// Output folder name for results
        #[arg(long, default_value = "synteny_check")]
        output: String,
        /// Method to compute false negatives: partition or blocks
        #[arg(long, value_enum, default_value = "partition")]
        fn_mode: FnMode,
        #[arg(long)]
        debug: bool,
        /// Ignore new breakpoints formed by the new blocks (default: check them)
        #[arg(long)]
        ignore_new_breakpoints: bool,
    },
}

fn main() -> Result<()> {
    let args = Args::parse();

    let num_cores = num_cpus::get();
    ThreadPoolBuilder::new()
        .num_threads(num_cores)
        .build_global()
        .unwrap();

    match args.command {
        Command::Compute { file_gff, no_dup, seqid2genome } => {
            let seqid_to_genome = match seqid2genome.as_ref() {
                Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
                None => None,
            };
            run_single_gff(&file_gff, seqid_to_genome.as_ref(), no_dup)
        }
        Command::Synteny {
            file_gff_original,
            file_gff_blocks,
            seqid2genome,
            breakpoints,
            output,
            fn_mode,
            debug,
            ignore_new_breakpoints,
        } => {
            let seqid_to_genome = match seqid2genome.as_ref() {
                Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
                None => None,
            };
            run_synteny(
                &file_gff_original,
                &file_gff_blocks,
                seqid_to_genome.as_ref(),
                breakpoints.as_ref(),
                &output,
                fn_mode,
                debug,
                ignore_new_breakpoints,
            )
        }
    }
}

fn run_single_gff(
    file_gff: &str,
    seqid_to_genome: Option<&HashMap<String, String>>,
    no_dup: bool,
) -> Result<()> {
    let (mut genomes, _seqid_to_genome_out) = load_genomes(file_gff, seqid_to_genome, true)
        .with_context(|| "Failed to parse the GFF")?;

    if no_dup {
        let duplicated_ids = duplicated_ids_by_genome(&genomes);
        remove_duplicates(&mut genomes, &duplicated_ids);
        eprintln!("duplicated: {}", duplicated_ids.len());
    }

    let seqs = genomes_to_sequences(&genomes);

    let breakpoints = compute_breakpoints(&seqs);

    let mut out = io::BufWriter::new(io::stdout().lock());
    for (a, b) in &breakpoints {
        writeln!(out, "{a} {b}")?;
    }

    Ok(())
}

fn run_synteny(
    file_gff_original: &str,
    file_gff_blocks: &str,
    seqid_to_genome: Option<&HashMap<String, String>>,
    breakpoints_file: Option<&PathBuf>,
    output_folder: &str,
    fn_mode: FnMode,
    debug: bool,
    ignore_new_breakpoints: bool,
) -> Result<()> {
    let (mut genomes_original, seqid_to_genome) =
        load_genomes(file_gff_original, seqid_to_genome, true)
            .with_context(|| "Failed to parse the original GFF")?;

    let (mut genomes_new, _seqid_to_genome) =
        load_genomes(file_gff_blocks, Some(&seqid_to_genome), true)
            .with_context(|| "Failed to parse the blocks GFF")?;

    let dup_orig = duplicated_ids_by_genome(&genomes_original);
    remove_duplicates(&mut genomes_original, &dup_orig);
    eprintln!("duplicated (original): {}", dup_orig.len());

    let dup_markers = duplicated_ids_by_genome(&genomes_new);
    remove_duplicates(&mut genomes_new, &dup_markers);
    eprintln!("duplicated (markers): {}", dup_markers.len());

    let breakpoints: HashSet<(usize, usize)> = match breakpoints_file {
        Some(path) => {
            let bp = load_breakpoints(path).with_context(|| "Failed to read breakpoints file")?;
            eprintln!("Loaded {} breakpoints from file", bp.len());
            bp
        }
        None => {
            let seqs = genomes_to_sequences(&genomes_original);
            let bp = compute_breakpoints(&seqs);
            eprintln!("Computed {} breakpoints from original GFF", bp.len());
            bp
        }
    };

    let genomes_new_blocks = compute_genome_blocks(&genomes_original, &genomes_new);

    eprintln!(
        "Created compute_genome_blocks: {} genome(s), {} sequence(s)",
        genomes_new_blocks.len(),
        genomes_new_blocks.values().map(|g| g.len()).sum::<usize>(),
    );

    std::fs::create_dir_all(output_folder)
        .with_context(|| format!("Failed to create output directory: {}", output_folder))?;

    //FP
    let (false_positives, debug_info) = compute_false_positive_breakpoints(
        &genomes_new,
        &genomes_new_blocks,
        &breakpoints,
        debug,
        !ignore_new_breakpoints,
    );

    let fp_path = PathBuf::from(output_folder).join("false_positive.txt");
    eprintln!("Writing false positives to: {:?}", fp_path);
    let fp_file = File::create(&fp_path)
        .with_context(|| format!("Failed to create file: {:?}", fp_path))?;
    let mut fp_out = io::BufWriter::new(fp_file);

    if let Some(debug_info) = debug_info {
        print_false_positives_debug(&debug_info, &mut fp_out)?;
    } else {
        for (a, b) in false_positives {
            writeln!(fp_out, "{a} {b}")?;
        }       
    }

    //FN + TP
    let false_negatives = compute_false_negative_breakpoints(
        &genomes_new,
        &genomes_new_blocks,
        &breakpoints,
        fn_mode,
    );

    let fn_path = PathBuf::from(output_folder).join("false_negative.txt");
    let fn_file = File::create(&fn_path)
        .with_context(|| format!("Failed to create file: {:?}", fn_path))?;
    let mut fn_out = io::BufWriter::new(fn_file);
    for (a, b) in false_negatives.iter() {
        writeln!(fn_out, "{a} {b}")?;
    }       

    let tp_path = PathBuf::from(output_folder).join("true_positive.txt");
    let tp_file = File::create(&tp_path)
        .with_context(|| format!("Failed to create file: {:?}", tp_path))?;
    let mut tp_out = io::BufWriter::new(tp_file);
    for (a, b) in breakpoints {
        if !false_negatives.contains(&(a,b)) {
            writeln!(tp_out, "{a} {b}")?;
        }
    }       

    Ok(())
}

fn print_false_positives_debug<W: Write>(
    debug_info: &[FalsePositiveDebugInfo],
    out: &mut W,
) -> Result<()> {
    for info in debug_info {
        writeln!(out, "# False positive breakpoint between markers {} and {}", info.marker_x.0, info.marker_y.0)?;
        writeln!(out, "# Genomes: {} (seqid: {}) <-> {} (seqid: {})",
            info.genome_a, info.seqid_a, info.genome_b, info.seqid_b)?;
        write!(out, "# Block X: [")?;
        for (i, (id, strand)) in info.block_x.iter().enumerate() {
            if i > 0 { write!(out, ", ")?; }
            write!(out, "{}{}", id, strand)?;
        }
        writeln!(out, "]")?;
        write!(out, "# Block Y: [")?;
        for (i, (id, strand)) in info.block_y.iter().enumerate() {
            if i > 0 { write!(out, ", ")?; }
            write!(out, "{}{}", id, strand)?;
        }
        writeln!(out, "]")?;
        writeln!(out, "{} {}", info.marker_x.0, info.marker_y.0)?;
        writeln!(out)?;
    }
    Ok(())
}

fn compute_false_negative_breakpoints(
    genomes_new: &Genomes,
    genomes_new_blocks: &GenomesBlocks,
    breakpoints: &HashSet<(usize, usize)>,
    fn_mode: FnMode,
) -> HashSet<(usize, usize)> {

    match fn_mode {
        FnMode::Partition => {
            compute_false_negative_breakpoints_from_partition(genomes_new, genomes_new_blocks, breakpoints)
        }
        FnMode::Blocks => {
            compute_false_negative_breakpoints_from_blocks(genomes_new_blocks, breakpoints)
        }
    }
}

fn compute_false_negative_breakpoints_from_partition(
    genomes_new: &Genomes,
    genomes_new_blocks: &GenomesBlocks,
    breakpoints: &HashSet<(usize, usize)>,
) -> HashSet<(usize, usize)> {

    let partition = compute_partition(genomes_new, genomes_new_blocks);

    let combined_false_negatives = partition
        .into_par_iter()
        .map(|part| {
            let mut local_false_negatives = HashSet::new();
            let part_vec: Vec<&usize> = part.iter().collect();
            for i in 0..part_vec.len() {
                  for j in (i + 1)..part_vec.len() {
                      let a = *part_vec[i];
                      let b = *part_vec[j];
                      let pair = if a < b { (a, b) } else { (b, a) };
                      if breakpoints.contains(&pair) {
                          local_false_negatives.insert(pair);
                      }
                  }
              }
            local_false_negatives
        })
        .reduce(HashSet::new, |mut acc, h| {
            acc.extend(h);
            acc
        });

    combined_false_negatives
}

fn compute_false_negative_breakpoints_from_blocks(
    genomes_new_blocks: &GenomesBlocks,
    breakpoints: &HashSet<(usize, usize)>,
) -> HashSet<(usize, usize)> {
    let combined_false_negatives = genomes_new_blocks
        .par_iter()
        .map(|(_, seqids_blocks)| {
            let mut local_false_negatives = HashSet::new();
            for blocks in seqids_blocks.values() {
                for block in blocks {
                    for i in 0..block.len() {
                          for j in (i + 1)..block.len() {
                              let (a, _) = block[i];
                              let (b, _) = block[j];
                              let pair = if a < b { (a, b) } else { (b, a) };
                              if breakpoints.contains(&pair) {
                                  local_false_negatives.insert(pair);
                              }
                          }
                    }
                }
            }
            local_false_negatives
        })
        .reduce(HashSet::new, |mut acc, h| {
            acc.extend(h);
            acc
        });

    combined_false_negatives
}

fn compute_false_positive_breakpoints(
    genomes_new: &Genomes,
    genomes_new_blocks: &GenomesBlocks,
    breakpoints: &HashSet<(usize, usize)>,
    debug: bool,
    add_new_breakpoints: bool,
) -> (HashSet<(usize, usize)>, Option<Vec<FalsePositiveDebugInfo>>) {
    let genome_names: Vec<&String> = genomes_new.keys().collect();
    let n_genomes = genome_names.len();

    let pairs: Vec<(usize, usize)> = (0..n_genomes)
        .flat_map(|i| ((i + 1)..n_genomes).map(move |j| (i, j)))
        .collect();

    let (false_positives, debug_info): (Vec<_>, Vec<_>) = pairs
        .into_par_iter()
        .map(|(i, j)| {
            let mut local_false_positives: HashSet<(usize, usize)> = HashSet::new();
            let mut local_debug_info = if debug {
                Some(Vec::new())
            } else {
                None
            };

            let genome_name_a = genome_names[i];
            let genome_name_b = genome_names[j];

            process_genome_pair(
                genome_name_a, genome_name_b,
                genomes_new, genomes_new_blocks, breakpoints,
                &mut local_false_positives,
                local_debug_info.as_mut(),
                add_new_breakpoints,
            );

            (local_false_positives, local_debug_info)
        })
        .unzip();

    let combined_false_positives = false_positives.into_iter()
        .fold(HashSet::new(), |mut acc, h| {
            acc.extend(h);
            acc
        });

    let combined_debug_info = if debug {
        Some(debug_info.into_iter().flatten().flatten().collect())
    } else {
        None
    };

    (combined_false_positives, combined_debug_info)
}

fn process_genome_pair(
    genome_name_a: &str,
    genome_name_b: &str,
    genomes_new: &Genomes,
    genomes_new_blocks: &GenomesBlocks,
    breakpoints: &HashSet<(usize, usize)>,
    false_positives: &mut HashSet<(usize, usize)>,
    mut debug_info: Option<&mut Vec<FalsePositiveDebugInfo>>,
    add_new_breakpoints: bool,
) {
    let seqs_a = &genomes_new[genome_name_a];
    let seqs_b = &genomes_new[genome_name_b];

    let blocks_a = &genomes_new_blocks[genome_name_a];
    let blocks_b = &genomes_new_blocks[genome_name_b];

    for (seqid_a, seq_a) in seqs_a {
        let blocks_seq_a = &blocks_a[seqid_a];

        for (seqid_b, seq_b) in seqs_b {
            let blocks_seq_b = &blocks_b[seqid_b];

            let breakpoints_signed = breakpoints_from_pair_of_seq(&seq_a.markers, &seq_b.markers);

            let additional_breakpoints = if add_new_breakpoints {
                let mut original_seq_a: Vec<(usize, char)> = Vec::new();
                for block in blocks_seq_a.iter() {
                    original_seq_a.extend(block);
                }
                let mut original_seq_b: Vec<(usize, char)> = Vec::new();
                for block in blocks_seq_b.iter() {
                    original_seq_b.extend(block);
                }
                let breakpoints_signed = breakpoints_from_pair_of_seq(&original_seq_a, &original_seq_b);
                let mut extra = HashSet::new();
                for ((a, _), (b, _)) in breakpoints_signed {
                    let pair = if a < b { (a, b) } else { (b, a) };
                    extra.insert(pair);
                }
                Some(extra)
            } else {
                None
            };

            let ctx_a = DebugContext {
                genome_a: genome_name_a,
                seqid_a,
                genome_b: genome_name_b,
                seqid_b,
            };
            let ctx_b = DebugContext {
                genome_a: genome_name_b,
                seqid_a: seqid_b,
                genome_b: genome_name_a,
                seqid_b: seqid_a,
            };

            let debug_a = debug_info.as_mut().map(|vec| (vec as &mut _, &ctx_a));
            collect_false_positives_with_breakpoints(
                &seq_a.markers, blocks_seq_a,
                &breakpoints_signed, breakpoints,
                additional_breakpoints.as_ref(),
                false_positives,
                debug_a,
            );

            let debug_b = debug_info.as_mut().map(|vec| (vec as &mut _, &ctx_b));
            collect_false_positives_with_breakpoints(
                &seq_b.markers, blocks_seq_b,
                &breakpoints_signed, breakpoints,
                additional_breakpoints.as_ref(),
                false_positives,
                debug_b,
            );
        }
    }
}

fn collect_false_positives_with_breakpoints(
    markers: &[(usize, char)],
    blocks: &Blocks,
    breakpoints_signed: &HashSet<((usize, char), (usize, char))>,
    breakpoints: &HashSet<(usize, usize)>,
    additional_breakpoints: Option<&HashSet<(usize, usize)>>,
    false_positives: &mut HashSet<(usize, usize)>,
    mut debug: Option<(&mut Vec<FalsePositiveDebugInfo>, &DebugContext)>,
) {
    for (i, window) in markers.windows(2).enumerate() {
        let marker_x = window[0];
        let marker_y = window[1];
        if breakpoints_signed.contains(&canonical_signed_adjacency((marker_x, marker_y))) {
            let block_x = &blocks[i];
            let block_y = &blocks[i + 1];
            let mut is_false_positive = true;
            'outer: for &(x, _) in block_x {
                for &(y, _) in block_y {
                    let pair = if x < y { (x, y) } else { (y, x) };
                    if breakpoints.contains(&pair)
                        || additional_breakpoints.map_or(false, |extra| extra.contains(&pair)) {
                        is_false_positive = false;
                        break 'outer;
                    }
                }
            }
            if is_false_positive {
                let pair = if marker_x.0 < marker_y.0 { (marker_x.0, marker_y.0) } else { (marker_y.0, marker_x.0) };
                false_positives.insert(pair);

                if let Some((info, ctx)) = debug.as_mut() {
                    info.push(FalsePositiveDebugInfo {
                        genome_a: ctx.genome_a.to_string(),
                        seqid_a: ctx.seqid_a.to_string(),
                        genome_b: ctx.genome_b.to_string(),
                        seqid_b: ctx.seqid_b.to_string(),
                        marker_x: marker_x,
                        marker_y: marker_y,
                        block_x: block_x.clone(),
                        block_y: block_y.clone(),
                    });
                }
            }
        }
    }
}

fn flip_element(t: (usize, char)) -> (usize, char) {
    let (id, s) = t;
    let s2 = match s {
        '+' => '-',
        '-' => '+',
        _ => s,
    };
    (id, s2)
}

fn canonical_signed_adjacency(
    adj: ((usize, char), (usize, char)),
) -> ((usize, char), (usize, char)) {
    let (x, y) = adj;
    let rf = (flip_element(y), flip_element(x));
    if (x, y) <= rf {
        (x, y)
    } else {
        rf
    }
}

fn canonical_adjacencies(path: &[(usize, char)]) -> HashSet<((usize, char), (usize, char))> {
    let mut h = HashSet::with_capacity(path.len());
    for i in 0..path.len().saturating_sub(1) {
        h.insert(canonical_signed_adjacency((path[i], path[i + 1])));
    }
    h
}

fn project_common(seq: &[(usize, char)], common: &HashSet<usize>) -> Vec<(usize, char)> {
    seq.iter()
        .copied()
        .filter(|(id, _)| common.contains(id))
        .collect()
}

fn breakpoints_from_pair_of_seq(seq_a: &[(usize, char)], seq_b: &[(usize, char)]) -> HashSet<((usize,char),(usize, char))> {

    let ids_a: HashSet<usize> = seq_a.iter().map(|(id, _)| *id).collect();
    let ids_b: HashSet<usize> = seq_b.iter().map(|(id, _)| *id).collect();

    let a_proj = project_common(seq_a, &ids_b);
    if a_proj.len() < 2 {
        return HashSet::new()
    }
    let b_proj = project_common(seq_b, &ids_a);
    if b_proj.len() < 2 {
        return HashSet::new()
    }

    let p1_canon = canonical_adjacencies(&a_proj);
    let p2_canon = canonical_adjacencies(&b_proj);

    let breakpoints_signed: HashSet<((usize, char), (usize, char))> = p1_canon
        .symmetric_difference(&p2_canon)
        .cloned()
        .collect();

    breakpoints_signed
}

fn compute_breakpoints(seqs: &[&Seq]) -> HashSet<(usize, usize)> {
    let n = seqs.len();

    (0..n)
        .into_par_iter()
        .map(|i| {
            let mut local_breakpoints: HashSet<(usize, usize)> = HashSet::new();

            for j in (i + 1)..n {
                let seq_a = seqs[i];
                let seq_b = seqs[j];

                let breakpoints_signed = breakpoints_from_pair_of_seq(&seq_a.markers, &seq_b.markers);

                for ((a, _), (b, _)) in breakpoints_signed {
                    let pair = if a < b { (a, b) } else { (b, a) };
                    local_breakpoints.insert(pair);
                }
            }

            local_breakpoints
        })
        .reduce(HashSet::new, |mut acc, h| {
            acc.extend(h);
            acc
        })
}

fn load_breakpoints(path: &Path) -> Result<HashSet<(usize, usize)>> {
    let fh = File::open(path).with_context(|| format!("Opening {path:?}"))?;
    let br = BufReader::new(fh);
    let mut breakpoints: HashSet<(usize, usize)> = HashSet::new();

    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {path:?}, line {i}"))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            return Err(anyhow!("Malformed breakpoints line (need 2 IDs): {line}"));
        }
        let a: usize = parts[0]
            .parse()
            .with_context(|| format!("Invalid ID in breakpoints line: {line}"))?;
        let b: usize = parts[1]
            .parse()
            .with_context(|| format!("Invalid ID in breakpoints line: {line}"))?;
        let pair = if a < b { (a, b) } else { (b, a) };
        breakpoints.insert(pair);
    }

    Ok(breakpoints)
}
