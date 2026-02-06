use anyhow::{anyhow, Context, Result};
use clap::Parser;
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
    marker_a: (usize, char),
    marker_b: (usize, char),
    block_a: Vec<(usize, char)>,
    block_b: Vec<(usize, char)>,
}

struct DebugContext<'a> {
    genome_a: &'a str,
    seqid_a: &'a str,
    genome_b: &'a str,
    seqid_b: &'a str,
}

#[derive(Parser, Debug)]
#[command(author, version, about, arg_required_else_help = true)]
struct Args {
    /// GFF file(s). One file: compute breakpoints. Two files: compare original vs new.
    #[arg(num_args = 1..=2)]
    file_gff: Vec<String>,
    #[arg(long)]
    no_dup: bool,
    #[arg(long)]
    seqid2genome: Option<PathBuf>,
    #[arg(long)]
    breakpoints: Option<PathBuf>,
    #[arg(long)]
    debug: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let num_cores = num_cpus::get();
    ThreadPoolBuilder::new()
        .num_threads(num_cores)
        .build_global()
        .unwrap();

    let seqid_to_genome = match args.seqid2genome.as_ref() {
        Some(path) => Some(load_seqid2genome(path).with_context(|| "Failed to read seqid2genome")?),
        None => None,
    };

    if args.file_gff.len() == 1 {
        run_single_gff(&args.file_gff[0], seqid_to_genome.as_ref(), args.no_dup)
    } else {
        if !args.no_dup {
            return Err(anyhow!("--no-dup is required when comparing two GFF files"));
        }
        run_two_gffs(
            &args.file_gff[0],
            &args.file_gff[1],
            seqid_to_genome.as_ref(),
            args.breakpoints.as_ref(),
            args.debug,
        )
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

fn run_two_gffs(
    file_gff_original: &str,
    file_gff_new: &str,
    seqid_to_genome: Option<&HashMap<String, String>>,
    breakpoints_file: Option<&PathBuf>,
    debug: bool,
) -> Result<()> {
    let (mut genomes_original, seqid_to_genome) =
        load_genomes(file_gff_original, seqid_to_genome, true)
            .with_context(|| "Failed to parse the original GFF")?;

    let (mut genomes_new, _seqid_to_genome) =
        load_genomes(file_gff_new, Some(&seqid_to_genome), true)
            .with_context(|| "Failed to parse the new GFF")?;

    let dup_orig = duplicated_ids_by_genome(&genomes_original);
    remove_duplicates(&mut genomes_original, &dup_orig);
    eprintln!("duplicated (original): {}", dup_orig.len());

    eprintln!(
        "Loaded original: {} genome(s), {} sequence(s)",
        genomes_original.len(),
        genomes_original.values().map(|g| g.len()).sum::<usize>(),
    );

    let dup_new = duplicated_ids_by_genome(&genomes_new);
    remove_duplicates(&mut genomes_new, &dup_new);
    eprintln!("duplicated (new): {}", dup_new.len());

    eprintln!(
        "Loaded new: {} genome(s), {} sequence(s)",
        genomes_new.len(),
        genomes_new.values().map(|g| g.len()).sum::<usize>(),
    );

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

    //print_genomes_new_blocks(&genomes_new_blocks)?;

    let (false_positives, debug_info) = compute_false_positive_breakpoints(&genomes_new, &genomes_new_blocks, &breakpoints, debug);

    let mut out = io::BufWriter::new(io::stdout().lock());

    if let Some(debug_info) = debug_info {
        print_false_positives_debug(&debug_info, &mut out)?;
    } else {
        print_false_positives(&false_positives, &mut out)?;
    }

    Ok(())
}

fn print_false_positives(
    false_positives: &HashSet<(usize, usize)>,
    out: &mut io::BufWriter<io::StdoutLock>,
) -> Result<()> {
    for (a, b) in false_positives {
        writeln!(out, "{a} {b}")?;
    }
    Ok(())
}

fn print_false_positives_debug(
    debug_info: &[FalsePositiveDebugInfo],
    out: &mut io::BufWriter<io::StdoutLock>,
) -> Result<()> {
    for info in debug_info {
        writeln!(out, "# False positive breakpoint between markers {} and {}", info.marker_a.0, info.marker_b.0)?;
        writeln!(out, "# Genomes: {} (seqid: {}) <-> {} (seqid: {})",
            info.genome_a, info.seqid_a, info.genome_b, info.seqid_b)?;
        write!(out, "# Block A: [")?;
        for (i, (id, strand)) in info.block_a.iter().enumerate() {
            if i > 0 { write!(out, ", ")?; }
            write!(out, "{}{}", id, strand)?;
        }
        writeln!(out, "]")?;
        write!(out, "# Block B: [")?;
        for (i, (id, strand)) in info.block_b.iter().enumerate() {
            if i > 0 { write!(out, ", ")?; }
            write!(out, "{}{}", id, strand)?;
        }
        writeln!(out, "]")?;
        writeln!(out, "{} {}", info.marker_a.0, info.marker_b.0)?;
        writeln!(out)?;
    }
    Ok(())
}

fn compute_false_positive_breakpoints(
    genomes_new: &Genomes,
    genomes_new_blocks: &GenomesBlocks,
    breakpoints: &HashSet<(usize, usize)>,
    debug: bool,
) -> (HashSet<(usize, usize)>, Option<Vec<FalsePositiveDebugInfo>>) {
    let genome_names: Vec<&String> = genomes_new.keys().collect();
    let n_genomes = genome_names.len();

    let pairs: Vec<(usize, usize)> = (0..n_genomes)
        .flat_map(|i| ((i + 1)..n_genomes).map(move |j| (i, j)))
        .collect();

    if debug {
        let (false_positives, debug_info): (Vec<_>, Vec<_>) = pairs
            .into_par_iter()
            .map(|(i, j)| {
                let mut local_false_positives: HashSet<(usize, usize)> = HashSet::new();
                let mut local_debug_info: Vec<FalsePositiveDebugInfo> = Vec::new();

                let genome_name_a = genome_names[i];
                let genome_name_b = genome_names[j];

                process_genome_pair(
                    genome_name_a, genome_name_b,
                    genomes_new, genomes_new_blocks, breakpoints,
                    &mut local_false_positives, Some(&mut local_debug_info),
                );

                (local_false_positives, local_debug_info)
            })
            .unzip();

        let combined_false_positives = false_positives.into_iter()
            .fold(HashSet::new(), |mut acc, h| {
                acc.extend(h);
                acc
            });

        let combined_debug_info = debug_info.into_iter()
            .flatten()
            .collect();

        (combined_false_positives, Some(combined_debug_info))
    } else {
        let false_positives = pairs
            .into_par_iter()
            .map(|(i, j)| {
                let mut local_false_positives: HashSet<(usize, usize)> = HashSet::new();

                let genome_name_a = genome_names[i];
                let genome_name_b = genome_names[j];

                process_genome_pair(
                    genome_name_a, genome_name_b,
                    genomes_new, genomes_new_blocks, breakpoints,
                    &mut local_false_positives, None,
                );

                local_false_positives
            })
            .reduce(HashSet::new, |mut acc, h| {
                acc.extend(h);
                acc
            });

        (false_positives, None)
    }
}

fn process_genome_pair(
    genome_name_a: &str,
    genome_name_b: &str,
    genomes_new: &Genomes,
    genomes_new_blocks: &GenomesBlocks,
    breakpoints: &HashSet<(usize, usize)>,
    false_positives: &mut HashSet<(usize, usize)>,
    mut debug_info: Option<&mut Vec<FalsePositiveDebugInfo>>,
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
            collect_false_positives(
                &seq_a.markers, blocks_seq_a,
                &breakpoints_signed, breakpoints,
                false_positives,
                debug_a,
            );

            let debug_b = debug_info.as_mut().map(|vec| (vec as &mut _, &ctx_b));
            collect_false_positives(
                &seq_b.markers, blocks_seq_b,
                &breakpoints_signed, breakpoints,
                false_positives,
                debug_b,
            );
        }
    }
}

fn collect_false_positives(
    markers: &[(usize, char)],
    blocks: &Blocks,
    breakpoints_signed: &HashSet<((usize, char), (usize, char))>,
    breakpoints: &HashSet<(usize, usize)>,
    false_positives: &mut HashSet<(usize, usize)>,
    mut debug: Option<(&mut Vec<FalsePositiveDebugInfo>, &DebugContext)>,
) {
    for (i, window) in markers.windows(2).enumerate() {
        let a = window[0];
        let b = window[1];
        if breakpoints_signed.contains(&canonical_signed_adjacency((a, b))) {
            let block_x = &blocks[i];
            let block_y = &blocks[i + 1];
            let mut is_false_positive = true;
            'outer: for &(x, _) in block_x {
                for &(y, _) in block_y {
                    let pair = if x < y { (x, y) } else { (y, x) };
                    if breakpoints.contains(&pair) {
                        is_false_positive = false;
                        break 'outer;
                    }
                }
            }
            if is_false_positive {
                let pair = if a.0 < b.0 { (a.0, b.0) } else { (b.0, a.0) };
                false_positives.insert(pair);

                if let Some((info, ctx)) = debug.as_mut() {
                    info.push(FalsePositiveDebugInfo {
                        genome_a: ctx.genome_a.to_string(),
                        seqid_a: ctx.seqid_a.to_string(),
                        genome_b: ctx.genome_b.to_string(),
                        seqid_b: ctx.seqid_b.to_string(),
                        marker_a: a,
                        marker_b: b,
                        block_a: block_x.clone(),
                        block_b: block_y.clone(),
                    });
                }
            }
        }
    }
}

#[allow(dead_code)]
fn print_genomes_new_blocks(genomes_new_blocks: &GenomesBlocks) -> Result<()>{
    let mut out = io::BufWriter::new(io::stdout().lock());
    let mut genome_names: Vec<String> = genomes_new_blocks.keys().cloned().collect();
    genome_names.sort();
    for genome in genome_names {
        if let Some(seqs) = genomes_new_blocks.get(&genome) {
            for (seqid, blocks) in seqs.iter() {
                write!(out, "{genome}\t{seqid}\t")?;
                for block in blocks {
                    let mut first = true;
                    write!(out, "[")?;
                    for (id, strand) in block.iter(){
                        if !first {
                            out.write(b",")?;
                        }
                        first = false;
                        write!(out, "{id}{strand}")?;
                    }
                    write!(out, "]")?;
                }
                out.write(b"\n")?;
            }
        }
    }
    Ok(())
}

fn compute_genome_blocks(
    genomes_original: &Genomes,
    genomes_new: &Genomes,
) -> GenomesBlocks {
    genomes_new
        .par_iter()
        .map(|(genome_name_new, seqid_table)| {
            let genome_ori = &genomes_original[genome_name_new];
            let seqid_blocks: HashMap<String, Vec<Vec<(usize, char)>>> = seqid_table
                .par_iter()
                .map(|(seqid_new, seq_new)| {
                    let seq_ori = &genome_ori[seqid_new];
                    let mut blocks = Vec::with_capacity(seq_new.len());
                    let mut j = 0; //ori
                    for i in 0..seq_new.len() {
                        let mut contained = Vec::new();
                        while j < seq_ori.len() && seq_new.starts[i] > seq_ori.starts[j] {
                            j += 1;
                        }
                        while j < seq_ori.len() && seq_new.ends[i] >= seq_ori.ends[j] {
                            contained.push(seq_ori.markers[j]);
                            j += 1;
                        }
                        blocks.push(contained);
                    }
                    (seqid_new.clone(), blocks)
                })
                .collect();
            (genome_name_new.clone(), seqid_blocks)
        })
        .collect()
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

    // Project sequences to common IDs
    let a_proj = project_common(seq_a, &ids_b);
    if a_proj.len() < 2 {
        return HashSet::new()
    }
    let b_proj = project_common(seq_b, &ids_a);
    if b_proj.len() < 2 {
        return HashSet::new()
    }

    // Canonical adjacencies
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

    // Parallel over i; each thread builds a local HashSet and we reduce at the end.
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
