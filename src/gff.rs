use anyhow::{anyhow, Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use rayon::prelude::*;

#[derive(Debug, Clone, Default)]
pub struct Seq {
    pub markers: Vec<(usize, char)>,
    pub starts: Vec<u64>,
    pub ends: Vec<u64>,
}

impl Seq {
    pub fn len(&self) -> usize {
        self.markers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.markers.is_empty()
    }

    pub fn push(&mut self, id: usize, strand: char, start: u64, end: u64) {
        self.markers.push((id, strand));
        self.starts.push(start);
        self.ends.push(end);
    }

    pub fn sort_by_start(&mut self) {
        let mut indices: Vec<usize> = (0..self.len()).collect();
        indices.sort_by(|&a, &b| self.starts[a].cmp(&self.starts[b]));
        self.reorder(&indices);
    }

    fn reorder(&mut self, order: &[usize]) {
        let markers: Vec<(usize, char)> = order.iter().map(|&i| self.markers[i]).collect();
        let starts: Vec<u64> = order.iter().map(|&i| self.starts[i]).collect();
        let ends: Vec<u64> = order.iter().map(|&i| self.ends[i]).collect();
        self.markers = markers;
        self.starts = starts;
        self.ends = ends;
    }

    /// Retain only entries where `keep[i]` is true.
    pub fn retain_by_mask(&mut self, keep: &[bool]) {
        let mut j = 0;
        for i in 0..self.len() {
            if keep[i] {
                self.markers[j] = self.markers[i];
                self.starts[j] = self.starts[i];
                self.ends[j] = self.ends[i];
                j += 1;
            }
        }
        self.markers.truncate(j);
        self.starts.truncate(j);
        self.ends.truncate(j);
    }
}

pub type Blocks = Vec<Vec<(usize,char)>>;
pub type Genomes = HashMap<String, HashMap<String, Seq>>;
pub type GenomesBlocks = HashMap<String, HashMap<String, Blocks>>;

#[derive(Debug, Clone)]
pub struct Record {
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub id: usize,          // ID= (1-based)
    pub strand: char,
    pub genome: Option<String>, // genome=
}

/// Parse GFF3 line into Option<Record> filtered to type == "SO:0000856".
fn parse_gff_line(line: &str) -> Result<Option<Record>> {
    if line.is_empty() || line.starts_with('#') {
        return Ok(None);
    }
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 9 {
        return Err(anyhow!("Malformed GFF line (need 9 columns): {line}"));
    }

    let seqid = cols[0].to_string();
    let ftype = cols[2];
    if ftype != "SO:0000856" {
        return Ok(None);
    }

    // GFF is 1-based inclusive
    let start: u64 = cols[3]
        .parse()
        .with_context(|| format!("Invalid start in line: {line}"))?;
    let end: u64 = cols[4]
        .parse()
        .with_context(|| format!("Invalid end in line: {line}"))?;

    let strand_field = cols[6];
    let strand = match strand_field.chars().next() {
        Some('+' | '-') => strand_field.chars().next().unwrap(),
        _ => return Err(anyhow!("Invalid strand in line: {line}")),
    };

    let attributes = cols[8];
    let mut id_opt: Option<usize> = None;
    let mut genome: Option<String> = None;

    for field in attributes.split(';') {
        if let Some(rest) = field.strip_prefix("ID=") {
            let num: usize = rest
                .parse()
                .with_context(|| format!("ID is not an integer in line: {line}"))?;
            id_opt = Some(num);
        } else if let Some(rest) = field.strip_prefix("genome=") {
            genome = Some(rest.to_string());
        }
    }
    let id = id_opt.ok_or_else(|| anyhow!("Missing ID= in attributes: {line}"))?;

    Ok(Some(Record {
        seqid,
        start,
        end,
        id,
        strand,
        genome,
    }))
}

pub fn load_genomes(
    path: &str,
    seqid_to_genome: Option<&HashMap<String, String>>,
    sort_by_seqid_start: bool,
) -> Result<(Genomes, HashMap<String, String>)> {
    let fh = File::open(path).with_context(|| format!("Opening {path}"))?;
    let br = BufReader::new(fh);
    let mut genomes: Genomes = HashMap::new();
    let mut built_map: HashMap<String, String> = match seqid_to_genome {
        Some(m) => m.clone(),
        None => HashMap::new(),
    };
    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {path}, line {i}"))?;
        if let Some(mut r) = parse_gff_line(&line)? {
            if r.genome.is_none() {
                if let Some(map) = seqid_to_genome {
                    if let Some(genome) = map.get(&r.seqid) {
                        r.genome = Some(genome.clone());
                    } else {
                        return Err(anyhow!(
                            "Missing genome= in GFF and no mapping for seqid {}",
                            r.seqid
                        ));
                    }
                } else {
                    return Err(anyhow!(
                        "Missing genome= in GFF for seqid {}. Please provide seqid2genome.",
                        r.seqid
                    ));
                }
            } else if let Some(map) = seqid_to_genome {
                if let Some(mapped) = map.get(&r.seqid) {
                    if r.genome.as_deref() != Some(mapped.as_str()) {
                        return Err(anyhow!(
                            "Genome mismatch for seqid {}: GFF has {:?}, mapping has {}",
                            r.seqid,
                            r.genome,
                            mapped
                        ));
                    }
                }
            }

            let genome = r
                .genome
                .as_ref()
                .expect("genome filled by mapping or present in GFF");

            if seqid_to_genome.is_none() {
                built_map
                    .entry(r.seqid.clone())
                    .or_insert_with(|| genome.clone());
            }

            genomes
                .entry(genome.clone())
                .or_default()
                .entry(r.seqid.clone())
                .or_default()
                .push(r.id, r.strand, r.start, r.end);
        }
    }

    if sort_by_seqid_start {
        for seqs in genomes.values_mut() {
            for seq in seqs.values_mut() {
                seq.sort_by_start();
            }
        }
    }
    Ok((genomes, built_map))
}

pub fn genomes_to_sequences(genomes: &Genomes) -> Vec<&Seq> {
    let mut new_seqs = Vec::new();

    for seqs in genomes.values() {
        for seq in seqs.values() {
            new_seqs.push(seq);
        }
    }
    new_seqs
}

//pub fn genomes_blocks_to_sequence_blocks(genomes: &GenomesBlocks) -> Vec<&Vec<Vec<(usize, char)>>> {
//    let mut new_seqs = Vec::new();
//
//    for seqs in genomes.values() {
//        for seq in seqs.values() {
//            new_seqs.push(seq);
//        }
//    }
//    new_seqs
//}

pub fn load_seqid2genome(path: &Path) -> Result<HashMap<String, String>> {
    let fh = File::open(path).with_context(|| format!("Opening {path:?}"))?;
    let br = BufReader::new(fh);
    let mut map: HashMap<String, String> = HashMap::new();

    for (i, line) in br.lines().enumerate() {
        let line = line.with_context(|| format!("Reading {path:?}, line {i}"))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (genome, rest) = match line.split_once(':') {
            Some((g, r)) => (g.trim(), r.trim()),
            None => {
                return Err(anyhow!(
                    "Malformed seqid2genome line (missing ':'): {line}"
                ))
            }
        };
        if genome.is_empty() {
            return Err(anyhow!("Malformed seqid2genome line (empty genome): {line}"));
        }
        for token in rest.split(|c: char| c == ',' || c.is_whitespace()) {
            let seqid = token.trim();
            if seqid.is_empty() {
                continue;
            }
            if let Some(existing) = map.get(seqid) {
                if existing != genome {
                    return Err(anyhow!(
                        "Seqid {seqid} is assigned to multiple genomes: {existing} and {genome}"
                    ));
                }
            } else {
                map.insert(seqid.to_string(), genome.to_string());
            }
        }
    }

    Ok(map)
}

pub fn remove_duplicates(genomes: &mut Genomes, duplicated: &HashSet<usize>) {
    for seqs in genomes.values_mut() {
        for seq in seqs.values_mut() {
            let keep: Vec<bool> = seq.markers.iter().map(|(id, _)| !duplicated.contains(id)).collect();
            seq.retain_by_mask(&keep);
        }
    }
}

pub fn duplicated_ids_by_genome(
    genomes: &Genomes,
) -> HashSet<usize> {
    let mut dup_ids: HashSet<usize> = HashSet::new();
    for seqs in genomes.values() {
        let mut seen: HashSet<usize> = HashSet::new();
        for seq in seqs.values() {
            for &(id, _) in &seq.markers {
                if !seen.insert(id) {
                    dup_ids.insert(id);
                }
            }
        }
    }
    dup_ids
}

pub fn compute_genome_blocks(
    genomes_original: &Genomes,
    genomes_new: &Genomes,
    extend: u64,
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
                        let extended_end = seq_new.ends[i] + extend;
                        while j < seq_ori.len() && extended_end >= seq_ori.ends[j] {
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

pub fn compute_partition(genomes_new: &Genomes, genomes_new_blocks: &GenomesBlocks) 
    -> Vec<HashSet<usize>> {
    let mut max_n = 0;
    for seqs in genomes_new.values() {
        for seq in seqs.values() {
            for (id, _) in seq.markers.iter() {
                if *id > max_n {
                    max_n = *id;
                }
            }
        }
    }
    max_n += 1;
    let mut partition: Vec<HashSet<usize>> = vec![HashSet::new(); max_n];

    for (genome, seqs) in genomes_new.iter() {
        let seqids_blocks = &genomes_new_blocks[genome];
        for (seqid, seq) in seqs.iter() {
            let blocks = &seqids_blocks[seqid];
            for (i, (marker, _)) in seq.markers.iter().enumerate() {
                let block = &blocks[i];
                for (element_id, _) in block.iter() {
                    partition[*marker].insert(*element_id);
                }
            }
        }
    }
    partition
}
