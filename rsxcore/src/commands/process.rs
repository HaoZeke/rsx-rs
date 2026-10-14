// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `process` command: create marker depth table from demultiplexed reads.
//!
//! Each file is counted on its own, then the per-individual maps are merged by
//! hash shard (no per-key lock). Sequences are stored as 2-bit packed DNA.

use crate::io::seq_reader::{
    PackedDnaKey, count_sequences_packed, get_input_files, unpack_2bit_into,
};
use std::io::Write;
use std::path::Path;

/// Parameters for the `process` command.
pub struct ProcessParams {
    pub input_dir_path: String,
    pub output_file_path: String,
    pub n_threads: u32,
    pub min_depth: u16,
    /// If set, group markers by min-hash of canonical k-mers of this size.
    /// Heuristic (not exact) collapse of sequencing error variants.
    /// Optional (default: disabled). See kmer.rs docs for limitations.
    pub kmer_dedup: Option<usize>,
}

#[cfg(feature = "parallel")]
const PARALLEL_MERGE_MIN_CAPACITY: usize = 1_024;
#[cfg(feature = "parallel")]
const PARALLEL_MERGE_MAX_CAPACITY: usize = 4_000_000;
// Tuned for typical RAD tag length plus the packed key and table-entry overhead.
#[cfg(feature = "parallel")]
const ESTIMATED_BYTES_PER_MARKER: u64 = 64;

#[cfg(feature = "parallel")]
fn parallel_merge_capacity_from_bytes(input_bytes: u64) -> usize {
    let estimated = input_bytes / ESTIMATED_BYTES_PER_MARKER;
    estimated.clamp(
        PARALLEL_MERGE_MIN_CAPACITY as u64,
        PARALLEL_MERGE_MAX_CAPACITY as u64,
    ) as usize
}

#[cfg(feature = "parallel")]
fn estimate_parallel_merge_capacity(input_files: &[crate::io::seq_reader::InputFile]) -> usize {
    let input_bytes = input_files
        .iter()
        .filter_map(|f| std::fs::metadata(&f.path).ok())
        .fold(0_u64, |acc, meta| acc.saturating_add(meta.len()));
    parallel_merge_capacity_from_bytes(input_bytes)
}

/// Run the `process` command.
pub fn run(params: &ProcessParams) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("RADSex process started");
    let start = std::time::Instant::now();

    let input_dir = Path::new(&params.input_dir_path);
    let input_files = get_input_files(input_dir)?;
    let n_individuals = input_files.len();

    #[cfg(feature = "parallel")]
    let mut global = {
        use rayon::prelude::*;

        rayon::ThreadPoolBuilder::new()
            .num_threads(params.n_threads as usize)
            .build_global()
            .ok();

        let per_file: Vec<ahash::AHashMap<PackedDnaKey, u16>> = input_files
            .par_iter()
            .map(|f| {
                let counts = count_sequences_packed(&f.path).map_err(|e| {
                    format!("process failed while reading {}: {e}", f.path.display())
                })?;
                log::debug!("Finished processing individual {}", f.individual_name);
                Ok(counts)
            })
            .collect::<Result<Vec<_>, String>>()?;

        if n_individuals >= 8 {
            let capacity = estimate_parallel_merge_capacity(&input_files);
            merge_counts_sharded(per_file, n_individuals, capacity)
        } else {
            merge_counts_serial(per_file, n_individuals)
        }
    };

    #[cfg(not(feature = "parallel"))]
    let mut global = {
        let mut per_file = Vec::with_capacity(n_individuals);
        for f in &input_files {
            match count_sequences_packed(&f.path) {
                Ok(counts) => {
                    log::info!("Finished processing individual {}", f.individual_name);
                    per_file.push(counts);
                }
                Err(e) => {
                    return Err(format!("Error processing {}: {e}", f.path.display()).into());
                }
            }
        }
        merge_counts_serial(per_file, n_individuals)
    };

    // Optional k-mer deduplication: group markers by canonical k-mer
    if let Some(k) = params.kmer_dedup {
        log::info!(
            "K-mer deduplication (k={}): grouping {} markers",
            k,
            global.len()
        );
        let keys: Vec<PackedDnaKey> = global.keys().cloned().collect();
        let sequences: Vec<Vec<u8>> = keys
            .iter()
            .map(|k| crate::io::seq_reader::unpack_2bit(k.as_slice()))
            .collect();
        let groups = crate::kmer::group_by_kmer(&sequences, k);
        let n_before = global.len();
        let n_groups = groups.len();
        log::info!(
            "K-mer dedup: {} markers -> {} groups ({:.1}% reduction)",
            n_before,
            n_groups,
            (1.0 - n_groups as f64 / n_before as f64) * 100.0
        );
        // For each group, keep the representative with highest total depth
        let mut keep: ahash::AHashSet<usize> = ahash::AHashSet::new();
        for (_hash, indices) in &groups {
            let best = indices.iter().copied().max_by_key(|&i| {
                global
                    .get(&keys[i])
                    .map(|d| d.iter().map(|&v| v as u64).sum::<u64>())
                    .unwrap_or(0)
            });
            if let Some(b) = best {
                keep.insert(b);
            }
        }
        let kept_keys: ahash::AHashSet<PackedDnaKey> =
            keep.iter().map(|&i| keys[i].clone()).collect();
        global.retain(|k, _| kept_keys.contains(k));
        log::info!("K-mer dedup: retained {} markers", global.len());
    }

    // Retain only markers that pass min_depth *before* writing the header so
    // `#Number of markers` matches the number of data rows (RADSex contract).
    let mut retained: Vec<(PackedDnaKey, Vec<u16>)> = global
        .into_iter()
        .filter(|(_, depths)| {
            params.min_depth <= 1 || depths.iter().any(|&d| d >= params.min_depth)
        })
        .collect();
    // Deterministic row order by packed sequence key. Keys are unique, so the
    // parallel sort agrees with a stable sort.
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        retained.par_sort_by(|(a, _), (b, _)| a.as_slice().cmp(b.as_slice()));
    }
    #[cfg(not(feature = "parallel"))]
    retained.sort_by(|(a, _), (b, _)| a.as_slice().cmp(b.as_slice()));

    // Write output (unpack 2-bit keys back to ASCII for TSV)
    let file = std::fs::File::create(&params.output_file_path)?;
    let mut output = std::io::BufWriter::with_capacity(1 << 20, file);
    writeln!(output, "#Number of markers : {}", retained.len())?;

    write!(output, "id\tsequence")?;
    for f in &input_files {
        write!(output, "\t{}", f.individual_name)?;
    }
    writeln!(output)?;

    log::info!(
        "Writing {} retained marker depths to output file",
        retained.len()
    );

    let mut sequence = Vec::with_capacity(256);
    let mut line = Vec::with_capacity(1024);
    for (id, (packed_seq, depths)) in (0_u64..).zip(retained.iter()) {
        unpack_2bit_into(packed_seq.as_slice(), &mut sequence);
        line.clear();
        push_u64(&mut line, id);
        line.push(b'\t');
        line.extend_from_slice(&sequence);
        for &depth in depths {
            line.push(b'\t');
            push_u16(&mut line, depth);
        }
        line.push(b'\n');
        output.write_all(&line)?;
    }

    let elapsed = start.elapsed();
    let secs = elapsed.as_secs();
    log::info!(
        "RADSex process ended (total runtime: {}h {}m {}s)",
        secs / 3600,
        (secs % 3600) / 60,
        secs % 60
    );

    Ok(())
}

fn merge_counts_serial(
    maps: Vec<ahash::AHashMap<PackedDnaKey, u16>>,
    n_individuals: usize,
) -> ahash::AHashMap<PackedDnaKey, Vec<u16>> {
    let mut global: ahash::AHashMap<PackedDnaKey, Vec<u16>> = ahash::AHashMap::new();
    for (index, counts) in maps.into_iter().enumerate() {
        for (packed_seq, count) in counts {
            let entry = global
                .entry(packed_seq)
                .or_insert_with(|| vec![0_u16; n_individuals]);
            entry[index] = count;
        }
    }
    global
}

#[cfg(feature = "parallel")]
fn merge_counts_sharded(
    maps: Vec<ahash::AHashMap<PackedDnaKey, u16>>,
    n_individuals: usize,
    capacity: usize,
) -> ahash::AHashMap<PackedDnaKey, Vec<u16>> {
    use rayon::prelude::*;
    use std::hash::{Hash, Hasher};

    let shards = rayon::current_num_threads().max(1);
    let mut buckets: Vec<Vec<(PackedDnaKey, u16, u32)>> = (0..shards).map(|_| Vec::new()).collect();
    for (index, counts) in maps.into_iter().enumerate() {
        for (key, count) in counts {
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            key.hash(&mut hasher);
            let shard = (hasher.finish() as usize) % shards;
            buckets[shard].push((key, count, index as u32));
        }
    }

    let parts: Vec<ahash::AHashMap<PackedDnaKey, Vec<u16>>> = buckets
        .into_par_iter()
        .map(|bucket| {
            let mut local: ahash::AHashMap<PackedDnaKey, Vec<u16>> =
                ahash::AHashMap::with_capacity(bucket.len());
            for (key, count, index) in bucket {
                let entry = local
                    .entry(key)
                    .or_insert_with(|| vec![0_u16; n_individuals]);
                entry[index as usize] = count;
            }
            local
        })
        .collect();

    let mut global: ahash::AHashMap<PackedDnaKey, Vec<u16>> =
        ahash::AHashMap::with_capacity(capacity);
    for part in parts {
        global.extend(part);
    }
    global
}

fn push_u16(buf: &mut Vec<u8>, mut value: u16) {
    let mut tmp = [0_u8; 5];
    let mut index = 5;
    loop {
        index -= 1;
        tmp[index] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[index..]);
}

fn push_u64(buf: &mut Vec<u8>, mut value: u64) {
    let mut tmp = [0_u8; 20];
    let mut index = 20;
    loop {
        index -= 1;
        tmp[index] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[index..]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_merge_capacity_scales_with_input_size() {
        assert_eq!(parallel_merge_capacity_from_bytes(0), 1_024);
        assert_eq!(parallel_merge_capacity_from_bytes(8 * 512), 1_024);
        assert_eq!(
            parallel_merge_capacity_from_bytes(5_000 * ESTIMATED_BYTES_PER_MARKER),
            5_000
        );
        assert_eq!(
            parallel_merge_capacity_from_bytes(10_000_000 * 128),
            4_000_000
        );
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_merge_capacity_uses_input_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let mut input_files = Vec::new();
        for i in 0..8 {
            let path = dir.path().join(format!("ind{i}.fq"));
            let file = std::fs::File::create(&path).unwrap();
            file.set_len(256 * 1_024).unwrap();
            input_files.push(crate::io::seq_reader::InputFile {
                path,
                individual_name: format!("ind{i}"),
            });
        }

        assert_eq!(estimate_parallel_merge_capacity(&input_files), 32_768);
    }

    /// Header `#Number of markers` must equal retained data rows when min_depth filters.
    #[test]
    fn process_header_matches_retained_rows_with_min_depth() {
        use std::io::Write;

        let dir = tempfile::tempdir().unwrap();
        // Two individuals; shared high-depth tag + per-ind low-depth unique tags.
        let mut f0 = std::fs::File::create(dir.path().join("ind0.fa")).unwrap();
        writeln!(f0, ">r1").unwrap();
        writeln!(f0, "AAAAAAAA").unwrap(); // high — will appear depth 2 in ind0 after two reads
        writeln!(f0, ">r2").unwrap();
        writeln!(f0, "AAAAAAAA").unwrap();
        writeln!(f0, ">r3").unwrap();
        writeln!(f0, "CCCCCCCC").unwrap(); // depth 1 only — filtered by min_depth=2
        drop(f0);

        let mut f1 = std::fs::File::create(dir.path().join("ind1.fa")).unwrap();
        writeln!(f1, ">r1").unwrap();
        writeln!(f1, "AAAAAAAA").unwrap();
        writeln!(f1, ">r2").unwrap();
        writeln!(f1, "AAAAAAAA").unwrap();
        writeln!(f1, ">r3").unwrap();
        writeln!(f1, "GGGGGGGG").unwrap(); // depth 1 only — filtered
        drop(f1);

        let out = dir.path().join("markers.tsv");
        run(&ProcessParams {
            input_dir_path: dir.path().to_str().unwrap().to_string(),
            output_file_path: out.to_str().unwrap().to_string(),
            n_threads: 1,
            min_depth: 2,
            kmer_dedup: None,
        })
        .expect("process must succeed");

        let body = std::fs::read_to_string(&out).unwrap();
        let mut lines = body.lines();
        let header_line = lines.next().expect("header comment");
        assert!(
            header_line.starts_with("#Number of markers :"),
            "unexpected header: {header_line}"
        );
        let claimed: u64 = header_line
            .split(':')
            .nth(1)
            .unwrap()
            .trim()
            .parse()
            .expect("parse marker count");
        let _col_header = lines.next().expect("column header");
        let data_rows = lines.filter(|l| !l.is_empty()).count() as u64;
        assert_eq!(
            claimed, data_rows,
            "header count must match data rows; body:\n{body}"
        );
        // Only AAAAAAAA is retained (depth >= 2 in at least one individual).
        assert_eq!(claimed, 1, "expected one retained marker; body:\n{body}");
    }

    #[test]
    fn process_errors_when_an_input_file_is_unreadable() {
        use std::io::Write;

        let dir = tempfile::tempdir().unwrap();
        let mut f0 = std::fs::File::create(dir.path().join("ind0.fa")).unwrap();
        writeln!(f0, ">r1").unwrap();
        writeln!(f0, "AAAAAAAA").unwrap();
        drop(f0);
        // Corrupt gzip forces needletail / flate2 to fail for this individual.
        let mut bad = std::fs::File::create(dir.path().join("ind1.fq.gz")).unwrap();
        bad.write_all(b"this is not valid gzip content").unwrap();
        drop(bad);

        let out = dir.path().join("markers.tsv");
        let err = run(&ProcessParams {
            input_dir_path: dir.path().to_str().unwrap().to_string(),
            output_file_path: out.to_str().unwrap().to_string(),
            n_threads: 1,
            min_depth: 1,
            kmer_dedup: None,
        });
        assert!(
            err.is_err(),
            "must fail when an individual cannot be read, got Ok"
        );
    }
}
