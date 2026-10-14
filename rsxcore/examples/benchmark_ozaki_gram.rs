// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Time the PCA Gram kernels: f64 rank-1, direct u64 GEMM, and Ozaki scheme II.
//!
//! Usage:
//!   cargo run -p rsxcore --release --example benchmark_ozaki_gram

use std::hint::black_box;
use std::time::Instant;

use rsx_core::ozaki_ii::{
    gram_upper_f32, gram_upper_f64, gram_upper_f64_simd, gram_upper_ozaki_ii, gram_upper_u64,
    u64_values_fit_f64,
};

fn fill(markers: usize, individuals: usize, seed: u64, depth_mod: u16) -> Vec<u16> {
    let mut state = seed;
    let mut rows = Vec::with_capacity(markers * individuals);
    for _ in 0..markers * individuals {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        rows.push((state >> 33) as u16 % depth_mod);
    }
    rows
}

fn time_ms(repeats: usize, mut body: impl FnMut()) -> (f64, f64) {
    body();
    let mut samples = Vec::with_capacity(repeats);
    for _ in 0..repeats {
        let start = Instant::now();
        body();
        samples.push(start.elapsed().as_secs_f64() * 1_000.0);
    }
    let mean = samples.iter().sum::<f64>() / repeats as f64;
    let var = samples
        .iter()
        .map(|sample| {
            let d = sample - mean;
            d * d
        })
        .sum::<f64>()
        / (repeats.saturating_sub(1).max(1) as f64);
    (mean, var.sqrt())
}

fn main() {
    let repeats = 5;
    let shapes = [(100_000, 40, 51_u16), (20_000, 128, 40), (4_096, 64, 30)];
    println!("kernel\tmarkers\tindividuals\tmean_ms\tstdev_ms\texact_f64");
    for (markers, individuals, depth_mod) in shapes {
        let rows = fill(markers, individuals, 2026, depth_mod);
        let (f64_gram, _) = gram_upper_f64(&rows, individuals);
        let (u64_gram, u64_mean) = gram_upper_u64(&rows, individuals);
        let (ozaki, ozaki_mean) = gram_upper_ozaki_ii(&rows, individuals).expect("ozaki range");
        assert_eq!(ozaki, u64_gram, "Ozaki II diverged from u64 GEMM");
        assert_eq!(ozaki_mean, u64_mean);
        let exact = u64_values_fit_f64(&u64_gram)
            && f64_gram
                .iter()
                .zip(u64_gram.iter())
                .all(|(left, right)| left.to_bits() == (*right as f64).to_bits());

        let (f64_mean, f64_sd) = time_ms(repeats, || {
            black_box(gram_upper_f64(black_box(&rows), individuals));
        });
        let (u64_mean_ms, u64_sd) = time_ms(repeats, || {
            black_box(gram_upper_u64(black_box(&rows), individuals));
        });
        let (oz_mean, oz_sd) = time_ms(repeats, || {
            black_box(gram_upper_ozaki_ii(black_box(&rows), individuals).unwrap());
        });
        let (f32_gram, _) = gram_upper_f32(&rows, individuals);
        let f32_exact = f32_gram.iter().zip(u64_gram.iter()).all(|(&low, &wide)| {
            wide <= (1_u64 << 24) && (wide as f64).to_bits() == f64::from(low).to_bits()
        });
        let f32_err = f32_gram
            .iter()
            .zip(u64_gram.iter())
            .map(|(&low, &wide)| (f64::from(low) - wide as f64).abs())
            .fold(0.0_f64, f64::max);
        let (simd_gram, simd_mean) = gram_upper_f64_simd(&rows, individuals);
        assert_eq!(simd_gram.len(), f64_gram.len(), "simd gram shape diverged");
        assert!(
            simd_gram
                .iter()
                .zip(f64_gram.iter())
                .all(|(left, right)| left.to_bits() == right.to_bits()),
            "simd f64 gram diverged from scalar f64"
        );
        assert!(
            simd_mean
                .iter()
                .zip(u64_mean.iter())
                .all(|(left, right)| { left.to_bits() == (*right as f64).to_bits() })
        );

        println!("f64_rank1\t{markers}\t{individuals}\t{f64_mean:.3}\t{f64_sd:.3}\t{exact}");
        println!("u64_gemm\t{markers}\t{individuals}\t{u64_mean_ms:.3}\t{u64_sd:.3}\t{exact}");
        println!("ozaki_ii\t{markers}\t{individuals}\t{oz_mean:.3}\t{oz_sd:.3}\t{exact}");
        let (f32_mean_ms, f32_sd) = time_ms(repeats, || {
            black_box(gram_upper_f32(black_box(&rows), individuals));
        });
        let (simd_mean_ms, simd_sd) = time_ms(repeats, || {
            black_box(gram_upper_f64_simd(black_box(&rows), individuals));
        });
        println!(
            "f32_rank1\t{markers}\t{individuals}\t{f32_mean_ms:.3}\t{f32_sd:.3}\t{f32_exact}\tmax_abs_err={f32_err:.1}"
        );
        println!("f64_simd\t{markers}\t{individuals}\t{simd_mean_ms:.3}\t{simd_sd:.3}\t{exact}");
    }

    println!("\n# bitset popcount, 500000 markers x 1 word (40 individuals)");
    let words: Vec<u64> = (0..500_000)
        .map(|i| 0x0000_FFFF_FFF0_u64.wrapping_mul(i as u64 + 1))
        .collect();
    let mask = (1_u64 << 40) - 1;
    let hw = popcount_hw(&words, mask);
    let naive = popcount_naive(&words, mask);
    assert_eq!(hw, naive);
    let (hw_ms, hw_sd) = time_ms(repeats, || {
        black_box(popcount_hw(black_box(&words), mask));
    });
    let (naive_ms, naive_sd) = time_ms(repeats, || {
        black_box(popcount_naive(black_box(&words), mask));
    });
    println!("popcount_hw\t{hw_ms:.3}\t{hw_sd:.3}");
    println!("popcount_naive\t{naive_ms:.3}\t{naive_sd:.3}");

    println!(
        "\n# canonical k-mer min-hash, 100000 tags x 94 bp, k=31 (not on the default process path)"
    );
    let tags = fill_tags(100_000, 94, 7);
    let naive_hashes: Vec<u64> = tags.iter().map(|tag| canonical_naive(tag, 31)).collect();
    let rolling_hashes: Vec<u64> = tags.iter().map(|tag| canonical_rolling(tag, 31)).collect();
    assert_eq!(naive_hashes, rolling_hashes);
    let (naive_ms, naive_sd) = time_ms(repeats, || {
        for tag in &tags {
            black_box(canonical_naive(black_box(tag), 31));
        }
    });
    let (roll_ms, roll_sd) = time_ms(repeats, || {
        for tag in &tags {
            black_box(canonical_rolling(black_box(tag), 31));
        }
    });
    println!("kmer_window\t{naive_ms:.3}\t{naive_sd:.3}");
    println!("kmer_rolling\t{roll_ms:.3}\t{roll_sd:.3}");
}

fn popcount_hw(words: &[u64], mask: u64) -> u32 {
    words.iter().map(|word| (word & mask).count_ones()).sum()
}

fn popcount_naive(words: &[u64], mask: u64) -> u32 {
    let mut total = 0_u32;
    for word in words {
        let mut bits = word & mask;
        while bits != 0 {
            total += (bits & 1) as u32;
            bits >>= 1;
        }
    }
    total
}

fn fill_tags(count: usize, length: usize, seed: u64) -> Vec<Vec<u8>> {
    let bases = *b"ACGT";
    let mut state = seed;
    let mut tags = Vec::with_capacity(count);
    for _ in 0..count {
        let mut tag = Vec::with_capacity(length);
        for _ in 0..length {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            tag.push(bases[(state >> 33) as usize & 3]);
        }
        tags.push(tag);
    }
    tags
}

fn base_code(base: u8) -> u64 {
    match base {
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

fn kmer_to_u64(kmer: &[u8]) -> u64 {
    let mut value = 0_u64;
    for &base in kmer {
        value = (value << 2) | base_code(base);
    }
    value
}

fn revcomp_u64(mut kmer: u64, k: usize) -> u64 {
    let mut rc = 0_u64;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (kmer & 3));
        kmer >>= 2;
    }
    rc
}

fn canonical_naive(seq: &[u8], k: usize) -> u64 {
    let mut best = u64::MAX;
    for window in seq.windows(k) {
        let forward = kmer_to_u64(window);
        best = best.min(forward.min(revcomp_u64(forward, k)));
    }
    best
}

fn canonical_rolling(seq: &[u8], k: usize) -> u64 {
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1_u64 << (2 * k)) - 1
    };
    let mut forward = 0_u64;
    for &base in &seq[..k] {
        forward = (forward << 2) | base_code(base);
    }
    let mut reverse = revcomp_u64(forward, k);
    let mut best = forward.min(reverse);
    let top = 2 * (k - 1);
    for &base in &seq[k..] {
        let code = base_code(base);
        forward = ((forward << 2) | code) & mask;
        reverse = (reverse >> 2) | ((3 - code) << top);
        best = best.min(forward.min(reverse));
    }
    best
}
