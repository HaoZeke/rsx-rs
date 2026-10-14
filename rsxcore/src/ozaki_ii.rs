// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Exact integer Gram products and an Ozaki scheme II emulation.
//!
//! The only GEMM-shaped kernel in rsx is the PCA Gram
//! \(G = X^{\mathsf T} X\) of the marker-by-individual depth matrix. Each
//! marker is one rank-1 update. Logistic regression also forms \(X^{\mathsf T} W X\),
//! but that design matrix is a handful of covariates, not the marker table.
//! Group contrasts are bitset popcounts, not matrix products.
//!
//! Depths are non-negative integers that fit in `u16`. A direct `u64` GEMM is
//! therefore exact. [Ozaki scheme II](https://arxiv.org/abs/2504.08009)
//! (Ozaki, Uchino, Imamura, 2025) reconstructs the same integer product from
//! several INT8-range residue GEMMs via the Chinese Remainder Theorem, which
//! is how FP64 GEMM is emulated on INT8 tensor cores. Both kernels live here
//! so the production PCA path can keep the direct integer product when the
//! entries still fit in the f64 mantissa, and so the CRT path can be measured
//! without changing command output.

/// Pairwise-coprime moduli used by Ozaki scheme II for INT8 tensor cores.
/// Each modulus is at most 256, so a residue fits in one byte. The list is the
/// \(s = 16\) set from Ozaki, Uchino, and Imamura (arXiv:2504.08009).
const OZAKI_II_MODULI: [u32; 16] = [
    256, 255, 253, 251, 247, 239, 233, 229, 227, 223, 217, 211, 199, 197, 193, 191,
];

/// Largest integer that is still an exact `f64` and whose every smaller
/// non-negative integer is exact as well.
pub const F64_EXACT_INTEGER_MAX: u64 = 1_u64 << 53;

/// Why an Ozaki II product could not be formed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OzakiError {
    /// `rows.len()` is not a multiple of the individual count.
    Shape,
    /// The product of the INT8 moduli does not cover the largest possible entry.
    Range,
}

impl std::fmt::Display for OzakiError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Shape => write!(formatter, "depth matrix length is not a multiple of n"),
            Self::Range => write!(
                formatter,
                "Ozaki II moduli cannot cover the integer Gram range"
            ),
        }
    }
}

impl std::error::Error for OzakiError {}

/// `true` when every entry converts to `f64` with no rounding.
#[inline]
pub fn u64_values_fit_f64(values: &[u64]) -> bool {
    values.iter().all(|&value| value <= F64_EXACT_INTEGER_MAX)
}

/// Serial f64 rank-1 Gram. Matches `GramAccumulator`'s CPU update, including
/// the zero lower triangle (the caller mirrors it).
pub fn gram_upper_f64(rows: &[u16], individuals: usize) -> (Vec<f64>, Vec<f64>) {
    let mut gram = vec![0.0; individuals * individuals];
    let mut mean = vec![0.0; individuals];
    for row in rows.chunks(individuals) {
        for (i, &di) in row.iter().enumerate() {
            let xi = f64::from(di);
            mean[i] += xi;
            let base = i * individuals;
            for (j, &dj) in row.iter().enumerate().skip(i) {
                gram[base + j] += xi * f64::from(dj);
            }
        }
    }
    (gram, mean)
}

/// f32 rank-1 Gram. Same update as [`gram_upper_f64`], narrowed to binary32.
///
/// Binary32 integers are exact only through `2^24`. A RAD-seq depth table
/// (about `10^5` markers, depths of a few tens) overflows that range, so this
/// kernel is a measurement of mixed precision, not a production accumulator.
pub fn gram_upper_f32(rows: &[u16], individuals: usize) -> (Vec<f32>, Vec<f32>) {
    let mut gram = vec![0.0_f32; individuals * individuals];
    let mut mean = vec![0.0_f32; individuals];
    for row in rows.chunks(individuals) {
        for (i, &di) in row.iter().enumerate() {
            let xi = f32::from(di);
            mean[i] += xi;
            let base = i * individuals;
            for (j, &dj) in row.iter().enumerate().skip(i) {
                gram[base + j] += xi * f32::from(dj);
            }
        }
    }
    (gram, mean)
}

/// Direct integer GEMM for \(X^{\mathsf T} X\). Upper triangle plus column sums.
/// Products of `u16` depths accumulate in `u64`.
pub fn gram_upper_u64(rows: &[u16], individuals: usize) -> (Vec<u64>, Vec<u64>) {
    let mut gram = vec![0_u64; individuals * individuals];
    let mut mean = vec![0_u64; individuals];
    for row in rows.chunks(individuals) {
        accumulate_row_u64(&mut gram, &mut mean, row, individuals);
    }
    (gram, mean)
}

#[inline]
pub(crate) fn accumulate_row_u64(
    gram: &mut [u64],
    mean: &mut [u64],
    row: &[u16],
    individuals: usize,
) {
    for (i, &di) in row.iter().enumerate() {
        let xi = u64::from(di);
        mean[i] += xi;
        let base = i * individuals;
        for (j, &dj) in row.iter().enumerate().skip(i) {
            gram[base + j] += xi * u64::from(dj);
        }
    }
}

/// f64 Gram using AVX-512, then AVX2, then the scalar rank-1 update.
///
/// Each Gram entry is its own accumulator and is updated once per marker, in
/// marker order, so when every partial sum is an exact `f64` integer the bits
/// match [`gram_upper_f64`].
pub fn gram_upper_f64_simd(rows: &[u16], individuals: usize) -> (Vec<f64>, Vec<f64>) {
    #[cfg(all(target_arch = "x86_64", rsx_avx512_intrinsics))]
    {
        if std::is_x86_feature_detected!("avx512f") {
            // Safety: the feature check above is the AVX-512F contract.
            return unsafe { gram_upper_f64_avx512(rows, individuals) };
        }
    }
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
            // Safety: the feature checks above are the AVX2+FMA contract.
            return unsafe { gram_upper_f64_avx2(rows, individuals) };
        }
    }
    gram_upper_f64(rows, individuals)
}

#[cfg(all(target_arch = "x86_64", rsx_avx512_intrinsics))]
#[allow(clippy::incompatible_msrv)]
#[target_feature(enable = "avx512f")]
unsafe fn gram_upper_f64_avx512(rows: &[u16], individuals: usize) -> (Vec<f64>, Vec<f64>) {
    use std::arch::x86_64::{_mm512_fmadd_pd, _mm512_loadu_pd, _mm512_set1_pd, _mm512_storeu_pd};
    let (mut gram, mut mean, mut row_f) = f64_buffers(individuals);
    for row in rows.chunks(individuals) {
        load_row_f64(row, &mut row_f, &mut mean);
        for i in 0..individuals {
            let xi = row_f[i];
            let broadcast = unsafe { _mm512_set1_pd(xi) };
            let base = i * individuals;
            let mut j = i;
            while j + 8 <= individuals {
                unsafe {
                    let acc = _mm512_loadu_pd(gram.as_ptr().add(base + j));
                    let rhs = _mm512_loadu_pd(row_f.as_ptr().add(j));
                    let updated = _mm512_fmadd_pd(broadcast, rhs, acc);
                    _mm512_storeu_pd(gram.as_mut_ptr().add(base + j), updated);
                }
                j += 8;
            }
            scalar_tail(&mut gram, &row_f, base, xi, j, individuals);
        }
    }
    (gram, mean)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn gram_upper_f64_avx2(rows: &[u16], individuals: usize) -> (Vec<f64>, Vec<f64>) {
    use std::arch::x86_64::{_mm256_fmadd_pd, _mm256_loadu_pd, _mm256_set1_pd, _mm256_storeu_pd};
    let (mut gram, mut mean, mut row_f) = f64_buffers(individuals);
    for row in rows.chunks(individuals) {
        load_row_f64(row, &mut row_f, &mut mean);
        for i in 0..individuals {
            let xi = row_f[i];
            let broadcast = unsafe { _mm256_set1_pd(xi) };
            let base = i * individuals;
            let mut j = i;
            while j + 4 <= individuals {
                unsafe {
                    let acc = _mm256_loadu_pd(gram.as_ptr().add(base + j));
                    let rhs = _mm256_loadu_pd(row_f.as_ptr().add(j));
                    let updated = _mm256_fmadd_pd(broadcast, rhs, acc);
                    _mm256_storeu_pd(gram.as_mut_ptr().add(base + j), updated);
                }
                j += 4;
            }
            scalar_tail(&mut gram, &row_f, base, xi, j, individuals);
        }
    }
    (gram, mean)
}

fn f64_buffers(individuals: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    (
        vec![0.0; individuals * individuals],
        vec![0.0; individuals],
        vec![0.0; individuals],
    )
}

fn load_row_f64(row: &[u16], row_f: &mut [f64], mean: &mut [f64]) {
    for (slot, &depth) in row_f.iter_mut().zip(row) {
        *slot = f64::from(depth);
    }
    for (sum, &value) in mean.iter_mut().zip(row_f.iter()) {
        *sum += value;
    }
}

fn scalar_tail(
    gram: &mut [f64],
    row_f: &[f64],
    base: usize,
    xi: f64,
    mut j: usize,
    individuals: usize,
) {
    while j < individuals {
        gram[base + j] += xi * row_f[j];
        j += 1;
    }
}

/// Ozaki scheme II integer Gram.
///
/// Each depth is reduced modulo a pairwise-coprime INT8 modulus. The residue
/// matrices are multiplied in blocks whose INT32 accumulators cannot overflow
/// (the same contract INT8 tensor-core GEMM provides for an inner dimension
/// below \(2^{17}\)). Block residues are reduced and CRT-combined into the
/// exact non-negative product.
pub fn gram_upper_ozaki_ii(
    rows: &[u16],
    individuals: usize,
) -> Result<(Vec<u64>, Vec<u64>), OzakiError> {
    if individuals == 0 || rows.len() % individuals != 0 {
        return Err(OzakiError::Shape);
    }
    let markers = rows.len() / individuals;
    let max_depth = rows.iter().copied().max().unwrap_or(0) as u128;
    let bound = (markers as u128) * max_depth * max_depth;
    let moduli = moduli_covering(bound)?;
    let modulus_product: u128 = moduli.iter().map(|&modulus| u128::from(modulus)).product();

    let mut gram = vec![0_u128; individuals * individuals];
    for &modulus in &moduli {
        let residue = modular_upper_gram(rows, individuals, modulus);
        let modulus_u = u128::from(modulus);
        let cofactor = modulus_product / modulus_u;
        let inverse = mod_inverse(cofactor % modulus_u, modulus_u);
        let scale = cofactor * inverse % modulus_product;
        for (entry, &residue_ij) in gram.iter_mut().zip(residue.iter()) {
            let term = (u128::from(residue_ij) * scale) % modulus_product;
            *entry = (*entry + term) % modulus_product;
        }
    }
    if gram.iter().any(|&entry| entry > u128::from(u64::MAX)) {
        return Err(OzakiError::Range);
    }
    let gram = gram.into_iter().map(|entry| entry as u64).collect();
    // Column sums are a vector reduction, not a product. They stay exact and
    // independent of the CRT so callers can center the Gram the same way.
    Ok((gram, column_sums(rows, individuals)))
}

fn column_sums(rows: &[u16], individuals: usize) -> Vec<u64> {
    let mut mean = vec![0_u64; individuals];
    for row in rows.chunks(individuals) {
        for (sum, &depth) in mean.iter_mut().zip(row) {
            *sum += u64::from(depth);
        }
    }
    mean
}

fn moduli_covering(bound: u128) -> Result<Vec<u32>, OzakiError> {
    let mut product = 1_u128;
    let mut chosen = Vec::new();
    for &modulus in &OZAKI_II_MODULI {
        product = product
            .checked_mul(u128::from(modulus))
            .ok_or(OzakiError::Range)?;
        chosen.push(modulus);
        if product > bound {
            return Ok(chosen);
        }
    }
    Err(OzakiError::Range)
}

/// INT8×INT8 → INT32 block GEMM, reduced modulo `modulus`.
///
/// A block of 8_192 markers satisfies `8192 * 255^2 < 2^31`, so each block
/// product is an error-free INT32 GEMM of byte residues. That is the CPU-side
/// stand-in for the INT8 tensor-core product in Ozaki scheme II.
fn modular_upper_gram(rows: &[u16], individuals: usize, modulus: u32) -> Vec<u32> {
    const BLOCK_MARKERS: usize = 8_192;
    let mut gram = vec![0_u32; individuals * individuals];
    let mut block = vec![0_u32; individuals * individuals];
    let mut residue = vec![0_u8; individuals];
    let mut filled = 0_usize;
    for row in rows.chunks(individuals) {
        for (slot, &depth) in residue.iter_mut().zip(row) {
            *slot = (u32::from(depth) % modulus) as u8;
        }
        for i in 0..individuals {
            let left = u32::from(residue[i]);
            let base = i * individuals;
            for j in i..individuals {
                block[base + j] += left * u32::from(residue[j]);
            }
        }
        filled += 1;
        if filled == BLOCK_MARKERS {
            fold_block(&mut gram, &mut block, modulus);
            filled = 0;
        }
    }
    if filled != 0 {
        fold_block(&mut gram, &mut block, modulus);
    }
    gram
}

fn fold_block(gram: &mut [u32], block: &mut [u32], modulus: u32) {
    for (entry, product) in gram.iter_mut().zip(block.iter()) {
        *entry = ((*entry as u64 + u64::from(*product)) % u64::from(modulus)) as u32;
    }
    block.fill(0);
}

fn mod_inverse(value: u128, modulus: u128) -> u128 {
    let mut t = 0_i128;
    let mut new_t = 1_i128;
    let mut r = modulus as i128;
    let mut new_r = value as i128;
    while new_r != 0 {
        let quotient = r / new_r;
        let next_t = t - quotient * new_t;
        t = new_t;
        new_t = next_t;
        let next_r = r - quotient * new_r;
        r = new_r;
        new_r = next_r;
    }
    if t < 0 {
        t += modulus as i128;
    }
    t as u128
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn u64_gram_matches_f64_rank1_inside_the_mantissa() {
        let rows = fill(80, 6, 11, 40);
        let (f64_gram, f64_mean) = gram_upper_f64(&rows, 6);
        let (u64_gram, u64_mean) = gram_upper_u64(&rows, 6);
        assert!(u64_values_fit_f64(&u64_gram));
        for (left, right) in f64_gram.iter().zip(u64_gram.iter()) {
            assert_eq!(left.to_bits(), (*right as f64).to_bits());
        }
        for (left, right) in f64_mean.iter().zip(u64_mean.iter()) {
            assert_eq!(left.to_bits(), (*right as f64).to_bits());
        }
    }

    #[test]
    fn ozaki_ii_matches_direct_integer_gemm() {
        let individuals = 5;
        let rows = fill(300, individuals, 99, 25);
        let (exact, exact_mean) = gram_upper_u64(&rows, individuals);
        let (ozaki, ozaki_mean) = gram_upper_ozaki_ii(&rows, individuals).unwrap();
        assert_eq!(ozaki, exact);
        assert_eq!(ozaki_mean, exact_mean);
    }

    #[test]
    fn ozaki_ii_rejects_a_ragged_matrix() {
        let rows = vec![1_u16, 2, 3];
        assert_eq!(gram_upper_ozaki_ii(&rows, 2), Err(OzakiError::Shape));
    }

    #[test]
    fn simd_f64_gram_matches_scalar_bits() {
        let rows = fill(60, 13, 3, 40);
        let (scalar, scalar_mean) = gram_upper_f64(&rows, 13);
        let (simd, simd_mean) = gram_upper_f64_simd(&rows, 13);
        assert_eq!(scalar, simd);
        assert_eq!(scalar_mean, simd_mean);
    }

    #[test]
    fn f32_gram_drops_bits_on_a_table_sized_product() {
        let individuals = 4;
        let rows = fill(100_000, individuals, 5, 49);
        let (exact, _) = gram_upper_u64(&rows, individuals);
        let (narrow, _) = gram_upper_f32(&rows, individuals);
        let diverged = exact.iter().zip(narrow.iter()).any(|(&wide, &low)| {
            wide <= F64_EXACT_INTEGER_MAX && (wide as f64).to_bits() != f64::from(low).to_bits()
        });
        assert!(diverged, "f32 stayed exact; the fixture is too small");
    }
}
