// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `pca` command: streaming PCA of the depth matrix.
//!
//! Computes Tucker mode-2 factors via streaming Gram eigendecomposition.
//! Memory: O(n_individuals^2) = ~320KB for 200 individuals.
//!
//! Algorithm:
//! 1. Stream markers, accumulate X^T X (n_ind x n_ind) and mean vector
//! 2. Center the Gram matrix: C = X^T X - n * mu * mu^T
//! 3. Eigendecompose C (Jacobi rotation, exact for symmetric matrices)
//! 4. Output eigenvalues + loadings
//!
//! Mathematical proof: scripts/sympy/tucker_covariance_proof.py

use crate::markers_table::{MarkersTableStream, ParserConfig};
use std::io::Write;
use std::path::Path;

pub struct PcaParams {
    pub markers_table_path: String,
    pub output_dir: String,
    pub min_depth: u16,
    pub n_components: Option<usize>, // None = all
}

pub fn run(params: &PcaParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);

    let config = ParserConfig {
        store_sequence: false,
        store_depths: true,
        compute_groups: false,
        min_depth: params.min_depth,
    };

    let stream = MarkersTableStream::open(table_path, None, config)?;
    let n = stream.header.n_individuals as usize;

    log::info!(
        "PCA: streaming {} individuals, building {}x{} Gram matrix",
        n,
        n,
        n
    );

    // Accumulate Gram matrix C = X^T X and mean vector
    let mut gram = vec![0.0f64; n * n]; // row-major n x n
    let mut mean = vec![0.0f64; n];
    let mut n_markers = 0u64;

    stream.for_each(|marker| {
        if marker.n_individuals == 0 {
            return;
        }
        n_markers += 1;

        // Accumulate outer product: gram += x^T * x
        // Only upper triangle (symmetric), fill lower later
        for i in 0..n {
            let xi = marker.individual_depths[i] as f64;
            mean[i] += xi;
            for j in i..n {
                let xj = marker.individual_depths[j] as f64;
                gram[i * n + j] += xi * xj;
            }
        }
    })?;

    if n_markers == 0 {
        return Err("No markers found".into());
    }

    log::info!("PCA: {} markers streamed, centering Gram matrix", n_markers);

    // Fill lower triangle
    for i in 0..n {
        for j in 0..i {
            gram[i * n + j] = gram[j * n + i];
        }
    }

    // Center: C = X^T X - n_markers * mu * mu^T
    let nm = n_markers as f64;
    for m in &mut mean {
        *m /= nm;
    }
    for i in 0..n {
        for j in 0..n {
            gram[i * n + j] -= nm * mean[i] * mean[j];
        }
    }

    // Eigendecompose the centered Gram matrix (Jacobi rotation)
    log::info!("PCA: eigendecomposing {}x{} matrix", n, n);
    let (eigenvalues, eigenvectors) = jacobi_eigen(&mut gram, n);

    // Sort by eigenvalue descending
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| eigenvalues[b].partial_cmp(&eigenvalues[a]).unwrap());

    let r = params.n_components.unwrap_or(n).min(n);
    let total_var: f64 = eigenvalues.iter().filter(|&&v| v > 0.0).sum();

    // Create output directory
    std::fs::create_dir_all(&params.output_dir)?;

    // Write eigenvalues
    let eigenval_path = Path::new(&params.output_dir).join("eigenvalues.tsv");
    let mut f = std::fs::File::create(&eigenval_path)?;
    writeln!(f, "component\teigenvalue\tvariance_fraction\tcumulative")?;
    let mut cumulative = 0.0;
    for (k, &idx) in indices.iter().take(r).enumerate() {
        let ev = eigenvalues[idx].max(0.0);
        let frac = if total_var > 0.0 { ev / total_var } else { 0.0 };
        cumulative += frac;
        writeln!(f, "PC{}\t{:.6}\t{:.6}\t{:.6}", k + 1, ev, frac, cumulative)?;
    }

    // Write loadings (individual x component)
    let loadings_path = Path::new(&params.output_dir).join("loadings.tsv");
    let mut f = std::fs::File::create(&loadings_path)?;
    write!(f, "individual")?;
    for k in 0..r {
        write!(f, "\tPC{}", k + 1)?;
    }
    writeln!(f)?;

    let header_cols = &stream.header.columns;
    for i in 0..n {
        let name = if i + 2 < header_cols.len() {
            &header_cols[i + 2]
        } else {
            "?"
        };
        write!(f, "{}", name)?;
        for &idx in indices.iter().take(r) {
            write!(f, "\t{:.6}", eigenvectors[i * n + idx])?;
        }
        writeln!(f)?;
    }

    // Write summary
    let summary_path = Path::new(&params.output_dir).join("summary.txt");
    let mut f = std::fs::File::create(&summary_path)?;
    writeln!(f, "Streaming PCA of depth matrix")?;
    writeln!(f, "Markers: {}", n_markers)?;
    writeln!(f, "Individuals: {}", n)?;
    writeln!(f, "Components: {}", r)?;
    writeln!(f, "Total variance: {:.2}", total_var)?;
    writeln!(f)?;
    writeln!(f, "Top components:")?;
    cumulative = 0.0;
    for (k, &idx) in indices.iter().take(r.min(10)).enumerate() {
        let ev = eigenvalues[idx].max(0.0);
        let frac = if total_var > 0.0 { ev / total_var } else { 0.0 };
        cumulative += frac;
        writeln!(
            f,
            "  PC{}: {:.4} variance ({:.1}% cumulative)",
            k + 1,
            frac,
            cumulative * 100.0
        )?;
    }

    log::info!(
        "PCA done: {} markers, {} individuals, {} components -> {}",
        n_markers,
        n,
        r,
        params.output_dir
    );

    Ok(())
}

/// Jacobi eigendecomposition for symmetric matrix.
/// Input: row-major n x n symmetric matrix (modified in place -> diagonal).
/// Returns: (eigenvalues, eigenvectors as row-major n x n).
fn jacobi_eigen(a: &mut [f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    // Initialize eigenvector matrix to identity
    let mut v = vec![0.0f64; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }

    let max_sweeps = 100;
    let tol = 1e-15;

    for _sweep in 0..max_sweeps {
        // Find max off-diagonal element
        let mut max_off = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                max_off = max_off.max(a[i * n + j].abs());
            }
        }
        if max_off < tol {
            break;
        }

        for i in 0..n {
            for j in (i + 1)..n {
                let aij = a[i * n + j];
                if aij.abs() < tol * 0.01 {
                    continue;
                }

                let aii = a[i * n + i];
                let ajj = a[j * n + j];
                let diff = ajj - aii;

                let t = if diff.abs() < tol {
                    1.0
                } else {
                    let tau = diff / (2.0 * aij);
                    let sign = if tau >= 0.0 { 1.0 } else { -1.0 };
                    sign / (tau.abs() + (1.0 + tau * tau).sqrt())
                };

                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;

                // Rotate matrix A
                a[i * n + i] = aii - t * aij;
                a[j * n + j] = ajj + t * aij;
                a[i * n + j] = 0.0;
                a[j * n + i] = 0.0;

                for k in 0..n {
                    if k == i || k == j {
                        continue;
                    }
                    let aki = a[k * n + i];
                    let akj = a[k * n + j];
                    a[k * n + i] = c * aki - s * akj;
                    a[i * n + k] = a[k * n + i];
                    a[k * n + j] = s * aki + c * akj;
                    a[j * n + k] = a[k * n + j];
                }

                // Rotate eigenvectors
                for k in 0..n {
                    let vki = v[k * n + i];
                    let vkj = v[k * n + j];
                    v[k * n + i] = c * vki - s * vkj;
                    v[k * n + j] = s * vki + c * vkj;
                }
            }
        }
    }

    let eigenvalues: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    (eigenvalues, v)
}
