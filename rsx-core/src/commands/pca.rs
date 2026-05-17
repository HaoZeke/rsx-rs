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

/// Internal result of the core PCA computation.
/// Shared by the file-based writer and the in-memory Arrow path.
struct ComputedPca {
    eigenvalues: Vec<f64>,
    eigenvectors: Vec<f64>,      // row-major n x n
    sorted_indices: Vec<usize>,  // descending by eigenvalue
    n_markers: u64,
    n_individuals: usize,
    n_components: usize,         // actual r returned
    total_variance: f64,
    individual_names: Vec<String>,
}

/// Core streaming PCA (Gram + centering + Jacobi eigendecomposition).
/// Used by both the classic file writer and the Arrow in-memory path.
fn compute_pca(
    markers_table_path: &str,
    min_depth: u16,
    n_components: Option<usize>,
) -> Result<ComputedPca, Box<dyn std::error::Error>> {
    let table_path = Path::new(markers_table_path);

    let config = ParserConfig {
        store_sequence: false,
        store_depths: true,
        compute_groups: false,
        min_depth,
    };

    let stream = MarkersTableStream::open(table_path, None, config)?;
    let n = stream.header.n_individuals as usize;

    log::info!(
        "PCA: streaming {} individuals, building {}x{} Gram matrix",
        n, n, n
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

    // Eigendecompose
    log::info!("PCA: eigendecomposing {}x{} matrix", n, n);
    let (eigenvalues, eigenvectors) = jacobi_eigen(&mut gram, n);

    // Sort descending
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| eigenvalues[b].partial_cmp(&eigenvalues[a]).unwrap());

    let r = n_components.unwrap_or(n).min(n);
    let total_var: f64 = eigenvalues.iter().filter(|&&v| v > 0.0).sum();

    // Individual names (same logic as before)
    let header_cols = &stream.header.columns;
    let individual_names: Vec<String> = (0..n)
        .map(|i| {
            if i + 2 < header_cols.len() {
                header_cols[i + 2].clone()
            } else {
                format!("ind{}", i + 1)
            }
        })
        .collect();

    Ok(ComputedPca {
        eigenvalues,
        eigenvectors,
        sorted_indices: indices,
        n_markers,
        n_individuals: n,
        n_components: r,
        total_variance: total_var,
        individual_names,
    })
}

pub fn run(params: &PcaParams) -> Result<(), Box<dyn std::error::Error>> {
    let c = compute_pca(&params.markers_table_path, params.min_depth, params.n_components)?;

    std::fs::create_dir_all(&params.output_dir)?;

    // Write eigenvalues.tsv
    let eigenval_path = Path::new(&params.output_dir).join("eigenvalues.tsv");
    let mut f = std::fs::File::create(&eigenval_path)?;
    writeln!(f, "component\teigenvalue\tvariance_fraction\tcumulative")?;
    let mut cumulative = 0.0;
    for (k, &idx) in c.sorted_indices.iter().take(c.n_components).enumerate() {
        let ev = c.eigenvalues[idx].max(0.0);
        let frac = if c.total_variance > 0.0 { ev / c.total_variance } else { 0.0 };
        cumulative += frac;
        writeln!(f, "PC{}\t{:.6}\t{:.6}\t{:.6}", k + 1, ev, frac, cumulative)?;
    }

    // Write loadings.tsv
    let loadings_path = Path::new(&params.output_dir).join("loadings.tsv");
    let mut f = std::fs::File::create(&loadings_path)?;
    write!(f, "individual")?;
    for k in 0..c.n_components {
        write!(f, "\tPC{}", k + 1)?;
    }
    writeln!(f)?;

    for (i, name) in c.individual_names.iter().enumerate() {
        write!(f, "{}", name)?;
        for &idx in c.sorted_indices.iter().take(c.n_components) {
            write!(f, "\t{:.6}", c.eigenvectors[i * c.n_individuals + idx])?;
        }
        writeln!(f)?;
    }

    // Write summary.txt
    let summary_path = Path::new(&params.output_dir).join("summary.txt");
    let mut f = std::fs::File::create(&summary_path)?;
    writeln!(f, "Streaming PCA of depth matrix")?;
    writeln!(f, "Markers: {}", c.n_markers)?;
    writeln!(f, "Individuals: {}", c.n_individuals)?;
    writeln!(f, "Components: {}", c.n_components)?;
    writeln!(f, "Total variance: {:.2}", c.total_variance)?;
    writeln!(f)?;
    writeln!(f, "Top components:")?;
    cumulative = 0.0;
    for (k, &idx) in c.sorted_indices.iter().take(c.n_components.min(10)).enumerate() {
        let ev = c.eigenvalues[idx].max(0.0);
        let frac = if c.total_variance > 0.0 { ev / c.total_variance } else { 0.0 };
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
        c.n_markers,
        c.n_individuals,
        c.n_components,
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

#[cfg(feature = "arrow-output")]
use arrow::array::builder::{Float64Builder, StringBuilder};
#[cfg(feature = "arrow-output")]
use arrow::datatypes::{DataType, Field, Schema};
#[cfg(feature = "arrow-output")]
use arrow::record_batch::RecordBatch;

/// Result of in-memory PCA Arrow emission.
/// Two RecordBatches that match the classic TSV outputs exactly
/// (eigenvalues + loadings) plus provenance scalars.
#[cfg(feature = "arrow-output")]
pub struct PcaArrowResult {
    pub eigenvalues: RecordBatch,
    pub loadings: RecordBatch,
    pub n_markers: u64,
    pub n_individuals: usize,
    pub n_components: usize,
    pub total_variance: f64,
}

/// Real in-memory Arrow emission for PCA (no temp files for the data path).
///
/// Computes the identical streaming Gram + Jacobi decomposition as `run`,
/// then materializes two RecordBatches:
///   - eigenvalues: component | eigenvalue | variance_fraction | cumulative
///   - loadings:    individual | PC1 | PC2 | ... | PCr
#[cfg(feature = "arrow-output")]
pub fn run_to_arrow(params: &PcaParams) -> Result<PcaArrowResult, Box<dyn std::error::Error>> {
    let c = compute_pca(&params.markers_table_path, params.min_depth, params.n_components)?;

    // Build eigenvalues RecordBatch
    let eigen_schema = Schema::new(vec![
        Field::new("component", DataType::Utf8, false),
        Field::new("eigenvalue", DataType::Float64, false),
        Field::new("variance_fraction", DataType::Float64, false),
        Field::new("cumulative", DataType::Float64, false),
    ]);

    let n_rows = c.n_components;
    let mut comp_b = StringBuilder::with_capacity(n_rows, n_rows * 4);
    let mut ev_b = Float64Builder::with_capacity(n_rows);
    let mut frac_b = Float64Builder::with_capacity(n_rows);
    let mut cum_b = Float64Builder::with_capacity(n_rows);

    let mut cumulative = 0.0;
    for (k, &idx) in c.sorted_indices.iter().take(c.n_components).enumerate() {
        let ev = c.eigenvalues[idx].max(0.0);
        let frac = if c.total_variance > 0.0 { ev / c.total_variance } else { 0.0 };
        cumulative += frac;

        comp_b.append_value(format!("PC{}", k + 1));
        ev_b.append_value(ev);
        frac_b.append_value(frac);
        cum_b.append_value(cumulative);
    }

    let eigen_batch = RecordBatch::try_new(
        std::sync::Arc::new(eigen_schema),
        vec![
            std::sync::Arc::new(comp_b.finish()),
            std::sync::Arc::new(ev_b.finish()),
            std::sync::Arc::new(frac_b.finish()),
            std::sync::Arc::new(cum_b.finish()),
        ],
    )?;

    // Build loadings RecordBatch
    let mut loading_fields = vec![Field::new("individual", DataType::Utf8, false)];
    for k in 0..c.n_components {
        loading_fields.push(Field::new(format!("PC{}", k + 1), DataType::Float64, false));
    }
    let loadings_schema = Schema::new(loading_fields);

    let n_ind = c.n_individuals;
    let n_pc = c.n_components;
    let mut ind_b = StringBuilder::with_capacity(n_ind, n_ind * 16);
    let mut pc_builders: Vec<Float64Builder> = (0..n_pc).map(|_| Float64Builder::with_capacity(n_ind)).collect();

    for (i, name) in c.individual_names.iter().enumerate() {
        ind_b.append_value(name);
        for (k, &idx) in c.sorted_indices.iter().take(c.n_components).enumerate() {
            pc_builders[k].append_value(c.eigenvectors[i * c.n_individuals + idx]);
        }
    }

    let mut loading_columns: Vec<std::sync::Arc<dyn arrow::array::Array>> =
        vec![std::sync::Arc::new(ind_b.finish())];
    for mut b in pc_builders {
        loading_columns.push(std::sync::Arc::new(b.finish()));
    }

    let loadings_batch = RecordBatch::try_new(
        std::sync::Arc::new(loadings_schema),
        loading_columns,
    )?;

    Ok(PcaArrowResult {
        eigenvalues: eigen_batch,
        loadings: loadings_batch,
        n_markers: c.n_markers,
        n_individuals: c.n_individuals,
        n_components: c.n_components,
        total_variance: c.total_variance,
    })
}

#[cfg(all(test, feature = "arrow-output"))]
mod tests {
    use super::*;

    fn make_pca_test_data(dir: &std::path::Path) -> std::path::PathBuf {
        let path = dir.join("markers.tsv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "#Number of markers : 4").unwrap();
        writeln!(f, "id\tsequence\tind1\tind2\tind3\tind4\tind5").unwrap();
        writeln!(f, "0\tATCGATCG\t10\t5\t8\t12\t7").unwrap();
        writeln!(f, "1\tGGGGAAAA\t15\t20\t10\t0\t0").unwrap();
        writeln!(f, "2\tCCCCTTTT\t0\t0\t0\t25\t30").unwrap();
        writeln!(f, "3\tAAAATTTT\t5\t0\t3\t8\t6").unwrap();
        path
    }

    #[test]
    fn run_to_arrow_matches_file_based_pca() {
        let dir = std::env::temp_dir().join("rsx_pca_arrow_test");
        std::fs::create_dir_all(&dir).unwrap();

        let table = make_pca_test_data(&dir);
        let out_dir = dir.join("pca_file");
        std::fs::create_dir_all(&out_dir).unwrap();

        let params = PcaParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            output_dir: out_dir.to_str().unwrap().to_string(),
            min_depth: 1,
            n_components: Some(3),
        };

        // File path (writes the classic TSVs + summary)
        run(&params).unwrap();

        // Read the file-based eigenvalues to compare counts
        let eigen_file = std::fs::read_to_string(out_dir.join("eigenvalues.tsv")).unwrap();
        let eigen_lines: Vec<&str> = eigen_file.lines().filter(|l| !l.starts_with("component")).collect();
        let n_file_rows = eigen_lines.len();

        // Real Arrow path
        let arrow = run_to_arrow(&params).expect("pca run_to_arrow must succeed");

        // Differential checks
        assert_eq!(
            arrow.eigenvalues.num_rows(),
            n_file_rows,
            "Eigenvalues batch must have same number of components as the TSV"
        );
        assert_eq!(
            arrow.loadings.num_rows(),
            5,
            "Loadings must have one row per individual (5 in the fixture)"
        );
        assert!(
            arrow.total_variance > 0.0,
            "Total variance must be positive"
        );
        assert_eq!(arrow.n_components, 3);

        // Schema sanity
        let ev_schema = arrow.eigenvalues.schema();
        let ev_names: Vec<_> = ev_schema.fields().iter().map(|f| f.name().as_str()).collect();
        assert!(ev_names.contains(&"eigenvalue"));
        assert!(ev_names.contains(&"variance_fraction"));

        let ld_schema = arrow.loadings.schema();
        let ld_names: Vec<_> = ld_schema.fields().iter().map(|f| f.name().as_str()).collect();
        assert!(ld_names.contains(&"individual"));
        assert!(ld_names.contains(&"PC1"));
        assert!(ld_names.contains(&"PC3"));
    }
}
