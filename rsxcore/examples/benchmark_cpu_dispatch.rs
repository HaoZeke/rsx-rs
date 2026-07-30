// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Compare the direct scalar call with the dispatched one on the CPU path.
//!
//! Generalising the batch entry point over the test method once cost about a
//! factor of two here, by moving the match inside the per-marker loop. This
//! example is how that is checked.

use std::time::Instant;

use rsx_core::compute_backend::{AssociationCounts, PValueBackend, compute_p_batch_with_metrics};
use rsx_core::test_method::TestMethod;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let markers: usize = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "10000000".into())
        .parse()?;
    let (t1, t2) = (24u32, 24u32);
    let counts: Vec<_> = (0..markers)
        .map(|i| AssociationCounts {
            group1: (i as u32).wrapping_mul(17) % (t1 + 1),
            group2: (i as u32).wrapping_mul(29) % (t2 + 1),
        })
        .collect();

    for round in 1..=3 {
        let started = Instant::now();
        let direct: Vec<f64> = counts
            .iter()
            .map(|c| rsx_core::stats::p_association(c.group1, c.group2, t1, t2))
            .collect();
        let direct_seconds = started.elapsed().as_secs_f64();

        let dispatched = compute_p_batch_with_metrics(
            PValueBackend::Cpu,
            TestMethod::ChiSquared,
            &counts,
            t1,
            t2,
        )?;
        let dispatched_seconds = dispatched.metrics.total_seconds;

        assert_eq!(direct.len(), dispatched.p_values.try_as_slice()?.len());
        println!(
            "round {round}: direct={direct_seconds:.6}s dispatched={dispatched_seconds:.6}s ratio={:.2}",
            dispatched_seconds / direct_seconds
        );
    }
    Ok(())
}
