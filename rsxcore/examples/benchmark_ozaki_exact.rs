// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Time the exact association kernels against the lgamma and f64 paths.
//!
//! Usage:
//!   cargo run -p rsxcore --release --example benchmark_ozaki_exact

use std::hint::black_box;
use std::time::Instant;

use rsx_core::ozaki_exact;
use rsx_core::stats::{self, BayesFactorModel};

fn time_us(iters: usize, mut body: impl FnMut()) -> f64 {
    body();
    let start = Instant::now();
    for _ in 0..iters {
        body();
    }
    start.elapsed().as_secs_f64() * 1e6 / iters as f64
}

fn main() {
    let model = BayesFactorModel::uniform_v1();
    println!("kernel\tgroup\tus_per_call\tratio_to_current");
    for &n in &[10u32, 20, 40, 80, 200] {
        let iters = match n {
            10 | 20 => 2_000,
            40 => 400,
            80 => 80,
            _ => 20,
        };
        let present1 = n / 2;
        let present2 = n / 5;
        let fisher_old = time_us(iters, || {
            black_box(stats::fisher_exact(present1, present2, n, n));
        });
        let fisher_new = time_us(iters, || {
            black_box(ozaki_exact::fisher_probability(present1, present2, n, n));
        });
        let bayes_old = time_us(iters, || {
            black_box(stats::bayes_factor_2x2(present1, present2, n, n));
        });
        let bayes_new = time_us(iters, || {
            black_box(
                ozaki_exact::bayes_factor_integer_prior(present1, present2, n, n, &model)
                    .expect("uniform prior"),
            );
        });
        let chi_old = time_us(iters, || {
            black_box(stats::chi_squared_yates(present1, present2, n, n));
        });
        let chi_new = time_us(iters, || {
            black_box(ozaki_exact::chi_squared_yates_rounded(
                present1, present2, n, n,
            ));
        });
        println!("fisher_lgamma\t{n}\t{fisher_old:.3}\t1");
        println!(
            "fisher_integer\t{n}\t{fisher_new:.3}\t{:.3}",
            fisher_new / fisher_old
        );
        println!("bayes_lgamma\t{n}\t{bayes_old:.3}\t1");
        println!(
            "bayes_integer\t{n}\t{bayes_new:.3}\t{:.3}",
            bayes_new / bayes_old
        );
        println!("chi_f64\t{n}\t{chi_old:.3}\t1");
        println!("chi_rational\t{n}\t{chi_new:.3}\t{:.3}", chi_new / chi_old);
    }

    let mut worst_rel = 0.0f64;
    let mut disagreed = 0u32;
    let mut tables = 0u32;
    for t in 1..41u32 {
        for n1 in 0..=t {
            for n2 in 0..=t {
                tables += 1;
                let old = stats::fisher_exact(n1, n2, t, t);
                let new = ozaki_exact::fisher_p_value(n1, n2, t, t);
                let rel = if old == 0.0 {
                    0.0
                } else {
                    (old - new).abs() / old
                };
                worst_rel = worst_rel.max(rel);
                if old.to_bits() != new.to_bits() {
                    disagreed += 1;
                }
            }
        }
    }
    println!(
        "fisher_compare\ttables\t{tables}\tdisagreed\t{disagreed}\tworst_rel\t{worst_rel:.3e}"
    );
}
