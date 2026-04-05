// GPL-3.0-or-later

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use radsex_core::stats;

fn bench_chi_squared(c: &mut Criterion) {
    c.bench_function("chi_squared_yates", |b| {
        b.iter(|| stats::chi_squared_yates(black_box(10), black_box(2), black_box(15), black_box(10)))
    });
}

fn bench_p_association(c: &mut Criterion) {
    c.bench_function("p_association", |b| {
        b.iter(|| stats::p_association(black_box(10), black_box(2), black_box(15), black_box(10)))
    });
}

criterion_group!(benches, bench_chi_squared, bench_p_association);
criterion_main!(benches);
