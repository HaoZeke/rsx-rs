// GPL-3.0-or-later
// Criterion microbenchmarks for rsx-rs hot paths.

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_chi_squared(c: &mut Criterion) {
    use rsx_core::stats;
    c.bench_function("chi_squared_yates", |b| {
        b.iter(|| stats::chi_squared_yates(black_box(10), black_box(2), black_box(15), black_box(10)))
    });
}

fn bench_p_association(c: &mut Criterion) {
    use rsx_core::stats;
    c.bench_function("p_association", |b| {
        b.iter(|| stats::p_association(black_box(10), black_box(2), black_box(15), black_box(10)))
    });
}

fn bench_cg_format(c: &mut Criterion) {
    use rsx_core::stats::Cg;
    c.bench_function("cg_format", |b| {
        b.iter(|| format!("{}", Cg(black_box(0.000456057))))
    });
}

fn bench_bitset_popcount_small(c: &mut Criterion) {
    use rsx_core::bitset::{BitsetRow, GroupMask};
    let mut row = BitsetRow::new(40);
    for i in (0..40).step_by(2) { row.set(i); }
    let groups: Vec<String> = (0..42).map(|i| {
        if i < 2 { String::new() } else if i < 22 { "M".into() } else { "F".into() }
    }).collect();
    let mask = GroupMask::from_columns(&groups, "M", 40);
    c.bench_function("bitset_popcount_40ind", |b| {
        b.iter(|| black_box(row.count_masked(&mask)))
    });
}

fn bench_bitset_popcount_large(c: &mut Criterion) {
    use rsx_core::bitset::{BitsetRow, GroupMask};
    let mut row = BitsetRow::new(200);
    for i in (0..200).step_by(3) { row.set(i); }
    let groups: Vec<String> = (0..202).map(|i| {
        if i < 2 { String::new() } else if i < 102 { "M".into() } else { "F".into() }
    }).collect();
    let mask = GroupMask::from_columns(&groups, "M", 200);
    c.bench_function("bitset_popcount_200ind", |b| {
        b.iter(|| black_box(row.count_masked(&mask)))
    });
}

fn bench_bitset_popcount_1000(c: &mut Criterion) {
    use rsx_core::bitset::{BitsetRow, GroupMask};
    let mut row = BitsetRow::new(1000);
    for i in (0..1000).step_by(2) { row.set(i); }
    let groups: Vec<String> = (0..1002).map(|i| {
        if i < 2 { String::new() } else if i < 502 { "M".into() } else { "F".into() }
    }).collect();
    let mask = GroupMask::from_columns(&groups, "M", 1000);
    c.bench_function("bitset_popcount_1000ind", |b| {
        b.iter(|| black_box(row.count_masked(&mask)))
    });
}

fn bench_fast_parse_u16(c: &mut Criterion) {
    use rsx_core::io::table_io::fast_parse_u16;
    c.bench_function("fast_parse_u16", |b| {
        b.iter(|| fast_parse_u16(black_box(b"12345")))
    });
}

criterion_group!(
    benches,
    bench_chi_squared,
    bench_p_association,
    bench_cg_format,
    bench_bitset_popcount_small,
    bench_bitset_popcount_large,
    bench_bitset_popcount_1000,
    bench_fast_parse_u16,
);
criterion_main!(benches);
