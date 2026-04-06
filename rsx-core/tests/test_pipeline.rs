// GPL-3.0-or-later
// Integration tests: end-to-end pipeline validation

use std::io::Write;
use std::path::PathBuf;

fn test_dir() -> PathBuf {
    let dir = std::env::temp_dir().join("radsex_rs_test");
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn create_test_markers_table(dir: &std::path::Path) -> PathBuf {
    let path = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    // 5 individuals, 4 markers
    writeln!(f, "#Number of markers : 4").unwrap();
    writeln!(f, "id\tsequence\tind1\tind2\tind3\tind4\tind5").unwrap();
    // Marker 0: present in all (depths 10, 5, 8, 12, 7)
    writeln!(f, "0\tATCGATCG\t10\t5\t8\t12\t7").unwrap();
    // Marker 1: present in M only (ind1, ind2, ind3 are M)
    writeln!(f, "1\tGGGGAAAA\t15\t20\t10\t0\t0").unwrap();
    // Marker 2: present in F only (ind4, ind5 are F)
    writeln!(f, "2\tCCCCTTTT\t0\t0\t0\t25\t30").unwrap();
    // Marker 3: present in most
    writeln!(f, "3\tAAAATTTT\t5\t0\t3\t8\t6").unwrap();
    path
}

fn create_test_popmap(dir: &std::path::Path) -> PathBuf {
    let path = dir.join("popmap.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "ind1\tM").unwrap();
    writeln!(f, "ind2\tM").unwrap();
    writeln!(f, "ind3\tM").unwrap();
    writeln!(f, "ind4\tF").unwrap();
    writeln!(f, "ind5\tF").unwrap();
    path
}

fn create_test_fastq_dir(dir: &std::path::Path) -> PathBuf {
    let input_dir = dir.join("reads");
    std::fs::create_dir_all(&input_dir).unwrap();

    // Create a small FASTQ file for each individual
    for name in &["ind1", "ind2", "ind3", "ind4", "ind5"] {
        let path = input_dir.join(format!("{name}.fq"));
        let mut f = std::fs::File::create(&path).unwrap();
        // Each individual gets 2 reads of the same sequence (depth 2)
        writeln!(f, "@read1").unwrap();
        writeln!(f, "ATCGATCGATCG").unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "IIIIIIIIIIII").unwrap();
        writeln!(f, "@read2").unwrap();
        writeln!(f, "ATCGATCGATCG").unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "IIIIIIIIIIII").unwrap();
        // ind1 and ind2 get an extra unique sequence
        if *name == "ind1" || *name == "ind2" {
            writeln!(f, "@read3").unwrap();
            writeln!(f, "GGGGAAAA").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "IIIIIIII").unwrap();
        }
    }

    input_dir
}

fn create_test_genome(dir: &std::path::Path) -> PathBuf {
    let path = dir.join("genome.fa");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, ">chr1").unwrap();
    writeln!(f, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();
    writeln!(f, "GGGGAAAACCCCTTTTAAAATTTTCCCCGGGGTTTTAAAA").unwrap();
    writeln!(f, ">chr2").unwrap();
    writeln!(f, "TTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCC").unwrap();
    path
}

#[test]
fn test_freq_command() {
    let dir = test_dir().join("freq");
    std::fs::create_dir_all(&dir).unwrap();
    let table = create_test_markers_table(&dir);
    let output = dir.join("freq_output.tsv");

    rsx_core::commands::freq::run(&rsx_core::commands::freq::FreqParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("#source:rsx-freq"));
    assert!(content.contains("Frequency\tCount"));
    // Marker 1 is present in 3 individuals, marker 2 in 2, marker 3 in 4, marker 0 in 5
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() > 2); // header + comment + data lines
}

#[test]
fn test_depth_command() {
    let dir = test_dir().join("depth");
    std::fs::create_dir_all(&dir).unwrap();
    let table = create_test_markers_table(&dir);
    let popmap = create_test_popmap(&dir);
    let output = dir.join("depth_output.tsv");

    rsx_core::commands::depth::run(&rsx_core::commands::depth::DepthParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_frequency: 0.5,
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("Sample\tGroup\tReads\tMarkers"));
    // Check that all 5 individuals appear
    assert!(content.contains("ind1\tM"));
    assert!(content.contains("ind5\tF"));
}

#[test]
fn test_distrib_command() {
    let dir = test_dir().join("distrib");
    std::fs::create_dir_all(&dir).unwrap();
    let table = create_test_markers_table(&dir);
    let popmap = create_test_popmap(&dir);
    let output = dir.join("distrib_output.tsv");

    rsx_core::commands::distrib::run(&rsx_core::commands::distrib::DistribParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        disable_correction: false,
        group1: String::new(),
        group2: String::new(),
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("#source:rsx-distrib"));
    assert!(content.contains("Markers\tP\tCorrectedP\tSignif\tBias"));
    // The sex-specific markers should show up in the distribution
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() > 2);
}

#[test]
fn test_signif_command() {
    let dir = test_dir().join("signif");
    std::fs::create_dir_all(&dir).unwrap();
    let table = create_test_markers_table(&dir);
    let popmap = create_test_popmap(&dir);
    let output = dir.join("signif_output.tsv");

    rsx_core::commands::signif::run(&rsx_core::commands::signif::SignifParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        disable_correction: true, // disable correction so small test data can have significant markers
        output_fasta: false,
        group1: String::new(),
        group2: String::new(),
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("#source:rsx-signif"));
}

#[test]
fn test_signif_fasta_output() {
    let dir = test_dir().join("signif_fasta");
    std::fs::create_dir_all(&dir).unwrap();
    let table = create_test_markers_table(&dir);
    let popmap = create_test_popmap(&dir);
    let output = dir.join("signif_output.fa");

    rsx_core::commands::signif::run(&rsx_core::commands::signif::SignifParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        disable_correction: true,
        output_fasta: true,
        group1: String::new(),
        group2: String::new(),
    })
    .unwrap();

    // FASTA output should exist (may be empty if no markers pass threshold)
    assert!(output.exists());
}

#[test]
fn test_subset_command() {
    let dir = test_dir().join("subset");
    std::fs::create_dir_all(&dir).unwrap();
    let table = create_test_markers_table(&dir);
    let popmap = create_test_popmap(&dir);
    let output = dir.join("subset_output.tsv");

    rsx_core::commands::subset::run(&rsx_core::commands::subset::SubsetParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        disable_correction: true,
        output_fasta: false,
        group1: String::new(),
        group2: String::new(),
        min_group1: 2,
        min_group2: 0,
        max_group1: 3,
        max_group2: 0,
        min_individuals: 2,
        max_individuals: 5,
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("#source:rsx-subset"));
}

#[test]
fn test_process_command() {
    let dir = test_dir().join("process");
    std::fs::create_dir_all(&dir).unwrap();
    let input_dir = create_test_fastq_dir(&dir);
    let output = dir.join("markers_table.tsv");

    rsx_core::commands::process::run(&rsx_core::commands::process::ProcessParams {
        input_dir_path: input_dir.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        n_threads: 2,
        min_depth: 1,
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("#Number of markers"));
    assert!(content.contains("id\tsequence"));
    // Should have the common sequence ATCGATCGATCG in all 5 individuals
    assert!(content.contains("ATCGATCGATCG"));
    // Should have the M-only sequence GGGGAAAA
    assert!(content.contains("GGGGAAAA"));
}

#[test]
fn test_map_command() {
    let dir = test_dir().join("map");
    std::fs::create_dir_all(&dir).unwrap();
    let table = create_test_markers_table(&dir);
    let popmap = create_test_popmap(&dir);
    let genome = create_test_genome(&dir);
    let output = dir.join("map_output.tsv");

    rsx_core::commands::map::run(&rsx_core::commands::map::MapParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        genome_file_path: genome.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        min_quality: 0, // low quality for test data (short sequences)
        min_frequency: 0.1,
        signif_threshold: 0.05,
        disable_correction: false,
        group1: String::new(),
        group2: String::new(),
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("#source:rsx-map"));
    assert!(content.contains("Contig\tPosition\tLength\tMarker_id\tBias\tP\tCorrectedP\tSignif"));
}

#[test]
fn test_full_pipeline() {
    // End-to-end: process -> freq -> distrib -> signif
    let dir = test_dir().join("full_pipeline");
    std::fs::create_dir_all(&dir).unwrap();

    // Step 1: Process reads into markers table
    let input_dir = create_test_fastq_dir(&dir);
    let markers_table = dir.join("markers.tsv");

    rsx_core::commands::process::run(&rsx_core::commands::process::ProcessParams {
        input_dir_path: input_dir.to_str().unwrap().to_string(),
        output_file_path: markers_table.to_str().unwrap().to_string(),
        n_threads: 2,
        min_depth: 1,
    })
    .unwrap();
    assert!(markers_table.exists());

    // Create popmap for our test individuals
    let popmap = create_test_popmap(&dir);

    // Step 2: Freq
    let freq_out = dir.join("freq.tsv");
    rsx_core::commands::freq::run(&rsx_core::commands::freq::FreqParams {
        markers_table_path: markers_table.to_str().unwrap().to_string(),
        output_file_path: freq_out.to_str().unwrap().to_string(),
        min_depth: 1,
    })
    .unwrap();
    assert!(freq_out.exists());

    // Step 3: Distrib
    let distrib_out = dir.join("distrib.tsv");
    rsx_core::commands::distrib::run(&rsx_core::commands::distrib::DistribParams {
        markers_table_path: markers_table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: distrib_out.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        disable_correction: false,
        group1: String::new(),
        group2: String::new(),
    })
    .unwrap();
    assert!(distrib_out.exists());

    // Step 4: Signif
    let signif_out = dir.join("signif.tsv");
    rsx_core::commands::signif::run(&rsx_core::commands::signif::SignifParams {
        markers_table_path: markers_table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: signif_out.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        disable_correction: true,
        output_fasta: false,
        group1: String::new(),
        group2: String::new(),
    })
    .unwrap();
    assert!(signif_out.exists());

    // Verify the processed markers table can be read back
    let content = std::fs::read_to_string(&markers_table).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() >= 3); // comment + header + at least 1 marker
}
