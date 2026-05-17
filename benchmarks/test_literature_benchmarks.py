import tempfile
import unittest
from pathlib import Path

from benchmarks.run_literature_benchmarks import (
    count_data_rows,
    marker_count_from_table,
    parse_sra_experiment_xml,
    prune_dataset_results,
    summarize_output,
    write_dataset_metadata,
)


SRA_XML = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <EXPERIMENT alias="RAD_danio_sample_1" />
    <SAMPLE>
      <SAMPLE_ATTRIBUTES>
        <SAMPLE_ATTRIBUTE><TAG>sex</TAG><VALUE>male</VALUE></SAMPLE_ATTRIBUTE>
      </SAMPLE_ATTRIBUTES>
    </SAMPLE>
    <RUN_SET>
      <RUN accession="SRR000001" total_spots="12" total_bases="1200" size="321" />
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
  <EXPERIMENT_PACKAGE>
    <EXPERIMENT alias="medaka_F_1" />
    <SAMPLE>
      <SAMPLE_ATTRIBUTES>
        <SAMPLE_ATTRIBUTE><TAG>sex</TAG><VALUE>female</VALUE></SAMPLE_ATTRIBUTE>
      </SAMPLE_ATTRIBUTES>
    </SAMPLE>
    <RUN_SET>
      <RUN accession="SRR000002" total_spots="34" total_bases="3400" size="654" />
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""


class LiteratureBenchmarkTests(unittest.TestCase):
    def test_parse_sra_experiment_xml_extracts_sample_metadata(self):
        samples = parse_sra_experiment_xml(SRA_XML)

        self.assertEqual([sample.name for sample in samples], ["danio_sample_1", "medaka_F_1"])
        self.assertEqual([sample.accession for sample in samples], ["SRR000001", "SRR000002"])
        self.assertEqual([sample.sex for sample in samples], ["male", "female"])
        self.assertEqual(sum(sample.spots for sample in samples), 46)
        self.assertEqual(sum(sample.bases for sample in samples), 4600)

    def test_write_dataset_metadata_creates_popmap_and_download_accessions(self):
        samples = parse_sra_experiment_xml(SRA_XML)
        with tempfile.TemporaryDirectory() as tmp:
            dataset_dir = Path(tmp)
            write_dataset_metadata(dataset_dir, samples)

            self.assertEqual(
                (dataset_dir / "popmap.tsv").read_text(),
                "danio_sample_1\tmale\nmedaka_F_1\tfemale\n",
            )
            self.assertEqual(
                (dataset_dir / "samples.tsv").read_text().splitlines()[0],
                "sample\taccession\tsex\tspots\tbases\tbytes",
            )
            self.assertEqual((dataset_dir / ".download" / "danio_sample_1.accession").read_text(), "SRR000001")

    def test_marker_count_from_table_prefers_radsex_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "markers.tsv"
            path.write_text("#Number of markers : 123\nid\tsequence\tind1\n")

            self.assertEqual(marker_count_from_table(path), 123)

    def test_count_data_rows_ignores_comments_and_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "freq.tsv"
            path.write_text("#source:rsx-freq\nFrequency\tCount\n1\t4\n2\t5\n")

            self.assertEqual(count_data_rows(path), 2)

    def test_summarize_output_reads_command_specific_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            distrib = root / "distrib.tsv"
            distrib.write_text(
                "#source:rsx-distrib;min_depth:10;signif_threshold:0.05;bonferroni:true;n_markers:7\n"
                "male\tfemale\tMarkers\tP\tCorrectedP\tSignif\tBias\n"
                "1\t0\t3\t1\t1\tFalse\t1\n"
                "2\t0\t4\t0.1\t0.2\tFalse\t1\n"
            )
            signif = root / "signif.fa"
            signif.write_text(">m1\nACGT\n>m2\nTGCA\n")

            self.assertEqual(summarize_output("distrib", distrib)["markers"], "7")
            self.assertEqual(summarize_output("distrib", distrib)["rows"], "2")
            self.assertEqual(summarize_output("signif", signif)["significant_markers"], "2")

    def test_prune_dataset_results_removes_selected_dataset_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "results.csv"
            path.write_text(
                "dataset,command\n"
                "oryzias_latipes,metadata\n"
                "danio_choprae,metadata\n"
            )

            prune_dataset_results(path, {"oryzias_latipes"})

            self.assertEqual(path.read_text(), "dataset,command\n" "danio_choprae,metadata\n")


if __name__ == "__main__":
    unittest.main()
