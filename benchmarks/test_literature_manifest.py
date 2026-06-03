import tempfile
import unittest
from pathlib import Path

from benchmarks.check_literature_manifest import load_rows


class LiteratureManifestTests(unittest.TestCase):
    def test_load_rows_skips_comment_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "literature_datasets.tsv"
            manifest.write_text(
                "dataset\tsource\taccession\tgenome\tbenchmark_role\tcommands\tnotes\n"
                "danio_albolineatus\tRADSex paper workflow\tPRJNA548074\tNA\tliterature panel\tprocess,depth\trow\n"
                "\n"
                "# note for operators\n"
                "# hippoglossus_stenolepis\tvalidation\tTBD\tNA\tvalidation\tprocess\tcommented\n"
            )

            rows = load_rows(manifest)

            self.assertEqual([row["dataset"] for row in rows], ["danio_albolineatus"])


if __name__ == "__main__":
    unittest.main()
