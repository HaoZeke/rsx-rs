import tempfile
import unittest
from pathlib import Path

from benchmarks.plot_cuda_benchmarks import load_benchmark, summarize_benchmark


class CudaBenchmarkPlotTests(unittest.TestCase):
    def test_repeated_headers_and_repetitions_are_summarized(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cuda.csv"
            header = (
                "markers,cpu_total_s,cuda_setup_s,cuda_h2d_s,cuda_kernel_s,cuda_d2h_s,"
                "cuda_total_s,h2d_bytes,d2h_bytes,h2d_gb_s,d2h_gb_s,kernel_speedup,"
                "total_speedup,output_buffer_reused,max_abs_error,device\n"
            )
            path.write_text(
                " WARN cache was redirected for this run\n"
                + "# repetition=1\n"
                + header
                + "100000,0.10,0.01,0.002,0.005,0.002,0.02,800000,800000,1,1,20,5,1,1e-16,A100\n"
                + "# repetition=2\n"
                + header
                + "100000,0.12,0.01,0.002,0.004,0.002,0.01,800000,800000,1,1,30,12,1,2e-16,A100\n"
            )

            raw = load_benchmark(path)
            summary = summarize_benchmark(raw)

        self.assertEqual(len(raw), 2)
        self.assertEqual(len(summary), 1)
        self.assertEqual(summary.iloc[0]["markers"], 100000)
        self.assertAlmostEqual(summary.iloc[0]["kernel_speedup"], 25.0)
        self.assertAlmostEqual(summary.iloc[0]["total_speedup"], 8.5)
        self.assertEqual(summary.iloc[0]["repetitions"], 2)
        self.assertEqual(summary.iloc[0]["device"], "A100")

    def test_numerical_disagreement_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cuda.csv"
            path.write_text(
                "markers,cpu_total_s,cuda_setup_s,cuda_h2d_s,cuda_kernel_s,cuda_d2h_s,"
                "cuda_total_s,h2d_bytes,d2h_bytes,h2d_gb_s,d2h_gb_s,kernel_speedup,"
                "total_speedup,output_buffer_reused,max_abs_error,device\n"
                "100,0.1,0.01,0.01,0.01,0.01,0.04,800,800,1,1,10,2.5,1,3e-15,A100\n"
            )

            with self.assertRaisesRegex(ValueError, "numerical agreement"):
                summarize_benchmark(load_benchmark(path))


if __name__ == "__main__":
    unittest.main()
