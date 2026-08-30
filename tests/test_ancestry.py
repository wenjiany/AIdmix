import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pysam

from aidmix.ancestry import (
    Marker,
    SiteObs,
    bootstrap,
    choose_fit_result,
    fit,
    read_log_likelihood,
    read_panel,
    score_ancestry_arrays,
    score_arrays,
    softmax,
)


class TestAIMAncestry(unittest.TestCase):
    def test_softmax_sums_to_one(self):
        self.assertAlmostEqual(float(softmax(np.array([1., 2., 3., 4.])).sum()), 1.0)

    def test_read_error_model_prefers_matching_genotype(self):
        self.assertGreater(read_log_likelihood("A", "G", "G", 30, 2),
                           read_log_likelihood("A", "G", "G", 30, 0))

    def test_fit_recovers_population_signal(self):
        markers = [Marker("chr1", i, "A", "G", (0.95, 0.05, 0.05, 0.05, 0.05)) for i in range(1, 30)]
        # Alt-homozygous observations should favor the population with high ALT AF.
        obs = [SiteObs(m, "background", (-16, -8, 0), 1) for m in markers]
        ancestry, likelihood, _converged = fit(obs, 0.25, "background")
        self.assertTrue(math.isfinite(likelihood))
        self.assertGreater(float(ancestry[0]), 0.5)

    def test_vectorized_score_analytic_gradient(self):
        frequencies = np.asarray([[0.1, 0.7, 0.3], [0.8, 0.2, 0.5]])
        log_likelihoods = np.asarray([[-3.0, -0.5, -2.0], [-4.0, -1.0, -0.2]])
        weights = np.asarray([1.0, 0.25])
        x = np.asarray([0.2, -0.4])
        _value, gradient = score_arrays(x, frequencies, log_likelihoods, weights)
        epsilon = 1e-6
        numeric = np.empty_like(x)
        for index in range(len(x)):
            delta = np.zeros_like(x)
            delta[index] = epsilon
            high = score_arrays(x + delta, frequencies, log_likelihoods, weights)[0]
            low = score_arrays(x - delta, frequencies, log_likelihoods, weights)[0]
            numeric[index] = (high - low) / (2 * epsilon)
        np.testing.assert_allclose(gradient, numeric, rtol=1e-5, atol=1e-7)

    def test_clipped_frequency_gradient_is_finite(self):
        frequencies = np.asarray([[0.0, 1.0], [0.0, 1.0]])
        log_likelihoods = np.asarray([[0.0, -1.0, -2.0], [0.0, -1.0, -2.0]])
        value, gradient = score_ancestry_arrays(
            np.asarray([1.0, 0.0]), frequencies, log_likelihoods, np.ones(2)
        )
        self.assertTrue(math.isfinite(value))
        self.assertTrue(np.all(np.isfinite(gradient)))

    def test_ancestry_score_analytic_gradient(self):
        frequencies = np.asarray([[0.1, 0.7, 0.3], [0.8, 0.2, 0.5]])
        log_likelihoods = np.asarray([[-3.0, -0.5, -2.0], [-4.0, -1.0, -0.2]])
        weights = np.asarray([1.0, 0.25])
        ancestry = np.asarray([0.2, 0.3, 0.5])
        _value, gradient = score_ancestry_arrays(
            ancestry, frequencies, log_likelihoods, weights
        )
        epsilon = 1e-6
        numeric = np.empty_like(ancestry)
        for index in range(len(ancestry)):
            delta = np.zeros_like(ancestry)
            delta[index] = epsilon
            high = score_ancestry_arrays(
                ancestry + delta, frequencies, log_likelihoods, weights
            )[0]
            low = score_ancestry_arrays(
                ancestry - delta, frequencies, log_likelihoods, weights
            )[0]
            numeric[index] = (high - low) / (2 * epsilon)
        np.testing.assert_allclose(gradient, numeric, rtol=1e-5, atol=1e-7)

    def test_fit_selection_prefers_equivalent_converged_solution(self):
        failed = SimpleNamespace(fun=100.0, success=False)
        converged = SimpleNamespace(fun=100.0 + 5e-7, success=True)
        self.assertIs(choose_fit_result([failed, converged]), converged)
        clearly_worse = SimpleNamespace(fun=100.01, success=True)
        self.assertIs(choose_fit_result([failed, clearly_worse]), failed)

    def test_panel_population_columns_are_dynamic(self):
        with tempfile.TemporaryDirectory() as directory:
            panel = Path(directory) / "panel.tsv"
            panel.write_text("chrom\tpos\tref\talt\tPOP_A\tPOP_B\nchr1\t10\tA\tG\t0.9\t0.1\n")
            markers, populations = read_panel(str(panel), 0.005, 0.995)
        self.assertEqual(populations, ("POP_A", "POP_B"))
        self.assertEqual(len(markers), 1)

    def test_native_iadmix_panel_is_oriented_to_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 20 + "\n")
            pysam.faidx(str(reference))
            panel = Path(directory) / "panel.ref.txt"
            panel.write_text(
                "#chrom position rsid A1 A2 POP_A POP_B\n"
                "chr1 10 rs1 A G 0.8 0.2\n"
            )
            markers, populations = read_panel(str(panel), 0.005, 0.995, str(reference))
        self.assertEqual(populations, ("POP_A", "POP_B"))
        self.assertEqual((markers[0].ref, markers[0].alt), ("A", "G"))
        np.testing.assert_allclose(markers[0].freq, (0.2, 0.8))

    def test_normalized_panel_is_oriented_to_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 20 + "\n")
            pysam.faidx(str(reference))
            panel = Path(directory) / "panel.tsv"
            panel.write_text(
                "chrom\tpos\tref\talt\tPOP_A\tPOP_B\n"
                "chr1\t10\tG\tA\t0.8\t0.2\n"
            )
            markers, populations = read_panel(
                str(panel), 0.005, 0.995, str(reference)
            )
        self.assertEqual(populations, ("POP_A", "POP_B"))
        self.assertEqual((markers[0].ref, markers[0].alt), ("A", "G"))
        np.testing.assert_allclose(markers[0].freq, (0.2, 0.8))

    def test_target_and_background_are_separate(self):
        m = Marker("chr1", 1, "A", "G", (0.5, 0.5, 0.5, 0.5, 0.5))
        obs = [SiteObs(m, "target", (0, 0, 0), 1), SiteObs(m, "background", (0, 0, 0), 1)]
        self.assertEqual(sum(x.region == "target" for x in obs), 1)
        self.assertEqual(sum(x.region == "background" for x in obs), 1)

    def test_bootstrap_empty_result_uses_dynamic_population_count(self):
        m = Marker("chr1", 1, "A", "G", (0.5, 0.5, 0.5))
        obs = [SiteObs(m, "background", (0, 0, 0), 1)]
        result = bootstrap(obs, 0.25, "background", 0, 1)
        self.assertEqual(result.shape, (0, 3))

    def test_fit_empty_observations_preserves_population_count(self):
        ancestry, likelihood, converged = fit([], 0.25, "target", n_pops=3)
        self.assertEqual(ancestry.shape, (3,))
        self.assertTrue(np.all(np.isnan(ancestry)))
        self.assertTrue(math.isnan(likelihood))
        self.assertFalse(converged)


if __name__ == "__main__":
    unittest.main()
