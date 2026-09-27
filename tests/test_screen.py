import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from receptor_screen.core import Settings, screen


class ScreenTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(91)
        n = 480
        r = rng.integers(0, 12, n)
        self.c = np.column_stack([r, rng.poisson(1 + 4*(r > 5), n),
                                  rng.poisson(3, n), rng.poisson(8, n),
                                  rng.poisson(2000, n)])
        self.genes = ["R", "signal", "null", "MT-X", "background"]
        self.obs = pd.DataFrame({"donor": np.repeat(np.arange(6), 80),
                                 "batch": np.repeat(["a", "b", "c"], 160)})
        self.s = Settings(folds=3, min_detected=10, max_iter=35, min_samples_leaf=15)

    def run_screen(self, c):
        return screen(c, self.genes, self.obs, ["R"], "donor",
                      [False, False, False, True, False], ["batch"], self.s)

    def test_signal_and_donor_separation(self):
        scores, donors, folds = self.run_screen(sparse.csc_matrix(self.c))
        scores = scores.set_index("gene")
        self.assertGreater(scores.loc["signal", "relative_mse_reduction"], 0.2)
        self.assertGreater(scores.loc["signal", "relative_mse_reduction"],
                           scores.loc["null", "relative_mse_reduction"])
        self.assertTrue((folds.groupby("donor").fold.nunique() == 1).all())
        self.assertEqual(len(donors), 24)
        self.assertNotIn("R", scores.index)

    def test_sparse_dense_agree(self):
        a, _, _ = self.run_screen(self.c)
        b, _, _ = self.run_screen(sparse.csr_matrix(self.c))
        pd.testing.assert_frame_equal(a, b)

    def test_invalid_counts(self):
        c = self.c.astype(float)
        c[0, 0] = -1
        with self.assertRaisesRegex(ValueError, "raw integer"):
            self.run_screen(c)

    def test_insufficient_donors(self):
        self.obs["donor"] = "one"
        with self.assertRaisesRegex(ValueError, "number of folds"):
            self.run_screen(self.c)


if __name__ == "__main__":
    unittest.main()
