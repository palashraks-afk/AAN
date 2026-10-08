import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import sldsc_lite  # noqa: E402


def test_merge_windows_joins_overlaps():
    s, e = sldsc_lite.merge_windows([10, 5, 40], [20, 12, 50])
    assert list(s) == [5, 40] and list(e) == [20, 50]


def test_in_windows_marks_points_inside():
    starts, ends = np.array([100, 300]), np.array([200, 400])
    pos = np.array([50, 100, 150, 250, 350, 450])
    assert list(sldsc_lite.in_windows(pos, starts, ends)) == [False, True, True, False, True, False]


def test_in_windows_empty():
    assert not sldsc_lite.in_windows(np.array([1, 2, 3]), np.array([]), np.array([])).any()


def test_genotype_lookup_decodes_plink_codes():
    # byte 0b11100100: individuals read right to left as 00,01,10,11 -> 2, missing, 1, 0 copies of A1
    row = sldsc_lite.GENO_LUT[0b11100100]
    assert row[0] == 2 and np.isnan(row[1]) and row[2] == 1 and row[3] == 0


def test_block_wls_recovers_a_known_coefficient():
    rng = np.random.default_rng(3)
    n = 40_000
    x = np.column_stack([rng.gamma(2, 10, n), rng.gamma(2, 5, n), np.ones(n)])
    y = x @ np.array([0.02, 0.05, 1.0]) + rng.normal(0, 1, n)
    est, se = sldsc_lite.block_wls(x, y, np.ones(n) / n, n_blocks=50)
    assert np.all(np.abs(est - np.array([0.02, 0.05, 1.0])) < 4 * se)
