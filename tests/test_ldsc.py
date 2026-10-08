import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import ldsc  # noqa: E402


def simulate(h2, intercept, n_snps=200_000, n=100_000, m=1_000_000, seed=1):
    rng = np.random.default_rng(seed)
    ld = rng.gamma(shape=2.0, scale=40.0, size=n_snps) + 1.0
    scale = intercept + n * h2 * ld / m
    chisq = scale * rng.chisquare(1, size=n_snps)
    return chisq, ld, np.full(n_snps, float(n)), float(m)


def test_recovers_heritability_and_intercept():
    chisq, ld, n, m = simulate(h2=0.30, intercept=1.0)
    res = ldsc.ld_score_regression(chisq, ld, ld, n, m)
    assert abs(res["h2_obs"] - 0.30) < 3 * res["h2_se"] + 0.01
    assert abs(res["intercept"] - 1.0) < 3 * res["intercept_se"] + 0.01


def test_picks_up_inflation_in_the_intercept():
    chisq, ld, n, m = simulate(h2=0.10, intercept=1.15, seed=2)
    res = ldsc.ld_score_regression(chisq, ld, ld, n, m)
    assert 1.10 < res["intercept"] < 1.20
    assert abs(res["h2_obs"] - 0.10) < 0.03


def test_no_signal_gives_h2_near_zero():
    chisq, ld, n, m = simulate(h2=0.0, intercept=1.0, seed=3)
    res = ldsc.ld_score_regression(chisq, ld, ld, n, m)
    assert abs(res["h2_obs"]) < 3 * res["h2_se"] + 0.01


def test_liability_scale_matches_hand_calculation():
    # sample prevalence equals population prevalence -> factor is K(1-K)/z^2
    from scipy.stats import norm
    k = 0.01
    h2, _ = ldsc.liability_scale(0.1, 0.01, cases=1000, controls=99000, prevalence=k)
    z = norm.pdf(norm.ppf(1 - k))
    assert abs(h2 - 0.1 * k * (1 - k) / z ** 2) < 1e-9
