"""Tests run on the labeled synthetic fixture. They check the *machinery*, never
clock accuracy (synthetic data can't tell you that — that's the whole ethos)."""
import warnings

import numpy as np
import pytest

from clockcheck import LinearClock, build_audit_card, card_to_markdown, synthetic_fixture
from clockcheck import audits


@pytest.fixture(scope="module")
def clocks():
    try:
        return [LinearClock.load("Horvath1"), LinearClock.load("Hannum")]
    except FileNotFoundError:
        pytest.skip("no coefficient source (install biolearn or set CLOCKCHECK_COEF_DIR)")


@pytest.fixture(scope="module")
def ds(clocks):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return synthetic_fixture(clocks, n_samples=40, seed=1)


def test_fixture_is_labeled(ds):
    assert ds.is_synthetic is True


def test_predict_shape(clocks, ds):
    pred = clocks[0].predict(ds.betas)
    assert len(pred) == ds.betas.shape[1]
    assert pred.notna().all()


def test_full_coverage_on_fixture(clocks, ds):
    # fixture is built from the union of clock CpGs -> everything present
    cov = audits.cpg_coverage(clocks[0], ds)
    assert cov["n_missing"] == 0
    assert cov["coverage_fraction"] == pytest.approx(1.0)
    assert cov["weight_coverage_fraction"] == pytest.approx(1.0)


def test_weight_coverage_drops_when_heavy_probe_removed(clocks, ds):
    clk = clocks[0]
    heavy = clk.weights.abs().idxmax()
    reduced = ds.betas.drop(index=heavy)
    from clockcheck.data import MethylationDataset
    ds2 = MethylationDataset(reduced, ds.meta, is_synthetic=True, name="reduced")
    cov = audits.cpg_coverage(clk, ds2)
    # only one probe gone, but it carries real weight -> weight coverage < naive coverage
    assert cov["n_missing"] == 1
    assert cov["weight_coverage_fraction"] < cov["coverage_fraction"]


def test_imputation_sensitivity_nonnegative(clocks, ds):
    imp = audits.imputation_sensitivity(clocks[0], ds)
    assert imp["median_spread_years"] >= 0
    assert (imp["per_sample_spread_years"] >= -1e-9).all()


def test_calibration_keys(clocks, ds):
    pred = clocks[0].predict(ds.betas)
    cal = audits.calibration_by_age(pred, ds.age)
    assert {"mae_years", "median_bias_years", "bias_vs_age_slope", "per_bin_bias"} <= set(cal)
    assert np.isfinite(cal["bias_vs_age_slope"])


def test_card_renders_and_stamps_synthetic(clocks, ds):
    card = build_audit_card(clocks[0], ds, compare_clocks=clocks[1:])
    md = card_to_markdown(card)
    assert "NON-EVIDENTIAL" in md
    assert card["is_synthetic"] is True
