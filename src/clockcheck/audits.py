"""The audits — the actual contribution of this project.

Each audit answers a question the original clock papers tend to leave out when a
clock is reused on a new dataset. Every function returns a plain dict/DataFrame so
the report layer can render it and so results are trivially diffable in git.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .clocks import LinearClock
from .data import MethylationDataset
from .utils import as_cpg_by_sample


def cpg_coverage(clock: LinearClock, ds: MethylationDataset) -> dict:
    """How much of the clock actually survives on this dataset?

    Naive coverage (fraction of CpGs present) hides the real risk: a clock can be
    "97% covered" yet be missing the single probe that carries most of its weight.
    We report weight-weighted coverage, which is the number that matters.
    """
    b = as_cpg_by_sample(ds.betas)
    absw = clock.weights.abs()
    present = absw.index.intersection(b.index)
    missing = absw.index.difference(present)
    total = float(absw.sum())
    covered = float(absw.loc[present].sum())
    heaviest_missing = absw.loc[missing].sort_values(ascending=False).head(10)
    return {
        "clock": clock.spec.name,
        "n_cpgs": int(len(absw)),
        "n_present": int(len(present)),
        "n_missing": int(len(missing)),
        "coverage_fraction": len(present) / len(absw),
        "weight_coverage_fraction": covered / total if total else float("nan"),
        "heaviest_missing_cpgs": heaviest_missing.to_dict(),
    }


def imputation_sensitivity(clock: LinearClock, ds: MethylationDataset,
                           strategies=("drop", "mean", "zero")) -> dict:
    """How far does a sample's predicted age move just by changing the (usually
    undisclosed) rule for handling absent CpGs?"""
    preds = {s: clock.predict(ds.betas, on_missing=s) for s in strategies}
    df = pd.DataFrame(preds)
    spread = df.max(axis=1) - df.min(axis=1)
    return {
        "clock": clock.spec.name,
        "strategies": list(strategies),
        "predictions": df,
        "per_sample_spread_years": spread,
        "median_spread_years": float(spread.median()),
        "p95_spread_years": float(spread.quantile(0.95)),
        "max_spread_years": float(spread.max()),
    }


def calibration_by_age(pred: pd.Series, chrono_age: pd.Series, n_bins: int = 5) -> dict:
    """Bias (predicted - chronological) as a function of age.

    A slope near -1 is the regression-to-the-mean signature: the clock drags old
    people young and young people old, which quietly inflates 'age acceleration'
    correlations at the cohort extremes.
    """
    s = pd.DataFrame({"pred": pred, "age": chrono_age}).dropna()
    bias = s["pred"] - s["age"]
    slope, intercept = np.polyfit(s["age"].to_numpy(), bias.to_numpy(), 1)
    edges = np.linspace(s["age"].min(), s["age"].max(), n_bins + 1)
    binned = (s.assign(bias=bias)
                .assign(bin=pd.cut(s["age"], edges, include_lowest=True))
                .groupby("bin", observed=True)["bias"].agg(["mean", "std", "count"]))
    return {
        "mae_years": float(bias.abs().mean()),
        "median_bias_years": float(bias.median()),
        "bias_vs_age_slope": float(slope),
        "per_bin_bias": binned,
    }


def cross_clock_agreement(preds: dict[str, pd.Series]) -> dict:
    """Do different clocks even agree on who is biologically old?

    Compares 'years'-unit clocks on the same samples. Large per-sample range or low
    Spearman correlation undercuts the idea of a single 'biological age'.
    """
    df = pd.DataFrame(preds).dropna()
    if df.shape[1] < 2:
        return {"note": "need >=2 comparable clocks", "n_clocks": int(df.shape[1])}
    rng = df.max(axis=1) - df.min(axis=1)
    return {
        "clocks": list(df.columns),
        "per_sample_range_years": rng,
        "median_range_years": float(rng.median()),
        "spearman_between_clocks": df.corr(method="spearman"),
    }


def age_acceleration_fragility(pred: pd.Series, chrono_age: pd.Series) -> dict:
    """'Age acceleration' is the residual everyone correlates with disease. Show how
    much that residual changes under two equally-defensible definitions."""
    s = pd.DataFrame({"pred": pred, "age": chrono_age}).dropna()
    coef = np.polyfit(s["age"].to_numpy(), s["pred"].to_numpy(), 1)
    resid = pd.DataFrame({
        "raw_diff": s["pred"] - s["age"],
        "regressed_out_age": s["pred"] - np.polyval(coef, s["age"].to_numpy()),
    })
    corr = float(resid["raw_diff"].corr(resid["regressed_out_age"]))
    return {
        "residuals": resid,
        "corr_between_definitions": corr,
        "note": ("Downstream disease associations rest on this residual. If two "
                 "reasonable definitions disagree, so will the associations built on them."),
    }
