"""Run published linear epigenetic clocks from their real coefficient files.

Design rule for this whole project: we never invent clock weights. Coefficients
are read from a published source on disk. By default we read the coefficient
CSVs bundled with `biolearn` (recommended: `pip install biolearn`), which packages
the weights straight from the original papers. You can also point at your own
directory of `CpG,weight` CSVs via `CLOCKCHECK_COEF_DIR`.

Only linear clocks are handled here on purpose: they run with numpy, need no GPU,
and are auditable line-by-line. Neural clocks (e.g. AltumAge) are out of scope for
v0.1 precisely because their failure modes need a different treatment.
"""
from __future__ import annotations
from dataclasses import dataclass
import importlib.util
import os

import numpy as np
import pandas as pd

from .utils import as_cpg_by_sample

ADULT_AGE = 20.0


def horvath_anti_transform(x):
    """Horvath's inverse age transformation (adult age = 20)."""
    x = np.asarray(x, dtype=float)
    return np.where(x < 0, (1 + ADULT_AGE) * np.exp(x) - 1, (1 + ADULT_AGE) * x + ADULT_AGE)


@dataclass(frozen=True)
class ClockSpec:
    name: str
    coef_file: str                 # filename in the coefficient source dir
    intercept: float = 0.0         # used only if the file has no intercept row
    transform: str = "identity"    # "identity" | "horvath"
    unit: str = "years"            # "years" | "pace" (pace clocks are not ages)
    citation: str = ""


# A deliberately small, well-known starter set. Extend by dropping specs here;
# the weights come from whatever coefficient dir is active.
BUILTIN: dict[str, ClockSpec] = {
    "Horvath1": ClockSpec("Horvath1", "Horvath1.csv", 0.695507258, "horvath", "years",
                          "Horvath 2013, Genome Biology"),
    "Hannum": ClockSpec("Hannum", "Hannum.csv", 0.0, "identity", "years",
                        "Hannum et al. 2013, Molecular Cell"),
    "PhenoAge": ClockSpec("PhenoAge", "PhenoAge.csv", 60.664, "identity", "years",
                          "Levine et al. 2018, Aging"),
    "DunedinPACE": ClockSpec("DunedinPACE", "DunedinPACE.csv", -1.949859, "identity", "pace",
                             "Belsky et al. 2022, eLife"),
}


def default_coef_dir() -> str:
    """Prefer biolearn's bundled published coefficient files, else env override."""
    env = os.environ.get("CLOCKCHECK_COEF_DIR")
    if env and os.path.isdir(env):
        return env
    spec = importlib.util.find_spec("biolearn")
    if spec and spec.submodule_search_locations:
        d = os.path.join(list(spec.submodule_search_locations)[0], "data")
        if os.path.isdir(d):
            return d
    return ""


class LinearClock:
    """A published linear clock: predicted value = weights . betas + intercept."""

    def __init__(self, spec: ClockSpec, coef_dir: str | None = None):
        self.spec = spec
        coef_dir = coef_dir or default_coef_dir()
        if not coef_dir:
            raise FileNotFoundError(
                "No coefficient source found. `pip install biolearn` or set "
                "CLOCKCHECK_COEF_DIR to a directory of published 'CpG,weight' CSVs."
            )
        path = os.path.join(coef_dir, spec.coef_file)
        df = pd.read_csv(path)
        df.columns = [c.strip() for c in df.columns]
        cpg_col, w_col = df.columns[0], df.columns[1]
        is_intercept = df[cpg_col].astype(str).str.lower().str.contains("intercept")
        self.intercept = spec.intercept
        if is_intercept.any():
            self.intercept = float(df.loc[is_intercept, w_col].iloc[0])
            df = df[~is_intercept]
        self.weights = pd.Series(
            df[w_col].to_numpy(dtype=float),
            index=df[cpg_col].astype(str).to_numpy(),
            name=spec.name,
        )

    @classmethod
    def load(cls, name: str, coef_dir: str | None = None) -> "LinearClock":
        if name not in BUILTIN:
            raise KeyError(f"Unknown clock {name!r}. Known: {list(BUILTIN)}")
        return cls(BUILTIN[name], coef_dir=coef_dir)

    @property
    def cpgs(self) -> list[str]:
        return list(self.weights.index)

    def predict(self, betas: pd.DataFrame, on_missing: str = "drop",
                reference_means: pd.Series | None = None) -> pd.Series:
        """Predict per sample.

        `on_missing` controls what happens to clock CpGs absent from the data —
        the choice papers rarely disclose:
          - "drop":      ignore absent CpGs (silently shrinks the model)
          - "mean":      fill absent CpGs with the dataset-wide mean beta
          - "zero":      fill absent CpGs with 0 (assume unmethylated)
          - "gold_mean": fill from an external reference-mean Series if given
        """
        b = as_cpg_by_sample(betas).astype(float)
        # within-CpG NaNs -> that CpG's row mean (the standard quiet default)
        b = b.apply(lambda r: r.fillna(r.mean()), axis=1)
        w = self.weights
        present = w.index.intersection(b.index)
        missing = w.index.difference(present)

        if on_missing == "drop" or len(missing) == 0:
            mat = b.loc[present]
            raw = w.loc[present].to_numpy() @ mat.to_numpy() + self.intercept
        else:
            mat = pd.DataFrame(0.0, index=w.index, columns=b.columns)
            mat.loc[present] = b.loc[present].to_numpy()
            if on_missing == "mean":
                mat.loc[missing] = float(np.nanmean(b.to_numpy()))
            elif on_missing == "zero":
                mat.loc[missing] = 0.0
            elif on_missing == "gold_mean":
                fallback = float(np.nanmean(b.to_numpy()))
                ref = reference_means if reference_means is not None else pd.Series(dtype=float)
                for cpg in missing:
                    mat.loc[cpg] = float(ref.get(cpg, fallback))
            else:
                raise ValueError(f"unknown on_missing={on_missing!r}")
            raw = w.to_numpy() @ mat.to_numpy() + self.intercept

        raw = np.asarray(raw, dtype=float).ravel()
        if self.spec.transform == "horvath":
            raw = horvath_anti_transform(raw)
        return pd.Series(raw, index=b.columns, name=self.spec.name)
