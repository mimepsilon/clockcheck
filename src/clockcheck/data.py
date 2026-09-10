"""Datasets in, honestly labeled.

The whole point of this project is to stop confusing synthetic plumbing with
evidence. So there is exactly one synthetic generator here, it is named loudly,
it emits a warning every time it runs, and nothing in the audit layer will let a
synthetic dataset masquerade as a validation result.
"""
from __future__ import annotations
from dataclasses import dataclass
import warnings

import numpy as np
import pandas as pd

from .utils import as_cpg_by_sample


@dataclass
class MethylationDataset:
    """CpG x sample beta matrix plus per-sample metadata (must include age)."""
    betas: pd.DataFrame              # index = CpG, columns = sample
    meta: pd.DataFrame               # index = sample
    age_col: str = "age"
    name: str = "dataset"
    is_synthetic: bool = False       # audits refuse to emit a validation verdict if True

    def __post_init__(self):
        self.betas = as_cpg_by_sample(self.betas)
        common = self.betas.columns.intersection(self.meta.index)
        if len(common) == 0:
            raise ValueError("No overlap between beta columns (samples) and meta index.")
        self.betas = self.betas[common]
        self.meta = self.meta.loc[common]
        if self.age_col not in self.meta.columns:
            raise ValueError(f"meta has no age column {self.age_col!r}; columns={list(self.meta.columns)}")

    @property
    def age(self) -> pd.Series:
        return self.meta[self.age_col].astype(float)


def load_csv(betas_path: str, meta_path: str, age_col: str = "age",
             name: str | None = None) -> MethylationDataset:
    """Load a real dataset from two CSVs: a beta matrix and a sample metadata table."""
    betas = pd.read_csv(betas_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)
    return MethylationDataset(betas, meta, age_col=age_col, name=name or betas_path)


def load_geo(gse_id: str, age_col: str = "age") -> MethylationDataset:
    """Load a real GEO methylation series via biolearn's loaders.

    This requires network access and `biolearn`. It is the recommended path for
    actual audits; see biolearn's GEO data library for series with age metadata.
    """
    from biolearn.data_library import DataLibrary  # imported lazily on purpose
    data = DataLibrary().get(gse_id).load()
    betas = data.dnam                                   # CpG x sample
    meta = data.metadata                                # sample x fields
    return MethylationDataset(betas, meta, age_col=age_col, name=gse_id)


def synthetic_fixture(clocks, n_samples: int = 60, seed: int = 0):
    """PLUMBING ONLY. Not data. Not evidence. Never cite an audit run on this.

    Builds a random beta matrix over the union of the given clocks' CpGs, with a
    faint age signal so the pipeline produces non-degenerate numbers. It exists so
    tests and the demo run offline with zero downloads. `is_synthetic=True` makes
    the report layer stamp every card as NON-EVIDENTIAL.
    """
    warnings.warn(
        "synthetic_fixture() produced SYNTHETIC betas for plumbing/CI only. "
        "Any age prediction from it is meaningless. Use load_geo/load_csv for real audits.",
        stacklevel=2,
    )
    rng = np.random.default_rng(seed)
    cpgs = sorted({c for clk in clocks for c in clk.cpgs})
    ages = rng.uniform(20, 85, size=n_samples)
    samples = [f"S{i:03d}" for i in range(n_samples)]
    # baseline betas + a small per-CpG age slope + noise, clipped to [0,1]
    base = rng.uniform(0.1, 0.9, size=(len(cpgs), 1))
    slope = rng.normal(0, 0.0015, size=(len(cpgs), 1))
    mat = base + slope * (ages[None, :] - 50) + rng.normal(0, 0.02, size=(len(cpgs), n_samples))
    betas = pd.DataFrame(np.clip(mat, 0, 1), index=cpgs, columns=samples)
    meta = pd.DataFrame({"age": ages}, index=samples)
    return MethylationDataset(betas, meta, name="SYNTHETIC_FIXTURE", is_synthetic=True)
