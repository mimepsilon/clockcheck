"""Small shared helpers. No fabricated data lives here."""
from __future__ import annotations
import numpy as np
import pandas as pd


def as_cpg_by_sample(df: pd.DataFrame) -> pd.DataFrame:
    """Return a beta matrix oriented as CpG (rows) x sample (cols).

    Illumina probe IDs start with 'cg'/'ch'/'rs'. We pick whichever axis looks
    more like probe IDs as the CpG axis, so the caller can pass either layout.
    """
    def frac_probe(labels) -> float:
        if len(labels) == 0:
            return 0.0
        return float(np.mean([str(x).startswith(("cg", "ch", "rs")) for x in labels]))

    if frac_probe(df.columns) > frac_probe(df.index):
        return df.T
    return df
