"""Offline demo. Runs the full pipeline on the LABELED SYNTHETIC fixture so you can
see the audit card shape without any downloads. The card is stamped NON-EVIDENTIAL
on purpose — synthetic data cannot tell you anything about a clock.

For a real audit, replace `synthetic_fixture(...)` with:
    from clockcheck import load_geo
    ds = load_geo("GSE40279")        # needs biolearn + network
or
    from clockcheck import load_csv
    ds = load_csv("betas.csv", "meta.csv", age_col="age")
"""
from clockcheck import LinearClock, build_audit_card, card_to_markdown, synthetic_fixture

clocks = [LinearClock.load("Horvath1"), LinearClock.load("Hannum"), LinearClock.load("PhenoAge")]
ds = synthetic_fixture(clocks, n_samples=60, seed=0)

card = build_audit_card(clocks[0], ds, compare_clocks=clocks[1:])
print(card_to_markdown(card))
