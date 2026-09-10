# Roadmap

One thesis runs through all of this: in a field flooded with models that predict
everything, **trustworthy interpretation is the scarce good.** Each tool is a
public, reproducible artifact that says clearly what a prediction does and does not
support. Order reflects where the unfair advantage is largest.

### Now — clockcheck (this repo)
Honest audits for epigenetic aging clocks: weight-weighted coverage, imputation
sensitivity, calibration, age-acceleration fragility, cross-clock disagreement.
Grows into a small, general biomarker-audit library.

- [ ] Add more published clocks (GrimAge, Zhang, Weidner, EpiTOC).
- [ ] Platform-shift audit (450K vs EPIC vs EPICv2 probe manifests).
- [ ] Batch/normalization sensitivity (BMIQ vs noob vs raw).
- [ ] A gallery of audit cards across public GEO cohorts, committed to the repo.
- [ ] Short methods write-up / preprint.

### Next — variant → regulatory-consequence explainer
Wrap AlphaGenome / the AVI score for one tissue with plain, caveat-inclusive
interpretation: what a non-coding variant is predicted to affect, the confidence,
and what experiment would validate it. Same honesty stance, new modality.

### Then — 3D-genome interpreter for AlphaGenome
The strongest personal wedge: AlphaGenome now emits **contact maps**. Take those
predictions and build the honest layer on top — predicted vs. observed local
architecture (boundary strength, loop presence) for a regulatory variant, flagging
where the model is confident and where it hand-waves. Very few people combine
polymer/extrusion intuition, real Hi-C experience, and pioneer-TF regulatory
grammar; this is where that combination pays off.

---
Each step is shippable on its own and feeds an "I read the paper" explainer for
the community. Small, rigorous, public — compounding beats boiling the ocean.
