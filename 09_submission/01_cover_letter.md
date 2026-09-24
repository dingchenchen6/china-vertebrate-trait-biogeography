# Cover letter

**To the Editors, *Nature Ecology & Evolution***

*(Written for* Nature Ecology & Evolution*. Figures, Extended Data figures and
Supplementary Tables are numbered in order of first citation, as all candidate journals
require. The section order (Methods before Results), author–year references and Discussion
subheadings follow ecology-journal convention and must be converted for a Nature-family
submission; see `07_submission_checklist.md` for the full list of format gaps.)*

---

Dear Editors,

Please consider the enclosed manuscript, **"Thermoregulatory mode sets the strength of
environmental filtering on vertebrate assemblage traits, and leaves a directional
conservation gap"**, for publication as an Article.

**What the paper reports.** Trait-based biogeography has a standard move: collapse many
traits into one composite axis, then ask which environmental gradient explains it. Using
2,837 terrestrial vertebrate species across 3,782 50-km cells of China, we show that both
steps discard information, and that what is discarded is exactly what conservation
planning needs.

Three results carry the paper.

1. **Assemblage trait structure is at least three separable dimensions, and compositing
   them costs sensitivity.** Body size, pace of life and niche breadth are near-independent
   (|*r*| ≤ 0.30), and combining traits into a composite axis *weakens* the
   environmental signal in six of eight comparisons. This is a methodological result with
   reach beyond our system: the conventional "PCA-then-regress" pipeline systematically
   understates trait–environment coupling.

2. **Climate dominates in every class, but thermoregulatory mode sets how strongly.**
   Ambient energy affects ectotherm assemblages about 2.3× more strongly than endotherm
   assemblages, stable across specifications (1.98–2.63; every 95% CI excluding 1). We
   test this with the fourth-corner max test rather than naive CWM regression, having
   first quantified the latter's Type I inflation on these data (1.27–1.72×, worst in
   ectotherms).

3. **The reserve network's functional shortfall has a direction, and the direction
   reverses between classes.** The network covers 10.36% of land but brings only 2.9% of
   species to the Kunming–Montreal 30%-of-range target; an optimal layout of the same area
   would reach 50.7%. Large mammals are the best-protected trait end (+0.283) while large
   reptiles are the worst (−0.308) — the signature of a network assembled around flagship
   mammals. A model-free institutional correlate agrees: of 1,028 reserves, 17.9% name
   mammals as a protection target against 1.1% for reptiles.

**Why we think it fits *Nature Ecology & Evolution*.** The paper moves a well-established
descriptive literature to a mechanistic contrast and then to an actionable spatial result,
in one system, on one data foundation. The conservation conclusion is not an appended
implication — it is a prediction derived from the mechanism and then tested.

**On how the evidence is reported.** We would draw the editors' attention to one feature
that is unusual enough to state plainly. We pre-specified refutation conditions for each
hypothesis, and one of them was nearly met: ectotherms respond more strongly to five of
six predictors, not only to thermal ones. We report the thermal-buffering hypothesis as
**partially supported** rather than confirmed, replace the failed single negative control
with a dose–response test across eight predictors (ρ = 0.74, *P* = 0.037), and carry that
qualification into the abstract and the conclusions. Where earlier analyses failed their own
checks, the paper reports the correction rather than the original result: an apparent
amphibian over-dispersion signal disappeared once the null model weighted species by
occupancy (+1.70 to +0.08), and a claim that range maps understate anthropogenic effects
held only for birds once identical cells were compared, and is stated that narrowly. We
believe this makes the surviving claims more, not less, credible, and the Limitations
section states what the design cannot settle.

**Data and code.** All code, documentation, figures and aggregated per-cell results are
public. IUCN Red List range polygons may not be redistributed, so species-by-cell matrices
and species-level trait tables are excluded; the scripts rebuild both from their original
public sources, with a fixed random seed.

The manuscript is original, is not under consideration elsewhere, and all data sources are
credited. The author declares no competing interests.

Thank you for your consideration.

Yours sincerely,

**Chenchen Ding**
School of Life Sciences, Peking University, Beijing 100871, China
jialeding1220@gmail.com

---

## 〔投稿前需作者填写 / to be completed before submission〕

- **Suggested reviewers (3–5).** Choose from the trait-biogeography, macroecological-inference
  and conservation-prioritisation communities the paper engages with — the authors of the
  fourth-corner and CWM-inference literature, of the global reptile/amphibian assessment
  work, and of the functional-prioritisation literature are the natural pool. Check each for
  conflicts (co-authorship within 5 years, shared institution, shared funding) before listing.
- **Opposed reviewers**, if any, with a one-line neutral reason.
- **Funding statement** — grant numbers.
- **ORCID**.
- **Archived-release DOI** for the code and aggregated results (see `06_statements/`); mint
  a Zenodo DOI from the GitHub release and paste it into both the cover letter and the
  Data Availability statement before submitting.
