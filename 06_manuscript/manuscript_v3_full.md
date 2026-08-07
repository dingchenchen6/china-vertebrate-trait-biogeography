# Thermoregulatory mode sets the strength of environmental filtering on vertebrate assemblage traits, and leaves a directional conservation gap

**Chenchen Ding**

School of Life Sciences, Peking University, Beijing 100871, China
Correspondence: jialeding1220@gmail.com

**Running title** Thermal mode and vertebrate trait biogeography

**Word count** Main text ~5,400 | Figures 8 | Extended Data figures 25 | References 48

---

## Abstract

Trait-based biogeography routinely collapses many traits into one composite axis
and asks which environmental gradient explains it. We show that both steps
discard information, and that what is discarded matters for conservation. Using
2,837 species of Chinese terrestrial vertebrates across 3,782 50-km cells — 2,595
of them with complete trait data — we
define three functional axes *a priori* — body size, slow–fast pace of life and
niche breadth — so each axis denotes the same quantity in birds, mammals,
reptiles and amphibians. The axes are near-independent (|r| ≤ 0.30) and,
contrary to common practice, combining traits into a composite axis *weakens*
the environmental signal in six of eight comparisons. Regressing each trait and
each axis on six named predictors, one per classical biodiversity hypothesis,
climate dominates in all four classes (unique adjusted *R*² = 0.24–0.40), but
its strength depends on thermoregulation: ambient energy affects ectotherms
about 2.3× more strongly than endotherms, robustly across specifications
(1.98–2.63; every 95% CI excluding 1). A single negative control proved
insufficient — heterogeneity proxies differ in how much temperature information
they carry — so we tested the prediction as a trend instead: the ectotherm
excess scales with a predictor's thermal content (Spearman ρ = 0.74, *P* =
0.037, eight predictors). Ectotherms nonetheless respond more strongly to five
of six predictors, so we report thermal buffering as partially supported rather
than established. All associations were re-tested with the fourth-corner
max test, which shows naive community-weighted-mean regression inflates
significance by 1.27–1.72×. Applying the framework to spatial planning, China's
nature reserves occupy 10.36% of the land yet leave only 2.9% of species and
1.5% of functional bins at the Kunming–Montreal target of 30% of range
protected; an optimal layout of the same area would reach 50.7% and 48.0%. The
shortfall is directional: the body-size effect reverses between mammals (+0.283)
and reptiles (−0.307), and amphibian gaps are nearly disjoint from those of
other classes (Jaccard 0.074–0.086 versus 0.315 between the two endotherm
classes). A model-free correlate points the same way — of 1,028 reserves, 17.9%
name mammals as a protection target against 1.1% for reptiles. We deliver the
priority cells that close these gaps.

**Keywords** functional biogeography; community-weighted moments; thermal
buffering; fourth-corner; systematic conservation planning; Kunming–Montreal
Target 3

---

## 1. Introduction

Three literatures converge on the question we ask, and each leaves a specific
gap.

**Trait-based community ecology** proposed replacing species identities with
comparable functional quantities (McGill et al. 2006; Violle et al. 2007), an
approach since scaled up to functional biogeography (Violle et al. 2014). For
plants it produced a widely accepted global spectrum of form and function (Díaz
et al. 2016). For terrestrial vertebrates the equivalent synthesis is still
forming, and studies typically summarise traits with a within-class principal
component analysis. That step carries an unremarked cost: the first principal
component denotes different quantities in different classes. In our data PC1 is
dominated by niche breadth in birds and amphibians, by body size and fecundity
in reptiles, and by vertical stratum use in mammals (Extended Data Fig. 21).
Comparing "PC1" across classes therefore compares four different things.

**Hypotheses for broad-scale diversity gradients** — ambient energy, water–energy
dynamics, productivity, habitat heterogeneity and climatic seasonality — have
been contrasted for three decades (Wright 1983; O'Brien 1998; Kerr & Packer
1997; Hawkins et al. 2003; Currie et al. 2004). Almost all of that work takes
**species richness** as the response. Whether the same hypotheses explain
**trait composition**, and whether they apply equally to endotherms and
ectotherms, has not been tested systematically. A parallel thermal-biology
literature supplies a reason to expect they do not: ectotherm body temperature
tracks environmental temperature directly (Deutsch et al. 2008; Buckley et al.
2012; Sunday et al. 2014), so ambient energy should filter ectotherm assemblages
more strongly. This yields a prediction that can fail — the asymmetry should be
specific to thermally loaded predictors.

**Systematic conservation planning** has used complementarity since Margules &
Pressey (2000) and has begun to incorporate phylogenetic and functional
dimensions (Brum et al. 2017; Pollock et al. 2017). The Kunming–Montreal Global
Biodiversity Framework's Target 3 has made the question concrete and urgent.
Yet the two closest antecedents to this study — Huang et al. (2023) and Sun et
al. (2025), which established the framework for Chinese vertebrate trait
biogeography — stop at pattern and driver. Neither asks whether a
species-currency network leaves a **predictably directional** functional gap,
nor delivers the spatial layer that would close it.

China suits these questions. Its environmental gradient spans tropical to boreal
and sea level to 8,000 m within one country, so energy and water can be
separated; four vertebrate classes are sympatric under sound taxonomy, so the
endotherm–ectotherm contrast can be drawn under one set of environments; and
recently published Chinese trait datasets give ectotherms a usable baseline for
the first time (Song et al. 2022; Zhong et al. 2022; Wang et al. 2023; Ding et
al. 2022).

We ask three questions and state three hypotheses with their refutation
conditions in advance (Supplementary Table S0).

**Q1 (Dimensions).** Is assemblage trait structure one dimension or several, and
does combining traits into composite axes help or hurt?
*H1:* body size, pace of life and niche breadth are separable dimensions that
respond to different gradients. *Refuted if* the axes are strongly
inter-correlated, or environmental effects are identical in sign and strength
across them.

**Q2 (Drivers).** Which biodiversity hypothesis governs assemblage trait
structure, and is its strength set by thermoregulatory mode?
*H2:* endothermy buffers assemblage trait composition against thermal filtering.
*Refuted if* ectotherms respond more strongly to *all* predictors — which would
indicate a signal-to-noise difference rather than a mechanism — or if the
asymmetry disappears under a different replication unit.

**Q3 (Design).** What does trait structure require of conservation planning?
*H3:* if H1 and H2 hold, a species-currency network leaves a functional gap with
a predictable direction. *Refuted if* axis coefficients scatter around zero, or
if significance vanishes once range size is controlled, or if different target
currencies select the same places.

---

## 2. Materials and Methods

### 2.1 Analysis units

We used 3,782 land-intersecting cells on a 50-km grid in China Albers
equal-area conic projection (`+proj=aea +lat_0=0 +lon_0=110 +lat_1=25
+lat_2=47`), with a 957-cell 100-km grid as a grain check. Thirty-two cells
falling inside the South China Sea inset frame of the standard national basemap
were flagged and excluded from all analyses: they pair correct species lists
with Pacific-Ocean climate. Cells were flagged rather than deleted so that cell
identifiers remain stable across all derived files.

### 2.2 Species distributions

Expert range maps (IUCN Red List; BirdLife International) were rasterised to the
grid, giving 1,177 birds, 631 mammals, 520 reptiles and 509 amphibians in the
assemblage matrix (2,837 species). Of these, 2,595 (1,086 birds, 599 mammals,
472 reptiles, 438 amphibians) have complete values for all six core traits and
therefore carry functional-axis scores; assemblage moments are computed from
this subset, and the two counts are kept distinct throughout. Range maps overstate local occupancy; we use them because
they are the basis on which the two antecedent studies are comparable, and we
validate the main conclusions against assemblages rebuilt from GBIF occurrence
records (Extended Data Figs 17–18). Assemblage moments were computed only for
cells holding ≥5 species of a class, so the maps are blank where fewer species
occur; in the far west only 2.6% of cells reach five amphibian species, and the
resulting blanks are not absences.

### 2.3 Traits and the three functional axes

Six core traits entered the analysis: body size, fecundity, nocturnality,
vertical stratum use, habitat breadth and range size. **One size measure per
class** — body mass for endotherms, body length for ectotherms — because mixing
the two makes "size" a different quantity in different classes and invalidates
cross-class comparison. Observed coverage is 100% for size, nocturnality,
vertical stratum use, habitat breadth and range size in all four classes.

Life-history coverage was audited before use (Supplementary Table 36).
Generation length exists only for endotherms (birds 80.8%, mammals 92.2%);
ectotherm age at maturity and longevity are 7–24% observed. **Fecundity is the
only fast–slow component usable across all four classes** (reptiles 66.9%,
amphibians 60.0%, birds and mammals 100%), integrating ReptTraits (Oskyrko et
al. 2024), Wang et al. (2023) for snakes, Zhong et al. (2022) for lizards,
AmphiBIO (Oliveira et al. 2017), Song et al. (2022), Etard et al. (2020) and
COMBINE (Soria et al. 2021). Remaining gaps were imputed with missForest
(Stekhoven & Bühlmann 2012) using phylogenetic eigenvectors and the other
traits; traits below 30% observed are reported but excluded from core analyses
by default.

Three axes were defined *a priori* so each denotes the same quantity in every
class:

| Axis | Definition | Interpretation |
|---|---|---|
| SIZE | body size | body size |
| PACE | −residual(fecundity ~ size) | fewer or more offspring than size predicts |
| NICHE | habitat breadth (+) and range size (+) | niche breadth, narrow → broad |

PACE uses the **residual** rather than "size minus fecundity". The latter embeds
size in PACE and correlated the two axes at 0.38–0.73, making a body-size effect
inseparable from a pace-of-life effect. The residual is orthogonal to SIZE by
construction. The data-driven PCA is retained as a check (Extended Data Fig. 21).

### 2.4 Hypothesis-based drivers

Six named predictors, one per hypothesis: mean annual temperature (ambient
energy), annual precipitation (water availability), NPP (productivity),
within-cell elevation range (habitat heterogeneity), temperature seasonality
(climatic seasonality) and the 2020 human footprint (anthropogenic filter).
Named variables rather than composite components, so each coefficient reads
directly as support for one hypothesis. All six enter every model together, so
each coefficient is that hypothesis's unique effect after the others are
controlled (VIF ≤ 5.4; Supplementary Table 38b).

Responses are assemblage means of each of the six traits and each of the three
axes, per class: 36 models. We used spatial simultaneous autoregressive error
models (*k* = 8 symmetrised nearest neighbours, row-standardised) rather than
OLS, whose residual Moran's *I* reaches 0.69; SAR reduces residual Moran's *I*
to between −0.055 and 0.016 across all models (λ = 0.935–0.994). The sparse
Cholesky ("Matrix") method was used because eigendecomposition is O(*n*³) at
*n* ≈ 3,700; we verified it against the exact eigen solution on the smallest
class, where coefficients agree to 8 × 10⁻¹⁰ and log-likelihoods to 3 × 10⁻¹².
Variance was partitioned among three hypothesis blocks (climate, habitat,
human).

### 2.5 Valid significance testing

Regressing community-weighted means on the environment inflates Type I error
because all sites share one trait vector (Peres-Neto et al. 2017; ter Braak et
al. 2017). Every trait–environment link was therefore re-tested with the
fourth-corner max test (`modeltype = 6`, 999 permutations, spatially stratified
subsample of 1,200 cells), which takes the larger of the site-permutation and
species-permutation *P*-values. SAR supplies effect size and direction; the max
test supplies significance.

### 2.6 Replication unit for the thermal contrast

Grid cells are not independent and classes number only four. The family-level
test therefore uses **families** as replicates (21 endothermic, 12 ectothermic
at the main threshold), in a mixed model nested as (1|class/order), with a
Wilcoxon test as a distribution-free check. Families enter if they hold ≥20
Chinese species and occupy ≥100 cells — both estimability requirements. This
retains 33 of 166 families but 87% of ectotherm and 63% of endotherm species.
We report the result at three thresholds (§3.4).

### 2.7 Conservation prioritisation

Targets were set at **30% of each feature's range area**, matching Target 3,
rather than "represented at least once" — representation targets saturate at
about 10% of land here and are unrelated to persistence (Rodrigues & Gaston
2001). A greedy algorithm maximises shortfall reduction per unit area; greedy
carries a (1 − 1/e) worst-case guarantee for maximal coverage. Features are
either species or functional bins (the three axes each cut into six within-class
quantiles, concatenated with class label; bin-number sensitivity in
Supplementary Table 41b). Three siting principles — species complementarity,
trait-space complementarity and richness ranking — were compared against 99
random selections and against the existing reserve network, all at equal area.

Per-species protection is protected area within the range divided by range area,
modelled as a quasi-binomial GLM on the protected/unprotected split with the
three axes and log range extent as predictors. Quasi-binomial rather than OLS on
logits because many species have exactly zero protection, where the required
clip drives the coefficients.

### 2.8 Land-use tolerance

Etard et al. (2020) record whether each species has been observed using
artificial habitat — the most direct trait-level measure of land-use tolerance
available. Coverage is 80.0% (birds), 84.5% (mammals) and 67.4% (amphibians);
reptiles at 21.4% are reported descriptively only. All models control range
size, because poorly studied species are both more likely to be uncoded and more
likely to be narrow-ranged.

### 2.9 Reproducibility

All analyses were run in R 4.5.1 with a fixed seed (20260801). Code, aggregated
per-cell results and figures are archived (§ Data and code availability).

---

## 3. Results

### 3.1 Assemblage trait structure is at least three separable dimensions

The three *a priori* axes are near-independent within every class: SIZE–PACE
−0.107 to 0.030, SIZE–NICHE −0.021 to 0.297, PACE–NICHE −0.197 to 0.097
(Fig. 2; Supplementary Table 37b). Their geographies do not coincide (Fig. 3),
and their trait–environment coefficient patterns differ (Fig. 4a). H1 is not
refuted.

### 3.2 Composite axes are not more sensitive than their component traits

Comparing each axis's mean |effect| across the six hypotheses against the mean of
its components, six of eight class-by-axis comparisons favour the **components**;
only amphibian PACE gains (+0.030), while reptile PACE loses most (−0.025)
(Fig. 4c; Supplementary Table 39b). The routine "PCA first, then regress"
workflow therefore **understates** trait–environment coupling. The value of
*a priori* axes lies in giving cross-class comparison a common language, not in
raising power. We accordingly report single traits and axes side by side
throughout.

### 3.3 Climate dominates, and thermoregulation sets its strength

Climate carries the great majority of unique explanatory power in every class
(unique adjusted *R*² 0.341 amphibians, 0.345 birds, 0.396 mammals, 0.240
reptiles), against 0.007–0.032 for habitat and 0.005–0.032 for human pressure
(Fig. 4b; Supplementary Table 39). Reptiles are the exception in both
directions: the lowest climate share and the highest human share.

Averaged over all responses, the ambient-energy effect is |β| = 0.169 in
ectotherms and 0.072 in endotherms — a ratio of **2.35** (Fig. 6; Supplementary
Table 38). Across the six hypotheses the ectotherm excess is 2.81 (human
pressure), 2.35 (ambient energy), 2.08 (productivity), 1.88 (seasonality), 1.30
(water) and **1.03** (habitat heterogeneity).

### 3.4 The thermal asymmetry survives specification, but the family-level layer is weaker

We first tested H2 with habitat heterogeneity as a negative control, on the
grounds that it carries no thermal information. That control is not
proxy-independent. Replacing within-cell elevation range with land-cover type
diversity gives an ectotherm excess of 1.81 (95% CI 1.07–2.87), and with
EarthEnv remote-sensing texture 0.67 (0.40–1.10), against 1.12 (0.64–1.95) for
elevation range (Extended Data Fig. 22; Supplementary Table 45b).

Diagnosing why rescues the test. The three proxies differ greatly in how much
temperature information they carry: |Spearman *r*| with mean annual temperature
is 0.05 (texture), 0.20 (elevation range) and **0.52** (land-cover diversity).
The one proxy showing an ectotherm excess is the one most correlated with
temperature — in China land cover is itself arranged along the thermal gradient
— so it is a thermal predictor in disguise, not a valid control. Reframed as a
**trend**, the prediction holds: the ectotherm excess scales with a predictor's
thermal content across eight predictors (Spearman ρ = 0.74, *P* = 0.037). A
trend prediction is harder to satisfy by chance than the single point it
replaces. The energy ratio itself is robust to specification (1.98–2.63; every
95% CI excluding 1).

The family-level test is directionally consistent but weaker than we first
reported. Ectothermic families show 1.79× the thermal coupling of endothermic
families at the main threshold (*P* = 0.017 mixed model; *P* = 0.026 Wilcoxon),
1.71× at ≥30 species (*P* = 0.049), but 1.33× and non-significant at ≥10 species
(*P* = 0.23), where 14 additional small families with noisier coupling estimates
enter (Extended Data Fig. 23; Supplementary Table 46b). An alternative
operationalisation — *R*² of the three axes on named thermal variables rather
than |β| on a composite thermal axis — gives 1.05. We therefore treat the
grid-cell analysis as the primary evidence and the family-level result as
directional corroboration.

**H2 comes close to its own refutation condition, and we report it as such.** We
stated in advance that H2 would be refuted if ectotherms responded more strongly
to *all* predictors, since that would indicate a signal-to-noise difference
rather than a thermal mechanism. In the event, five of six predictors give
ectotherm/endotherm ratios above 1.2 — human pressure 2.80, ambient energy 2.35,
productivity 2.08, seasonality 1.88, water 1.30 — and only habitat heterogeneity
is flat at 1.03. Ectotherms do respond more strongly to nearly everything.

Three considerations keep the thermal interpretation alive without rescuing it
entirely. (i) The ratio is not constant across predictors, as a pure
signal-to-noise account predicts, but ordered: it scales with thermal content
(ρ = 0.74, *P* = 0.037). (ii) The two predictors with the least thermal content
are the two flattest (elevation range 0.198 → 1.12; remote-sensing texture 0.054
→ 0.67). (iii) The ordering is nonetheless imperfect — human pressure has the
highest ratio on only intermediate thermal content (0.609), and water
availability has substantial thermal content (0.574) but a low ratio (1.23) —
so the dose-response explains the gradient in aggregate rather than each point.

We therefore state H2 as **partially supported**: the energy asymmetry is large
and specification-robust, and its magnitude tracks thermal content, but we
cannot exclude that part of the ectotherm excess reflects a general difference
in how tightly ectotherm assemblages track any environmental gradient. A design
that could separate these — for example contrasting elevational with latitudinal
temperature gradients at matched heterogeneity — is beyond the present data.

### 3.5 Most links survive a valid test; naive testing inflates by up to 1.72×

Across 216 hypothesis × response × class tests, naive CWM regression declares
89–98% significant; the fourth-corner max test declares 54–76% — an inflation of
1.27× (birds), 1.36× (mammals), 1.55× (amphibians) and 1.72× (reptiles)
(Fig. 5b; Supplementary Table 40b). Inflation is larger in ectotherms. Of the
216 tests, 140 survive, led by ambient energy (29), climatic seasonality (28)
and water availability (26), with human pressure last (16) — the same ordering
the SAR effect sizes give, so two statistically independent routes agree
(Fig. 5; Supplementary Table 40c). Under SAR, 49.1% of the 216 coefficients are
significant after BH correction.

### 3.6 The same conclusions hold under other indices and unconstrained ordination

Three hull-free functional-diversity indices (FEve, FDis, RaoQ) against the same
six hypotheses agree in sign 83–95% of the time within each class, so the
conclusions are not an artefact of one index. Non-metric multidimensional
scaling of assemblage trait composition converges well in all four classes
(stress 0.087–0.100), and post hoc environmental fitting places the strongest
vector at *R*² = 0.72–0.87 (Extended Data Fig. 24; Supplementary Tables 47b,
48). An unconstrained ordination therefore recovers the same gradients the
constrained CWM models identify.

### 3.7 No class substitutes for another

Cross-class congruence of functional regionalisations is moderate at best
(adjusted Rand index 0.24–0.59), and the least congruent pair is
**birds–mammals** (0.241), not birds–reptiles (0.588) (Fig. 7; Supplementary
Table 13). The conservation-relevant asymmetry is instead functional redundancy:
functional β relative to taxonomic β gives redundancy 0.810 in birds, 0.659 in
mammals and **0.338** in reptiles (Supplementary Table 20). In reptiles each
species carries the most unique function, so species-count protection converts
least readily into functional protection.

Conservation **gaps** separate the classes even more sharply. Jaccard overlap of
gap cells is 0.315 for the two endotherm classes, 0.238 and 0.182 for reptiles
against mammals and birds, but **0.074 and 0.086** for amphibians against
mammals and birds (Supplementary Table 28). Looking at birds or mammals will not
find the amphibian gaps.

A model-free institutional correlate points the same way. Of China's 1,028
nature reserves, 17.9% name mammals among their protection targets and 12.8%
name birds, against 3.2% for amphibians and **1.1%** for reptiles
(Supplementary Table 30). Consistent with this, only for mammals is assemblage
trait volume positively associated with reserve coverage (+0.11); for reptiles
it is negative (−0.05), and birds and amphibians are indistinguishable from
zero.

### 3.8 The reserve network's shortfall is large, and the gap is directional

China's nature reserves occupy 10.36% of the land, but only **2.9%** of species
and **1.5%** of functional bins reach 30% of range protected; an optimal layout
of the same area reaches 50.7% and 48.0% (Fig. 8a,b; Supplementary Table 42).
Part of this gap reflects concentration — reserves are thinly spread across many
cells whereas the greedy solution takes whole cells — and we do not read it as
evidence that reserves are simply misplaced.

Target currency matters most when the budget is tight. At 10% of land a
trait-space target covers 46.7% of functional bins against 35.0% for a species
target, and at 17% of land 62.8% against 52.9%; by 30% the two converge (92.2%
versus 92.3%) (Supplementary Table 41). What matters at the 30% budget is
complementarity rather than richness ranking, worth about 9 percentage points
(86.4% versus 77.1% of species). The two currencies nevertheless place
priorities differently: of the 30% solutions, 764 cells are shared but 377 are
selected only by the species target and 379 only by the trait target — about a
third of the priority network depends on the currency (Fig. 8e).

Under-protection is not random. Of twelve axis coefficients, seven are
significant after BH correction, and their directions are interpretable
(Fig. 8c; Supplementary Table 43). Slow-paced species are under-protected in
mammals (−0.123, *P* = 5 × 10⁻⁴) and amphibians (−0.182, *P* = 6 × 10⁻⁵). The
niche axis runs narrow to broad, so its negative coefficients in mammals
(−0.177) and amphibians (−0.240, *P* = 5 × 10⁻¹⁶) identify **wide-ranging
generalists**, not specialists, as the under-protected end; reptiles reverse
(+0.211). Body size also **reverses**: large mammals are better protected
(+0.283, *P* = 7 × 10⁻¹⁵) while large reptiles are under-protected (−0.307,
*P* = 2 × 10⁻¹⁰). No bird axis is significant; birds have the widest ranges and
are the easiest to cover incidentally. H3 is not refuted.

### 3.9 The generalist gap is not the alarming direction it appears

Using recorded use of artificial habitat as a direct trait-level measure of
land-use tolerance, niche breadth dominates tolerance (log-odds +2.34
amphibians, +1.86 birds, +2.51 mammals; all *P* < 10⁻¹²), large-bodied species
are less tolerant (−0.30 to −0.45), and pace of life is irrelevant
(Extended Data Fig. 25; Supplementary Table 49). At the assemblage level the
tolerant fraction rises with human footprint (ρ = 0.37 amphibians, 0.30 birds,
0.55 mammals) — the anthropogenic filter read off a trait rather than inferred
from a CWM.

The key test reverses the expectation. Species that tolerate artificial habitat
are the **less** protected group in all three classes with adequate coverage
(quasi-binomial coefficient −0.445 amphibians, −0.217 birds, −0.424 mammals;
all *P* < 2 × 10⁻⁵). Median range protected is 4.34%, 5.68% and 5.65% for
intolerant species against 2.74%, 4.89% and 5.29% for tolerant ones
(Supplementary Table 50). Reserves sit in remote cold high country, which is
where specialists are; generalists live in the human-dominated east where
reserves are scarce. The generalist gap is therefore real but reassuring in
direction: the species least able to persist in modified habitat are, at
present, the better protected.

---

## 4. Discussion

The three results form one argument. Trait structure is not one thing, so it
must be examined dimension by dimension (Q1). Each dimension is governed by
climate, but the strength of that governance is set by thermoregulatory mode
(Q2). Because classes are filtered at different strengths and are poor
surrogates for one another, a species-currency network leaves a functional gap
whose direction is predictable — and, being predictable, locatable and closable
(Q3).

**Thermal buffering as a first-order factor.** The ectotherm excess for ambient
energy is large and specification-robust, and it scales with how much thermal
information a predictor carries. That dose-response is the strongest form the
evidence takes, and it emerged only because a single negative control failed
under an alternative operationalisation. We report the failure and the
replacement rather than the version that first worked, because the failure is
informative: "habitat heterogeneity" is not one thing either, and proxies for it
differ by an order of magnitude in thermal content.

**Two methodological conclusions travel beyond this system.** First, composite
trait axes dilute rather than concentrate environmental signal at the assemblage
level; investigators who summarise traits with PC1 before regression should
expect to understate their effects, and the practice of reporting only the
composite should be reconsidered. Second, equiprobable null models manufacture
trait clustering in taxa with many narrow-ranged species. Switching to occupancy
weighting moved amphibian SES from +1.70 to +0.08 and eliminated an apparent
over-dispersion signal that we had initially reported and have withdrawn
(Extended Data Fig. 8).

**Conservation implications are specific rather than general.** The problem is
not principally that China protects too little land, but where and how
concentratedly: the same 10.36% could deliver an order of magnitude more target
attainment. Expansion toward 30% should be assessed in trait as well as species
currency, because the two disagree on about a third of the priority network, and
the disagreement is largest exactly at the tight budgets that describe the
present situation. The directional findings translate into concrete screening
criteria: large reptiles, slow-paced mammals and amphibians, and — with the
caveat above — wide-ranging generalists are the ends a reserve system built
around large flagship mammals systematically misses. That the designation record
itself names mammals sixteen times as often as reptiles suggests the mechanism
is institutional as much as ecological, and therefore addressable.

**Amphibians as the sharpest case.** Amphibian gaps overlap those of other
classes by 0.074–0.086, against 0.315 between birds and mammals. This was
invisible while amphibians were excluded from the gap analysis on a
trait-coverage argument that no longer applies once one size measure per class
is adopted. It is the strongest form the surrogacy result takes, and it means
amphibian priorities must be derived from amphibian data.

### 4.1 Limitations

1. **Range maps overstate local occupancy.** We use them for comparability with
   the antecedent studies and validate against GBIF-derived assemblages, where
   conclusions hold in direction. One earlier claim — that range maps understate
   anthropogenic effects — survives only for birds once the same cells are
   compared, and has been narrowed accordingly.
2. **Protection fraction is an upper bound**: it does not deduct unsuitable
   habitat inside reserves. The bias runs the same way for all species, so it
   does not affect comparisons *between traits*, which is what we use it for.
3. **The efficiency gap partly reflects concentration** rather than siting; we
   have not separated the two components, which would require a constrained
   optimisation that preserves the observed size distribution of reserves.
4. **Only nature reserves are included**, not national parks or other
   categories, so coverage is a lower bound on China's protected-area system.
5. **Ectotherm life-history values are substantially imputed.** Supplementary
   Table 36 gives observed fractions trait by trait; traits below 30% observed
   are excluded from core analyses. Generation length remains unavailable for
   both ectotherm classes and is the single largest remaining data gap.
6. **The family-level thermal contrast is threshold-sensitive** (§3.4) and
   should not be cited independently of the grid-cell result.

### 4.2 Outlook

Three extensions follow directly. Compiling ectotherm generation length from the
Chinese-language literature would allow a fast–slow axis defined by both ends
rather than by a residual. Because ectotherm trait structure is energy-dominated,
warming scenarios predict faster compositional change in ectotherm assemblages —
a falsifiable prediction the present framework can be projected onto. And the
50-km priority layer is a starting point, not a plan: within the identified
cells, land tenure, land use and connectivity would determine implementable
parcels.

---

## Data and code availability

All code, documentation, figures and aggregated per-cell results are available
at https://github.com/dingchenchen6/china-vertebrate-trait-biogeography
(seed 20260801; results fully reproducible). IUCN Red List spatial data may not
be redistributed, so species-by-cell matrices and species-level trait tables are
excluded and their provenance is documented in `DATA.md`; both are rebuilt by
the scripts from their public sources.

## Author contributions

C.D. designed the study, assembled the data, performed all analyses, produced
the figures and wrote the manuscript.

## Competing interests

The author declares no competing interests.

---

## Figures

| No. | Content |
|---|---|
| Fig. 1 | Study design, trait coverage and the assemblage framework |
| Fig. 2 | Three *a priori* functional axes and their independence |
| Fig. 3 | Geography of assemblage trait means |
| Fig. 4 | Which biodiversity hypothesis drives assemblage trait structure |
| Fig. 5 | Fourth-corner max test of trait–environment links |
| Fig. 6 | Endotherms and ectotherms differ in trait–environment coupling |
| Fig. 7 | No single class is a reliable surrogate |
| Fig. 8 | Conservation design: complementarity, under-protected traits and 30×30 |

Extended Data Figs 1–25 and the full manifest are listed in
`04_results/tables/table_00_figure_manifest.csv`.

---

## References

Brum, F.T. et al. (2017) Global priorities for conservation across multiple
dimensions of mammalian diversity. *PNAS* 114, 7641–7646.

Buckley, L.B. et al. (2012) Broad-scale ecological implications of ectothermy
and endothermy in changing environments. *Global Ecology and Biogeography* 21,
873–885.

Currie, D.J. et al. (2004) Predictions and tests of climate-based hypotheses of
broad-scale variation in taxonomic richness. *Ecology Letters* 7, 1121–1134.

Deutsch, C.A. et al. (2008) Impacts of climate warming on terrestrial ectotherms
across latitude. *PNAS* 105, 6668–6672.

Díaz, S. et al. (2016) The global spectrum of plant form and function. *Nature*
529, 167–171.

Ding, C. et al. (2022) A dataset on the traits of terrestrial mammals in China.
*Biodiversity Science* 30, 21520. [丁晨晨等, 中国陆生哺乳动物特征数据集]

Dray, S. & Legendre, P. (2008) Testing the species traits–environment
relationships: the fourth-corner problem revisited. *Ecology* 89, 3400–3412.

Etard, A. et al. (2020) Global gaps in trait data for terrestrial vertebrates.
*Global Ecology and Biogeography* 29, 2143–2158.

Etard, A. & Newbold, T. (2024) Species-level correlates of land-use responses
and climate-change sensitivity in terrestrial vertebrates. *Conservation
Biology* 38, e14208.

Hawkins, B.A. et al. (2003) Energy, water, and broad-scale geographic patterns
of species richness. *Ecology* 84, 3105–3117.

Huang, M. et al. (2023) [Chinese terrestrial vertebrate trait biogeography]
*The Innovation* 4, 100371.

Kerr, J.T. & Packer, L. (1997) Habitat heterogeneity as a determinant of mammal
species richness in high-energy regions. *Nature* 385, 252–254.

Margules, C.R. & Pressey, R.L. (2000) Systematic conservation planning. *Nature*
405, 243–253.

McGill, B.J. et al. (2006) Rebuilding community ecology from functional traits.
*Trends in Ecology & Evolution* 21, 178–185.

Moura, M.R. & Jetz, W. (2021) Shortfalls and opportunities in terrestrial
vertebrate species discovery. *Nature Ecology & Evolution* 5, 631–639.

O'Brien, E.M. (1998) Water–energy dynamics, climate, and prediction of woody
plant species richness. *Journal of Biogeography* 25, 379–398.

Oliveira, B.F. et al. (2017) AmphiBIO, a global database for amphibian
ecological traits. *Scientific Data* 4, 170123.

Oskyrko, O. et al. (2024) ReptTraits: a comprehensive dataset of ecological
traits in reptiles. *Scientific Data* 11, 243.

Peres-Neto, P.R. et al. (2017) Linking trait variation to the environment:
critical issues with community-weighted mean correlation resolved by the
fourth-corner approach. *Ecography* 40, 806–816.

Pollock, L.J. et al. (2017) Large conservation gains possible for global
biodiversity facets. *Nature* 546, 141–144.

Rodrigues, A.S.L. & Gaston, K.J. (2001) How large do reserve networks need to
be? *Ecology Letters* 4, 602–609.

Song, Y. et al. (2022) [A dataset on the traits of amphibians in China]
*Biodiversity Science* 30. [宋云枫等]

Soria, C.D. et al. (2021) COMBINE: a coalesced mammal database of intrinsic and
extrinsic traits. *Ecology* 102, e03344.

Stekhoven, D.J. & Bühlmann, P. (2012) MissForest — non-parametric missing value
imputation for mixed-type data. *Bioinformatics* 28, 112–118.

Sun, Y. et al. (2025) Environmental and evolutionary factors jointly shape
life-history trait diversity of terrestrial vertebrates across China.
*Zoological Research* 46.

Sunday, J.M. et al. (2014) Thermal-safety margins and the necessity of
thermoregulatory behaviour across latitude and elevation. *PNAS* 111,
5610–5615.

ter Braak, C.J.F. et al. (2017) Fourth-corner correlation is a score test
statistic in a log-linear trait–environment model. *Environmental and Ecological
Statistics* 24, 1–17.

Tobias, J.A. et al. (2022) AVONET: morphological, ecological and geographical
data for all birds. *Ecology Letters* 25, 581–597.

Tuanmu, M.-N. & Jetz, W. (2015) A global, remote sensing-based characterization
of terrestrial habitat heterogeneity. *Global Ecology and Biogeography* 24,
1329–1339.

Violle, C. et al. (2007) Let the concept of trait be functional! *Oikos* 116,
882–892.

Violle, C. et al. (2014) The emergence and promise of functional biogeography.
*PNAS* 111, 13690–13696.

Wang, J. et al. (2023) [A dataset on the traits of snakes in China]
*Biodiversity Science*. [王江等]

Wright, D.H. (1983) Species-energy theory: an extension of species-area theory.
*Oikos* 41, 496–506.

Zhong, Y. et al. (2022) [A dataset on the traits of lizards in China]
*Biodiversity Science* 30. [钟雨茜等]
