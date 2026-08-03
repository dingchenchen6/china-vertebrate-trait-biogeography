# Thermoregulatory mode sets the strength of environmental filtering on vertebrate assemblage traits, and leaves a directional gap in conservation

**Chenchen Ding**
Peking University, Beijing, China

---

## Abstract

Trait-based biogeography usually collapses many traits into a single composite
axis and asks which environmental gradient explains it. We show that both steps
lose information, and that what is lost matters for conservation. Using 2,595
species of Chinese terrestrial vertebrates across 3,782 50-km grid cells, we
define three functional axes *a priori* — body size, slow–fast pace of life, and
niche breadth — so that each axis means the same thing in birds, mammals,
reptiles and amphibians. The three axes are near-independent (|r| ≤ 0.30), and,
contrary to standard practice, combining traits into a composite axis *weakens*
the environmental signal in six of eight comparisons. We then regress each trait
and each axis on six named predictors, one per classical biodiversity
hypothesis, entering all six simultaneously so every coefficient is that
hypothesis's unique effect. Climate dominates in all four classes (unique
adjusted R² = 0.24–0.40), but its strength depends on thermoregulation: the
effect of ambient energy is 2.35 times larger in ectotherms than in endotherms,
whereas habitat heterogeneity — a predictor carrying no thermal information — is
indistinguishable between them (1.03 times). That negative control separates
thermal buffering from a generic difference in signal-to-noise, and a
family-level test (n = 33 families) reproduces the asymmetry independently
(1.79×, P = 0.017). All associations are re-tested with the fourth-corner max
test, which shows that naive community-weighted-mean regression inflates
significance by 1.27–1.72×. Finally, we ask what this implies for spatial
planning. China's nature reserve network occupies 10.36% of the land but leaves
only 2.9% of species and 1.5% of functional bins at the Kunming–Montreal target
of 30% of range protected; an optimal layout of the same area would reach 50.7%
and 48.0%. The shortfall is directional rather than random: slow-paced and
narrow-niched species are missed in three classes, and the body-size effect
reverses between mammals (+0.283) and reptiles (−0.308) — a network built around
large flagship mammals gives large reptiles no equivalent shelter. We provide
the priority grid cells that close these gaps.

**Keywords** functional biogeography · community-weighted moments · thermal
buffering · fourth-corner · systematic conservation planning · 30×30

---

## 1. Introduction

Three literatures converge on this study, and each has left a specific gap.

**Trait-based community ecology** proposed replacing species identities with
comparable functional quantities (McGill et al. 2006; Violle et al. 2007), and
the approach has since scaled up to functional biogeography (Violle et al.
2014). For plants it produced a widely accepted global spectrum of form and
function (Díaz et al. 2016). For terrestrial vertebrates the equivalent
synthesis is still forming, and studies typically summarise traits with a
within-class PCA. That step has an unremarked cost: the first principal
component means different things in different classes. In our data PC1 is
dominated by niche breadth in birds and amphibians, by body size and fecundity
in reptiles, and by vertical habitat use in mammals (Extended Data Fig. 21).
Comparing "PC1" across classes therefore compares four different quantities.

**Hypotheses for broad-scale diversity gradients** — ambient energy, water–energy
dynamics, productivity, habitat heterogeneity, climatic seasonality — have been
tested against one another for three decades (Wright 1983; O'Brien 1998; Kerr &
Packer 1997; Hawkins et al. 2003; Currie et al. 2004). Almost all of that work
takes **species richness** as the response. Whether the same hypotheses explain
**trait composition**, and whether they apply equally to endotherms and
ectotherms, has not been tested systematically. A parallel thermal-biology
literature gives a reason to expect they do not: ectotherm body temperature
tracks environmental temperature directly (Deutsch et al. 2008; Buckley et al.
2012; Sunday et al. 2014), so ambient energy should filter ectotherm
assemblages more strongly. This yields a prediction that can fail: the asymmetry
should be **specific to thermally loaded predictors**, and absent for a purely
structural one.

**Systematic conservation planning** has used complementarity since Margules &
Pressey (2000), and has begun to incorporate phylogenetic and functional
dimensions (Brum et al. 2017; Pollock et al. 2017). The Kunming–Montreal Global
Biodiversity Framework's Target 3 has made the question urgent and concrete. Yet
the two closest antecedents to this study — Huang et al. (2023) and Sun et al.
(2025), which established the framework for Chinese vertebrate trait
biogeography — stop at pattern and driver. Neither asks whether a
species-currency network leaves a **predictably directional** functional gap,
nor delivers the spatial layer that would close it.

China is a good place to close these gaps. Its environmental gradient spans
tropical to boreal and sea level to 8,000 m within one country, so energy and
water can be separated; four vertebrate classes are sympatric with sound
taxonomy, so the endotherm–ectotherm contrast can be drawn under one set of
environments; and recently published Chinese trait datasets give ectotherms a
usable trait baseline for the first time.

We ask three questions.

**Q1.** Is assemblage trait structure one dimension or several, and does
combining traits into composite axes help or hurt?

**Q2.** Which biodiversity hypothesis governs assemblage trait structure, and is
its strength set by thermoregulatory mode?

**Q3.** What does the answer require of conservation planning — where are the
gaps, are they directional, and how much land would close them?

We state three hypotheses with their refutation conditions in advance
(Supplementary Table S0), so that a negative result would be visible as such.

---

## 2. Methods

### 2.1 Analysis units

3,782 land-intersecting cells on a 50-km grid in China Albers equal-area conic
(`+proj=aea +lat_0=0 +lon_0=110 +lat_1=25 +lat_2=47`), with a 957-cell 100-km
grid as a grain check. Thirty-two cells falling inside the South China Sea inset
frame of the standard basemap are flagged and excluded from all analyses; they
pair correct species lists with Pacific-Ocean climate. Cells are flagged rather
than deleted so that cell identifiers stay stable across all derived files.

### 2.2 Species distributions

Expert range maps (IUCN Red List; BirdLife) rasterised to the grid, giving
1,177 birds, 631 mammals, 520 reptiles and 509 amphibians in the assemblage
matrix. Range maps overstate local occupancy; we use them because they are the
basis on which the two antecedent studies are comparable, and we validate the
main conclusions against assemblages rebuilt from GBIF occurrence records
(Extended Data Figs 17–18).

### 2.3 Traits and functional axes

Six core traits: body size, fecundity, nocturnality, verticality, habitat
breadth and range size. **One size measure per class** — body mass for
endotherms, body length for ectotherms — because mixing the two makes "size" a
different quantity in different classes and invalidates cross-class comparison.
Observed coverage is 100% for size, nocturnality, verticality, habitat breadth
and range size in all four classes.

Life-history coverage was audited before use (Table 36). Generation length
exists only for endotherms (birds 80.8%, mammals 92.2%); ectotherm age at
maturity and longevity are 7–24% observed. **Fecundity is the only fast–slow
component usable across all four classes** (reptiles 66.9%, amphibians 60.0%,
birds and mammals 100%), integrating ReptTraits (Oskyrko et al. 2024), Wang et
al. (2023) for snakes, Zhong et al. (2022) for lizards, AmphiBIO (Oliveira et
al. 2017), Song et al. (2022), Etard et al. (2020) and COMBINE. Remaining gaps
were imputed with missForest using phylogenetic eigenvectors and the other
traits; traits below 30% observed are written out but excluded from the core
analyses by default.

Three axes are defined *a priori* so each means the same thing in every class:

| Axis | Definition | Meaning |
|---|---|---|
| SIZE | body size | body size |
| PACE | −residual(fecundity ~ size) | fewer or more offspring than size predicts |
| NICHE | habitat breadth (+) and range size (+) | niche breadth |

PACE uses the **residual** rather than "size minus fecundity". The latter embeds
size in PACE and correlated the two axes at 0.38–0.73, making it impossible to
separate a body-size effect from a pace-of-life effect. The residual is
orthogonal to SIZE by construction. The data-driven PCA is retained as a check
(Extended Data Fig. 21).

### 2.4 Hypothesis-based drivers

Six named predictors, one per hypothesis: mean annual temperature (ambient
energy), annual precipitation (water availability), NPP (productivity),
within-cell elevation range (habitat heterogeneity), temperature seasonality
(climatic seasonality) and the 2020 human footprint (anthropogenic filter).
Named variables rather than composite components, so each coefficient reads
directly as support for one hypothesis. All six enter every model together, so
each coefficient is that hypothesis's unique effect after the others are
controlled; VIF ≤ 5.4 (Table 38b).

Responses are assemblage means of each of the six traits and each of the three
axes, per class: 36 models. Spatial simultaneous autoregressive error models
(k = 8 symmetrised nearest neighbours, row-standardised) rather than OLS, whose
residual Moran's I reaches 0.69; SAR reduces it to −0.001. The sparse Cholesky
("Matrix") method is used because eigendecomposition is O(n³) at n ≈ 3,700; we
verified against the exact eigen solution on the smallest class, where
coefficients agree to 8 × 10⁻¹⁰ and the log-likelihoods to 3 × 10⁻¹². Variance
is partitioned among three hypothesis blocks (climate, habitat, human).

### 2.5 Valid significance testing

Regressing community-weighted means on the environment inflates Type I error
because all sites share one trait vector (Peres-Neto et al. 2017; ter Braak et
al. 2017). Every trait–environment link is therefore re-tested with the
fourth-corner max test (`modeltype = 6`, 999 permutations, spatially stratified
subsample of 1,200 cells), which takes the larger of the site-permutation and
species-permutation *P*-values. SAR supplies effect size and direction; the max
test supplies significance.

### 2.6 Replication unit for the thermal contrast

Grid cells are not independent, and classes number only four. The headline
endotherm–ectotherm test therefore uses **families** as replicates (21
endothermic, 12 ectothermic), in a mixed model nested as (1|class/order), with a
Wilcoxon test as a distribution-free check.

### 2.7 Conservation prioritisation

Targets are set at **30% of each feature's range area**, matching Target 3,
rather than "represented at least once" — representation targets saturate at
about 10% of land here and are unrelated to persistence (Rodrigues & Gaston
2001). A greedy algorithm maximises shortfall reduction per unit area; greedy
carries a (1 − 1/e) worst-case guarantee for maximal coverage. Features are
either species or functional bins (the three axes each cut into six
within-class quantiles, concatenated with class label; bin-number sensitivity in
Table 41b). Three siting principles — species complementarity, trait-space
complementarity, and richness ranking — are compared against 99 random
selections and against the existing reserve network, all at equal area.

Per-species protection is the protected area within the range divided by range
area, modelled as a quasi-binomial GLM on the protected/unprotected split with
the three axes and log range extent as predictors. Quasi-binomial rather than
OLS on logits because many species have exactly zero protection, where the
required clip drives the coefficients.

---

## 3. Results

### 3.1 Trait structure is at least three separable dimensions (Q1)

The three *a priori* axes are near-independent within every class: SIZE–PACE
−0.107 to 0.030, SIZE–NICHE −0.021 to 0.297, PACE–NICHE −0.197 to 0.097
(Fig. 2; Table 37b). Their geographies do not coincide (Fig. 3), and the
trait–environment coefficient patterns differ among them (Fig. 4a).

A result we did not expect: **composite axes are not more sensitive than their
component traits**. Comparing each axis's mean |effect| across the six
hypotheses against the mean of its components, six of eight class-by-axis
comparisons favour the components; only amphibian PACE gains (+0.030), while
reptile PACE loses most (−0.025) (Fig. 4c; Table 39b). The routine "PCA first,
then regress" workflow therefore **understates** trait–environment coupling. The
value of *a priori* axes lies in giving cross-class comparison a common
language, not in raising power.

### 3.2 Climate dominates, and thermoregulation sets its strength (Q2)

Climate carries the great majority of unique explanatory power in every class
(unique adjusted R² 0.24–0.40), against 0.007–0.032 for habitat and 0.005–0.032
for human pressure (Fig. 4b; Table 39). Reptiles are the exception in both
directions: the lowest climate share (0.240) and the highest human share
(0.032).

The strength of the climate effect depends on thermoregulatory mode. Averaged
over all responses, the ambient-energy effect is |β| = 0.169 in ectotherms and
0.072 in endotherms — **2.35×** (Fig. 6; Table 38). The critical comparison is
the negative control: habitat heterogeneity, which carries no thermal
information, gives 0.043 versus 0.042 — **1.03×**. The asymmetry is therefore
specific to thermally loaded hypotheses, not a general tendency of ectotherms to
respond more strongly to everything.

The family-level test reproduces this independently: ectothermic families show
1.79× the thermal coupling of endothermic families (P = 0.017 mixed model,
P = 0.026 Wilcoxon), and 1.95× for human coupling (Table 6b).

### 3.3 Most links survive a valid test, but naive testing inflates by up to 1.72× (Q2)

Across 216 hypothesis × response × class tests, naive CWM regression declares
89–98% significant; the fourth-corner max test declares 54–76% — an inflation of
1.27× (birds) to 1.72× (reptiles) (Fig. 5b; Table 40b). Inflation is larger in
ectotherms. The surviving links are led by climatic seasonality (19), ambient
energy (18) and water availability (16), matching the SAR effect-size ordering:
two statistically independent routes give the same answer (Fig. 5; Table 40c).

### 3.4 No class substitutes for another (Q3)

Cross-class congruence of functional regionalisations is moderate at best
(ARI 0.24–0.59), and the least congruent pair is **birds–mammals** (0.241), not
birds–reptiles (0.588) (Fig. 7; Table 13). The relevant asymmetry is instead
functional redundancy: functional β relative to taxonomic β gives redundancy
0.810 in birds, 0.659 in mammals and **0.338** in reptiles (Table 20). In
reptiles each species carries the most unique function, so species-count
protection converts least readily into functional protection.

### 3.5 The reserve network's shortfall is large, and the gap is directional (Q3)

China's nature reserves occupy 10.36% of the land, but only **2.9%** of species
and **1.5%** of functional bins reach 30% of range protected; an optimal layout
of the same area reaches 50.7% and 48.0% (Fig. 8a,b; Table 42). Part of this gap
reflects concentration — reserves are thinly spread across many cells whereas
the greedy solution takes whole cells — and we do not read it as evidence that
reserves are simply misplaced.

Target currency matters most when the budget is tight. At 10–17% of land, a
trait-space target covers 10–12 percentage points more functional space than a
species target; by 30% the two converge (92.2% versus 92.3%) (Table 41). What
matters at the 30% budget is complementarity rather than richness ranking,
which is worth about 9 percentage points. The two currencies nevertheless place
priorities differently: of the 30% solutions, 764 cells are shared but 377 are
selected only by the species target and 379 only by the trait target — about a
third of the priority network depends on the currency (Fig. 8e).

Under-protection is not random. Of twelve axis coefficients, seven are
significant after BH correction, and their directions are interpretable
(Fig. 8c; Table 43). Slow-paced species are under-protected in mammals
(−0.123) and amphibians (−0.182); narrow-niched species in mammals (−0.177)
and, most strongly of all, amphibians (−0.240, P = 5 × 10⁻¹⁶). Body size
**reverses**: large mammals are better protected (+0.283) while large reptiles
are under-protected (−0.308). This points to a specific historical cause — a
reserve system sited around large flagship mammals — and it gives large reptiles
no equivalent shelter. No bird axis is significant; birds have the widest ranges
and are the easiest to cover incidentally.

---

## 4. Discussion

The three results form one argument. Trait structure is not one thing, so it
must be examined dimension by dimension (Q1). Each dimension is governed by
climate, but the strength of that governance is set by thermoregulatory mode,
with a negative control that rules out the obvious alternative explanation (Q2).
Because classes are filtered at different strengths and are poor surrogates for
one another, a species-currency network leaves a functional gap whose direction
is predictable — and, being predictable, locatable and closable (Q3).

Two methodological conclusions travel beyond this system. First, composite trait
axes dilute rather than concentrate the environmental signal at the assemblage
level; investigators who summarise traits with PC1 before regression should
expect to understate their effects. Second, equiprobable null models manufacture
trait clustering in taxa with many narrow-ranged species. In our data, switching
to occupancy weighting moved amphibian SES from +1.70 to +0.08 and eliminated an
apparent over-dispersion signal that we had initially reported and have
withdrawn.

The conservation implication is specific rather than general. The problem is not
principally that China protects too little land, but where and how
concentratedly that land is protected: the same 10.36% could deliver an order of
magnitude more target attainment. And expansion toward 30% should be assessed in
trait as well as species currency, because the two disagree on about a third of
the priority network.

### Limitations

1. Range maps overstate local occupancy. We use them for comparability with the
   antecedent studies and validate against GBIF-derived assemblages, where the
   conclusions hold in direction; one earlier claim — that range maps understate
   anthropogenic effects — survives only for birds once the same cells are
   compared, and has been narrowed accordingly.
2. Protection fraction is an upper bound: it does not deduct unsuitable habitat
   inside reserves. The bias runs the same way for all species, so it does not
   affect comparisons **between traits**, which is what we use it for.
3. The efficiency gap partly reflects concentration rather than siting.
4. Only nature reserves are included, not national parks or other categories, so
   coverage is a lower bound on China's protected-area system.
5. Ectotherm life-history values are substantially imputed; Table 36 gives
   observed fractions trait by trait, and traits below 30% observed are excluded
   from the core analyses.

---

## Data and code availability

All code, documentation, figures and aggregated per-cell results are at
https://github.com/dingchenchen6/china-vertebrate-trait-biogeography
(SEED = 20260801; results fully reproducible). IUCN Red List spatial data may
not be redistributed, so species-by-cell matrices and species-level trait tables
are excluded and their provenance is documented in `DATA.md`; both are rebuilt
by the scripts from their public sources.

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
| Fig. 8 | Conservation design: complementarity, under-protected traits, 30×30 |

Extended Data Figs 1–21 and the full figure manifest are listed in
`04_results/tables/table_00_figure_manifest.csv`.
