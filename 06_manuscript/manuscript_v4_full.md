# Thermoregulatory mode sets the strength of environmental filtering on vertebrate assemblage traits, and leaves a directional conservation gap

**Chenchen Ding**

School of Life Sciences, Peking University, Beijing 100871, China
Correspondence: jialeding1220@gmail.com

**Running title** Thermal mode and vertebrate trait biogeography

**Word count** Main text ~8,000 | Figures 8 | Extended Data figures 25 | References 100

---

## Abstract

Trait-based biogeography routinely collapses many traits into one composite axis
and asks which environmental gradient explains it. We show that both steps
discard information, and that what is discarded matters for conservation. Using
2,837 species of Chinese terrestrial vertebrates across 3,782 50-km cells —
2,595 of them with complete trait data — we define three functional axes
*a priori* — body size, slow–fast pace of life and niche breadth — so each axis
denotes the same quantity in birds, mammals, reptiles and amphibians. The axes
are near-independent (|r| ≤ 0.30) and, contrary to common practice, combining
traits into a composite axis *weakens* the environmental signal in six of eight
comparisons. Regressing each trait and each axis on six named predictors, one
per classical biodiversity hypothesis, climate dominates in all four classes
(unique adjusted *R*² = 0.24–0.40), but its strength depends on
thermoregulation: ambient energy affects ectotherms about 2.3× more strongly
than endotherms, robustly across specifications (1.98–2.63; every 95% CI
excluding 1). A single negative control proved insufficient — heterogeneity
proxies differ in how much temperature information they carry — so we tested the
prediction as a trend instead: the ectotherm excess scales with a predictor's
thermal content (Spearman ρ = 0.74, *P* = 0.037, eight predictors). Ectotherms
nonetheless respond more strongly to five of six predictors, so we report
thermal buffering as partially supported rather than established. All
associations were re-tested with the fourth-corner max test, which shows naive
community-weighted-mean regression inflates significance by 1.27–1.72×.
Applying the framework to spatial planning, China's nature reserves occupy
10.36% of the land yet leave only 2.9% of species and 1.5% of functional bins
at the Kunming–Montreal target of 30% of range protected; an optimal layout of
the same area would reach 50.7% and 48.0%. The shortfall is directional: the
body-size effect reverses between mammals (+0.283) and reptiles (−0.308), and
amphibian gaps are nearly disjoint from those of the endotherm classes
(Jaccard 0.074 and 0.086, versus 0.315 between birds and mammals) and only
weakly aligned with reptile gaps (0.160). A model-free
correlate points the same way — of 1,028 reserves, 17.9% name mammals as a
protection target against 1.1% for reptiles. We deliver the priority cells that
close these gaps.

**Keywords** functional biogeography; community-weighted mean; thermal
buffering; fourth-corner; systematic conservation planning; Kunming–Montreal
Target 3; protected areas; China

---

## 1. Introduction

The functional composition of assemblages — which trait values co-occur, where —
links biodiversity pattern to ecological process more directly than species
counts do. Traits govern how organisms acquire resources, tolerate abiotic
stress and respond to human pressure (McGill et al. 2006; Violle et al. 2007),
and shifts in assemblage-level trait composition foreshadow losses of ecological
strategies that species richness alone cannot register (Cooke et al. 2019;
Carmona et al. 2021; Toussaint et al. 2021). Scaling this logic up from plots to
maps founded functional biogeography (Violle et al. 2014; Lavorel & Garnier
2002). For plants the programme has matured into a widely accepted global
spectrum of form and function (Díaz et al. 2016) and continental syntheses of
community trait–environment relationships (Swenson et al. 2012; Lamanna et al.
2014; Bruelheide et al. 2018). For terrestrial vertebrates the equivalent
synthesis is still assembling, propelled by a rapid accumulation of trait
databases (Wilman-era compilations superseded by AVONET, COMBINE, AmphiBIO,
ReptTraits and the phylogeny-informed TetrapodTraits backbone; Tobias et al.
2022; Soria et al. 2021; Oliveira et al. 2017; Oskyrko et al. 2024; Moura et al.
2024). Yet three specific gaps remain between what these data could answer and
what has been asked of them.

**Gap 1 — Dimensionality and cross-class comparability.** Multi-trait studies
routinely compress traits into a leading principal component or a multivariate
index before asking what drives them. Two problems follow. First, a within-class
PCA is free to rotate differently in every class: in our data PC1 is dominated
by niche breadth in birds and amphibians, by body size and fecundity in
reptiles, and by vertical stratum use in mammals (Extended Data Fig. 21), so
comparing "PC1" across classes compares four different quantities. Whether a
fixed, *a priori* set of axes — the assemblage-level analogue of the plant
strategy schemes that made cross-taxon comparison possible in botany (Díaz et
al. 2016) — can be defined for vertebrates, and whether such composites gain or
lose statistical sensitivity relative to their component traits, has not been
tested. Second, the workhorse statistic of the field, the community-weighted
mean (CWM), carries a known inferential defect: because every site shares one
trait vector, CWM–environment regressions inflate Type I error severely
(Dray & Legendre 2008; Peres-Neto et al. 2017; Zelený 2018; Hawkins et al.
2017). The corrected procedure — the fourth-corner max test, which permutes
sites and species simultaneously (ter Braak et al. 2012) — exists but is rarely
applied in broad-scale trait biogeography, and the magnitude of the inflation
it corrects has seldom been quantified on real assemblage data.

**Gap 2 — Do richness hypotheses govern trait composition, and equally for
endotherms and ectotherms?** The classical hypotheses for broad-scale diversity
gradients — ambient energy (Wright 1983), water–energy dynamics (O'Brien 1998),
productivity, habitat heterogeneity (Kerr & Packer 1997) and climatic
seasonality (Janzen 1967) — have been contrasted for three decades, almost
always with **species richness** as the response (Hawkins et al. 2003; Currie
et al. 2004). Whether the same hypotheses explain **trait composition**, and
whether they apply with equal force to endotherms and ectotherms, remains
untested at the assemblage scale. Thermal biology supplies a directional
prediction. Ectotherm body temperature, performance and activity budgets track
environmental temperature directly (Huey & Kingsolver 1989; Angilletta 2009;
Deutsch et al. 2008; Sunday et al. 2011, 2014; Kearney & Porter 2009), whereas
endotherm thermal tolerances vary far less across climate (Khaliq et al. 2014;
Araújo et al. 2013). Endothermy should therefore *buffer* assemblage trait
composition against thermal filtering, making the ambient-energy hypothesis
stronger for ectotherms. Crucially, this prediction carries its own falsifier:
if ectotherm assemblages simply track *every* environmental gradient more
tightly — a signal-to-noise difference, not a thermal mechanism — the asymmetry
should appear for thermally neutral predictors too. Single-class studies cannot
make this contrast; a four-class design under one set of environments can.

**Gap 3 — From pattern to plan.** Systematic conservation planning rests on
complementarity (Margules & Pressey 2000; Kukkala & Moilanen 2013), and global
analyses show that adding phylogenetic and functional dimensions changes where
priorities fall (Brum et al. 2017; Pollock et al. 2017). The Kunming–Montreal
Global Biodiversity Framework's Target 3 — 30% of land effectively conserved by
2030 (CBD 2022; Maxwell et al. 2020) — has made the question of *where* urgent
and concrete. Three strands of evidence suggest that a species-currency network
will miss functionally distinctive biota in a *predictable direction*:
surrogacy among taxa is notoriously imperfect (Rodrigues & Brooks 2007);
reptiles and amphibians, the classes with the most distinctive distributions,
are globally under-served by protected areas and by assessment effort (Roll et
al. 2017; Cox et al. 2022; Böhm et al. 2013; Stuart et al. 2004; Luedtke et al.
2023); and traits themselves predict sensitivity to land use and extinction
risk (Cardillo et al. 2005; Newbold et al. 2015; Chichorro et al. 2019; Etard &
Newbold 2024). China's protected-area system — over 2,700 nature reserves
covering roughly 15% of the land, with documented representation biases (Wu et
al. 2011; Xu et al. 2017) — is a consequential place to test this. The two
closest antecedents to the present study, which established Chinese vertebrate
trait biogeography at exactly this grain (Huang et al. 2023; Sun et al. 2025),
stop at pattern and driver: neither asks whether the reserve network leaves a
directional functional gap, nor delivers the spatial layer that would close it.

China suits all three questions. Its environmental gradient spans tropical to
boreal and sea level to 8,000 m within one jurisdiction, so energy and water can
be statistically separated; four vertebrate classes are sympatric under sound
taxonomy, so the endotherm–ectotherm contrast can be drawn under one set of
environments; and recently published Chinese trait datasets (Song et al. 2022;
Zhong et al. 2022; Wang et al. 2023; Ding et al. 2022), combined with the
global compilations above, give ectotherms a usable trait baseline for the
first time.

We therefore ask three questions, stating hypotheses and their refutation
conditions in advance (Supplementary Table S0).

**Q1 (Dimensions).** Is assemblage trait structure one dimension or several, and
does combining traits into composite axes help or hurt?
*H1:* body size, pace of life and niche breadth are separable dimensions that
respond to different gradients. *Refuted if* the axes are strongly
inter-correlated, or environmental effects are identical in sign and strength
across them.

**Q2 (Drivers).** Which biodiversity hypothesis governs assemblage trait
structure, and is its strength set by thermoregulatory mode?
*H2:* endothermy buffers assemblage trait composition against thermal
filtering. *Refuted if* ectotherms respond more strongly to *all* predictors —
indicating a signal-to-noise difference rather than a mechanism — or if the
asymmetry disappears under a different replication unit.

**Q3 (Design).** What does trait structure require of conservation planning?
*H3:* if H1 and H2 hold, a species-currency network leaves a functional gap with
a predictable direction. *Refuted if* axis coefficients scatter around zero, or
if significance vanishes once range size is controlled, or if different target
currencies select the same places.

This study makes four advances. (i) It introduces *a priori*, cross-class
functional axes for vertebrate assemblages and shows — against expectation —
that composites usually dilute rather than concentrate environmental signal.
(ii) It tests the classical richness hypotheses on trait composition with
hypothesis-labelled predictors and inference that survives the CWM problem,
quantifying the Type I inflation that naive practice incurs. (iii) It evaluates
thermal buffering against its own pre-stated falsifier, reports it as partially
supported, and replaces a failed single negative control with a stronger
dose–response test. (iv) It converts the pattern into planning: a target-based
functional prioritisation aligned with Target 3, a diagnosis of *which* trait
ends the existing network misses, and the priority cells that close the gap.

---

## 2. Materials and Methods

### 2.1 Study area and analysis units

We used 3,782 land-intersecting cells on a 50-km grid in the China Albers
equal-area conic projection (`+proj=aea +lat_0=0 +lon_0=110 +lat_1=25
+lat_2=47`), with a 957-cell 100-km grid as a grain check. Thirty-two cells
falling inside the South China Sea inset frame of the standard national basemap
were flagged and excluded from all analyses: they pair correct species lists
with open-ocean climate. Cells were flagged rather than deleted so that cell
identifiers remain stable across all derived files. All maps use the standard
national basemap, including the inset.

### 2.2 Species distributions and assemblage matrices

Expert range maps (IUCN 2024; BirdLife International & HBW 2022) were
rasterised to the grid, giving 1,177 birds, 631 mammals, 520 reptiles and 509
amphibians in the assemblage matrix (2,837 species). Of these, 2,595 species
(1,086 birds, 599 mammals, 472 reptiles, 438 amphibians) have complete values
for all six core traits and therefore carry functional-axis scores; assemblage
moments are computed from this subset, and the two counts are kept distinct
throughout. Scientific names were harmonised across trait and range sources by
a normalisation routine (trimming, capitalisation, subspecies removal) followed
by synonym matching; unmatched names are listed in Supplementary Table S3b–c.

Assemblage metrics were computed only for cells holding ≥5 species of a class.
The resulting blanks — most visibly for amphibians in arid north-west China,
where only 2.6% of far-western cells reach five species — are estimability
limits, not absences. Range maps overstate local occupancy (Hurlbert-type
commission error); we use them because they are the basis on which the
antecedent studies are comparable, and we validate the main conclusions against
assemblages rebuilt from GBIF occurrence records (Extended Data Figs 17–18),
where conclusions hold in direction.

### 2.3 Trait compilation: one size measure per class

Six core traits entered the analysis: body size, fecundity, nocturnality,
vertical stratum use (0 = fossorial to 1 = aerial), habitat breadth and range
size, assembled on the phylogeny-informed TetrapodTraits backbone (Moura et al.
2024) and overridden wherever a measured value exists in a specialist source:
AVONET for birds (Tobias et al. 2022), COMBINE and a national mammal dataset
for mammals (Soria et al. 2021; Ding et al. 2022), ReptTraits plus national
snake and lizard datasets for reptiles (Oskyrko et al. 2024; Wang et al. 2023;
Zhong et al. 2022), and AmphiBIO plus a national amphibian dataset for
amphibians (Oliveira et al. 2017; Song et al. 2022).

We used **one size measure per class** — body mass for endotherms, body length
for ectotherms — because mixing mass and length makes "size" a different
quantity in different classes and invalidates cross-class comparison. This
choice also removes a large imputation artefact: amphibian body *mass* in the
global backbone is 74.7% allometrically inferred, whereas amphibian body
*length* is 0.7% imputed. Observed coverage is 100% for size, nocturnality,
vertical stratum use, habitat breadth and range size in all four classes.

Life-history coverage was audited before use (Supplementary Table 36).
Generation length exists essentially only for endotherms (birds 80.8%, mammals
92.2%; both ectotherm classes 0%). Age at maturity is observed for 7.1% of
amphibians and for no reptile in the dataset; longevity for 23.5% of reptiles
and 8.4% of amphibians. **Fecundity is the only fast–slow component usable across
all four classes** (reptiles 66.9%, amphibians 60.0%, birds and mammals 100%).
Sex-specific lengths in the Chinese amphibian compilation are mostly written
as ranges ("350.0-417.0"). Parsing them to interval midpoints rather than
coercing them directly raised usable male–female pairs from 75 to 472 species;
stripping non-numeric characters instead, as is often done, would have
concatenated "145-191" into 145191 mm.

### 2.4 Imputation and its limits

Remaining gaps were imputed with missForest (Stekhoven & Bühlmann 2012), a
random-forest approach that performs well on ecological trait data (Penone et
al. 2014), using the ten leading phylogenetic eigenvectors and the other traits
as predictors, per class (out-of-bag NRMSE: birds 0.033, mammals 0.058,
amphibians 0.270, reptiles 0.413). Imputed entries are flagged; traits below
30% observed are reported but excluded from core analyses by default, because
after imputation more than two-thirds of their values would be model prediction
rather than measurement. Cross-validation against independent global values
gave *r* = 0.96–0.99 for body size.

### 2.5 Three a priori functional axes

Axes were defined *a priori* so each denotes the same quantity in every class:

| Axis | Definition | Interpretation |
|---|---|---|
| SIZE | body size (one measure per class) | body size |
| PACE | −residual(fecundity ~ size) | fewer or more offspring than size predicts (slow → fast) |
| NICHE | habitat breadth (+) and range size (+), z-mean | niche breadth, narrow → broad |

SIZE and PACE follow the two dominant dimensions of vertebrate life-history
variation (Stearns 1992; Promislow & Harvey 1990; Bielby et al. 2007; Healy et
al. 2019); the fast–slow continuum structures strategy variation across the
tree of life (Salguero-Gómez et al. 2016). PACE uses the **residual** rather
than "size minus fecundity": the latter embeds size in PACE and correlated the
two axes at 0.38–0.73, making a body-size effect inseparable from a
pace-of-life effect, whereas the residual is orthogonal to SIZE by
construction. NICHE combines habitat breadth with range size, exploiting their
general association (Slatyer et al. 2013); its geographic behaviour connects to
classical range-size gradients (Stevens 1989). Nocturnality and vertical
stratum use were kept as single traits — they are not homologous components of
any of the three axes. Traits were z-standardised within class before axis
construction, and the data-driven PCA is retained as a check (Extended Data
Fig. 21).

### 2.6 Environmental predictors and hypothesis mapping

Six named predictors, one per hypothesis: mean annual temperature (bio1;
ambient energy), annual precipitation (bio12; water availability), net primary
productivity (productivity), within-cell elevation range (habitat
heterogeneity), temperature seasonality (bio4; climatic seasonality) and the
2020 human footprint (anthropogenic filter). Climate and elevation come from
WorldClim 2.1 (Fick & Hijmans 2017), NPP from MODIS MOD17A3, land cover from
ESA WorldCover, habitat texture from EarthEnv (Tuanmu & Jetz 2015), and the
human footprint from the 2000–2020 annual series (Mu et al. 2022; lineage
Venter et al. 2016). We deliberately use named variables rather than composite
components, so each coefficient reads directly as support for one hypothesis;
all six enter every model together, making each coefficient that hypothesis's
unique effect after the others are controlled (VIF ≤ 5.4; Supplementary Table
38b).

### 2.7 Assemblage metrics and null models

For each cell × class we computed community-weighted moments (mean, variance,
skewness, kurtosis) of each trait and axis, and functional diversity as
functional richness, evenness, divergence, dispersion and Rao's quadratic
entropy (Mason et al. 2005; Villéger et al. 2008; Laliberté & Legendre 2010;
Botta-Dukát 2005) on Gower distances. Standardised effect sizes used 999
draws from an **occupancy-weighted** null model: species are drawn with
probability proportional to their observed occupancy, holding richness fixed
(cf. Gotelli 2000). The conventional equiprobable null treats narrow endemics
and widespread species as exchangeable and thereby manufactures apparent trait
clustering in endemic-rich taxa; switching to occupancy weighting moved
amphibian SES from +1.70 to +0.08 and eliminated a spurious over-dispersion
signal (Extended Data Fig. 8). The equiprobable null is retained as a
sensitivity analysis. Functional β-diversity was decomposed into turnover and
nestedness (Baselga 2010; Villéger et al. 2013).

### 2.8 Driver models

Responses are assemblage means of each of the six traits and each of the three
axes, per class: 36 models. Because assemblage data are strongly spatially
autocorrelated (Legendre 1993) and OLS residual Moran's *I* reached 0.69, we
used simultaneous autoregressive error models (Dormann et al. 2007; Kissling &
Carl 2008) with *k* = 8 symmetrised nearest neighbours and row-standardised
weights; SAR reduces residual Moran's *I* to between −0.055 and 0.016 across
all models (λ = 0.935–0.994). The sparse Cholesky ("Matrix") estimator was used
because eigendecomposition is O(*n*³) at *n* ≈ 3,700; we verified it against
the exact eigen solution on the smallest class, where coefficients agree to
8 × 10⁻¹⁰ and log-likelihoods to 3 × 10⁻¹². Variance was partitioned among
three hypothesis blocks (climate, habitat, human) with adjusted *R*².

### 2.9 Valid significance: the fourth-corner max test

Every trait–environment link was re-tested with the fourth-corner max test
(Dray & Legendre 2008; ter Braak et al. 2012; Peres-Neto et al. 2017):
`modeltype = 6`, 999 permutations, on a spatially stratified subsample of 1,200
cells, taking the larger of the site-permutation and species-permutation
*P*-values, with Benjamini–Hochberg correction. SAR supplies effect size and
direction; the max test supplies significance. Naive CWM regressions were run
on the identical matrices to quantify their Type I inflation (Zelený 2018).

### 2.10 The thermal contrast and its robustness

**Grid-cell contrast.** The ectotherm/endotherm ratio of mean |β| was computed
per predictor, with 95% bootstrap CIs over class × response fits (2,000
resamples), and the difference of log-ratios between predictors tested by joint
bootstrap.

**Heterogeneity operationalisations.** Because "habitat heterogeneity" is not
one variable, the negative control was repeated with three proxies:
within-cell elevation range, land-cover type diversity (Shannon diversity of
the six WorldCover class shares, row-normalised; cells with <5% mapped cover
excluded) and EarthEnv remote-sensing texture. Their thermal content —
|Spearman *r*| with mean annual temperature — spans 0.05 to 0.52, and the
ectotherm excess was then tested as a **dose–response** across all eight
predictors (Spearman correlation of log-ratio with thermal content).

**Family-level replication.** Grid cells are not independent and classes number
four, so the class contrast was replicated with **families** as units (mixed
model with (1|class/order) nesting; Wilcoxon as a distribution-free check).
Families enter at ≥20 Chinese species and ≥100 occupied cells — both
estimability requirements — retaining 33 of 166 families but 87% of ectotherm
and 63% of endotherm species. The full pipeline was re-run at ≥10/≥50 and
≥30/≥150 to test threshold sensitivity, and coupling was re-operationalised
(axis *R*² on named thermal variables) to test definition sensitivity.

The family-level responses are drawn from the fixed six-trait core computed in
the community-metrics step, in which body size enters as body mass for every
class so that the cross-class contrast concerns a single named quantity.
Amphibian mass, however, is 74.7% derived allometrically from length
(Extended Data Fig. 1), and amphibian assemblage mass and length correlate only
*r* = 0.76 (reptiles 0.87; birds 0.98; mammals 1.00), so the choice is not
cosmetic. We therefore re-ran the whole family-level pipeline with one size
measure per class — length for amphibians, matching the trait protocol used for
the functional axes and the grid-cell driver models — and report both
(Supplementary Table 19). Only the size response changes: the Gower distance
underlying functional dispersion is computed within family, so the measures it
contains affect all species of a family alike and do not enter the cross-class
comparison.

### 2.11 Robustness to metric and method

Conclusions were re-checked with three hull-free functional-diversity indices
(FEve, FDis, RaoQ) against the same six predictors, and with non-metric
multidimensional scaling of assemblage trait composition (*k* = 2, Euclidean on
z-scored CWM profiles, 1,200-cell stratified subsample) followed by post hoc
environmental vector fitting (999 permutations). The 100-km grid replicates the
main models at a coarser grain.

### 2.12 Conservation analyses

**Coverage and gaps.** Reserve boundaries (1,028 national and provincial nature
reserves, 2021 vintage) were intersected with the grid to give per-cell
protected fractions. Gap cells per class are those in the top richness or
functional-diversity quartile with below-median protection; cross-class gap
overlap used Jaccard similarity. As a model-free institutional correlate, we
tallied which classes each reserve's stated protection targets mention.

**Target-based prioritisation.** Representation targets ("at least one cell")
saturate at ~10% of land and are unrelated to persistence (Rodrigues & Gaston
2001; Rodrigues et al. 2004), so targets were set at **30% of each feature's
range area**, matching Target 3. Features are species, or functional bins (each
axis cut into six within-class quantiles, concatenated with class; 717 occupied
bins; sensitivity at 4–8 bins in Supplementary Table 41b). A greedy algorithm
maximising shortfall reduction per km² was used; greedy selection carries a
(1 − 1/e) worst-case guarantee for problems of this form (Nemhauser et al.
1978). Species complementarity, trait-space complementarity and richness
ranking were compared against 99 random selections and against the existing
network, all at equal area.

**Directional under-protection.** Per-species protected fraction (protected ÷
total range area) was modelled as a quasi-binomial GLM on the
protected/unprotected area split with the three axes and log range extent as
predictors — quasi-binomial because many species have exactly zero protection,
where the clipping required by logit-OLS drives the coefficients.

**Land-use tolerance.** Recorded use of artificial habitat (Etard et al. 2020)
served as a direct trait-level measure of land-use tolerance (coded for 80.0%
of birds, 84.5% of mammals, 67.4% of amphibians; reptiles 21.4%, descriptive
only). Tolerance was modelled on the three axes (logistic), its assemblage
fraction related to human footprint, and per-species protection related to
tolerance (quasi-binomial), always controlling range size because poorly
studied species are both more likely uncoded and more likely narrow-ranged.

### 2.13 Software and reproducibility

All analyses used R 4.5.1 (R Core Team 2025) with sf, terra, spdep, spatialreg,
vegan, FD, ade4 (Dray & Dufour 2007), lme4 and missForest, under a fixed seed
(20260801). Code, aggregated per-cell results and figures are archived
(§ Data and code availability).

---

## 3. Results

### 3.1 Assemblage trait structure is at least three separable dimensions

The three *a priori* axes are near-independent within every class: SIZE–PACE
−0.107 to 0.030, SIZE–NICHE −0.021 to 0.297, PACE–NICHE −0.197 to 0.097
(Fig. 2; Supplementary Table 37b). Their geographies do not coincide (Fig. 3),
and their trait–environment coefficient patterns differ (Fig. 4a). H1 is not
refuted.

### 3.2 Composite axes are not more sensitive than their component traits

Comparing each axis's mean |effect| across the six hypotheses against the mean
of its components, six of eight class-by-axis comparisons favour the
**components**; only amphibian PACE gains (+0.030), while reptile PACE loses
most (−0.025) (Fig. 4c; Supplementary Table 39b). The routine "PCA first, then
regress" workflow therefore **understates** trait–environment coupling. The
value of *a priori* axes lies in giving cross-class comparison a common
language, not in raising power. We accordingly report single traits and axes
side by side throughout.

### 3.3 Climate dominates, and thermoregulation sets its strength

Climate carries the great majority of unique explanatory power in every class
(unique adjusted *R*² 0.341 amphibians, 0.345 birds, 0.396 mammals, 0.240
reptiles), against 0.007–0.032 for habitat and 0.005–0.032 for human pressure
(Fig. 4b; Supplementary Table 39). Reptiles are the exception in both
directions: the lowest climate share and the highest human share.

Averaged over all responses, the ambient-energy effect is |β| = 0.169 in
ectotherms and 0.072 in endotherms — a ratio of **2.35** (Fig. 6; Supplementary
Table 38). Across the six hypotheses the ectotherm excess is 2.80 (human
pressure), 2.35 (ambient energy), 2.08 (productivity), 1.88 (seasonality), 1.30
(water) and **1.03** (habitat heterogeneity).

### 3.4 The thermal asymmetry survives specification, but comes close to its own refutation condition

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

The family-level replication is directionally consistent but weaker than the
grid-cell result. Ectothermic families show 1.79× the thermal coupling of
endothermic families at the main threshold (*P* = 0.017 mixed model; *P* =
0.026 Wilcoxon), 1.71× at ≥30 species (*P* = 0.049), but 1.33× and
non-significant at ≥10 species (*P* = 0.23), where 14 additional small families
with noisier coupling estimates enter (Extended Data Fig. 23; Supplementary
Table 46b). An alternative operationalisation — *R*² of the three axes on named
thermal variables rather than |β| on a composite thermal axis — gives 1.05. The
contrast is, by contrast, insensitive to how body size is measured: substituting
length for mass in amphibians — the one class where mass is mostly inferred —
moves the thermal ratio from 1.79 to 1.77 and the overall ratio from 1.72 to
1.67, with the Wilcoxon test significant either way (Supplementary Table 19). We
therefore treat the grid-cell analysis as the primary evidence and the
family-level result as directional corroboration, threshold-sensitive in the
number of families admitted but not in how their traits are measured.

**H2 comes close to its own refutation condition, and we report it as such.**
We stated in advance that H2 would be refuted if ectotherms responded more
strongly to *all* predictors, since that would indicate a signal-to-noise
difference rather than a thermal mechanism. In the event, five of six
predictors give ectotherm/endotherm ratios above 1.2 — human pressure 2.80,
ambient energy 2.35, productivity 2.08, seasonality 1.88, water 1.30 — and only
habitat heterogeneity is flat at 1.03. Ectotherms do respond more strongly to
nearly everything.

Three considerations keep the thermal interpretation alive without rescuing it
entirely. (i) The ratio is not constant across predictors, as a pure
signal-to-noise account predicts, but ordered: it scales with thermal content
(ρ = 0.74, *P* = 0.037). (ii) The two predictors with the least thermal content
are the two flattest (elevation range 0.198 → 1.12; remote-sensing texture
0.054 → 0.67). (iii) The ordering is nonetheless imperfect — human pressure has
the highest ratio on only intermediate thermal content (0.609), and water
availability has substantial thermal content (0.574) but a low ratio (1.23) —
so the dose–response explains the gradient in aggregate rather than each point.

We therefore state H2 as **partially supported**: the energy asymmetry is large
and specification-robust, and its magnitude tracks thermal content, but we
cannot exclude that part of the ectotherm excess reflects a general difference
in how tightly ectotherm assemblages track any environmental gradient. A design
that could separate these — for example contrasting elevational with
latitudinal temperature gradients at matched heterogeneity — is beyond the
present data.

### 3.5 Most links survive a valid test; naive testing inflates by up to 1.72×

Across 216 hypothesis × response × class tests, naive CWM regression declares
89–98% significant; the fourth-corner max test declares 54–76% — an inflation
of 1.27× (birds), 1.36× (mammals), 1.55× (amphibians) and 1.72× (reptiles)
(Fig. 5b; Supplementary Table 40b). Inflation is larger in ectotherms. Of the
216 tests, 140 survive, led by ambient energy (29), climatic seasonality (28)
and water availability (26), with human pressure last (16) — the same ordering
the SAR effect sizes give, so two statistically independent routes agree
(Fig. 5; Supplementary Table 40c). Under SAR, 49.1% of the 216 coefficients are
significant after BH correction.

### 3.6 The same conclusions hold under other indices and unconstrained ordination

Three hull-free functional-diversity indices (FEve, FDis, RaoQ) against the
same six hypotheses agree in sign 83–95% of the time within each class, so the
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
Table 13). The conservation-relevant asymmetry is instead functional
redundancy: functional β relative to taxonomic β gives redundancy 0.810 in
birds, 0.659 in mammals and **0.338** in reptiles (Supplementary Table 20). In
reptiles each species carries the most unique function, so species-count
protection converts least readily into functional protection.

Conservation **gaps** separate the classes even more sharply. Jaccard overlap
of gap cells is 0.315 for the two endotherm classes, 0.238 and 0.182 for
reptiles against mammals and birds, **0.074 and 0.086** for amphibians
against mammals and birds, and 0.160 for amphibians against reptiles
(Supplementary Table 28). The ordering follows thermoregulatory mode: amphibian
gaps overlap the other ectotherm class about twice as much as they overlap
either endotherm class — but even that is half the bird–mammal value. Looking
at birds or mammals will not find the amphibian gaps.

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
Part of this gap reflects concentration — reserves are thinly spread across
many cells whereas the greedy solution takes whole cells — and we do not read
it as evidence that reserves are simply misplaced.

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
(+0.283, *P* = 7 × 10⁻¹⁵) while large reptiles are under-protected (−0.308,
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

### 4.1 A multidimensional trait geography, and why composites dilute

That body size, pace of life and niche breadth vary near-independently across
Chinese vertebrate assemblages parallels what plant syntheses found for size
and economics (Díaz et al. 2016; Bruelheide et al. 2018), and it extends the
classical single-trait gradients — Bergmann-type size clines (Meiri & Dayan
2003; Olalla-Tárraga et al. 2006; Olson et al. 2009) and Rapoport-type
range-size gradients (Stevens 1989) — into a joint, cross-class framework. The
practical lesson is sharper than the parallel: composites usually *dilute*.
Averaging traits whose environmental responses differ in sign or strength
cancels signal, which is why six of eight composites underperform their
components. Studies that summarise assemblage traits with PC1 before
regression — still the modal workflow — should expect to understate
trait–environment coupling, and cross-class comparisons of free-rotating
components risk comparing different quantities altogether. Fixed, *a priori*
axes solve the comparability problem; they do not, and on this evidence cannot,
solve the power problem. Reporting single traits alongside axes costs little
and guards against both failure modes.

### 4.2 Thermal buffering: partially supported, honestly bounded

The prediction that endothermy buffers assemblage composition against thermal
filtering (Huey & Kingsolver 1989; Angilletta 2009; Deutsch et al. 2008; Khaliq
et al. 2014) survives in its strongest quantitative form — a ~2.3× energy
asymmetry, stable across specifications, scaling with a predictor's thermal
content — yet our own pre-stated falsifier was nearly met, because ectotherms
track almost every gradient more tightly. Candidate explanations for that
general excess are not mutually exclusive: genuinely tighter environmental
control of ectotherm distributions (Sunday et al. 2011; Kearney & Porter 2009);
smaller ranges and lower occupancy producing steeper assemblage turnover;
seasonality acting as a second thermal axis for organisms that cannot buffer it
(Janzen 1967) — consistent with seasonality's high ratio here (1.88); and
residual imputation noise attenuating endotherm coefficients less than
ectotherm ones. The dose–response argues that thermal content is the organising
variable, but its exceptions (human pressure high, water low) mark the boundary
of what these data can attribute. Marine comparisons show ectotherms living
closer to their thermal limits than terrestrial ones (Pinsky et al. 2019); an
analogous within-terrestrial design — elevational versus latitudinal
temperature gradients at matched heterogeneity — is the experiment our result
sets up. Meanwhile the applied implication is not hostage to mechanism:
whatever its ultimate cause, ectotherm assemblage composition is the more
climate-exposed, and warming should reorganise it faster.

### 4.3 Inference practice: what this study changes

Three methodological results travel beyond this system. First, the CWM Type I
problem (Dray & Legendre 2008; Peres-Neto et al. 2017; Zelený 2018; Hawkins et
al. 2017) is not an abstract worry: on real national-scale data the naive test
declares 89–98% of links significant where the valid max test (ter Braak et al.
2012) sustains 54–76%, and the inflation is largest exactly in the endemic-rich
ectotherm classes where conclusions matter most for conservation. Second,
equiprobable null models manufacture trait clustering wherever narrow-ranged
species are common; occupancy weighting removed an amphibian "over-dispersion"
signal we had initially believed. Third, environmental-filtering language
(Keddy 1992) needs the discipline urged by its critics (Kraft et al. 2015;
Cadotte & Tucker 2017): our "filtering" is assemblage-level association,
validated against permutation nulls, not demonstrated exclusion — and one of
our headline lessons is how easily an invalid control (a "structural" variable
that secretly carries temperature) can masquerade as mechanism.

### 4.4 Directional gaps, their institutional signature, and 30×30

The reserve network's functional gap has a sign, and the sign has a history.
Large mammals are the best-protected trait end (+0.283) while large reptiles
are the worst (−0.308) — precisely what a network assembled around flagship
mammal conservation would produce, and consistent with the global neglect of
reptiles in assessment and protection (Böhm et al. 2013; Roll et al. 2017; Cox
et al. 2022). The designation record makes the mechanism legible without any
model: mammals are named as protection targets sixteen times as often as
reptiles. This institutional signature matters because it is addressable by
policy, not just by land.

Amphibians sharpen the surrogacy warning to its limit. Their gap cells overlap
those of birds and mammals by 0.074 and 0.086, and those of reptiles by 0.160 —
versus 0.315 between birds and mammals — so no
amount of bird- or mammal-based planning will find them (Rodrigues & Brooks
2007). Given the global amphibian decline (Stuart et al. 2004; Luedtke et al.
2023), amphibian priorities must be derived from amphibian data; in China those
data now exist (Song et al. 2022).

For Target 3 (CBD 2022; Maxwell et al. 2020), the analysis yields three
planning rules. (i) The binding constraint today is siting and concentration,
not currency: the same 10.36% of land, optimally arranged, would raise
target attainment by an order of magnitude. (ii) At tight budgets (10–17% of
land) a functional currency buys 10–12 extra percentage points of trait-space
coverage; by 30% the currencies converge, and what matters is complementarity
over richness ranking (Kukkala & Moilanen 2013) — worth ~9 points of species
attainment. (iii) Because a third of the 30% network depends on the currency,
the choice should be explicit, not incidental. The apparent paradox that
under-protected generalists are also the land-use-tolerant species (Etard &
Newbold 2024; Newbold et al. 2015) resolves the alarm but not the obligation:
tolerance of modified habitat is persistence in degraded states, not a
conservation outcome equivalent to protection.

### 4.5 Limitations

1. **Range maps overstate local occupancy.** We use them for comparability with
   the antecedent studies and validate against GBIF-derived assemblages, where
   conclusions hold in direction. One earlier claim — that range maps
   understate anthropogenic effects — survives only for birds once identical
   cells are compared, and has been narrowed accordingly.
2. **Protection fraction is an upper bound**: it does not deduct unsuitable
   habitat inside reserves. The bias runs the same way for all species, so it
   does not affect comparisons *between traits*, which is what we use it for.
3. **The efficiency gap partly reflects concentration** rather than siting; we
   have not separated the two, which would require a constrained optimisation
   preserving the observed reserve size distribution.
4. **Only nature reserves are included**, not national parks or other
   categories, so coverage is a lower bound on China's protected-area system
   (cf. Xu et al. 2017).
5. **Ectotherm life-history values are substantially imputed.** Supplementary
   Table 36 gives observed fractions trait by trait; traits below 30% observed
   are excluded from core analyses. Generation length remains unavailable for
   both ectotherm classes and is the single largest remaining data gap.
6. **The family-level thermal contrast is threshold-sensitive** (§3.4) and
   should not be cited independently of the grid-cell result.

### 4.6 Conclusions and outlook

Assemblage trait structure in Chinese terrestrial vertebrates is
multidimensional, climate-governed with a strength set largely by
thermoregulatory mode, and mirrored by a conservation network whose functional
blind spots are directional and institutionally legible. Three extensions
follow directly. Compiling ectotherm generation length from the
Chinese-language natural-history literature would let the fast–slow axis be
defined by both ends rather than a residual. Because ectotherm composition is
the more energy-exposed, projecting the fitted models through warming scenarios
yields a falsifiable forecast: faster compositional reorganisation in reptile
and amphibian assemblages than in birds and mammals. And the 50-km priority
layer is a beginning, not a plan — within the identified cells, land use,
tenure and connectivity will determine implementable parcels.

---

## Data and code availability

The aggregated data underlying every figure and table — per-cell assemblage
metrics, community-weighted trait moments, functional-axis scores, model
coefficients and all result tables — together with the complete analysis code
are available at https://github.com/dingchenchen6/china-vertebrate-trait-biogeography
and archived at [Zenodo DOI]. The pipeline runs end to end from the raw public
sources with the random seed fixed at 20260801.

Species range polygons come from the IUCN Red List and BirdLife International,
whose licences prohibit redistribution. The species-by-cell occurrence matrices
and species-level trait tables derived from them are therefore not deposited.
They are not available on request either: both providers supply the source data
free of charge to any user under their standard terms, and the scripts rebuild
the derived files from those downloads without manual intervention. `DATA.md`
lists every source file, its version, its access date and the script that
consumes it. Source data for each main-text figure are indexed by
`table_00_figure_manifest.csv`.

## Author contributions

C.D. designed the study, assembled the data, performed all analyses, produced
the figures and wrote the manuscript.

## Competing interests

The author declares no competing interests.

---

## Figure legends

### Fig. 1 | Study design, trait coverage and the assemblage framework

**a**, The thermal-buffering hypothesis stated as a prediction about slopes. Because
endothermy decouples body temperature from ambient conditions, the set of trait values that
can persist in an endothermic assemblage should be less tightly constrained by the
environment than in an ectothermic one; coupling strength is the slope of assemblage trait
composition on the environmental gradient (bars at right). **b**, Species analysed per
class on the 50-km equal-area grid: 1,177 birds, 631 mammals, 520 reptiles and 509
amphibians, 2,837 in total across 3,782 analysable cells. **c–f**, Species richness of each
class. Grey cells hold no species of that class.

### Fig. 2 | Three *a priori* functional axes and their independence

**a**, Definition of the three axes. Rather than rotating a principal-component solution
freely within each class — which would make the leading axis denote a different quantity in
each class — each axis is a fixed combination of named traits, so that "body size", "pace of
life" and "niche breadth" mean the same thing in birds, mammals, reptiles and amphibians.
Pace of life is the negative residual of fecundity on body size, and is therefore orthogonal
to size by construction. **b**, Spearman correlations between axis pairs within each class.
All |*r*| ≤ 0.30 (maximum 0.297, size–niche in reptiles), so the axes carry separate
information. **c–f**, The fecundity–size relationship from which the pace axis is built,
one panel per class, with the Spearman correlation of the *residuals* against size shown.

### Fig. 3 | Geography of assemblage trait means across China

Community-weighted means of three traits — body size (**a–d**), habitat breadth (**e–h**)
and nocturnality (**i–l**) — for birds, mammals, reptiles and amphibians on the 50-km grid.
Each row shares one colour scale. Body size uses one measure per class: body mass for birds,
mammals and reptiles, body length for amphibians, whose mass is 74.7% derived allometrically
from length (Extended Data Fig. 1). Because traits are standardised within class before the
moments are computed, panels remain comparable despite the different underlying measures.
Blank cells hold fewer than five species of that class, so assemblage moments are not
estimable there; **they are not absences**. Cells mapped: birds 100%, mammals 100%, reptiles
87%, amphibians 53% of analysable land cells.

### Fig. 4 | Which biodiversity hypothesis drives assemblage trait structure

**a**, Standardised effect of each of six named predictors — one per classical biodiversity
hypothesis (ambient energy, water availability, productivity, habitat heterogeneity, climatic
seasonality, human pressure) — on every trait and every functional axis, by class, from
spatial simultaneous autoregressive error models. Every model controls for all six predictors
at once, so each cell is that hypothesis's unique effect. Dots mark BH-adjusted *P* < 0.05.
**b**, Variance uniquely attributable to the climate, habitat and human predictor groups.
Climate dominates in all four classes (class means 0.240–0.396 of adjusted *R*²) while the
habitat and human group means each stay below 0.04. **c**, Mean absolute standardised effect
of the component traits against that of the composite axis they form. In six of eight
comparisons the composite is *weaker* than its components, so combining traits into an axis
dilutes rather than concentrates the environmental signal.

### Fig. 5 | Fourth-corner max test of trait–environment links

**a**, Signed fourth-corner statistic for every trait × predictor combination, by class. Dots
mark links that survive the max test (`modeltype = 6`, which permutes sites and species
independently and takes the larger of the two *P* values) at BH-adjusted *P* < 0.05.
**b**, Percentage of tests declared significant by naive community-weighted-mean regression
against the max test. The naive test inflates significance because all sites share a single
trait vector, so residuals are not independent: 88.9–98.1% of links are called significant
naively against 53.7–75.9% under the valid test, an inflation of 1.27× (birds), 1.36×
(mammals), 1.55× (amphibians) and 1.72× (reptiles) — worst in the ectotherm classes.
**c**, Links surviving the max test per predictor, separated into single traits and composite
axes and expressed as proportions because the two have different denominators (24 tests per
trait bar, 12 per axis bar).

### Fig. 6 | Endotherms and ectotherms differ in trait–environment coupling

**a**, Coupling strength — the mean absolute standardised environmental effect — for each of
33 families that meet the inclusion thresholds (≥20 Chinese species and ≥100 occupied cells),
coloured by class, with boxes summarising the two thermoregulatory modes. Ectotherm families
are 1.72× more tightly coupled (Wilcoxon *P* = 0.019). **b**, The same contrast resolved by
environmental axis; ectotherms exceed endotherms on every axis. Greyscale distinguishes the
thermoregulatory modes so that no colour means both a class and a mode within this figure.
**c**, Family-level effect of thermal energy on the community-weighted mean of body mass
(point estimate ± 95% CI), the single strongest contrast. Body mass is used for all four
classes here so that the cross-class comparison concerns one named quantity; substituting
length for mass in amphibians — the only class where mass is mostly inferred — moves the
overall ratio from 1.72 to 1.67 and the thermal ratio from 1.79 to 1.77 (Supplementary
Table 19).

*Families are the replicate unit because classes provide only two replicates per
thermoregulatory mode and grid cells are pseudoreplicated. Estimates come from
family-specific spatial autoregressive error models on the 50-km grid.*

### Fig. 7 | No single class is a reliable surrogate for another

**a–c**, *k*-means partition of grid cells on their profile of trait-moment standardised
effect sizes, with the number of regions chosen by average silhouette width. The dashed line
is the Hu Huanyong (Heihe–Tengchong) line, an independent population-geographic divide that
was **not** an input to the clustering; each class independently recovers a north-west versus
south-east division agreeing with it in 76–78% of cells. Mammals additionally separate the
Tibetan Plateau and the north-east. **d**, Adjusted Rand index between class-specific
partitions. The bracketed value forces every class to two regions, which separates
"a different boundary" from "a different number of regions". **e**, Cluster centres, giving
each region an interpretable trait signature. Amphibians are omitted because too few cells
pass the five-species threshold for a stable partition (Fig. 3).

### Fig. 8 | Conservation design: complementarity, under-protected traits and 30 × 30

**a**, **b**, Accumulation curves for target-based greedy complementarity, in which a feature
counts as covered once 30% of its range is protected — the Kunming–Montreal Target 3
criterion — for functional bins (**a**) and species (**b**). Grey band: 95% of random
selections. The existing reserve network (orange point) covers 10.36% of land and brings only
1.5% of functional bins and 2.9% of species to target; an optimal layout of the same area
reaches 48.0% and 50.7%. **c**, Effect of each functional axis on the protected fraction of a
species' range (quasi-binomial GLM weighted by range area, controlling range size). Negative
means the high end of that axis is under-protected. Body size reverses between mammals
(+0.283) and reptiles (−0.308); slow pace and broad niche are under-protected consistently.
Note that the niche axis runs narrow → broad, so its negative coefficients identify
**wide-ranging generalists**, not specialists. Asterisks mark BH-adjusted *P* < 0.05.
**d**, Percentage of species reaching the 30%-of-range target under the current network and
under three 30%-of-land strategies. **e**, Priority cells selected at a 30% budget under the
species and trait-space targets, and their overlap.

---

## Extended Data figure legends

**ED 1 | Trait coverage and imputation.** **a**, Observed coverage of each core and extended
trait before imputation. **b**, Assignment of traits to the core set (complete in all four
classes), the extended set (birds, mammals and reptiles only) and the excluded set
(data-deficient). **c**, Values already imputed upstream by TetrapodTraits, using the flags
supplied there — including the 74.7% of amphibian body masses derived allometrically.
**d**, Out-of-bag error of the missForest imputation, which included the ten leading
phylogenetic eigenvectors as predictors.

**ED 2 | Trait space and its principal coordinates.** Principal-coordinate ordination of the
Gower distance among the core traits, one panel per class. Arrows are trait–axis correlations;
percentages are variance explained.

**ED 3 | Environmental predictors and their axes.** Loadings of individual variables on each
ecological axis, the spatial pattern of each axis on the 50-km grid, and the correlations
among axes confirming that axis construction removed the collinearity present among the raw
variables.

**ED 4 | Spatial model diagnostics.** Moran's *I* of ordinary-least-squares residuals (up to
0.69, significant in every model, which is why spatial models are required); the fitted
autoregressive parameter λ; and the AIC improvement over ordinary least squares. Residual
Moran's *I* after correction falls to −0.055 to 0.016.

**ED 5 | Full trait-by-predictor effect matrix.** Standardised effects of every predictor on
every moment (mean, variance, skewness, kurtosis) of every analysed trait, by class. Dots mark
BH-adjusted *P* < 0.05.

**ED 6 | Type I error inflation of naive CWM regression.** Fourth-corner statistics by class;
naive against max-test *P* values, where points in the lower-right quadrant are links the
naive test calls significant and the valid test does not; and the resulting inflation ratio.

**ED 7 | Robustness to grain: 100 km results.** The main assemblage results replicated on the
100-km grid, with a cell-wise comparison between grains exploiting the exact nesting of the
50-km grid within the 100-km grid.

**ED 8 | Null-model behaviour under occupancy weighting.** Standardised effect size against
species richness — a flat relationship indicates the richness-constrained null has removed
the metric's mechanical dependence on richness — with observed values for comparison. Species
are drawn with probability proportional to occupancy rather than equiprobably; the
equiprobable null treats narrow endemics and widespread species as exchangeable and
manufactures apparent trait clustering in endemic-rich taxa. Moving to occupancy weighting
changed the amphibian standardised effect from +1.70 to +0.08 and removed an apparent
over-dispersion signal.

**ED 9 | Additional robustness checks.** Ectotherm-to-endotherm coupling ratio when each class
is dropped in turn; mean standardised effect computed from all species against only those with
fully observed core traits; and the percentage of species per class with fully observed traits.

**ED 10 | Standardised effect sizes of functional diversity.** Functional dispersion relative
to the occupancy-weighted null for each class, and its response to mean annual temperature and
to human footprint.

**ED 11 | The anthropogenic filter on assemblage traits.** Human Footprint on the 50-km grid,
and the community-weighted mean, variance and skewness of body size along that gradient, by
class.

**ED 12 | Grain sensitivity and cross-class congruence.** Distribution of standardised effect
sizes at 50 and 100 km, and the cross-class correlations of functional dispersion. Congruence
is highest between the two endotherm classes.

**ED 13 | Geography of assemblage trait variance.** Standardised effect of community-weighted
trait *variance*, the correlation between the mean's and the variance's effect sizes for each
trait, and a filtering-regime map classifying cells by the signs of both.

**ED 14 | Which trait is most strongly filtered, mapped.** For each cell, the core trait with
the largest absolute standardised effect of its community-weighted mean; the composition of
dominant traits overall; and the same composition within deciles of the thermal-energy axis.

**ED 15 | Geographically weighted trait–environment coupling.** Local coupling estimated in a
moving window of the 200 nearest cells, its difference between reptiles and the endotherm
mean, and a diagnostic of local coupling against cell richness.

> These local estimates must not be read against Fig. 6. The family-level analysis is free of
> the attenuation shown in the diagnostic panel; the local estimates are not. The maps are
> reported as a description of spatial heterogeneity, the diagnostic as a methodological
> caution, not as a biological contrast.

**ED 16 | Functional beta diversity and its turnover component.** Sørensen beta diversity
between adjacent cell pairs, decomposed following Baselga, computed for species and for
functional groups; functional against taxonomic beta diversity, whose vertical distance below
the 1:1 line is functional redundancy; and the share of beta diversity attributable to
turnover rather than nestedness.

**ED 17 | Occurrence records versus expert range maps.** Chao sample coverage per cell from
GBIF/eBird records; richness from occurrences against richness from range maps; and the
standardised effect of human pressure on assemblage metrics fitted to range-map assemblages,
occurrence assemblages, and occurrence assemblages with coverage as a covariate.

**ED 18 | Occurrence-based assemblages for all four classes.** Chao coverage per class, and
the human-pressure effect under four settings, including the decisive "range map, same cells"
control that separates a data-source effect from a cell-selection effect.

**ED 19 | Conservation gaps and reserve representation.** Based on China's 1,028 designated
nature reserves; other protected-area categories are excluded, so coverage is a lower bound.
Reserve coverage and distance to the nearest reserve; conservation-gap cells for each of the
four classes; the effect of assemblage trait volume on reserve coverage from spatial models
controlling the environmental axes; Jaccard overlap of gap cells between classes; and the
percentage of reserves whose designation record names each class as a protection target
(mammals 17.9%, birds 12.8%, amphibians 3.2%, reptiles 1.1%).

**ED 20 | Sexual size dimorphism across assemblages.** Species-level dimorphism from the
Chinese trait compilations, and its assemblage-level geography. Reptile dimorphism comes
entirely from the snake dataset — the lizard compilation does not separate the sexes — and is
labelled as snakes rather than reptiles throughout.

**ED 21 | Data-driven trait axes retained as a check on the *a priori* axes.** Free
within-class principal components, shown to document how far the leading data-driven axis
differs in meaning between classes, which is the reason the main analysis fixes the axes
*a priori*.

**ED 22 | Three ways of measuring habitat heterogeneity, and the dose–response with thermal
content.** **a**, Ectotherm/endotherm effect ratio for ambient energy and for habitat
heterogeneity under three heterogeneity proxies (topographic elevation range, land-cover type
diversity, EarthEnv remote-sensing texture). The energy ratio is stable at 1.98–2.63; the
heterogeneity ratio is unstable and straddles 1. Bars are 95% bootstrap intervals over
class × response fits. **b**, Correlations among the three proxies. **c**, The replacement
test: the ectotherm excess plotted against how much temperature information each predictor
carries (|Spearman *r*| with mean annual temperature) across eight predictors.
Spearman ρ = 0.74, *P* = 0.037. The single heterogeneity proxy that shows an ectotherm excess
is the one most strongly correlated with temperature — it is a disguised thermal variable,
which is why the single negative control failed and was replaced by a trend test.

**ED 23 | Which families enter the clade-level test, and how the threshold changes it.**
Family inclusion under three threshold pairs and the resulting ectotherm/endotherm ratio:
1.79× at the main threshold (*P* = 0.017), 1.71× at ≥30 species (*P* = 0.049), and 1.33× and
non-significant at ≥10 species (*P* = 0.23), where 14 additional small families with noisier
coupling estimates enter.

**ED 24 | Functional-diversity indices and unconstrained NMDS ordination.** Functional
evenness, dispersion and Rao's quadratic entropy against the same six hypotheses, showing that
the sign of the effect is consistent across indices; and non-metric multidimensional scaling
of assemblage trait composition with environmental vectors fitted afterwards. Functional
richness is not shown: the convex-hull computation does not scale to the species numbers here.

**ED 25 | Land-use tolerance: which traits tolerate modification, and are they the protected
ones.** Species-level tolerance of artificial habitat against the three functional axes;
assemblage-mean tolerance against human footprint; and protected fraction against tolerance.

---


## References

Angilletta, M.J. (2009) *Thermal Adaptation: A Theoretical and Empirical
Synthesis*. Oxford University Press, Oxford.

Araújo, M.B., Ferri-Yáñez, F., Bozinovic, F., Marquet, P.A., Valladares, F. &
Chown, S.L. (2013) Heat freezes niche evolution. *Ecology Letters* 16,
1206–1219.

Baselga, A. (2010) Partitioning the turnover and nestedness components of beta
diversity. *Global Ecology and Biogeography* 19, 134–143.

Bielby, J., Mace, G.M., Bininda-Emonds, O.R.P., Cardillo, M., Gittleman, J.L.,
Jones, K.E., Orme, C.D.L. & Purvis, A. (2007) The fast–slow continuum in
mammalian life history: an empirical reevaluation. *The American Naturalist*
169, 748–757.

BirdLife International & Handbook of the Birds of the World (2022) Bird
species distribution maps of the world, version 2022.2. BirdLife
International, Cambridge.

Böhm, M. et al. (2013) The conservation status of the world's reptiles.
*Biological Conservation* 157, 372–385.

Botta-Dukát, Z. (2005) Rao's quadratic entropy as a measure of functional
diversity based on multiple traits. *Journal of Vegetation Science* 16,
533–540.

Bruelheide, H. et al. (2018) Global trait–environment relationships of plant
communities. *Nature Ecology & Evolution* 2, 1906–1917.

Brum, F.T., Graham, C.H., Costa, G.C., Hedges, S.B., Penone, C., Radeloff,
V.C., Rondinini, C., Loyola, R. & Davidson, A.D. (2017) Global priorities for
conservation across multiple dimensions of mammalian diversity. *Proceedings
of the National Academy of Sciences USA* 114, 7641–7646.

Buckley, L.B., Hurlbert, A.H. & Jetz, W. (2012) Broad-scale ecological
implications of ectothermy and endothermy in changing environments. *Global
Ecology and Biogeography* 21, 873–885.

Cadotte, M.W. & Tucker, C.M. (2017) Should environmental filtering be
abandoned? *Trends in Ecology & Evolution* 32, 429–437.

Cardillo, M., Mace, G.M., Jones, K.E., Bielby, J., Bininda-Emonds, O.R.P.,
Sechrest, W., Orme, C.D.L. & Purvis, A. (2005) Multiple causes of high
extinction risk in large mammal species. *Science* 309, 1239–1241.

Carmona, C.P. et al. (2021) Erosion of global functional diversity across the
tree of life. *Science Advances* 7, eabf2675.

CBD (2022) Kunming–Montreal Global Biodiversity Framework. Decision
CBD/COP/DEC/15/4. Convention on Biological Diversity, Montreal.

Chichorro, F., Juslén, A. & Cardoso, P. (2019) A review of the relation
between species traits and extinction risk. *Biological Conservation* 237,
220–229.

Cooke, R.S.C., Eigenbrod, F. & Bates, A.E. (2019) Projected losses of global
mammal and bird ecological strategies. *Nature Communications* 10, 2279.

Cox, N. et al. (2022) A global reptile assessment highlights shared
conservation needs of tetrapods. *Nature* 605, 285–290.

Currie, D.J. et al. (2004) Predictions and tests of climate-based hypotheses
of broad-scale variation in taxonomic richness. *Ecology Letters* 7,
1121–1134.

Deutsch, C.A., Tewksbury, J.J., Huey, R.B., Sheldon, K.S., Ghalambor, C.K.,
Haak, D.C. & Martin, P.R. (2008) Impacts of climate warming on terrestrial
ectotherms across latitude. *Proceedings of the National Academy of Sciences
USA* 105, 6668–6672.

Díaz, S. et al. (2016) The global spectrum of plant form and function.
*Nature* 529, 167–171.

Ding, C., Hu, Y., Li, C. & Jiang, Z. (2022) A dataset on the morphological,
life-history and ecological traits of the mammals in China. *Biodiversity
Science* 30, 21520. [丁晨晨等]

Dormann, C.F. et al. (2007) Methods to account for spatial autocorrelation in
the analysis of species distributional data: a review. *Ecography* 30,
609–628.

Dray, S. & Dufour, A.-B. (2007) The ade4 package: implementing the duality
diagram for ecologists. *Journal of Statistical Software* 22(4), 1–20.

Dray, S. & Legendre, P. (2008) Testing the species traits–environment
relationships: the fourth-corner problem revisited. *Ecology* 89, 3400–3412.

Etard, A., Morrill, S. & Newbold, T. (2020) Global gaps in trait data for
terrestrial vertebrates. *Global Ecology and Biogeography* 29, 2143–2158.

Etard, A. & Newbold, T. (2024) Species-level correlates of land-use responses
and climate-change sensitivity in terrestrial vertebrates. *Conservation
Biology* 38, e14208.

Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial resolution
climate surfaces for global land areas. *International Journal of Climatology*
37, 4302–4315.

Gotelli, N.J. (2000) Null model analysis of species co-occurrence patterns.
*Ecology* 81, 2606–2621.

Hawkins, B.A. et al. (2003) Energy, water, and broad-scale geographic patterns
of species richness. *Ecology* 84, 3105–3117.

Hawkins, B.A., Leroy, B., Rodríguez, M.Á., Singer, A., Vilela, B., Villalobos,
F., Wang, X. & Zelený, D. (2017) Structural bias in aggregated species-level
variables driven by repeated species co-occurrences: a pervasive problem in
community and assemblage data. *Journal of Biogeography* 44, 1199–1211.

Healy, K., Ezard, T.H.G., Jones, O.R., Salguero-Gómez, R. & Buckley, Y.M.
(2019) Animal life history is shaped by the pace of life and the distribution
of age-specific mortality and reproduction. *Nature Ecology & Evolution* 3,
1217–1224.

Huang, M., Huang, G., Fan, H. & Wei, F. (2023) Influence of Last Glacial
Maximum legacies on functional diversity and community assembly of extant
Chinese terrestrial vertebrates. *The Innovation* 4, 100379.

Huey, R.B. & Kingsolver, J.G. (1989) Evolution of thermal sensitivity of
ectotherm performance. *Trends in Ecology & Evolution* 4, 131–135.

IUCN (2024) The IUCN Red List of Threatened Species, spatial data. IUCN,
Gland. https://www.iucnredlist.org

Janzen, D.H. (1967) Why mountain passes are higher in the tropics. *The
American Naturalist* 101, 233–249.

Jenkins, C.N., Pimm, S.L. & Joppa, L.N. (2013) Global patterns of terrestrial
vertebrate diversity and conservation. *Proceedings of the National Academy of
Sciences USA* 110, E2602–E2610.

Kearney, M. & Porter, W. (2009) Mechanistic niche modelling: combining
physiological and spatial data to predict species' ranges. *Ecology Letters*
12, 334–350.

Keddy, P.A. (1992) Assembly and response rules: two goals for predictive
community ecology. *Journal of Vegetation Science* 3, 157–164.

Kerr, J.T. & Packer, L. (1997) Habitat heterogeneity as a determinant of
mammal species richness in high-energy regions. *Nature* 385, 252–254.

Khaliq, I., Hof, C., Prinzinger, R., Böhning-Gaese, K. & Pfenninger, M. (2014)
Global variation in thermal tolerances and vulnerability of endotherms to
climate change. *Proceedings of the Royal Society B* 281, 20141097.

Kissling, W.D. & Carl, G. (2008) Spatial autocorrelation and the selection of
simultaneous autoregressive models. *Global Ecology and Biogeography* 17,
59–71.

Kraft, N.J.B., Adler, P.B., Godoy, O., James, E.C., Fuller, S. & Levine, J.M.
(2015) Community assembly, coexistence and the environmental filtering
metaphor. *Functional Ecology* 29, 592–599.

Kukkala, A.S. & Moilanen, A. (2013) Core concepts of spatial prioritisation in
systematic conservation planning. *Biological Reviews* 88, 443–464.

Laliberté, E. & Legendre, P. (2010) A distance-based framework for measuring
functional diversity from multiple traits. *Ecology* 91, 299–305.

Lamanna, C. et al. (2014) Functional trait space and the latitudinal diversity
gradient. *Proceedings of the National Academy of Sciences USA* 111,
13745–13750.

Lavorel, S. & Garnier, E. (2002) Predicting changes in community composition
and ecosystem functioning from plant traits: revisiting the Holy Grail.
*Functional Ecology* 16, 545–556.

Legendre, P. (1993) Spatial autocorrelation: trouble or new paradigm?
*Ecology* 74, 1659–1673.

Luedtke, J.A. et al. (2023) Ongoing declines for the world's amphibians in the
face of emerging threats. *Nature* 622, 308–314.

Margules, C.R. & Pressey, R.L. (2000) Systematic conservation planning.
*Nature* 405, 243–253.

Mason, N.W.H., Mouillot, D., Lee, W.G. & Wilson, J.B. (2005) Functional
richness, functional evenness and functional divergence: the primary
components of functional diversity. *Oikos* 111, 112–118.

Maxwell, S.L. et al. (2020) Area-based conservation in the twenty-first
century. *Nature* 586, 217–227.

McGill, B.J., Enquist, B.J., Weiher, E. & Westoby, M. (2006) Rebuilding
community ecology from functional traits. *Trends in Ecology & Evolution* 21,
178–185.

Meiri, S. & Dayan, T. (2003) On the validity of Bergmann's rule. *Journal of
Biogeography* 30, 331–351.

Moura, M.R. & Jetz, W. (2021) Shortfalls and opportunities in terrestrial
vertebrate species discovery. *Nature Ecology & Evolution* 5, 631–639.

Moura, M.R. et al. (2024) A phylogeny-informed characterisation of global
tetrapod traits addresses data gaps and biases. *PLoS Biology* 22, e3002658.

Mu, H. et al. (2022) A global record of annual terrestrial Human Footprint
dataset from 2000 to 2018. *Scientific Data* 9, 176.

Nemhauser, G.L., Wolsey, L.A. & Fisher, M.L. (1978) An analysis of
approximations for maximizing submodular set functions—I. *Mathematical
Programming* 14, 265–294.

Newbold, T. et al. (2015) Global effects of land use on local terrestrial
biodiversity. *Nature* 520, 45–50.

O'Brien, E.M. (1998) Water–energy dynamics, climate, and prediction of woody
plant species richness: an interim general model. *Journal of Biogeography*
25, 379–398.

Olalla-Tárraga, M.Á., Rodríguez, M.Á. & Hawkins, B.A. (2006) Broad-scale
patterns of body size in squamate reptiles of Europe and North America.
*Journal of Biogeography* 33, 781–793.

Oliveira, B.F., São-Pedro, V.A., Santos-Barrera, G., Penone, C. & Costa, G.C.
(2017) AmphiBIO, a global database for amphibian ecological traits.
*Scientific Data* 4, 170123.

Olson, V.A., Davies, R.G., Orme, C.D.L., Thomas, G.H., Meiri, S., Blackburn,
T.M., Gaston, K.J., Owens, I.P.F. & Bennett, P.M. (2009) Global biogeography
and ecology of body size in birds. *Ecology Letters* 12, 249–259.

Oskyrko, O., Mi, C., Meiri, S. & Du, W. (2024) ReptTraits: a comprehensive
dataset of ecological traits in reptiles. *Scientific Data* 11, 243.

Penone, C., Davidson, A.D., Shoemaker, K.T., Di Marco, M., Rondinini, C.,
Brooks, T.M., Young, B.E., Graham, C.H. & Costa, G.C. (2014) Imputation of
missing data in life-history trait datasets: which approach performs the best?
*Methods in Ecology and Evolution* 5, 961–970.

Peres-Neto, P.R., Dray, S. & ter Braak, C.J.F. (2017) Linking trait variation
to the environment: critical issues with community-weighted mean correlation
resolved by the fourth-corner approach. *Ecography* 40, 806–816.

Pinsky, M.L., Eikeset, A.M., McCauley, D.J., Payne, J.L. & Sunday, J.M. (2019)
Greater vulnerability to warming of marine versus terrestrial ectotherms.
*Nature* 569, 108–111.

Pollock, L.J., Thuiller, W. & Jetz, W. (2017) Large conservation gains
possible for global biodiversity facets. *Nature* 546, 141–144.

Promislow, D.E.L. & Harvey, P.H. (1990) Living fast and dying young: a
comparative analysis of life-history variation among mammals. *Journal of
Zoology* 220, 417–437.

R Core Team (2025) *R: A Language and Environment for Statistical Computing*.
R Foundation for Statistical Computing, Vienna.

Rodrigues, A.S.L. & Brooks, T.M. (2007) Shortcuts for biodiversity
conservation planning: the effectiveness of surrogates. *Annual Review of
Ecology, Evolution, and Systematics* 38, 713–737.

Rodrigues, A.S.L. & Gaston, K.J. (2001) How large do reserve networks need to
be? *Ecology Letters* 4, 602–609.

Rodrigues, A.S.L. et al. (2004) Effectiveness of the global protected area
network in representing species diversity. *Nature* 428, 640–643.

Roll, U. et al. (2017) The global distribution of tetrapods reveals a need for
targeted reptile conservation. *Nature Ecology & Evolution* 1, 1677–1682.

Salguero-Gómez, R., Jones, O.R., Jongejans, E., Blomberg, S.P., Hodgson, D.J.,
Mbeau-Ache, C., Zuidema, P.A., de Kroon, H. & Buckley, Y.M. (2016) Fast–slow
continuum and reproductive strategies structure plant life-history variation
worldwide. *Proceedings of the National Academy of Sciences USA* 113, 230–235.

Slatyer, R.A., Hirst, M. & Sexton, J.P. (2013) Niche breadth predicts
geographical range size: a general ecological pattern. *Ecology Letters* 16,
1104–1114.

Song, Y. et al. (2022) A dataset on the life-history and ecological traits
of Chinese amphibians. *Biodiversity Science* 30, 22053. [宋云枫等]

Soria, C.D., Pacifici, M., Di Marco, M., Stephen, S.M. & Rondinini, C. (2021)
COMBINE: a coalesced mammal database of intrinsic and extrinsic traits.
*Ecology* 102, e03344.

Stearns, S.C. (1992) *The Evolution of Life Histories*. Oxford University
Press, Oxford.

Stekhoven, D.J. & Bühlmann, P. (2012) MissForest — non-parametric missing
value imputation for mixed-type data. *Bioinformatics* 28, 112–118.

Stevens, G.C. (1989) The latitudinal gradient in geographical range: how so
many species coexist in the tropics. *The American Naturalist* 133, 240–256.

Stuart, S.N., Chanson, J.S., Cox, N.A., Young, B.E., Rodrigues, A.S.L.,
Fischman, D.L. & Waller, R.W. (2004) Status and trends of amphibian declines
and extinctions worldwide. *Science* 306, 1783–1786.

Sun, Z.-J., Sun, B.-J., Wang, Y.-P. et al. (2025) Environmental and
evolutionary factors jointly shape life-history trait diversity of terrestrial
vertebrates across China. *Zoological Research* 46, 983–995.

Sunday, J.M., Bates, A.E. & Dulvy, N.K. (2011) Global analysis of thermal
tolerance and latitude in ectotherms. *Proceedings of the Royal Society B*
278, 1823–1830.

Sunday, J.M., Bates, A.E., Kearney, M.R., Colwell, R.K., Dulvy, N.K.,
Longino, J.T. & Huey, R.B. (2014) Thermal-safety margins and the necessity of
thermoregulatory behavior across latitude and elevation. *Proceedings of the
National Academy of Sciences USA* 111, 5610–5615.

Swenson, N.G., Enquist, B.J., Pither, J. et al. (2012) The biogeography and
filtering of woody plant functional diversity in North and South America.
*Global Ecology and Biogeography* 21, 798–808.

ter Braak, C.J.F., Cormont, A. & Dray, S. (2012) Improved testing of species
traits–environment relationships in the fourth-corner problem. *Ecology* 93,
1525–1534.

Tobias, J.A. et al. (2022) AVONET: morphological, ecological and geographical
data for all birds. *Ecology Letters* 25, 581–597.

Toussaint, A., Brosse, S., Bueno, C.G., Pärtel, M., Tamme, R. & Carmona, C.P.
(2021) Extinction of threatened vertebrates will lead to idiosyncratic spatial
and taxonomic changes in functional diversity. *Nature Communications* 12,
5162.

Tuanmu, M.-N. & Jetz, W. (2015) A global, remote sensing-based
characterization of terrestrial habitat heterogeneity for biodiversity and
ecosystem modelling. *Global Ecology and Biogeography* 24, 1329–1339.

Venter, O. et al. (2016) Sixteen years of change in the global terrestrial
human footprint and implications for biodiversity conservation. *Nature
Communications* 7, 12558.

Villéger, S., Mason, N.W.H. & Mouillot, D. (2008) New multidimensional
functional diversity indices for a multifaceted framework in functional
ecology. *Ecology* 89, 2290–2301.

Villéger, S., Grenouillet, G. & Brosse, S. (2013) Decomposing functional
β-diversity reveals that low functional β-diversity is driven by low
functional turnover in European fish assemblages. *Global Ecology and
Biogeography* 22, 671–681.

Violle, C., Navas, M.-L., Vile, D., Kazakou, E., Fortunel, C., Hummel, I. &
Garnier, E. (2007) Let the concept of trait be functional! *Oikos* 116,
882–892.

Violle, C., Reich, P.B., Pacala, S.W., Enquist, B.J. & Kattge, J. (2014) The
emergence and promise of functional biogeography. *Proceedings of the National
Academy of Sciences USA* 111, 13690–13696.

Wang, J. et al. (2023) A dataset of the morphological, life-history, and
ecological traits of snakes in China. *Biodiversity Science* 31, 23126.
[王江等]

Wright, D.H. (1983) Species-energy theory: an extension of species-area
theory. *Oikos* 41, 496–506.

Wu, R., Zhang, S., Yu, D.W., Zhao, P., Li, X., Wang, L., Yu, Q., Ma, J., Chen,
A. & Long, Y. (2011) Effectiveness of China's nature reserves in representing
ecological diversity. *Frontiers in Ecology and the Environment* 9, 383–389.

Xu, W. et al. (2017) Strengthening protected areas for biodiversity and
ecosystem services in China. *Proceedings of the National Academy of Sciences
USA* 114, 1601–1606.

Zelený, D. (2018) Which results of the standard test for community-weighted
mean approach are too optimistic? *Journal of Vegetation Science* 29, 953–966.

Zhong, Y. et al. (2022) A dataset on the life-history and ecological traits
of Chinese lizards. *Biodiversity Science* 30, 22071. [钟雨茜等]
