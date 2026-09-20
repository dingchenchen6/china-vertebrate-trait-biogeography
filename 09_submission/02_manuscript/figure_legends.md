# Figure legends

Numbering follows `table_00_figure_manifest.csv`. All figures are supplied as vector PDF
and 600 dpi PNG at 183 mm width (Nature double column), Helvetica throughout.

**Colour convention.** Hue encodes taxonomic class (birds red, mammals orange, reptiles
blue, amphibians green); the endotherm/ectotherm contrast is shown in greyscale so that no
colour carries two meanings within a figure. Diverging quantities use `vik` and sequential
quantities `batlow` (Crameri's perceptually uniform, colour-vision-deficiency-safe scales).

**Maps.** All maps of China use a standard provincial basemap in China Albers equal-area
conic projection with the South China Sea islands drawn at their true geographic positions
and the ten-segment line as a separate layer.

**Trait standardisation.** Traits are z-standardised *within each class* before assemblage
moments are computed, so every community-weighted mean is in within-class standard-deviation
units. Values are comparable across traits and classes but not to raw trait units.

---

## Main figures

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

## Extended Data figures

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

## Superseded

`06_manuscript/figure_captions.md` predates the null-model correction and the final figure
numbering. It is retained only as a record of the earlier draft and must not be used for
submission.
