# 图件说明 / Figure captions

所有图件同时输出矢量 PDF 与 600 dpi PNG，尺寸按 Nature 规格（双栏 183 mm，字号 7 pt，Helvetica）。
所有中国地图基于自然资源部标准地图 **GS(2023)2767**，国界与南海诸岛附图按官方底图原样呈现。

All figures are supplied as vector PDF and 600 dpi PNG at Nature specifications (183 mm double column, 7 pt Helvetica). All maps of China are based on the officially approved standard map GS(2023)2767; national boundaries and the South China Sea inset are reproduced without modification.

---

## 正文图 / Main figures

### Fig. 1 | Conceptual framework and data overview
**a**, The thermal-buffering hypothesis. Because endothermy decouples body temperature from ambient conditions, the set of trait values that can persist in an endothermic assemblage should be less tightly constrained by the environment than in an ectothermic assemblage; coupling strength is the slope of assemblage trait composition on the environmental gradient (bars at right). **b**, Number of species analysed per class (2,837 in total) on the 50-km equal-area grid (3,814 cells). **c–f**, Species richness of birds, mammals, reptiles and amphibians. Grey indicates cells with no species of that class; the notch in the north-east is the Bohai Sea.

*Data: IUCN Red List spatial data (amphibians, reptiles, mammals) and the BirdLife/IUCN-derived grid matrix of Huang et al. (2023).*

---

### Fig. 2 | Community-weighted trait means across China
Community-weighted means of three core traits — body mass (**a–d**), habitat breadth (**e–h**) and nocturnality (**i–l**) — for birds, mammals, reptiles and amphibians on the 50-km grid. Traits were z-standardised **within each class** before the moments were computed, so values are in within-class standard-deviation units and are comparable across traits and classes but not to raw trait units. Each row shares one colour scale.

---

### Fig. 3 | Functional structure relative to null expectation
**a–d**, Standardised effect size (SES) of functional dispersion for each class, all four panels on a single shared colour scale. Blue indicates functional clustering (assemblages more similar than random draws of the same richness from the class species pool), brown indicates over-dispersion. Birds, mammals and reptiles are clustered (mean SES −1.68, −1.26, −1.88) and remain so when only species with fully observed traits are used. The amphibian panel is shown for completeness only: 84% of amphibian core trait values are imputed upstream and the apparent over-dispersion reverses sign in the observed-only subset (Extended Data Fig. 9), so no assembly signal is claimed for that class. **e**, **f**, Response of SES functional dispersion to mean annual temperature and to human footprint (generalised additive smooths with 95% confidence bands).

*Null model: 999 randomisations holding richness constant, computed once per unique richness value.*

---

### Fig. 4 | Trait–environment coupling is stronger in ectotherms
**a**, Coupling strength (mean absolute standardised environmental effect) for each of 33 families, coloured by class; boxes summarise the two thermoregulatory modes. Ectotherm families are 1.72× more tightly coupled (Wilcoxon *P* = 0.010). **b**, The same contrast resolved by environmental axis; ectotherms exceed endotherms on every axis. **c**, Family-level effect of thermal energy on community-weighted body mass (point estimate ± 95% CI), the single strongest contrast (ectotherm 0.348 versus endotherm 0.148; *t* = 3.91, BH-adjusted *P* = 0.032).

*Families are the replicate unit because classes provide only two replicates per thermoregulatory mode and grid cells are pseudoreplicated. Estimates come from family-specific spatial simultaneous autoregressive error models on the 50-km grid.*

---

### Fig. 5 | Assemblage trait moments along the human-pressure gradient
**a**, Human Footprint in 2020 on the 50-km grid. **b–d**, Community-weighted mean, variance and skewness of body mass as a function of human footprint, by class (generalised additive smooths with 95% confidence bands). Mammals lose large-bodied species (declining mean); birds show a contraction of body-mass variance, i.e. functional homogenisation; reptiles show increasing right skew, indicating truncation of large-bodied species (β = 0.046, SE = 0.018, BH-adjusted *P* = 0.030).

*Because expert range maps are insensitive to fine-scale land-use change, the anthropogenic effects shown here are conservative.*

---

### Fig. 6 | Grain invariance and cross-taxon congruence
**a**, Distribution of SES functional dispersion at 50 km and 100 km; the assembly signal is essentially unchanged, so the results are not an artefact of the modifiable areal unit problem. **b**, Spearman correlations of SES functional dispersion among classes across grid cells. Congruence is substantial only between the two endothermic classes (birds–mammals *r* = 0.63); all other pairs, including reptiles–amphibians (*r* = −0.07), are negligible. Avian functional diversity is therefore an informative conservation surrogate for mammals but not for reptiles or amphibians.

---

### Fig. 7 | Geography of assemblage trait variance
**a–d**, Standardised effect size of community-weighted body-mass **variance** on the 50-km grid, all four panels on one shared scale. Blue indicates a narrower spread of body masses than random draws of the same richness; brown a wider spread. **e**, Spearman correlation between the SES of the mean and the SES of the variance for each core trait. Range size is the only trait whose mean and variance move in opposite directions in every class (*r* = −0.50 to −0.90): cells dominated by wide-ranging species also show the most convergent range sizes. **f**, Filtering-regime map for avian body mass, classifying each cell by the signs of its mean-SES and variance-SES; this combines the direction of filtering with its strength, which neither a mean map nor a variance map conveys alone. The amphibian panel is shown for comparison only (see Extended Data Fig. 9).

---

### Fig. 8 | The identity of the filtered trait varies geographically
**a–c**, For each cell, the core trait with the largest absolute SES of its community-weighted mean. Range size dominates almost everywhere in birds (93.7% of cells) and mammals (78.0%), whereas reptiles split between habitat breadth in the north and north-east (43.9%) and range size in the south (42.3%). **d**, Percentage of cells in which each trait dominates. **e**, The same composition within deciles of the thermal-energy axis: mammals switch to body length at the cold extreme (the Tibetan Plateau) and birds switch to body mass and verticality at the warm extreme, so the identity of the filter itself tracks climate. Amphibians are omitted because their traits are largely imputed.

---

### Fig. 9 | Functional regionalisation from trait structure
**a–c**, k-means partition of grid cells on their 12-dimensional profile of trait-moment SES (six core traits × mean and variance), with the number of regions chosen by average silhouette width. Dashed line: the Hu Huanyong (Heihe–Tengchong) line, an independent population-geographic divide that was **not** used as an input. All three classes independently recover a north-west versus south-east division that agrees with it in 74.5–79.7% of cells (Cohen's kappa 0.50–0.56, i.e. moderate rather than exact correspondence — the functional boundary lies south-east of the line in the south-west); mammals additionally separate the Tibetan Plateau. **d**, Adjusted Rand index between class-specific partitions; the value in brackets forces every class to two regions, separating a genuinely different boundary from a different number of regions. **e**, Cluster centres, giving each region an interpretable trait signature.

---

### Fig. 10 | Geographically weighted trait–environment coupling
**a–c**, Local strength of the coupling between community-weighted body mass and the thermal-energy axis, estimated within a moving window of the 200 nearest cells. A fixed neighbour count rather than a fixed bandwidth keeps estimates comparable between the sparsely and densely gridded parts of the country. Coupling is strongly heterogeneous in space: strongest in the Hengduan Mountains for birds, on the Tibetan Plateau and in central China for mammals, and in the south-west and south-east coastal zone for reptiles. **d**, Difference in local coupling between reptiles and the mean of birds and mammals; the difference is mapped rather than the ratio because a ratio of two noisy estimates is unstable when the denominator approaches zero. **e**, Diagnostic: local coupling plotted against the species richness of the cell. Reptile coupling rises steeply with richness (ρ = +0.41) whereas the endotherm classes do not, the signature of regression attenuation in species-poor assemblages (reptiles hold a median of 13 species per cell against 129 for birds).

> **These local estimates must not be read against Fig. 4.** The family-level analysis in Fig. 4 is free of this attenuation — coupling strength there is uncorrelated with the number of species in a family (ρ = 0.012, *P* = 0.95) — whereas the local estimates are not. Panels **a–c** are reported as a description of spatial heterogeneity; panels **d–e** are reported as a methodological diagnostic, not as a biological contrast.

---

### Fig. 11 | Functional beta diversity and functional redundancy
Sorensen beta diversity between the 7,344 rook-adjacent pairs of 50-km cells, decomposed following Baselga (2010) and computed both for species and for functional groups (k-means on the trait principal-coordinate space, K = 20). **a-c**, Species turnover. **d-f**, Functional beta diversity on a common scale. **g**, Functional against taxonomic beta diversity for every adjacent pair; the vertical distance below the 1:1 line is functional redundancy. Reptile assemblages track the 1:1 line whereas bird assemblages flatten, meaning their species change without their function changing. **h**, Ratio of functional to taxonomic beta diversity with 95% bootstrap intervals over cell pairs. The ratio rises monotonically from birds (0.19) through mammals (0.34) to reptiles (0.66); the ordering is unchanged at K = 10 and K = 30. **i**, Share of beta diversity attributable to turnover rather than nestedness. Reptile beta diversity is nestedness-dominated (28% turnover), the signature of a richness gradient produced by strong environmental filtering.

*This is evidence for the thermal-buffering hypothesis that does not depend on any regression model or standardised effect size, and is therefore immune to the attenuation described for Fig. 10.*

---

### Fig. 12 | Range-map versus occurrence-derived assemblages
**a**, Chao sample coverage per 50-km cell, computed from 4.58 million GBIF/eBird bird records for China (2000-2025). Grey cells hold no records or fall below the retention threshold of coverage 0.7 and 20 species; 940 of 3,814 cells were retained. **b**, Richness from occurrences against richness from range maps; occurrences detect 62% of mapped species and the two measures correlate only weakly (Spearman rho = 0.29), whereas the corresponding community-weighted trait means correlate far more strongly (0.57-0.70). Assemblage trait composition is therefore more robust to the choice of data source than species richness is. **c**, Standardised effect of the human-pressure axis on five assemblage metrics, fitted with the same spatial autoregressive framework to range-map assemblages, occurrence assemblages, and occurrence assemblages with Chao coverage as an additional covariate.

*Range maps understate anthropogenic filtering: the effect on community-weighted range size is 5.6x larger in occurrence data (0.100 versus 0.018), and the effect on habitat-breadth variance is present in occurrence data (0.115, P = 0.034) but absent from range maps (0.0004, P = 0.98). Controlling for sampling coverage strengthens rather than weakens these effects, which argues against a sampling-bias explanation. The wider intervals for occurrence models reflect the smaller number of usable cells, not weaker effects.*

---

### Fig. 13 | Occurrence assemblages across all four classes
**a-d**, Chao sample coverage per 50-km cell for each class, from GBIF records for China (birds 4.58 million, amphibians 102,386, mammals 70,450, reptiles 35,446). Subtitles give the retention threshold selected for each class and the resulting number of usable cells. **e**, Standardised effect of the human-pressure axis on four assemblage metrics under four settings: range-map assemblages on all cells; range-map assemblages restricted to exactly the cells the occurrence models use; occurrence assemblages; and occurrence assemblages with Chao coverage as a covariate.

*The "range map, same cells" setting is the decisive control. Cells holding occurrence records carry substantially higher human pressure than cells without (for example +0.31 versus −0.34 for birds, P < 2 × 10⁻¹⁶), so a larger effect in occurrence data could in principle reflect which cells enter the analysis rather than which data source is used. For birds — the only class with enough usable cells (804) — the range-map effect is unchanged by that restriction (0.017 versus 0.018 for range size), so the difference is genuinely attributable to the data source. For mammals, reptiles and amphibians only 163–170 cells pass even the most permissive threshold, estimates are unstable under both sources, and no per-class inference is drawn.*

---

### Fig. 14 | Conservation gaps in vertebrate functional diversity
Based on China's 1,028 designated nature reserves (2021 boundaries); other protected-area categories such as national parks are not included, so coverage is a lower bound. **a**, Fraction of each 50-km cell inside a nature reserve. Nationwide coverage is 10.31% (8.38% at national level) and 62% of cells hold no reserve at all; coverage is concentrated in the vast desert and plateau reserves of western China. **b**, Distance from each cell to the nearest reserve. **c-e**, Conservation-gap cells for each class, defined as the top quartile of assemblage trait volume combined with the bottom quartile of reserve coverage. **f**, Effect of assemblage trait volume on reserve coverage from spatial autoregressive models that also control the five environmental axes. Coverage tracks mammal trait volume (+0.106, BH-adjusted P = 0.012) but is unrelated to birds and tends negative for reptiles. **g**, Jaccard overlap of gap cells between classes; bird and reptile gaps share only 18% of cells. **h**, Percentage of reserves whose official designation record names each class as a protection target.

*Panel h supplies the mechanism for panel f: China's reserve system was built around large mammals (17.9% of reserves) and birds (12.8%), with reptiles named in only 1.1%. The representation shortfall for ectotherms is a direct consequence of designation history, and it confirms the corollary of the surrogacy result in Fig. 6 — planning from avian hotspots would systematically miss where ectotherm functional diversity is unprotected.*

---

## 扩展数据图 / Extended Data figures

### Extended Data Fig. 1 | Trait coverage and imputation
**a**, Percentage of species with observed values for each core and extended trait, before imputation. **b**, Assignment of each trait to the core set (complete in all four classes), the extended set (birds, mammals and reptiles only) or the excluded set (data-deficient). **c**, Percentage of values already imputed upstream by TetrapodTraits v2.0.1, using the flags supplied there. **d**, Out-of-bag normalised root-mean-square error of the missForest imputation, which included the first ten phylogenetic eigenvectors as predictors.

### Extended Data Fig. 2 | Class-specific trait space
Principal-coordinate ordination of the Gower distance among the six core traits, one panel per class. Arrows are trait–axis correlations. Percentages are variance explained by each axis.

### Extended Data Fig. 3 | Environmental predictor axes
**a**, Loadings of the individual variables on the first principal component of each of the five ecological axes. **b–f**, Spatial pattern of each axis score on the 50-km grid. **g**, Pearson correlations among the five axes, confirming that the axis construction removed the collinearity present among the raw variables.

### Extended Data Fig. 4 | Spatial model diagnostics
**a**, Moran's I of ordinary-least-squares residuals; spatial autocorrelation is significant in every model, which is why spatial autoregressive models are required. **b**, Fitted spatial autoregressive parameter λ. **c**, AIC improvement of the spatial model over ordinary least squares on a signed-log scale.

### Extended Data Fig. 5 | Complete standardised-effect matrix
Standardised effects of the five environmental axes on all four moments (mean, variance, skewness, kurtosis) of every analysed trait, by class, on the 50-km grid. Dots mark effects significant at BH-adjusted *P* < 0.05. Amphibian panels are blank for diet and litter size because those traits were too incomplete to retain.

### Extended Data Fig. 6 | Fourth-corner analysis and Type I error inflation
**a**, Fourth-corner statistics for every trait × environment combination, by class; dots mark associations significant under the max test. **b**, Naive community-weighted-mean regression *P* values against max-test *P* values; points in the lower-right quadrant are associations that the naive test declares significant and the correct test does not. **c**, Ratio of the number of associations declared significant by the naive test to that declared by the max test (1.27× to 2.72×).

### Extended Data Fig. 7 | Robustness to grain
**a–d**, SES functional dispersion on the 100-km grid, replicating Fig. 3a–d. **e**, Cell-wise comparison of SES functional dispersion between grains, using the exact nesting of the 50-km grid within the 100-km grid. **f**, The same comparison for species richness.

### Extended Data Fig. 8 | Null-model behaviour
**a**, SES functional dispersion against species richness; a flat relationship indicates that the richness-constrained null model has removed the mechanical dependence of the metric on richness. **b**, Observed functional dispersion against richness, for comparison. **c**, Distribution of standardised effect sizes for three multivariate metrics, with the conventional |SES| = 2 thresholds marked.


### Extended Data Fig. 9 | Robustness checks
**a**, Ectotherm-to-endotherm coupling ratio when each class is dropped in turn; filled symbols are Wilcoxon *P* < 0.05. Ratios exceed 1 in all 15 scenario × metric combinations. **b**, Mean SES of functional dispersion computed from all species versus only species whose core traits were fully observed rather than imputed upstream; only the amphibian result changes sign. **c**, Percentage of species per class with fully observed core traits.

---

## 表格 / Tables

正文表 Table 1–3 与附表 Table S1–S10 的完整标题与脚注见 [`table_captions.md`](table_captions.md)；
全部表格另汇编为单一 Excel 工作簿 `04_results/vertebrate_trait_analysis_tables.xlsx`（含目录页与逐表说明）。

Full captions and footnotes for main Tables 1–3 and Supplementary Tables S1–S10 are in `table_captions.md`; all tables are also bundled into one Excel workbook with a contents sheet.