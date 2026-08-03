# 表格说明 / Table captions

**Table 1 | Data summary for the four terrestrial vertebrate classes analysed on the 50-km equal-area grid across China.**

Species counts are those retained after IUCN filtering (presence 1-2, origin 1/2/6, seasonal 1-2, terrestrial only); bird assemblages derive from the equivalent BirdLife/IUCN-based matrix of Huang et al. (2023). Traits analysed is class-specific: maximum longevity was 19-55% incomplete in every class and was excluded throughout, and amphibian diet and litter size were additionally too incomplete to retain. Occupied cells is the number of cells holding at least one species; cells analysed is the subset with at least two species, for which assemblage metrics are defined. SES = standardised effect size against a richness-constrained null model with 999 randomisations; negative values indicate functional clustering, positive values over-dispersion.

**Table 2 | Trait-environment coupling strength in endothermic versus ectothermic families of Chinese terrestrial vertebrates.**

Coupling strength is the mean absolute standardised effect of the environmental axes on assemblage trait metrics, estimated per family from spatial simultaneous autoregressive error models on the 50-km grid (n = 33 families; 21 endothermic, 12 ectothermic). The mixed model is y ~ thermoregulatory mode + (1 | class/order); because thermoregulatory mode is nearly collinear with class, this test is deliberately conservative and the Wilcoxon rank-sum test is reported alongside. Directional consistency across the 30 individual response x predictor contrasts: 30 of 30 favoured ectotherms (binomial P = 1.863e-09); those contrasts are mutually correlated and provide descriptive support only. Significance marks: *** P < 0.001, ** P < 0.01, * P < 0.05, dagger P < 0.1.

**Table 3 | Trait-environment coupling strength for each of the 33 families used as replicates.**

Families were retained when they contained at least 20 Chinese species and occupied at least 100 grid cells. Coupling is the mean absolute standardised effect from family-specific spatial autoregressive models.

**Table S1 | Properties of the two equal-area analysis grids.**

Grids are in the China Albers equal-area conic projection matching the officially approved basemap GS(2023)2767. Cells with less than 50% land were excluded.

**Table S2 | Assemblage composition summary by class and grain.**

**Table S3 | Trait coverage before imputation, and assignment of traits to the core, extended or excluded set.**

The core set comprises the six traits complete in all four classes and underpins every cross-class contrast. Upstream imputation refers to values already imputed by TetrapodTraits v2.0.1 and flagged there.

**Table S4 | Outcome of the three-pass taxonomic harmonisation against the TetrapodTraits backbone.**

Pass 1 matched scientific names directly; pass 2 used the backbone's IUCN_Binomial field; pass 3 resolved the remainder through the GBIF backbone taxonomy.

**Table S5 | Trait imputation diagnostics.**

Imputation used missForest with the first ten phylogenetic eigenvectors and collapsed order/family factors as auxiliary predictors. Traits missing in more than 40% of species within a class were excluded rather than imputed.

**Table S6 | Summary statistics for all 47 environmental variables extracted per grid cell.**

**Table S7 | Loadings of the individual environmental variables on the five ecological predictor axes.**

Each axis is the first principal component of its z-standardised member variables, sign-aligned so that the first listed variable loads positively.

**Table S8 | Type I error inflation of naive community-weighted-mean regression relative to the fourth-corner max test.**

The max test permutes sites and species simultaneously and takes the larger of the two P values, controlling the error rate that naive CWM regression inflates because all sites share one trait vector. Conclusions reversed is the percentage of trait x environment combinations that the naive test declares significant and the max test does not. Trait sets match the main analysis exactly, so amphibians contribute 60 rather than 100 combinations.

**Table S9 | Correlations of each core trait with the first two principal-coordinate axes of the class-specific trait space.**

**Table S10 | All 30 individual response x environmental-axis contrasts between endothermic and ectothermic families.**

Ectotherm coupling exceeded endotherm coupling in all 30 contrasts. These contrasts are mutually correlated, so the headline inference is the single family-level test in Table 2 rather than these individual comparisons.

**Table S11 | Leave-one-class-out robustness of the thermoregulatory-mode contrast.**

Each class was dropped in turn and the family-level contrast refitted. The ectotherm-to-endotherm ratio exceeded 1 in all 15 scenario x metric combinations (range 1.39-2.69). Because dropping a class leaves only 5-13 families per group, loss of significance reflects loss of power rather than reversal of the effect.

**Table S12 | Sensitivity of assemblage functional structure to upstream trait imputation.**

Species were retained only if none of their core trait values carried a TetrapodTraits imputation flag. Birds, mammals and reptiles retain 81-84% of species and their mean SES of functional dispersion shifts by at most 0.21. Amphibians retain 16% and the mean SES reverses from +1.85 to -0.35, so no assembly signal is reported for that class. The observed-only subset is itself biased towards well-studied species, so this test can only show that a result may be imputation-driven.

**Table S13 | Observed trait coverage for Chinese amphibians, taken from the primary databases rather than from imputed values.**

Values were read directly from AmphiBIO v1 and Etard et al. (2020); range size is computed from the IUCN polygons and is therefore a measurement. Body mass and diel activity are observed for only 8% and 11% of species, which reflects a genuine gap in the primary literature rather than an integration failure.

**Table S14 | Quasi-factorial test of what drives the apparent over-dispersion of Chinese amphibian assemblages.**

Comparing A with D holds the species set constant and changes the trait set; comparing D with C holds the trait set constant and changes imputed for observed values. Each contrast removes about half of the apparent signal, and the cleanest estimate (C) is indistinguishable from zero.

**Table S15 | Functional and taxonomic beta diversity between adjacent 50-km cells.**

Sorensen dissimilarity decomposed following Baselga (2010). Functional groups are k-means clusters (K = 20) in the trait principal-coordinate space. Confidence intervals are 999 bootstrap resamples of cell pairs. The ordering of the ratio across classes is unchanged at K = 10 and K = 30 (Table S16).

**Table S16 | Sensitivity of the functional-to-taxonomic beta ratio to the number of functional groups.**

**Table S17 | Agreement between bird assemblages derived from range maps and from 4.58 million GBIF/eBird occurrence records.**

Occurrence assemblages were retained only for cells with Chao sample coverage of at least 0.7 and at least 20 recorded species. Trait composition agrees far better between the two sources than species richness does.

**Table S18 | Effect of the human-pressure axis on bird assemblage metrics under three assemblage definitions.**

All models are spatial simultaneous autoregressive error models with the same five environmental axes on the 50-km grid. The effect on community-weighted range size is 5.6 times larger in occurrence assemblages than in range-map assemblages, and the effect on habitat-breadth variance is absent from range maps entirely. Adding Chao coverage as a covariate strengthens rather than weakens these effects, which argues against a sampling-bias explanation.

