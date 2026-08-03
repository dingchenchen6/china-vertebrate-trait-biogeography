# Endothermy buffers assemblage trait composition against environmental filtering across Chinese terrestrial vertebrates

**Running title**: Thermoregulatory mode and trait–environment coupling

---

## Abstract

Environmental filtering is thought to shape which functional traits persist in local assemblages, yet whether the *strength* of that filtering is itself predictable from organismal physiology remains untested at biogeographic scales. Endothermy decouples body temperature from ambient conditions, and should therefore weaken the environmental control of assemblage trait composition — but this prediction has never been evaluated with a design that provides genuine replication across thermoregulatory modes. Here we assembled distributions of **2,837 terrestrial vertebrate species** (1,177 birds, 631 mammals, 520 reptiles, 509 amphibians) across China on nested **50-km and 100-km equal-area grids**, integrated a common set of six functional traits defined identically in all four classes, and quantified the first four moments of each assemblage's trait distribution together with null-model standardised functional diversity. Using **families as replicates** (n = 33) — rather than the four classes, which afford only two replicates per thermoregulatory mode, or grid cells, which are pseudoreplicated — we fitted spatial autoregressive models per family and formally contrasted endothermic and ectothermic clades. Ectotherm assemblages showed **1.72-fold stronger overall trait–environment coupling** than endotherms (Wilcoxon *P* = 0.010) and **1.68-fold stronger coupling to thermal energy** (mixed model *P* = 0.023); the contrast was directionally consistent in **30 of 30** response × predictor combinations and in all five leave-one-class-out scenarios (ratios 1.39–2.69), although individual significance depends on retaining the full family set. Contrary to the prediction that large-bodied, slow-lived endotherms should be the more sensitive to human pressure, coupling to contemporary human footprint was also **2.05-fold stronger in ectotherms** (*P* = 0.0016). Assemblages of birds, mammals and reptiles were functionally clustered relative to null expectation (SES of functional dispersion = −1.68, −1.26, −1.88), and this held when the analysis was restricted to species with fully observed traits. Amphibians could not be assessed reliably: 84% of their core trait values are imputed upstream, and their apparent over-dispersion reversed sign when restricted to the 16% of species with observed traits. Spatial congruence in functional structure was substantial only between the two endothermic classes (birds–mammals *r* = 0.63) and negligible for every other pair (|*r*| ≤ 0.18), so avian functional diversity is an informative conservation surrogate for mammals but not for reptiles or amphibians. Patterns were essentially identical at 50 and 100 km, indicating scale-invariance over this grain range. Our results indicate that endothermy buffers assemblage trait composition not only against thermal gradients but against environmental filtering in general, and that ectotherm assemblages are correspondingly the more exposed to both climatic and anthropogenic change.

**Keywords**: functional biogeography; community assembly; thermoregulation; trait–environment coupling; community-weighted moments; China; terrestrial vertebrates

---

## 1. Introduction

*(to be expanded)*

Functional-biogeographic studies routinely map community-weighted trait means and functional diversity across environmental gradients and relate them to climate, productivity and topography. Two recent syntheses have mapped the functional and phylogenetic structure of Chinese terrestrial vertebrates and related it to Quaternary climate legacies (Huang et al. 2023, *The Innovation*) and to contemporary environmental and evolutionary factors (Sun et al. 2025, *Zoological Research*). Both concluded that the environmental signal "varies among taxonomic groups", but neither derived that variation from first principles, tested it formally, or considered contemporary anthropogenic pressure.

We address three gaps.

**(i) A mechanistic, testable expectation.** Thermoregulatory mode offers a first-principles prediction: endothermy renders body temperature largely independent of ambient conditions, so the traits that can persist in an endothermic assemblage should be less tightly constrained by the local environment than those in an ectothermic assemblage. We formalise this as the **thermal-buffering hypothesis** and test it directly.

**(ii) A design with genuine replication.** Contrasting endotherms with ectotherms using classes as units gives n = 4 (two per mode); using grid cells is pseudoreplication, because cells within a class share one species pool and one trait matrix. We therefore use **families** as replicates, with class/order nested random effects to absorb phylogenetic non-independence.

**(iii) Separating directional filtering from filtering strength.** Community-weighted means capture *where* an assemblage sits in trait space; variance captures *how broadly* it is spread; skewness identifies *which tail* is truncated; kurtosis indexes niche packing. Analysing all four moments per trait separates processes that conventional multivariate functional-diversity indices confound.

---

## 2. Methods

### 2.1 Analysis units
Grid cells were generated in the China Albers equal-area conic projection (`+proj=aea +lat_0=0 +lon_0=110 +lat_1=25 +lat_2=47`) identical to the officially approved national basemap GS(2023)2767, at 50 km (n = 3,814) and 100 km (n = 957); the 50-km grid nests exactly within the 100-km grid. Cells with < 50% land were excluded, leaving 9.38 × 10⁶ km² and 9.26 × 10⁶ km² of China respectively.

### 2.2 Species distributions
Range polygons for amphibians, reptiles and mammals were taken from the IUCN Red List spatial database (version 2025-1), retaining polygons coded Extant or Probably Extant (`presence` 1–2), Native, Reintroduced or Assisted Colonisation (`origin` 1, 2, 6), and Resident or Breeding Season (`seasonal` 1–2), and restricted to terrestrial species. Bird assemblages were derived from the published grid-level presence matrix of Huang et al. (2023), itself compiled from BirdLife/IUCN range polygons under equivalent criteria, and resampled onto our grids by area-weighted overlap (species retained where occupied source cells covered ≥ 50% of a target cell).

Range polygons were projected, simplified (2-km tolerance — negligible at a 50-km grain) and rasterised onto a 5-km sub-grid; a species was scored present in a cell when its range covered ≥ 1% of that cell. This yielded 1,177 bird, 631 mammal, 520 reptile and 509 amphibian species.

### 2.3 Traits
`TetrapodTraits` v2.0.1 (Moura et al. 2024) provided a taxonomically unified backbone across all four classes, supplemented by AVONET, COMBINE, AmphiBIO, ReptTraits, Meiri (2018), GlobTherm, Etard et al. (2020), the Amniote Life History Database and AnAge. Names were harmonised in three passes (direct match → the backbone's `IUCN_Binomial` field → the GBIF backbone), matching 2,595 of 2,837 species (91.5%).

Six traits were complete in **all four classes** and form the **core set** used for every cross-class contrast: body mass, body length, nocturnality, verticality (fossorial → aerial), IUCN habitat breadth and range size. Four further traits (diet breadth, vertebrate and plant diet fractions, litter/clutch size) were adequately covered in birds, mammals and reptiles only and were used in class-specific analyses. Maximum longevity was 19–55% incomplete in every class and was **excluded** from the main analyses and reported as data-deficient.

Remaining gaps were imputed with `missForest` including the first ten phylogenetic eigenvectors as predictors (out-of-bag NRMSE 0.005–0.009). Traits were z-standardised within class before analysis, so all moments are expressed in within-class standard-deviation units.

### 2.4 Assemblage trait moments and functional diversity
For each cell × class × grain we computed, per trait, the community-weighted **mean, variance, skewness and excess kurtosis** (assemblages are presence–absence, so weights are equal), plus multivariate functional dispersion, Rao's Q, functional richness (convex-hull volume) and the trait volume and density metrics of Sun et al. (2025), from a PCoA of the Gower distance among the six core traits.

**Null model.** For every observed richness *S*, 999 assemblages of *S* species were drawn at random from the class-wide (or family-wide) pool and all metrics recomputed; standardised effect sizes are SES = (obs − mean_null)/sd_null. Because the null distribution depends only on *S*, it was computed once per unique richness value — an exact optimisation that made the 999-replicate design tractable.

### 2.5 Environmental predictors
Forty-seven variables were extracted per cell and reduced to five *a priori* ecological axes by within-group PCA (first axis; loadings in Table S7): **thermal energy** (bio1, bio6, bio4, bio3), **water availability** (bio12, bio15, bio17), **productivity** (MODIS NPP, tree cover), **habitat structure** (elevation range, EarthEnv habitat-texture Shannon, within-cell temperature SD) and **human pressure** (Human Footprint 2020, cropland and built-up cover).

### 2.6 Statistical analysis
Per class and per family, each SES metric was modelled against the five environmental axes with **spatial simultaneous autoregressive error models** (`spatialreg::errorsarlm`; symmetrised 8-nearest-neighbour weights). Spatial dependence was strong (λ ≈ 0.98), confirming the need for the spatial term.

The **headline test** collapsed each family to a single coupling-strength value — the mean absolute standardised effect across responses and predictors — and contrasted thermoregulatory modes in a linear mixed model with class/order nested random effects, with a Wilcoxon rank-sum test as a distribution-free back-up. Because thermoregulatory mode is almost collinear with class (two classes per mode), this mixed model is deliberately conservative; we therefore also report the directional consistency of the 30 individual response × predictor contrasts, noting explicitly that these are correlated and provide descriptive rather than independent support.

---

## 3. Results

### 3.1 Contrasting assembly signals among classes
Bird, mammal and reptile assemblages were functionally **clustered** relative to null expectation (mean SES of functional dispersion −1.68, −1.26 and −1.88), most strongly in the cold, arid west and north (Fig. 3a–c, e). Amphibian assemblages were, uniquely, **over-dispersed** (+1.70), most strongly in the eastern monsoon region (Fig. 3d). Values at 100 km differed from those at 50 km by ≤ 0.15 in every class, so the assembly signal is scale-invariant over this grain range (Fig. 6a).

### 3.2 Ectotherm trait composition is more tightly coupled to the environment
Across 33 families, mean absolute standardised environmental effects were **1.72× larger in ectotherms** than endotherms (0.133 vs 0.078; Wilcoxon *P* = 0.010; Fig. 4a). Restricting to the thermal-energy axis gave a 1.68-fold difference that was significant even under the conservative nested mixed model (0.252 vs 0.150; *P* = 0.023). The contrast held for every environmental axis (Fig. 4b), and in **30 of 30** response × predictor combinations the ectotherm estimate exceeded the endotherm estimate.

The single strongest contrast, and the only one to survive Benjamini–Hochberg correction across the 30 individual tests, was the effect of thermal energy on community-weighted **body mass**: 0.348 in ectotherm families versus 0.148 in endotherm families (contrast +0.201, SE 0.051, *t* = 3.91, *P*_BH = 0.032) — a 2.4-fold difference (Fig. 4c).

### 3.3 Human pressure filters ectotherms more strongly, not less
Coupling to the human-pressure axis was **2.05× stronger in ectotherms** (0.062 vs 0.030; Wilcoxon *P* = 0.0016), the opposite of the expectation that large-bodied, slow-lived endotherms would be the more exposed. Consistent with this, the only significant skewness response to human pressure was in reptiles, whose community-weighted body-mass distribution became increasingly right-skewed with human footprint (β = 0.046, SE = 0.018, *P*_BH = 0.030), indicating truncation of large-bodied species. Among endotherms, mammals lost large-bodied species (declining CWM body mass) and birds showed a contraction of body-mass variance — functional homogenisation — but neither translated into stronger overall coupling than in ectotherms (Fig. 5).

### 3.4 Functional structure is congruent only between the two endothermic classes
Spatial congruence in SES functional dispersion was substantial for the bird–mammal pair (Spearman *r* = 0.63 at 50 km, 0.59 at 100 km) but negligible for all five remaining pairs (|*r*| ≤ 0.18), including the reptile–amphibian pair, which was slightly negative (*r* = −0.07; Fig. 6b). Congruence is therefore not a property of shared thermoregulatory mode as such: the two endothermic classes converge on a common, weakly filtered spatial response, whereas the two ectothermic classes are each governed by different limiting factors — thermal energy and aridity in reptiles, moisture in amphibians. Practically, avian functional diversity is an informative surrogate only for mammals; it carries almost no information about reptile or amphibian functional structure (*r* = 0.02 and 0.16).

### 3.5 Naive CWM regression inflates trait–environment inference by 1.25–4.0×
Because every site shares the same trait vector, regressing community-weighted means on environment inflates Type I error (Peres-Neto et al. 2017; ter Braak et al. 2018). We therefore repeated all trait–environment tests with fourth-corner analysis under the max test, which permutes sites and species simultaneously and takes the larger of the two *P* values. Across every trait × environment combination, naive regression declared 1.27× (birds) to 2.72× (amphibians at 100 km) more associations significant than the max test, and 21–54% of individual conclusions were reversed (Table S8). Reversals were extreme: the association between reptile community-weighted verticality and annual precipitation had a naive *P* = 4.6 × 10⁻¹⁷⁵ but a max-test *P* = 0.054; amphibian body mass versus mean annual temperature fell from *P* = 3.2 × 10⁻⁴² to *P* = 0.375; and avian body mass versus human footprint fell from *P* = 1.5 × 10⁻¹⁸ to *P* = 0.154. Inflation was largest in the two ectothermic classes, whose smaller species pools and stronger trait autocorrelation make the naive test most misleading.

Note that the *number* of max-test-significant associations was nonetheless higher in endotherms (birds 74%, mammals 60%) than ectotherms (reptiles 38%, amphibians 32%). This reflects statistical power, which scales with species-pool size (1,056 birds versus 425 amphibians), and is not in conflict with §3.2: our headline conclusion rests on effect *magnitude* (mean absolute standardised effect per family), which is independent of pool size, not on counts of significant tests.

### 3.6 Robustness: leave-one-class-out and trait imputation
Because thermoregulatory mode is nearly collinear with class, we repeated the family-level contrast four times, each time dropping one class (Extended Data Fig. 9a; Table S11). The ectotherm-to-endotherm coupling ratio exceeded 1 in **all 15 scenario × metric combinations** (range 1.39–2.69). Individual significance depended on retaining the full family set: dropping amphibians left only seven ectotherm families and the overall coupling contrast became non-significant (ratio 1.39, *P* = 0.185), whereas dropping reptiles strengthened it (ratio 2.18, *P* = 0.005). The human-pressure asymmetry was the most robust single result, significant in **all five** scenarios.

We also tested every assemblage metric against the subset of species whose core traits were fully observed rather than imputed upstream by TetrapodTraits (Extended Data Fig. 9b, c; Table S12). For birds, mammals and reptiles, 81–84% of species were retained and mean SES of functional dispersion changed by at most 0.21 (birds −1.67 → −1.73; mammals −1.26 → −1.06; reptiles −1.92 → −1.77), so the functional clustering of these three classes is not an artefact of imputation. Amphibians were different: only 69 of 438 species (16%) had fully observed core traits — 75% of body masses and 71% of nocturnality values are imputed upstream — and mean SES reversed from +1.85 to −0.35. **We therefore do not report an assembly signal for Chinese amphibians**; their functional structure cannot be determined reliably from currently available trait data. Because the observed-only subset is itself biased towards well-studied, larger species, this test shows only that the original result may be imputation-driven; it does not establish clustering as the true pattern.

### 3.7 Trait moments are not redundant
Means, variances and skewness of the same trait responded to different environmental axes and with different frequencies of significance (Table S: 20–56% of tests significant depending on moment and class), supporting the treatment of moments as complementary rather than interchangeable descriptors. In the two ectothermic classes, skewness yielded a higher proportion of significant environmental effects than variance, indicating that tail truncation is a more sensitive signature of filtering than overall narrowing.

---

## 4. Discussion

*(to be expanded)*

**Endothermy buffers more than temperature.** The expectation that ectotherms track thermal gradients more closely was supported, but the effect generalised: ectotherm assemblage traits were more tightly coupled to *every* environmental axis we examined, including contemporary human pressure. This suggests that the physiological independence conferred by endothermy — homeothermy, higher mobility, greater behavioural and dietary flexibility — buffers assemblage composition against environmental variation broadly, rather than against thermal variation specifically.

**A reversal of the expected anthropogenic asymmetry.** Body size and life-history pace are the classic predictors of vertebrate sensitivity to human pressure, which led us to predict stronger anthropogenic filtering in endotherms. The data reverse this. A plausible mechanism is that land-use change acts primarily by altering **microhabitat structure and microclimate**, on which ectotherms depend directly for thermoregulation and water balance, whereas endotherms can traverse and tolerate modified matrices. This should be tested with microclimate-explicit data.

**Amphibian over-dispersion.** Amphibians were the only class whose assemblages were more functionally dispersed than random. Two non-exclusive explanations require separation in future work: (i) limiting similarity or competitive exclusion among ecologically similar congeners, and (ii) the strong allopatric structure of the Chinese amphibian fauna, in which geographically replacing sister species inflate the national pool with ecologically redundant forms that never co-occur. Phylogenetic beta-diversity and congener co-occurrence null models would discriminate these.

### Limitations
Thermoregulatory mode is nearly collinear with class, imposing a ceiling on the power of any within-China test; our family-level design mitigates but cannot eliminate this, and the strength of the conclusion rests on directional consistency as much as on individual *P* values. Expert range maps overestimate fine-scale occupancy and are insensitive to land-use change, so the anthropogenic effects reported here are conservative; repeating the analysis with occurrence records and explicit sampling-completeness correction is a priority. Longevity and amphibian diet were too incomplete to include, so the life-history "pace" axis is not represented in the core trait set. Most seriously, 84% of amphibian core trait values are imputed upstream, which is why we withhold any assembly conclusion for that class; resolving it requires primary morphological and ecological measurements for Chinese amphibians rather than better statistics.

---

## Data and code availability
All analysis code is in `01_code/` (`00`–`07`), fully scripted from raw data to figures. Derived data, result tables (Tables 2, 5, 6, 6b, 6c, 7; S1–S8) and figures are in `03_data_derived/`, `04_results/` and `05_figures/`. Primary data are third-party: IUCN Red List spatial data (iucnredlist.org; redistribution restricted by the IUCN Terms of Use), TetrapodTraits v2.0.1 (Zenodo 10.5281/zenodo.18926700), AVONET, COMBINE, AmphiBIO, ReptTraits, Etard et al. (2020), GlobTherm, the Amniote Life History Database, AnAge, WorldClim 2.1, MODIS NPP, ESA WorldCover, EarthEnv and the Human Footprint dataset. Phylogenies: Jetz et al. (2012), Upham et al. (2019), Tonini et al. (2016), Jetz & Pyron (2018).

All maps are based on the standard map GS(2023)2767 issued by the Ministry of Natural Resources of China, with no modification to national boundaries.
