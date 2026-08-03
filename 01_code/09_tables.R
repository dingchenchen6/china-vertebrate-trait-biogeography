# ============================================================
# 本文件 / This file:
#   生成可直接排版的正文表与附表，并汇编为单一 Excel 工作簿
#   Publication-ready main and supplementary tables, plus one Excel workbook
#
# 设计意图 / Rationale:
#   分析脚本输出的是"机器可读"的长表；投稿需要的是"人可读"的窄表——
#   列名规范、单位明确、有效数字统一、显著性标注一致、每张表自带标题
#   与脚注。本脚本完成这一转换，不重新计算任何统计量。
#   The analysis scripts emit machine-readable long tables; submission needs
#   human-readable ones with clean headers, explicit units, consistent
#   rounding and significance marks. This script only reformats — it never
#   recomputes a statistic.
#
# 输出 / Output:
#   04_results/tables_formatted/Table1 .. Table3, TableS1 .. TableS9 (.csv)
#   04_results/vertebrate_trait_analysis_tables.xlsx  （全部表 + 说明页）
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
suppressPackageStartupMessages({ library(tidyr); library(openxlsx) })

log_msg("=== 09 规范化表格 / Formatted tables ===")

TBL  <- PATH$tables
OUT  <- file.path(PATH$results, "tables_formatted")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
rd <- function(f) data.table::fread(file.path(TBL, f))

TAXA_EN <- c(Aves = "Birds", Mammalia = "Mammals",
             Reptilia = "Reptiles", Amphibia = "Amphibians")
TRAIT_EN <- c(body_mass = "Body mass", body_length = "Body length",
              nocturnality = "Nocturnality", verticality = "Verticality",
              habitat_breadth = "Habitat breadth", range_size = "Range size",
              diet_breadth = "Diet breadth", diet_vert = "Vertebrate diet",
              diet_plant = "Plant diet", litter_size = "Litter/clutch size",
              max_longevity = "Maximum longevity")
PRED_EN <- c(ax_thermal = "Thermal energy", ax_water = "Water availability",
             ax_productivity = "Productivity", ax_structure = "Habitat structure",
             ax_human = "Human pressure")

#' 统一的 P 值格式：小于 0.001 用科学计数，其余三位小数
#' Consistent P formatting: scientific below 0.001, else three decimals
fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
    ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
           formatC(p, format = "f", digits = 3)))
}
#' 显著性星号 / significance marks
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*",
      ifelse(p < 0.1, "†", "")))))
}
#' "估计值 ± 标准误" 的紧凑写法 / compact estimate with SE
est_se <- function(e, s, d = 3) sprintf(paste0("%.", d, "f ± %.", d, "f"), e, s)

sheets <- list(); caps <- list()
put <- function(name, df, caption, note = "") {
  data.table::fwrite(df, file.path(OUT, paste0(name, ".csv")))
  sheets[[name]] <<- df
  caps[[name]]   <<- list(caption = caption, note = note)
  log_msg("  -> ", name, " (", nrow(df), " rows)")
}

# ---------------------------------------------------------------
# Table 1 | 数据概览 / Data summary
# ---------------------------------------------------------------
{
  s2 <- rd("table_S2b_iucn_assemblage_summary.csv")
  s2a <- rd("table_S2_assemblage_summary.csv")
  asm <- bind_rows(s2, s2a |> filter(class == "Aves")) |> filter(grain == "50km")
  # 每个纲实际进入分析的性状数（两栖类的食性与窝仔数因缺失过高被排除）
  # Traits actually analysed per class; amphibian diet and litter size were
  # too incomplete and are excluded, so the count is class-specific
  ntr <- rd("table_S6_imputation_diagnostics.csv") |>
    filter(used_in_main) |>
    group_by(class) |> summarise(n_traits = n_distinct(trait), .groups = "drop")
  m5 <- rd("table_S5_community_metrics_summary.csv") |> filter(grain == "50km")

  t1 <- asm |>
    left_join(ntr, by = "class") |>
    left_join(m5[, c("class", "n_cells", "mean_SES_FDis")], by = "class") |>
    left_join(TAXA[, c("class", "thermal")], by = "class") |>
    transmute(
      Class = TAXA_EN[class],
      `Thermoregulatory mode` = thermal,
      `Species (n)` = n_species,
      `Traits analysed (n)` = n_traits,
      `Occupied cells (n)` = n_cells_occupied,
      `Cells analysed (n)` = n_cells,
      `Mean richness per cell` = sprintf("%.1f", mean_richness),
      `Maximum richness` = max_richness,
      `Mean SES of functional dispersion` = sprintf("%+.2f", mean_SES_FDis)) |>
    arrange(match(Class, TAXA_EN))
  put("Table1_data_summary", t1,
      "Table 1 | Data summary for the four terrestrial vertebrate classes analysed on the 50-km equal-area grid across China.",
      "Species counts are those retained after IUCN filtering (presence 1-2, origin 1/2/6, seasonal 1-2, terrestrial only); bird assemblages derive from the equivalent BirdLife/IUCN-based matrix of Huang et al. (2023). Traits analysed is class-specific: maximum longevity was 19-55% incomplete in every class and was excluded throughout, and amphibian diet and litter size were additionally too incomplete to retain. Occupied cells is the number of cells holding at least one species; cells analysed is the subset with at least two species, for which assemblage metrics are defined. SES = standardised effect size against a richness-constrained null model with 999 randomisations; negative values indicate functional clustering, positive values over-dispersion.")
}

# ---------------------------------------------------------------
# Table 2 | H1 核心检验 / Headline thermoregulatory-mode test
# ---------------------------------------------------------------
{
  h <- rd("table_6b_headline_thermal_test.csv")
  sc <- rd("table_6c_sign_consistency.csv")
  METRIC_EN <- c(coupling = "All environmental axes",
                 coupling_thermal = "Thermal-energy axis only",
                 coupling_human = "Human-pressure axis only")
  t2 <- h |> transmute(
    `Coupling measure` = METRIC_EN[metric],
    `Endotherm families (n)` = n_endo,
    `Ectotherm families (n)` = n_ecto,
    `Endotherm mean` = sprintf("%.3f", endo_mean),
    `Ectotherm mean` = sprintf("%.3f", ecto_mean),
    `Ratio (ecto / endo)` = sprintf("%.2f", ratio_ecto_endo),
    `Contrast ± SE` = est_se(contrast, se),
    `P (mixed model)` = paste0(fmt_p(p_mixed), stars(p_mixed)),
    `P (Wilcoxon)` = paste0(fmt_p(p_wilcoxon), stars(p_wilcoxon)))
  put("Table2_headline_thermal_test", t2,
      "Table 2 | Trait-environment coupling strength in endothermic versus ectothermic families of Chinese terrestrial vertebrates.",
      sprintf("Coupling strength is the mean absolute standardised effect of the environmental axes on assemblage trait metrics, estimated per family from spatial simultaneous autoregressive error models on the 50-km grid (n = 33 families; 21 endothermic, 12 ectothermic). The mixed model is y ~ thermoregulatory mode + (1 | class/order); because thermoregulatory mode is nearly collinear with class, this test is deliberately conservative and the Wilcoxon rank-sum test is reported alongside. Directional consistency across the 30 individual response x predictor contrasts: %d of %d favoured ectotherms (binomial P = %s); those contrasts are mutually correlated and provide descriptive support only. Significance marks: *** P < 0.001, ** P < 0.01, * P < 0.05, dagger P < 0.1.",
              sc$n_ectotherm_stronger[1], sc$n_contrasts[1], format(sc$binomial_p[1])))
}

# ---------------------------------------------------------------
# Table 3 | 各科耦合强度 / Family-level coupling strength
# ---------------------------------------------------------------
{
  fam <- rd("table_7_family_coupling_strength.csv")
  t3 <- fam |> transmute(
    Class = TAXA_EN[class], Order = order, Family = clade,
    `Thermoregulatory mode` = thermal,
    `Species (n)` = n_species,
    `Coupling: all axes` = sprintf("%.3f", coupling),
    `Coupling: thermal energy` = sprintf("%.3f", coupling_thermal),
    `Coupling: human pressure` = sprintf("%.3f", coupling_human)) |>
    arrange(`Thermoregulatory mode`, Class, desc(`Species (n)`))
  put("Table3_family_coupling", t3,
      "Table 3 | Trait-environment coupling strength for each of the 33 families used as replicates.",
      "Families were retained when they contained at least 20 Chinese species and occupied at least 100 grid cells. Coupling is the mean absolute standardised effect from family-specific spatial autoregressive models.")
}

# ---------------------------------------------------------------
# 附表 / Supplementary tables
# ---------------------------------------------------------------
put("TableS1_grid", rd("table_S1_grid_summary.csv") |>
      transmute(Grain = grain, `Cells (n)` = n_cells,
                `Cell area (km2)` = cell_area_km2,
                `Land covered (km2)` = format(total_land_km2, big.mark = ","),
                `Mean land fraction` = mean_land_frac,
                `Provinces represented` = n_provinces),
    "Table S1 | Properties of the two equal-area analysis grids.",
    "Grids are in the China Albers equal-area conic projection matching the officially approved basemap GS(2023)2767. Cells with less than 50% land were excluded.")

put("TableS2_assemblages", bind_rows(
      rd("table_S2b_iucn_assemblage_summary.csv"),
      rd("table_S2_assemblage_summary.csv") |> filter(class == "Aves")) |>
      transmute(Class = TAXA_EN[class], Grain = grain, `Species (n)` = n_species,
                `Occupied cells (n)` = n_cells_occupied,
                `Total cells (n)` = n_cells_total,
                `Mean richness` = mean_richness, `Maximum richness` = max_richness,
                `Species-cell records (n)` = n_records) |>
      arrange(Class, desc(Grain)),
    "Table S2 | Assemblage composition summary by class and grain.")

put("TableS3_trait_coverage", rd("table_S3_trait_coverage_before_imputation.csv") |>
      transmute(Class = TAXA_EN[class], Trait = TRAIT_EN[trait],
                `Species (n)` = n, `With data (n)` = n_obs,
                `Coverage before imputation (%)` = coverage,
                `Imputed upstream by TetrapodTraits (%)` = upstream_imputed_pct,
                `Trait set` = trait_set) |>
      arrange(Class, Trait),
    "Table S3 | Trait coverage before imputation, and assignment of traits to the core, extended or excluded set.",
    "The core set comprises the six traits complete in all four classes and underpins every cross-class contrast. Upstream imputation refers to values already imputed by TetrapodTraits v2.0.1 and flagged there.")

put("TableS4_taxonomy", rd("table_S3c_taxonomy_harmonisation.csv") |>
      # fread 把缺失的匹配路径读成空字符串而非 NA，两者都要归为 "unmatched"
      # fread reads the missing route as "" rather than NA; catch both
      mutate(route = ifelse(is.na(match_route) | !nzchar(match_route),
                            "unmatched", match_route)) |>
      count(class, route, name = "n_species") |>
      mutate(Class = TAXA_EN[class]) |>
      select(Class, route, n_species) |>
      pivot_wider(names_from = route, values_from = n_species, values_fill = 0) |>
      mutate(Total = rowSums(across(where(is.numeric)))),
    "Table S4 | Outcome of the three-pass taxonomic harmonisation against the TetrapodTraits backbone.",
    "Pass 1 matched scientific names directly; pass 2 used the backbone's IUCN_Binomial field; pass 3 resolved the remainder through the GBIF backbone taxonomy.")

put("TableS5_imputation", rd("table_S6_imputation_diagnostics.csv") |>
      transmute(Class = TAXA_EN[class], Trait = TRAIT_EN[trait],
                `Species (n)` = n_species,
                `Missing before imputation (%)` = missing_pct_before,
                `Used in main analysis` = ifelse(used_in_main, "yes", "no"),
                `missForest OOB NRMSE` = OOB_NRMSE) |>
      arrange(Class, Trait),
    "Table S5 | Trait imputation diagnostics.",
    "Imputation used missForest with the first ten phylogenetic eigenvectors and collapsed order/family factors as auxiliary predictors. Traits missing in more than 40% of species within a class were excluded rather than imputed.")

put("TableS6_environment", rd("table_S4_env_summary.csv") |>
      transmute(Grain = grain, Variable = variable, `Cells with data (n)` = n_valid,
                Mean = mean, SD = sd, Minimum = min, Maximum = max),
    "Table S6 | Summary statistics for all 47 environmental variables extracted per grid cell.")

put("TableS7_predictor_axes", bind_rows(
      rd("table_S7_predictor_axis_loadings_50km.csv") |> mutate(Grain = "50km"),
      rd("table_S7_predictor_axis_loadings_100km.csv") |> mutate(Grain = "100km")) |>
      transmute(Grain, `Ecological axis` = PRED_EN[group], Variable = variable,
                `Loading on PC1` = loading,
                `Variance explained by PC1 (%)` = var_explained),
    "Table S7 | Loadings of the individual environmental variables on the five ecological predictor axes.",
    "Each axis is the first principal component of its z-standardised member variables, sign-aligned so that the first listed variable loads positively.")

# 「结论翻转」= 朴素检验显著而 max test 不显著；直接从原始结果计算
# "Reversed" = significant under naive regression but not under the max test
flip <- rd("table_4_fourthcorner_maxtest.csv") |>
  group_by(class, grain) |>
  summarise(pct_reversed = round(100 * mean(p_naive_adj < 0.05 & p_maxtest >= 0.05,
                                            na.rm = TRUE), 1), .groups = "drop")
put("TableS8_typeI", rd("table_S8_typeI_inflation_diagnostic.csv") |>
      left_join(flip, by = c("class", "grain")) |>
      transmute(Class = TAXA_EN[class], Grain = grain, `Tests (n)` = n_tests,
                `Significant, naive CWM regression` = sig_naive,
                `Significant, fourth-corner max test` = sig_maxtest,
                `Naive (%)` = pct_naive, `Max test (%)` = pct_maxtest,
                `Inflation ratio` = inflation_ratio,
                `Conclusions reversed (%)` = pct_reversed) |>
      arrange(Class, desc(Grain)),
    "Table S8 | Type I error inflation of naive community-weighted-mean regression relative to the fourth-corner max test.",
    "The max test permutes sites and species simultaneously and takes the larger of the two P values, controlling the error rate that naive CWM regression inflates because all sites share one trait vector. Conclusions reversed is the percentage of trait x environment combinations that the naive test declares significant and the max test does not. Trait sets match the main analysis exactly, so amphibians contribute 60 rather than 100 combinations.")

f9 <- file.path(TBL, "table_S9_trait_space_loadings.csv")
if (file.exists(f9)) {
  put("TableS9_trait_space", rd("table_S9_trait_space_loadings.csv") |>
        transmute(Class = TAXA_EN[class], Trait = TRAIT_EN[trait],
                  `Correlation with PCo1` = round(r1, 3),
                  `Correlation with PCo2` = round(r2, 3)),
      "Table S9 | Correlations of each core trait with the first two principal-coordinate axes of the class-specific trait space.")
}

put("TableS10_all_contrasts", rd("table_6_thermal_mode_contrast.csv") |>
      filter(quantity == "strength") |>
      transmute(Response = response, `Environmental axis` = PRED_EN[predictor],
                `Endotherm mean` = endo_mean, `Ectotherm mean` = ecto_mean,
                `Contrast ± SE` = est_se(contrast, se),
                t = t, P = fmt_p(p), `P (BH-adjusted)` = fmt_p(p_adj)) |>
      arrange(P),
    "Table S10 | All 30 individual response x environmental-axis contrasts between endothermic and ectothermic families.",
    "Ectotherm coupling exceeded endotherm coupling in all 30 contrasts. These contrasts are mutually correlated, so the headline inference is the single family-level test in Table 2 rather than these individual comparisons.")

f11 <- file.path(TBL, "table_8_leave_one_class_out.csv")
if (file.exists(f11)) {
  put("TableS11_leave_one_class_out", rd("table_8_leave_one_class_out.csv") |>
        transmute(Scenario = scenario,
                  `Coupling measure` = c(coupling = "All environmental axes",
                                         coupling_thermal = "Thermal-energy axis",
                                         coupling_human = "Human-pressure axis")[metric],
                  `Endotherm families (n)` = n_endo,
                  `Ectotherm families (n)` = n_ecto,
                  `Endotherm mean` = endo_mean, `Ectotherm mean` = ecto_mean,
                  `Ratio (ecto / endo)` = ratio,
                  `P (Wilcoxon)` = paste0(fmt_p(p_wilcoxon), stars(p_wilcoxon)),
                  `P (mixed model)` = paste0(fmt_p(p_mixed), stars(p_mixed))),
      "Table S11 | Leave-one-class-out robustness of the thermoregulatory-mode contrast.",
      "Each class was dropped in turn and the family-level contrast refitted. The ectotherm-to-endotherm ratio exceeded 1 in all 15 scenario x metric combinations (range 1.39-2.69). Because dropping a class leaves only 5-13 families per group, loss of significance reflects loss of power rather than reversal of the effect.")
}

f12 <- file.path(TBL, "table_9_imputation_sensitivity.csv")
if (file.exists(f12)) {
  put("TableS12_imputation_sensitivity", rd("table_9_imputation_sensitivity.csv") |>
        transmute(Class = TAXA_EN[class],
                  `Species, all (n)` = n_species_all,
                  `Species with fully observed core traits (n)` = n_species_observed_only,
                  `Retained (%)` = pct_species_retained,
                  `Mean SES, all species` = sprintf("%+.2f", mean_SES_all),
                  `Mean SES, observed traits only` = sprintf("%+.2f", mean_SES_observed_only),
                  `Cells, observed-only analysis (n)` = n_cells_observed_only),
      "Table S12 | Sensitivity of assemblage functional structure to upstream trait imputation.",
      "Species were retained only if none of their core trait values carried a TetrapodTraits imputation flag. Birds, mammals and reptiles retain 81-84% of species and their mean SES of functional dispersion shifts by at most 0.21. Amphibians retain 16% and the mean SES reverses from +1.85 to -0.35, so no assembly signal is reported for that class. The observed-only subset is itself biased towards well-studied species, so this test can only show that a result may be imputation-driven.")
}

f13 <- file.path(TBL, "table_17_amphibian_trait_rescue.csv")
if (file.exists(f13)) {
  put("TableS13_amphibian_trait_rescue", rd("table_17_amphibian_trait_rescue.csv") |>
        transmute(Trait = TRAIT_EN[trait], `Species with observed values (n)` = n_observed,
                  `Species total (n)` = n_total, `Coverage (%)` = pct,
                  `Retained in the well-measured set` = ifelse(retained, "yes", "no")),
      "Table S13 | Observed trait coverage for Chinese amphibians, taken from the primary databases rather than from imputed values.",
      "Values were read directly from AmphiBIO v1 and Etard et al. (2020); range size is computed from the IUCN polygons and is therefore a measurement. Body mass and diel activity are observed for only 8% and 11% of species, which reflects a genuine gap in the primary literature rather than an integration failure.")
}

f14 <- file.path(TBL, "table_18_amphibian_ses_comparison.csv")
if (file.exists(f14)) {
  put("TableS14_amphibian_ses_comparison", rd("table_18_amphibian_ses_comparison.csv") |>
        transmute(Scenario = scenario, `Traits (n)` = n_traits,
                  `Species (n)` = n_species, `Cells (n)` = n_cells,
                  `Mean SES of functional dispersion` = sprintf("%+.2f", mean_SES),
                  `Cells with positive SES (%)` = pct_positive),
      "Table S14 | Quasi-factorial test of what drives the apparent over-dispersion of Chinese amphibian assemblages.",
      "Comparing A with D holds the species set constant and changes the trait set; comparing D with C holds the trait set constant and changes imputed for observed values. Each contrast removes about half of the apparent signal, and the cleanest estimate (C) is indistinguishable from zero.")
}

f15 <- file.path(TBL, "table_19_functional_beta.csv")
if (file.exists(f15)) {
  put("TableS15_functional_beta", rd("table_19_functional_beta.csv") |>
        transmute(Class = TAXA_EN[class], `Thermoregulatory mode` = thermal,
                  `Adjacent cell pairs (n)` = n_pairs,
                  `Taxonomic beta` = tax_sor, `Functional beta` = fun_sor,
                  `Functional / taxonomic` = round(beta_ratio, 3),
                  `95% CI` = paste0(ratio_lo, "-", ratio_hi),
                  `Turnover share, species (%)` = round(pct_turnover_tax, 1),
                  `Turnover share, functional (%)` = round(pct_turnover_fun, 1)),
      "Table S15 | Functional and taxonomic beta diversity between adjacent 50-km cells.",
      "Sorensen dissimilarity decomposed following Baselga (2010). Functional groups are k-means clusters (K = 20) in the trait principal-coordinate space. Confidence intervals are 999 bootstrap resamples of cell pairs. The ordering of the ratio across classes is unchanged at K = 10 and K = 30 (Table S16).")
}

f16 <- file.path(TBL, "table_20_redundancy_contrast.csv")
if (file.exists(f16)) {
  put("TableS16_redundancy_K_sensitivity", rd("table_20_redundancy_contrast.csv") |>
        transmute(Class = TAXA_EN[class], `Functional groups (K)` = K,
                  `Functional / taxonomic beta` = beta_ratio),
      "Table S16 | Sensitivity of the functional-to-taxonomic beta ratio to the number of functional groups.")
}

f17 <- file.path(TBL, "table_21_occurrence_vs_rangemap.csv")
if (file.exists(f17)) {
  put("TableS17_occurrence_vs_rangemap", rd("table_21_occurrence_vs_rangemap.csv") |>
        transmute(`Cells compared (n)` = n_cells,
                  `Mean richness, occurrences` = mean_richness_occ,
                  `Mean richness, range maps` = mean_richness_map,
                  `Detection ratio` = detection_ratio,
                  `Spearman r, richness` = r_richness,
                  `Spearman r, body-mass CWM` = r_body_mass_cwm,
                  `Spearman r, habitat-breadth CWM` = r_habitat_breadth_cwm),
      "Table S17 | Agreement between bird assemblages derived from range maps and from 4.58 million GBIF/eBird occurrence records.",
      "Occurrence assemblages were retained only for cells with Chao sample coverage of at least 0.7 and at least 20 recorded species. Trait composition agrees far better between the two sources than species richness does.")
}

f18 <- file.path(TBL, "table_22_anthropogenic_effect_comparison.csv")
if (file.exists(f18)) {
  put("TableS18_anthropogenic_effect_comparison",
      rd("table_22_anthropogenic_effect_comparison.csv") |>
        filter(term == "ax_human") |>
        transmute(Response = response, `Assemblage source` = source,
                  `Estimate ± SE` = est_se(estimate, se, 4),
                  `Cells (n)` = n, `P (BH-adjusted)` = fmt_p(p_adj)),
      "Table S18 | Effect of the human-pressure axis on bird assemblage metrics under three assemblage definitions.",
      "All models are spatial simultaneous autoregressive error models with the same five environmental axes on the 50-km grid. The effect on community-weighted range size is 5.6 times larger in occurrence assemblages than in range-map assemblages, and the effect on habitat-breadth variance is absent from range maps entirely. Adding Chao coverage as a covariate strengthens rather than weakens these effects, which argues against a sampling-bias explanation.")
}

# ---------------------------------------------------------------
# Excel 汇编 / Excel workbook
# ---------------------------------------------------------------
wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#EFEFEF",
                   border = "bottom", halign = "left", wrapText = TRUE)
cap_style <- createStyle(textDecoration = "italic", wrapText = TRUE,
                         valign = "top", fontSize = 10)

addWorksheet(wb, "Contents")
contents <- data.frame(
  Sheet = names(sheets),
  Caption = vapply(caps, function(x) x$caption, character(1)),
  Note = vapply(caps, function(x) x$note, character(1)), row.names = NULL)
writeData(wb, "Contents", contents, headerStyle = hdr)
setColWidths(wb, "Contents", 1:3, widths = c(28, 70, 90))
addStyle(wb, "Contents", cap_style, rows = 2:(nrow(contents) + 1), cols = 2:3,
         gridExpand = TRUE)

for (nm in names(sheets)) {
  sn <- substr(nm, 1, 31)
  addWorksheet(wb, sn)
  writeData(wb, sn, caps[[nm]]$caption, startRow = 1)
  addStyle(wb, sn, cap_style, rows = 1, cols = 1)
  writeData(wb, sn, sheets[[nm]], startRow = 3, headerStyle = hdr)
  setColWidths(wb, sn, seq_len(ncol(sheets[[nm]])), widths = "auto")
  freezePane(wb, sn, firstActiveRow = 4)
}
xlsx <- file.path(PATH$results, "vertebrate_trait_analysis_tables.xlsx")
saveWorkbook(wb, xlsx, overwrite = TRUE)
log_msg("Excel 汇编 / workbook: ", basename(xlsx), " (", length(sheets), " 张表)")

# 表格说明清单 / caption index for the manuscript
cap_md <- c("# 表格说明 / Table captions", "",
  unlist(lapply(names(caps), function(nm) c(
    paste0("**", caps[[nm]]$caption, "**"), "",
    if (nzchar(caps[[nm]]$note)) c(caps[[nm]]$note, "") else NULL))))
writeLines(cap_md, file.path(PATH$ms, "table_captions.md"))
log_msg("=== 09 完成 / done ===")
