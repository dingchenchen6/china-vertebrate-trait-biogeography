# ============================================================
# 目的 / Purpose:
#   把各分析脚本各自命名的图件，按最终论文的编号装配为正文图与扩展数据图。
#   Assemble the figures produced by the individual analysis scripts into the
#   final main-text and Extended Data numbering used by the manuscript.
#
# 为什么单独做这一步 / Why a separate step:
#   分析脚本以「分析内容」命名图件（如 Fig4_hypothesis_drivers），
#   便于重跑与追溯；论文编号则随行文结构变动。把两者解耦，改编号时
#   不必去改十几个分析脚本，也不会因改名而破坏脚本间的引用。
#   Analysis scripts name figures by content so they stay traceable; the
#   manuscript numbering changes with the narrative. Decoupling the two means
#   renumbering never requires touching an analysis script.
#
# 输入 / Input:  05_figures/*.pdf|png
# 输出 / Output: 05_figures/main/Fig1..Fig8, 05_figures/extended/ED1..
#                04_results/tables/table_00_figure_manifest.csv
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))

FIG <- PATH$figures
DIR_MAIN <- file.path(FIG, "main"); DIR_ED <- file.path(FIG, "extended")
dir.create(DIR_MAIN, showWarnings = FALSE); dir.create(DIR_ED, showWarnings = FALSE)

# ---------------------------------------------------------------
# 论文的三个关键科学问题与图件归属
# The three key questions and which figure answers which
# ---------------------------------------------------------------
Q1 <- "Q1 Dimensions: is assemblage trait structure one axis or several?"
Q2 <- "Q2 Drivers: which hypothesis governs it, and does thermoregulation decide?"
Q3 <- "Q3 Design: what does trait structure require of conservation planning?"

MAIN <- rbind(
  data.frame(num = "Fig1", src = "Fig1_framework_and_data", question = "Framework",
             title = "Study design, trait coverage and the assemblage framework"),
  data.frame(num = "Fig2", src = "Fig2_functional_axes", question = Q1,
             title = "Three a priori functional axes and their independence"),
  data.frame(num = "Fig3", src = "Fig2_CWM_maps", question = Q1,
             title = "Geography of assemblage trait means across China"),
  data.frame(num = "Fig4", src = "Fig4_hypothesis_drivers", question = Q2,
             title = "Which biodiversity hypothesis drives assemblage trait structure"),
  data.frame(num = "Fig5", src = "Fig5_fourth_corner", question = Q2,
             title = "Fourth-corner max test of trait-environment links"),
  data.frame(num = "Fig6", src = "Fig4_thermal_mode_asymmetry", question = Q2,
             title = "Endotherms and ectotherms differ in trait-environment coupling"),
  data.frame(num = "Fig7", src = "Fig9_functional_regionalisation", question = Q3,
             title = "No single class is a reliable surrogate: functional regions and redundancy differ"),
  data.frame(num = "Fig8", src = "Fig6_prioritisation", question = Q3,
             title = "Conservation design: complementarity, under-protected traits and 30x30"))

ED <- rbind(
  data.frame(num = "ED1", src = "ED1_trait_coverage_imputation",
             title = "Trait coverage before and after imputation"),
  data.frame(num = "ED2", src = "ED2_trait_space_pcoa",
             title = "Trait space and its principal coordinates"),
  data.frame(num = "ED3", src = "ED3_environmental_axes",
             title = "Environmental predictors and their axes"),
  data.frame(num = "ED4", src = "ED4_spatial_model_diagnostics",
             title = "Spatial model diagnostics"),
  data.frame(num = "ED5", src = "ED5_effect_heatmap_all",
             title = "Full trait-by-predictor effect matrix"),
  data.frame(num = "ED6", src = "ED6_fourthcorner_typeI",
             title = "Type I error inflation of naive CWM regression"),
  data.frame(num = "ED7", src = "ED7_grain_robustness_100km",
             title = "Robustness to grain: 100 km results"),
  data.frame(num = "ED8", src = "ED8_null_model_behaviour",
             title = "Null-model behaviour under occupancy weighting"),
  data.frame(num = "ED9", src = "ED9_robustness_checks",
             title = "Additional robustness checks"),
  data.frame(num = "ED10", src = "Fig3_SES_functional_diversity",
             title = "Standardised effect sizes of functional diversity"),
  data.frame(num = "ED11", src = "Fig5_anthropogenic_filter",
             title = "The anthropogenic filter on assemblage traits"),
  data.frame(num = "ED12", src = "Fig6_scale_and_congruence",
             title = "Grain sensitivity and cross-class congruence"),
  data.frame(num = "ED13", src = "Fig7_trait_variance_geography",
             title = "Geography of assemblage trait variance"),
  data.frame(num = "ED14", src = "Fig8_dominant_filtered_trait",
             title = "Which trait is most strongly filtered, mapped"),
  data.frame(num = "ED15", src = "Fig10_geographically_weighted_coupling",
             title = "Geographically weighted trait-environment coupling"),
  data.frame(num = "ED16", src = "Fig11_functional_beta",
             title = "Functional beta diversity and its turnover component"),
  data.frame(num = "ED17", src = "Fig12_occurrence_vs_rangemap",
             title = "Occurrence records versus expert range maps"),
  data.frame(num = "ED18", src = "Fig13_occurrence_four_classes",
             title = "Occurrence-based assemblages for all four classes"),
  data.frame(num = "ED19", src = "Fig14_conservation_gaps",
             title = "Conservation gaps and reserve representation"),
  data.frame(num = "ED20", src = "Fig15_sexual_dimorphism",
             title = "Sexual size dimorphism across assemblages"),
  data.frame(num = "ED21", src = "Fig2_trait_axes",
             title = "Data-driven trait axes retained as a check on the a priori axes"),
  data.frame(num = "ED22", src = "Fig6b_negative_control",
             title = "Three ways of measuring habitat heterogeneity, and the dose-response with thermal content"),
  data.frame(num = "ED23", src = "Fig6c_family_threshold",
             title = "Which families enter the clade-level test, and how the threshold changes it"),
  data.frame(num = "ED24", src = "Fig9_fd_and_nmds",
             title = "Functional-diversity indices and unconstrained NMDS ordination"),
  data.frame(num = "ED25", src = "Fig10_landuse_tolerance",
             title = "Land-use tolerance: which traits tolerate modification, and are they the protected ones"))

copy_set <- function(tab, dest) {
  ok <- logical(nrow(tab))
  for (i in seq_len(nrow(tab))) {
    hit <- FALSE
    for (ext in c("pdf", "png")) {
      src <- file.path(FIG, paste0(tab$src[i], ".", ext))
      if (file.exists(src)) {
        file.copy(src, file.path(dest, paste0(tab$num[i], ".", ext)),
                  overwrite = TRUE)
        hit <- TRUE
      }
    }
    ok[i] <- hit
    if (!hit) log_msg("  [缺失 missing] ", tab$src[i])
  }
  ok
}

MAIN$found <- copy_set(MAIN, DIR_MAIN)
ED$found   <- copy_set(ED, DIR_ED)

man <- rbind(cbind(panel = "Main", MAIN[, c("num", "src", "title", "found")],
                   question = MAIN$question),
             cbind(panel = "Extended Data", ED[, c("num", "src", "title", "found")],
                   question = NA_character_))
write_table(man, "table_00_figure_manifest")
log_msg("正文图 / main figures: ", sum(MAIN$found), "/", nrow(MAIN),
        " | 扩展数据图 / extended: ", sum(ED$found), "/", nrow(ED))
print(man[, c("panel", "num", "src", "found")], row.names = FALSE)
log_msg("=== 95 完成 / done ===")
