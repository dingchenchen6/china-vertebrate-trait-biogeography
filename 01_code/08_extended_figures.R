# ============================================================
# 本文件 / This file:
#   扩展数据图（Extended Data Figures）——方法诊断、数据质量与稳健性
#   Extended Data Figures: method diagnostics, data quality, robustness
#
# 设计意图 / Rationale:
#   顶刊正文图只承载核心论证；所有"审稿人会问"的证据——性状覆盖度、
#   插补误差、性状空间、预测变量共线性、空间自相关残差、零模型行为、
#   另一粒度的重复、第一类错误膨胀——放在扩展数据图中，逐条对应
#   研究设计文件 §7 的审稿人雷区清单。
#   Main figures carry the argument; every "a reviewer will ask" item goes
#   here, mapped one-to-one onto the reviewer-landmine checklist.
#
# 输出 / Output:
#   05_figures/ED1 .. ED9 (.pdf + .png)
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(ape); library(cluster)
})

log_msg("=== 08 扩展数据图 / Extended Data Figures ===")
BM  <- load_basemap(PATH$derived)
FIG <- PATH$figures
TBL <- PATH$tables

rd <- function(f) data.table::fread(file.path(TBL, f))
g50  <- st_read(file.path(PATH$derived, "grid_50km.gpkg"),  quiet = TRUE)
g100 <- st_read(file.path(PATH$derived, "grid_100km.gpkg"), quiet = TRUE)
read_metrics <- function(lab) {
  fs <- list.files(PATH$derived, pattern = paste0("^metrics_.*_", lab, "\\.rds$"),
                   full.names = TRUE)
  lapply(fs, readRDS) |> bind_rows() |>
    left_join(TAXA[, c("class", "thermal")], by = "class")
}
m50  <- read_metrics("50km"); m100 <- read_metrics("100km")

TRAIT_LAB <- c(body_mass = "Body mass", body_length = "Body length",
               nocturnality = "Nocturnality", verticality = "Vertical stratum use",
               habitat_breadth = "Habitat breadth", range_size = "Range size",
               diet_breadth = "Diet breadth", diet_vert = "Vertebrate diet",
               diet_plant = "Plant diet", litter_size = "Litter/clutch size",
               max_longevity = "Maximum longevity")
PRED_LAB <- c(ax_thermal = "Thermal energy", ax_water = "Water availability",
              ax_productivity = "Productivity", ax_structure = "Habitat structure",
              ax_human = "Human pressure")

# ===============================================================
# ED1 | 性状覆盖度与插补诊断
# ED1 | Trait coverage and imputation diagnostics
# ===============================================================
ed1 <- function() {
  cov <- rd("table_S3_trait_coverage_before_imputation.csv")
  ORD <- c("body_mass", "body_length", "nocturnality", "verticality",
           "habitat_breadth", "range_size", "diet_breadth", "diet_vert",
           "diet_plant", "litter_size", "max_longevity")
  cov$trait_f <- factor(TRAIT_LAB[cov$trait], levels = rev(TRAIT_LAB[ORD]))
  cov$class_f <- factor(TAXA_LAB[cov$class], levels = unname(TAXA_LAB))

  # (a) 插补前覆盖率。用单调递增的深浅表示"越深越完整"，并按背景亮度切换字色
  # (a) Coverage before imputation; darker = more complete, text colour flips
  #     with background luminance so every number stays legible
  pa <- ggplot(cov, aes(class_f, trait_f, fill = coverage)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f", coverage),
                  colour = coverage > 60), size = 1.9) +
    scale_fill_gradient(low = "#F3F7F5", high = "#1F5B4E", limits = c(0, 100),
                        name = "Coverage before\nimputation (%)") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                        guide = "none") +
    labs(x = NULL, y = NULL, tag = "a") +
    theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                        legend.position = "right",
                        legend.key.width = unit(6, "pt"))

  # (b) 性状集归属，与 (a) 行对齐 / trait-set membership, aligned to panel a
  st <- cov |> distinct(trait, trait_set) |>
    mutate(trait_f = factor(TRAIT_LAB[trait], levels = rev(TRAIT_LAB[ORD])),
           trait_set = factor(trait_set,
             levels = c("core (all 4 classes)",
                        "extended (birds/mammals/reptiles)",
                        "excluded (data-deficient)")))
  pb <- ggplot(st, aes(x = 1, y = trait_f, fill = trait_set)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("core (all 4 classes)" = "#2E6E8E",
                                 "extended (birds/mammals/reptiles)" = "#8FBF9F",
                                 "excluded (data-deficient)" = "grey78"),
                      name = "Trait set") +
    labs(x = NULL, y = NULL, tag = "b") +
    theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "right", legend.text = element_text(size = 5),
          legend.key.size = unit(6, "pt"))

  # (c) TetrapodTraits 上游已插补的比例（仅这 5 个性状带标记）
  # (c) Share imputed upstream; only these five traits carry a flag
  up <- cov |> filter(is.finite(upstream_imputed_pct))
  pc <- ggplot(up, aes(class_f, trait_f, fill = upstream_imputed_pct)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f", upstream_imputed_pct),
                  colour = upstream_imputed_pct > 45), size = 1.9) +
    scale_fill_gradient(low = "#FFF8EC", high = "#B5541F", limits = c(0, 100),
                        name = "Imputed upstream by\nTetrapodTraits (%)") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                        guide = "none") +
    labs(x = NULL, y = NULL, tag = "c") +
    theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                        axis.text.x = element_text(angle = 20, hjust = 1),
                        legend.position = "right",
                        legend.key.width = unit(6, "pt"))

  # (d) missForest 袋外误差。两栖类为 0，因为剔除高缺失性状后其余性状已完整，
  #     无需插补——必须标注，否则会被误读为异常值。
  # (d) OOB error; amphibians are exactly zero because, once the data-deficient
  #     traits were dropped, nothing remained to impute. Label it explicitly.
  dg <- rd("table_S6_imputation_diagnostics.csv") |>
    distinct(class, OOB_NRMSE) |> filter(is.finite(OOB_NRMSE))
  dg$class_f <- factor(TAXA_LAB[dg$class], levels = unname(TAXA_LAB))
  dg$lab <- ifelse(dg$OOB_NRMSE == 0, "0 (no imputation\nrequired)",
                   sprintf("%.4f", dg$OOB_NRMSE))
  pd <- ggplot(dg, aes(class_f, OOB_NRMSE, fill = class)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = lab), vjust = -0.35, size = 1.8, lineheight = 0.9) +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(x = NULL, y = "missForest OOB NRMSE", tag = "d") +
    theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))

  p <- (pa | pb) / (pc | pd) +
    plot_layout(widths = c(3.1, 1), heights = c(1.45, 1))
  save_fig(p, "ED1_trait_coverage_imputation", W2, 145, FIG)
}

# ===============================================================
# ED2 | 性状空间：各类群 PCoA 与性状载荷
# ED2 | Trait space: per-class PCoA and trait loadings
# ===============================================================
ed2 <- function() {
  traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))
  CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
            "habitat_breadth", "range_size")
  panels <- list(); loads <- list()
  for (cl in TAXA$class) {
    comm <- readRDS(file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl)))
    tr <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
    sp <- intersect(unique(comm$species), tr$species)
    tr <- tr[match(sp, tr$species), ]
    TM <- apply(as.matrix(sapply(tr[, CORE], as.numeric)), 2, zscore)
    pc <- ape::pcoa(as.dist(cluster::daisy(as.data.frame(TM), metric = "gower")))
    ax <- as.data.frame(pc$vectors[, 1:2]); names(ax) <- c("PCo1", "PCo2")
    ve <- round(100 * pc$values$Relative_eig[1:2], 1)
    # 性状与前两轴的相关 = 载荷 / trait-axis correlations as loadings
    ld <- data.frame(trait = CORE,
                     r1 = cor(TM, ax$PCo1), r2 = cor(TM, ax$PCo2))
    loads[[cl]] <- ld |> mutate(class = cl)
    sc <- 0.9 * max(abs(c(ax$PCo1, ax$PCo2)))
    panels[[length(panels) + 1]] <-
      ggplot(ax, aes(PCo1, PCo2)) +
      geom_point(size = 0.25, alpha = 0.25, colour = PAL$taxa[[cl]]) +
      geom_segment(data = ld, inherit.aes = FALSE,
                   aes(x = 0, y = 0, xend = r1 * sc, yend = r2 * sc),
                   arrow = arrow(length = unit(1.6, "pt")),
                   linewidth = 0.25, colour = "grey25") +
      ggrepel::geom_text_repel(data = ld, inherit.aes = FALSE,
                               aes(r1 * sc, r2 * sc, label = TRAIT_LAB[trait]),
                               size = 1.7, colour = "grey15", segment.size = 0.15,
                               min.segment.length = 0, max.overlaps = 20) +
      labs(title = sprintf("%s (n = %d)", TAXA_LAB[cl], nrow(TM)),
           x = sprintf("PCo1 (%.1f%%)", ve[1]),
           y = sprintf("PCo2 (%.1f%%)", ve[2])) +
      theme_pub()
  }
  write_table(bind_rows(loads), "table_S9_trait_space_loadings")
  p <- wrap_plots(panels, nrow = 2) + plot_annotation(tag_levels = "a")
  save_fig(p, "ED2_trait_space_pcoa", W2, 130, FIG)
}

# ===============================================================
# ED3 | 环境预测轴：载荷、空间格局与相关结构
# ED3 | Environmental axes: loadings, maps and correlation structure
# ===============================================================
ed3 <- function() {
  ld <- rd("table_S7_predictor_axis_loadings_50km.csv")
  ld$group_f <- factor(PRED_LAB[ld$group], levels = unname(PRED_LAB))
  pa <- ggplot(ld, aes(loading, reorder(variable, loading), fill = group_f)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    facet_wrap(~group_f, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = "Loading on first PCA axis", y = NULL, tag = "a") +
    theme_pub() + theme(strip.text = element_text(size = 5.5))

  d <- readRDS(file.path(PATH$derived, "model_data.rds"))[["50km"]] |>
    distinct(cell_id, .keep_all = TRUE)
  maps <- lapply(names(PRED_LAB), function(v) {
    dd <- d[, c("cell_id", v)]; names(dd)[2] <- "val"
    map_china(g50 |> inner_join(dd, by = "cell_id"), "val", BM, "div",
              title = PRED_LAB[[v]], legend = "Axis score (z)", inset = FALSE) +
      theme(legend.position = "bottom", legend.key.width = unit(18, "pt"),
            legend.key.height = unit(3, "pt"))
  })

  # 五轴之间的相关 / correlations among the five axes
  cm <- cor(d[, names(PRED_LAB)], use = "pairwise.complete.obs")
  cd <- as.data.frame(as.table(cm)); names(cd) <- c("x", "y", "r")
  pc <- ggplot(cd, aes(PRED_LAB[as.character(x)], PRED_LAB[as.character(y)], fill = r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.9) +
    scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-1, 1),
                     name = "Pearson r") +
    labs(x = NULL, y = NULL, tag = "g") +
    theme_pub() + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                        axis.line = element_blank(), axis.ticks = element_blank())

  right <- wrap_plots(maps, nrow = 2) + plot_annotation(tag_levels = "b")
  p <- (pa | wrap_plots(maps, nrow = 3)) / pc +
    plot_layout(widths = c(1, 2), heights = c(2.6, 1))
  save_fig(p, "ED3_environmental_axes", W2, 190, FIG)
}

# ===============================================================
# ED4 | 空间自回归诊断：Moran's I、lambda、AIC 改善
# ED4 | SAR diagnostics: Moran's I, lambda, AIC improvement
# ===============================================================
ed4 <- function() {
  co <- rd("table_2_sar_coefficients_all.csv") |>
    distinct(class, grain, response, .keep_all = TRUE)
  co$class_f <- factor(TAXA_LAB[co$class], levels = unname(TAXA_LAB))

  pa <- ggplot(co, aes(class_f, -log10(pmax(moran_ols_p, 1e-300)), fill = class)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25,
               colour = "grey45") +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    facet_wrap(~grain) +
    labs(x = NULL, y = expression(-log[10]~italic(P)*", Moran's I of OLS residuals"),
         tag = "a") +
    theme_pub()

  pb <- ggplot(co, aes(class_f, lambda, fill = class)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.6) +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    facet_wrap(~grain) +
    labs(x = NULL, y = expression("SAR spatial parameter "*lambda), tag = "b") +
    theme_pub()

  # ΔAIC 用双向对数刻度：SAR 几乎总是更优，但极少数模型 ΔAIC <= 0，
  # 直接取对数会产生 -Inf，故用 signed-log 变换保留全部数据点。
  # Signed-log keeps the rare non-positive ΔAIC values instead of dropping them.
  co$daic <- co$aic_ols - co$aic_sar
  slog <- function(x) sign(x) * log10(1 + abs(x))
  pc <- ggplot(co, aes(class_f, slog(daic), fill = class)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.6) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, colour = "grey45") +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    scale_y_continuous(breaks = slog(c(0, 10, 100, 1000, 10000)),
                       labels = c("0", "10", "100", "1000", "10000")) +
    facet_wrap(~grain) +
    labs(x = NULL, y = expression(Delta*"AIC (OLS - SAR), signed-log scale"),
         tag = "c",
         caption = sprintf("SAR outperformed OLS in %.1f%% of models",
                           100 * mean(co$daic > 0, na.rm = TRUE))) +
    theme_pub() + theme(plot.caption = element_text(size = 5, colour = "grey35"))

  p <- pa / pb / pc
  save_fig(p, "ED4_spatial_model_diagnostics", W2, 165, FIG)
}

# ===============================================================
# ED5 | 全部标准化效应量热图（矩 x 环境轴 x 类群）
# ED5 | Full standardised-effect heatmap (moment x axis x class)
# ===============================================================
ed5 <- function() {
  co <- rd("table_2_sar_coefficients_all.csv") |>
    filter(grain == "50km", term %in% names(PRED_LAB),
           moment %in% c("mean", "var", "skew", "kurt"))
  co$trait_f  <- factor(TRAIT_LAB[co$trait], levels = rev(unname(TRAIT_LAB)))
  co$term_f   <- factor(PRED_LAB[co$term], levels = unname(PRED_LAB))
  co$class_f  <- factor(TAXA_LAB[co$class], levels = unname(TAXA_LAB))
  co$moment_f <- factor(co$moment, levels = c("mean", "var", "skew", "kurt"),
                        labels = c("Mean (CWM)", "Variance (CWV)",
                                   "Skewness (CWS)", "Kurtosis (CWK)"))
  lim <- stats::quantile(abs(co$estimate), 0.98, na.rm = TRUE)

  p <- ggplot(co, aes(term_f, trait_f, fill = estimate)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_point(data = co |> filter(p_adj < 0.05), size = 0.35, colour = "grey10") +
    scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim),
                     oob = scales::squish, name = "Standardised effect") +
    facet_grid(moment_f ~ class_f) +
    labs(x = NULL, y = NULL,
         caption = "Dots mark effects significant at BH-adjusted P < 0.05") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 5),
          axis.text.y = element_text(size = 5),
          axis.line = element_blank(), axis.ticks = element_blank(),
          strip.text = element_text(size = 6),
          plot.caption = element_text(size = 5, colour = "grey35"),
          legend.position = "bottom", legend.key.width = unit(30, "pt"),
          legend.key.height = unit(3.5, "pt"))
  save_fig(p, "ED5_effect_heatmap_all", W2, 175, FIG)
}

# ===============================================================
# ED6 | 第四角 max test 与第一类错误膨胀
# ED6 | Fourth-corner max test and Type I error inflation
# ===============================================================
ed6 <- function() {
  f <- file.path(TBL, "table_4_fourthcorner_maxtest.csv")
  if (!file.exists(f)) { log_msg("  [skip] ED6"); return(invisible(NULL)) }
  fc <- rd("table_4_fourthcorner_maxtest.csv") |> filter(grain == "50km")
  ENV_LAB <- c(bio1_mean = "Mean annual temp.", bio6_mean = "Min. temp. coldest month",
               bio4_mean = "Temp. seasonality", bio12_mean = "Annual precip.",
               bio15_mean = "Precip. seasonality", npp = "NPP",
               elev_range = "Elevation range", het_shannon = "Habitat heterogeneity",
               hfp2020 = "Human footprint", lc_cropland = "Cropland cover")
  fc$env_f   <- factor(ENV_LAB[fc$env], levels = unname(ENV_LAB))
  fc$trait_f <- factor(TRAIT_LAB[fc$trait], levels = rev(unname(TRAIT_LAB)))
  fc$class_f <- factor(TAXA_LAB[fc$class], levels = unname(TAXA_LAB))
  lim <- stats::quantile(abs(fc$stat), 0.98, na.rm = TRUE)

  pa <- ggplot(fc, aes(env_f, trait_f, fill = stat)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_point(data = fc |> filter(p_maxtest < 0.05), size = 0.35, colour = "grey10") +
    scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim),
                     oob = scales::squish, name = "Fourth-corner statistic") +
    facet_wrap(~class_f, nrow = 1) +
    labs(x = NULL, y = NULL, tag = "a",
         caption = "Dots mark associations significant under the max test (BH-adjusted P < 0.05)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 4.6),
          axis.text.y = element_text(size = 5),
          axis.line = element_blank(), axis.ticks = element_blank(),
          plot.caption = element_text(size = 5, colour = "grey35"),
          legend.position = "bottom", legend.key.width = unit(28, "pt"),
          legend.key.height = unit(3.5, "pt"))

  # 朴素 p vs max test p / naive versus max-test P values
  fc$flip <- fc$p_naive_adj < 0.05 & fc$p_maxtest >= 0.05
  pb <- ggplot(fc, aes(-log10(pmax(p_naive_adj, 1e-300)),
                       -log10(pmax(p_maxtest, 1e-4)), colour = class)) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = "grey45") +
    geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = "grey45") +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    labs(x = expression("Naive CWM regression, "*-log[10]~italic(P)),
         y = expression("Fourth-corner max test, "*-log[10]~italic(P)), tag = "b") +
    theme_pub() + theme(legend.position = c(0.98, 0.98),
                        legend.justification = c(1, 1))

  inf <- rd("table_S8_typeI_inflation_diagnostic.csv")
  inf$class_f <- factor(TAXA_LAB[inf$class], levels = unname(TAXA_LAB))
  pc <- ggplot(inf, aes(class_f, inflation_ratio, fill = grain)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey45") +
    geom_text(aes(label = sprintf("%.2f", inflation_ratio)),
              position = position_dodge(width = 0.75), vjust = -0.4, size = 1.8) +
    scale_fill_manual(values = c(`50km` = "grey70", `100km` = "grey35"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Type I inflation\n(naive / max-test significant)", tag = "c") +
    theme_pub()

  p <- pa / (pb | pc) + plot_layout(heights = c(1.5, 1))
  save_fig(p, "ED6_fourthcorner_typeI", W2, 150, FIG)
}

# ===============================================================
# ED7 | 100 km 粒度的重复（尺度稳健性）
# ED7 | 100-km replicate maps (grain robustness)
# ===============================================================
ed7 <- function() {
  qs <- stats::quantile(m100$SES.FDis, c(0.02, 0.98), na.rm = TRUE)
  LIM <- c(-1, 1) * max(abs(qs))
  maps <- lapply(TAXA$class, function(cl) {
    dd <- m100 |> filter(class == cl) |> select(cell_id, val = SES.FDis)
    if (!nrow(dd)) return(NULL)
    map_china(g100 |> inner_join(dd, by = "cell_id"), "val", BM, "div",
              title = TAXA_LAB[cl], legend = "SES functional dispersion",
              limits = LIM, inset = FALSE)
  }) |> Filter(f = Negate(is.null))
  for (i in seq_along(maps))
    maps[[i]] <- maps[[i]] + labs(tag = letters[i]) +
      theme(legend.position = "bottom", legend.key.width = unit(34, "pt"),
            legend.key.height = unit(3.5, "pt"))
  top <- wrap_plots(maps, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom")

  # 50 km 与 100 km 的逐格对应（用嵌套关系）/ cell-wise 50-vs-100 km comparison
  nest <- readRDS(file.path(PATH$derived, "grid_nesting_50_in_100.rds"))
  cmp <- m50 |> select(cell_50km = cell_id, class, s50 = SES.FDis) |>
    left_join(nest, by = "cell_50km") |>
    left_join(m100 |> select(cell_100km = cell_id, class, s100 = SES.FDis),
              by = c("cell_100km", "class")) |>
    filter(is.finite(s50), is.finite(s100))
  rr <- cmp |> group_by(class) |>
    summarise(r = cor(s50, s100, method = "spearman"), .groups = "drop")
  pe <- ggplot(cmp, aes(s50, s100, colour = class)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25,
                colour = "grey45") +
    geom_point(size = 0.25, alpha = 0.15) +
    geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.4, se = FALSE) +
    geom_text(data = rr, inherit.aes = FALSE,
              aes(x = -Inf, y = Inf, label = sprintf("%s: r = %.2f",
                                                     TAXA_LAB[class], r),
                  colour = class),
              hjust = -0.08, vjust = seq(1.6, by = 1.35, length.out = nrow(rr)),
              size = 1.9, show.legend = FALSE) +
    scale_colour_manual(values = PAL$taxa, guide = "none") +
    labs(x = "SES functional dispersion, 50 km",
         y = "SES functional dispersion, 100 km", tag = "e") +
    theme_pub()

  rich <- m50 |> select(cell_50km = cell_id, class, r50 = richness) |>
    left_join(nest, by = "cell_50km") |>
    left_join(m100 |> select(cell_100km = cell_id, class, r100 = richness),
              by = c("cell_100km", "class")) |>
    filter(is.finite(r50), is.finite(r100))
  pf <- ggplot(rich, aes(r50, r100, colour = class)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25,
                colour = "grey45") +
    geom_point(size = 0.25, alpha = 0.15) +
    scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    labs(x = "Species richness, 50 km", y = "Species richness, 100 km", tag = "f") +
    theme_pub() + theme(legend.position = c(0.98, 0.02),
                        legend.justification = c(1, 0))

  p <- top / (pe | pf) + plot_layout(heights = c(1.3, 1))
  save_fig(p, "ED7_grain_robustness_100km", W2, 115, FIG)
}

# ===============================================================
# ED8 | 零模型行为：SES 对丰富度的依赖与观测-零分布对照
# ED8 | Null-model behaviour: SES versus richness, observed vs null
# ===============================================================
ed8 <- function() {
  pa <- ggplot(m50[is.finite(m50$SES.FDis), ],
               aes(richness, SES.FDis, colour = class)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, colour = "grey45") +
    geom_point(size = 0.2, alpha = 0.12) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.5, se = FALSE) +
    scale_x_continuous(trans = "log10") +
    scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    labs(x = expression("Species richness ("*log[10]*" scale)"),
         y = "SES functional dispersion", tag = "a",
         subtitle = "A flat relationship indicates the null model has removed the richness effect") +
    theme_pub() + theme(legend.position = c(0.98, 0.02),
                        legend.justification = c(1, 0),
                        plot.subtitle = element_text(size = 5))

  # 观测 FDis 与零期望的对照 / observed FDis against null expectation
  pb <- ggplot(m50[is.finite(m50$FDis), ], aes(richness, FDis, colour = class)) +
    geom_point(size = 0.2, alpha = 0.12) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.5, se = FALSE) +
    scale_x_continuous(trans = "log10") +
    scale_colour_manual(values = PAL$taxa, guide = "none") +
    labs(x = expression("Species richness ("*log[10]*" scale)"),
         y = "Observed functional dispersion", tag = "b") +
    theme_pub()

  # SES 的分布：|SES|>2 的比例 / distribution of SES with the |SES|>2 threshold
  ses_long <- m50 |>
    select(class, SES.FDis, SES.FRic, SES.RaoQ) |>
    pivot_longer(-class, names_to = "metric", values_to = "v") |>
    filter(is.finite(v)) |>
    mutate(metric = recode(metric, SES.FDis = "Functional dispersion",
                           SES.FRic = "Functional richness", SES.RaoQ = "Rao's Q"))
  pc <- ggplot(ses_long, aes(v, fill = class)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, colour = "grey45") +
    geom_vline(xintercept = c(-2, 2), linetype = 3, linewidth = 0.25, colour = "grey45") +
    geom_density(alpha = 0.45, linewidth = 0.2) +
    scale_fill_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    labs(x = "Standardised effect size", y = "Density", tag = "c") +
    theme_pub() + theme(legend.position = "bottom")

  p <- (pa | pb) / pc + plot_layout(heights = c(1, 1))
  save_fig(p, "ED8_null_model_behaviour", W2, 125, FIG)
}

# ===============================================================
# ED9 | 稳健性：留一纲检验与上游插补敏感性
# ED9 | Robustness: leave-one-class-out and imputation sensitivity
# ===============================================================
ed9 <- function() {
  f8 <- file.path(TBL, "table_8_leave_one_class_out.csv")
  f9 <- file.path(TBL, "table_9_imputation_sensitivity.csv")
  if (!all(file.exists(c(f8, f9)))) { log_msg("  [skip] ED9"); return(invisible(NULL)) }
  loo <- rd("table_8_leave_one_class_out.csv")
  sen <- rd("table_9_imputation_sensitivity.csv")

  MET <- c(coupling = "All axes", coupling_thermal = "Thermal energy",
           coupling_human = "Human pressure")
  loo$metric_f <- factor(MET[loo$metric], levels = unname(MET))
  loo$scen_f <- factor(loo$scenario, levels = rev(unique(loo$scenario)))
  loo$sig <- ifelse(loo$p_wilcoxon < 0.05, "P < 0.05", "n.s.")

  pa <- ggplot(loo, aes(ratio, scen_f, colour = metric_f, shape = sig)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey45") +
    geom_point(size = 1.6, position = position_dodge(width = 0.55)) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    scale_shape_manual(values = c(`P < 0.05` = 16, `n.s.` = 1), name = NULL) +
    scale_x_continuous(limits = c(0.9, 2.9)) +
    labs(x = "Coupling ratio (ectotherm / endotherm)", y = NULL, tag = "a",
         subtitle = "Ratios exceed 1 in every leave-one-class-out scenario") +
    theme_pub() + theme(plot.subtitle = element_text(size = 5, colour = "grey30"),
                        legend.position = "bottom", legend.box = "vertical",
                        legend.spacing.y = unit(0, "pt"))

  sen$class_f <- factor(TAXA_LAB[sen$class], levels = unname(TAXA_LAB))
  sl <- sen |>
    select(class_f, class, `All species` = mean_SES_all,
           `Observed traits only` = mean_SES_observed_only) |>
    tidyr::pivot_longer(c(`All species`, `Observed traits only`),
                        names_to = "subset", values_to = "ses")
  pb <- ggplot(sl, aes(class_f, ses, fill = subset)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, colour = "grey45") +
    geom_col(position = position_dodge(width = 0.75), width = 0.66) +
    geom_text(aes(label = sprintf("%+.2f", ses),
                  vjust = ifelse(ses >= 0, -0.4, 1.3)),
              position = position_dodge(width = 0.75), size = 1.8) +
    scale_fill_manual(values = c(`All species` = "grey65",
                                 `Observed traits only` = "#2E6E8E"), name = NULL) +
    labs(x = NULL, y = "Mean SES of functional dispersion", tag = "b",
         subtitle = "Only the amphibian result changes sign") +
    theme_pub() + theme(legend.position = "bottom",
                        plot.subtitle = element_text(size = 5, colour = "grey30"))

  pc <- ggplot(sen, aes(class_f, pct_species_retained, fill = class)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.0f%%\n(%d spp.)", pct_species_retained,
                                  n_species_observed_only)),
              vjust = -0.3, size = 1.7, lineheight = 0.9) +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30)), limits = c(0, 100)) +
    labs(x = NULL, y = "Species with fully observed\ncore traits (%)", tag = "c",
         subtitle = "Amphibian traits are overwhelmingly imputed upstream") +
    theme_pub() + theme(plot.subtitle = element_text(size = 5, colour = "grey30"))

  p <- pa / (pb | pc) + plot_layout(heights = c(1.15, 1))
  save_fig(p, "ED9_robustness_checks", W2, 135, FIG)
}

for (fn in list(ed1, ed2, ed3, ed4, ed5, ed6, ed7, ed8, ed9)) {
  r <- try(fn(), silent = FALSE)
  if (inherits(r, "try-error")) log_msg("  [warn] 一张扩展图失败，继续 / one ED figure failed")
}
log_msg("=== 08 完成 / done ===")
