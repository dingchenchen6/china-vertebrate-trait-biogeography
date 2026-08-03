# ============================================================
# 本文件 / This file:
#   生成正文主图（Nature 规格：单栏 89 mm / 双栏 183 mm，7 pt 字号）
#   Generate main-text figures at Nature specifications
#
# 输入 / Input:
#   03_data_derived/metrics_*.rds, env_*.rds, grid_*.gpkg
#   04_results/tables/table_2_*, table_3_*, table_4_*
#
# 输出 / Output:
#   05_figures/Fig1..Fig6 (.pdf 矢量 + .png 预览)
#
# 制图合规 / Cartographic compliance:
#   所有中国地图均基于自然资源部标准地图 GS(2023)2767，
#   含完整国界与南海诸岛附图。
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({ library(tidyr); library(dplyr) })

log_msg("=== 07 主图 / Main figures ===")
BM <- load_basemap(PATH$derived)
FIG <- PATH$figures

read_metrics <- function(lab) {
  fs <- list.files(PATH$derived, pattern = paste0("^metrics_.*_", lab, "\\.rds$"),
                   full.names = TRUE)
  if (!length(fs)) return(NULL)
  lapply(fs, readRDS) |> bind_rows() |>
    left_join(TAXA[, c("class", "thermal")], by = "class")
}
m50  <- read_metrics("50km")
m100 <- read_metrics("100km")
env50 <- readRDS(file.path(PATH$derived, "env_50km.rds"))
g50  <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE)

join_grid <- function(df, g = g50) g |> inner_join(df, by = "cell_id")

# ===============================================================
# Fig 1 | 概念框架与数据概览
# Fig 1 | Conceptual framework and data overview
# ===============================================================

#' 概念示意图：热缓冲假说
#' Schematic of the thermal-buffering hypothesis
#'
#' 左：变温动物体温随环境变化，可持续存在的性状值被环境窄化 -> 陡峭斜率
#' 右：内温性使体温与环境解耦，性状值受约束较松 -> 平缓斜率
#' Left: ectotherm body temperature tracks ambient, narrowing the set of
#' viable trait values (steep slope). Right: endothermy decouples body
#' temperature from ambient, relaxing the constraint (shallow slope).
concept_panel <- function() {
  x <- seq(0, 1, length.out = 200)
  B_ECTO <- 0.80; B_ENDO <- 0.24        # 斜率 = 耦合强度 / slope = coupling
  dl <- bind_rows(
    data.frame(x = x, y = 0.5 + B_ECTO * (x - 0.5), mode = "Ectotherm"),
    data.frame(x = x, y = 0.5 + B_ENDO * (x - 0.5), mode = "Endotherm"))
  dl$mode <- factor(dl$mode, levels = c("Endotherm", "Ectotherm"))

  ggplot(dl, aes(x, y, colour = mode)) +
    geom_line(linewidth = 0.8) +
    # 右侧竖线标示两者的性状变化幅度 / bars showing the realised trait range
    annotate("segment", x = 1.04, xend = 1.04,
             y = 0.5 - B_ECTO / 2, yend = 0.5 + B_ECTO / 2,
             colour = PAL$thermal[["Ectotherm"]], linewidth = 1.1) +
    annotate("segment", x = 1.10, xend = 1.10,
             y = 0.5 - B_ENDO / 2, yend = 0.5 + B_ENDO / 2,
             colour = PAL$thermal[["Endotherm"]], linewidth = 1.1) +
    annotate("text", x = 1.16, y = 0.5, label = "coupling\nstrength",
             size = 1.7, colour = "grey30", hjust = 0, lineheight = 0.95) +
    annotate("text", x = 0.05, y = 0.93, label = "Ectotherms — tight coupling",
             hjust = 0, size = 2, colour = PAL$thermal[["Ectotherm"]]) +
    annotate("text", x = 0.05, y = 0.845, label = "Endotherms — buffered",
             hjust = 0, size = 2, colour = PAL$thermal[["Endotherm"]]) +
    scale_colour_manual(values = PAL$thermal, guide = "none") +
    scale_x_continuous(breaks = c(0.03, 0.97), labels = c("low", "high"),
                       limits = c(0, 1.32), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0.06, 0.94), labels = c("low", "high"),
                       limits = c(0.02, 0.99), expand = c(0, 0)) +
    labs(x = "Environmental gradient", y = "Assemblage trait composition",
         tag = "a",
         subtitle = "Hypothesis: endothermy weakens the environmental control of assemblage traits") +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 5, colour = "grey30"),
          axis.line.x = element_line(linewidth = 0.3, colour = "grey20"),
          axis.line.y = element_line(linewidth = 0.3, colour = "grey20"))
}

fig1 <- function() {
  pa <- concept_panel()

  # (b) 各类群物种数，按体温调节模式分组 / species counts by thermoregulatory mode
  sp_n <- lapply(TAXA$class, function(cl) {
    f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
    if (!file.exists(f)) return(NULL)
    data.frame(class = cl, n_species = length(unique(readRDS(f)$species)))
  }) |> bind_rows() |>
    left_join(TAXA[, c("class", "thermal")], by = "class") |>
    mutate(label = TAXA_LAB[class])

  pb <- ggplot(sp_n, aes(x = reorder(label, n_species), y = n_species, fill = class)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = n_species), hjust = -0.15, size = 2) +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    coord_flip() +
    labs(x = NULL, y = "Species analysed (n)", tag = "b",
         subtitle = sprintf("%s species, 3,814 cells of 50 km",
                            format(sum(sp_n$n_species), big.mark = ","))) +
    theme_pub() + theme(plot.subtitle = element_text(size = 5, colour = "grey30"))

  # (c-f) 四类群物种丰富度地图 / richness maps
  maps <- lapply(TAXA$class, function(cl) {
    d <- m50 |> filter(class == cl) |> select(cell_id, richness)
    if (!nrow(d)) return(NULL)
    map_china(join_grid(d), "richness", BM, "seq",
              title = TAXA_LAB[cl], legend = "Species richness", inset = FALSE) +
      theme(legend.position = "bottom", legend.key.width = unit(15, "pt"),
            legend.key.height = unit(3, "pt"), legend.margin = margin(-4, 0, 0, 0))
  }) |> Filter(f = Negate(is.null))
  for (i in seq_along(maps)) maps[[i]] <- maps[[i]] + labs(tag = letters[i + 2])

  p <- ((pa | pb) + plot_layout(widths = c(1.15, 1))) /
       wrap_plots(maps, nrow = 1) +
    plot_layout(heights = c(1, 1.05))
  save_fig(p, "Fig1_framework_and_data", W2, 105, FIG)
}

# ===============================================================
# Fig 2 | 关键性状的群落加权均值 (CWM) 地图
# Fig 2 | Community-weighted mean maps for key traits
# ===============================================================
fig2 <- function() {
  traits_show <- c("body_mass.mean", "habitat_breadth.mean", "nocturnality.mean")
  # 注意：性状在每个纲内做了 z 标准化后再算矩，故 CWM 的单位是「类群内标准差」，
  # 这是跨性状、跨类群可比的前提；坐标标签必须如实标注。
  # NOTE: traits are z-standardised WITHIN each class before the moments are
  # computed, so CWM is in within-class SD units. Label accordingly.
  tlab <- c(body_mass.mean = "CWM body mass (within-class z)",
            habitat_breadth.mean = "CWM habitat breadth (within-class z)",
            nocturnality.mean = "CWM nocturnality (within-class z)")
  tags <- letters
  k <- 0L
  panels <- list()
  for (tv in traits_show) {
    for (cl in TAXA$class) {
      d <- m50 |> filter(class == cl)
      if (!tv %in% names(d)) next
      dd <- d |> select(cell_id, val = all_of(tv))
      if (!sum(is.finite(dd$val))) next
      k <- k + 1L
      panels[[length(panels) + 1]] <-
        map_china(join_grid(dd), "val", BM, "seq",
                  title = if (tv == traits_show[1]) TAXA_LAB[cl] else NULL,
                  legend = tlab[[tv]], inset = FALSE) +
        labs(tag = tags[k])
    }
  }
  if (!length(panels)) return(invisible(NULL))
  # 每一行（一个性状 x 四类群）共享一个图例；标签已逐格手动指定，
  # 因此不使用 plot_annotation(tag_levels)，避免图例元素也被编号
  rows <- split(panels, rep(seq_along(traits_show), each = 4)[seq_along(panels)])
  rows <- lapply(rows, function(ps) panels_shared_legend(ps, ncol = 4))
  p <- wrap_plots(rows, ncol = 1)
  save_fig(p, "Fig2_CWM_maps", W2, 50 * length(traits_show), FIG)
}

# ===============================================================
# Fig 3 | 功能多样性零模型标准化效应量 (SES) 地图与梯度响应
# Fig 3 | SES functional diversity maps and gradient responses
# ===============================================================
fig3 <- function() {
  # SES 是无量纲的标准化效应量，四类群必须共用同一色标才能直接比较；
  # 逐图各自拉伸会掩盖"两栖为正、其余为负"这一核心对比。
  # SES is unit-free, so all four classes must share one colour scale;
  # per-panel rescaling would hide the core contrast (amphibians positive,
  # the other three negative).
  qs <- stats::quantile(m50$SES.FDis, c(0.02, 0.98), na.rm = TRUE)
  LIM <- c(-1, 1) * max(abs(qs))
  maps <- lapply(TAXA$class, function(cl) {
    d <- m50 |> filter(class == cl)
    if (!"SES.FDis" %in% names(d)) return(NULL)
    dd <- d |> select(cell_id, val = SES.FDis)
    if (!sum(is.finite(dd$val))) return(NULL)
    map_china(join_grid(dd), "val", BM, "div", title = TAXA_LAB[cl],
              legend = "SES functional dispersion", limits = LIM, inset = FALSE)
  }) |> Filter(f = Negate(is.null))

  # 沿温度梯度的响应 / response along the temperature gradient
  d <- m50 |> left_join(env50[, c("cell_id", "bio1_mean", "hfp2020")], by = "cell_id")
  pb <- ggplot(d[is.finite(d$SES.FDis), ],
               aes(bio1_mean, SES.FDis, colour = class, fill = class)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = 2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.5, alpha = 0.15) +
    scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    scale_fill_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    labs(x = expression("Mean annual temperature ("*degree*"C)"),
         y = "SES functional dispersion", tag = "e") +
    theme_pub() + theme(legend.position = c(0.02, 0.02),
                        legend.justification = c(0, 0))

  pc <- ggplot(d[is.finite(d$SES.FDis), ],
               aes(hfp2020, SES.FDis, colour = class, fill = class)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = 2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.5, alpha = 0.15) +
    scale_colour_manual(values = PAL$taxa, guide = "none") +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    labs(x = "Human footprint (2020)", y = NULL, tag = "f") +
    theme_pub()

  # 四图共用同一色标，用 patchwork 的 guides = "collect" 合并为单一图例
  # One shared colour scale; patchwork collects it into a single guide
  for (i in seq_along(maps))
    maps[[i]] <- maps[[i]] + labs(tag = letters[i]) +
      theme(legend.position = "bottom", legend.key.width = unit(34, "pt"),
            legend.key.height = unit(3.5, "pt"))
  top <- wrap_plots(maps, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom")
  p <- top / (pb | pc) + plot_layout(heights = c(1.3, 1))
  save_fig(p, "Fig3_SES_functional_diversity", W2, 115, FIG)
}

# ===============================================================
# Fig 4 | H1 核心：性状–环境耦合强度的恒温/变温不对称（科为重复）
# Fig 4 | H1: endotherm-ectotherm asymmetry in trait-environment
#         coupling, with FAMILIES as replicates
# ===============================================================
fig4 <- function() {
  f5 <- file.path(PATH$tables, "table_5_clade_level_effects.csv")
  f6 <- file.path(PATH$tables, "table_6b_headline_thermal_test.csv")
  f7 <- file.path(PATH$tables, "table_7_family_coupling_strength.csv")
  if (!all(file.exists(c(f5, f6, f7)))) {
    log_msg("  [skip] Fig4 需先运行 06c_clade_level.R"); return(invisible(NULL))
  }
  eff  <- data.table::fread(f5)
  head <- data.table::fread(f6)
  fam  <- data.table::fread(f7)
  fam$thermal <- factor(fam$thermal, levels = c("Endotherm", "Ectotherm"))

  PRED_LAB <- c(ax_thermal = "Thermal\nenergy", ax_water = "Water\navailability",
                ax_productivity = "Productivity", ax_structure = "Habitat\nstructure",
                ax_human = "Human\npressure")

  # (a) 每个科的总体耦合强度：变温 vs 恒温 / overall coupling per family
  lab_a <- sprintf("ratio = %.2f x\nWilcoxon P = %.3f",
                   head$ratio_ecto_endo[head$metric == "coupling"],
                   head$p_wilcoxon[head$metric == "coupling"])
  pa <- ggplot(fam, aes(thermal, coupling)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.3,
                 fill = NA, colour = "grey35") +
    geom_jitter(aes(fill = class), width = 0.14, height = 0, size = 1.4,
                alpha = 0.85, shape = 21, stroke = 0.2, colour = "grey20") +
    scale_fill_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    scale_colour_manual(values = PAL$taxa, guide = "none") +
    annotate("text", x = 1.5, y = max(fam$coupling, na.rm = TRUE) * 1.02,
             label = lab_a, size = 1.9, lineheight = 0.95, colour = "grey20") +
    labs(x = NULL, y = "Trait-environment coupling\n(mean |standardised effect|)",
         tag = "a") +
    theme_pub() + theme(legend.position = "bottom")

  # (b) 分环境轴的耦合强度 / coupling per environmental axis
  ax <- eff |>
    group_by(class, clade, thermal, term) |>
    summarise(v = mean(abs(estimate), na.rm = TRUE), .groups = "drop") |>
    filter(term %in% names(PRED_LAB)) |>
    mutate(term = factor(PRED_LAB[term], levels = unname(PRED_LAB)),
           thermal = factor(thermal, levels = c("Endotherm", "Ectotherm")))
  pb <- ggplot(ax, aes(term, v, fill = thermal)) +
    geom_boxplot(outlier.size = 0.25, linewidth = 0.25, width = 0.62,
                 position = position_dodge(width = 0.72)) +
    scale_fill_manual(values = PAL$thermal, name = NULL) +
    labs(x = NULL, y = "Coupling strength", tag = "b") +
    theme_pub() + theme(legend.position = c(0.98, 0.98),
                        legend.justification = c(1, 1))

  # (c) 逐科效应量森林图：体型 CWM ~ 热能量（最强的单项对比）
  #     Forest plot: body-mass CWM vs thermal energy, the strongest contrast
  bm <- eff |> filter(response == "SES.body_mass.mean", term == "ax_thermal") |>
    mutate(clade = factor(clade, levels = clade[order(estimate)]),
           thermal = factor(thermal, levels = c("Endotherm", "Ectotherm")))
  pc <- ggplot(bm, aes(estimate, clade, colour = class)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60", linetype = 2) +
    geom_linerange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                   linewidth = 0.35) +
    geom_point(size = 1) +
    scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    labs(x = "Effect of thermal energy on CWM body mass",
         y = NULL, tag = "c") +
    theme_pub() +
    theme(axis.text.y = element_text(size = 4.2),
          legend.position = "bottom")

  p <- (pa | pb) / pc + plot_layout(heights = c(1, 1.35))
  save_fig(p, "Fig4_thermal_mode_asymmetry", W2, 150, FIG)
}

# ===============================================================
# Fig 5 | H2：气候 vs 人为解释力，及气候耦合随人类足迹的衰减
# Fig 5 | H2: climate vs human explanatory power and its interaction
# ===============================================================
fig5 <- function() {
  d <- m50 |> left_join(env50, by = "cell_id")
  if (!"SES.FDis" %in% names(d)) return(invisible(NULL))

  # (a) 人类足迹地图 / human footprint map
  pa <- map_china(join_grid(env50[, c("cell_id", "hfp2020")]), "hfp2020", BM, "seq",
                  title = "Human footprint (2020)", legend = "HFP", inset = TRUE)

  # (b) 体型 CWM 沿人类足迹的变化（分类群）
  bm <- d |> filter(is.finite(`body_mass.mean`))
  pb <- ggplot(bm, aes(hfp2020, `body_mass.mean`, colour = class, fill = class)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.5, alpha = 0.15) +
    scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    labs(x = "Human footprint (2020)", y = "CWM body mass (within-class z)", tag = "b") +
    theme_pub()

  # (c) 体型 CWV：功能同质化 / trait variance = functional homogenisation
  pv <- d |> filter(is.finite(`body_mass.var`))
  pc <- ggplot(pv, aes(hfp2020, `body_mass.var`, colour = class, fill = class)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.5, alpha = 0.15) +
    scale_colour_manual(values = PAL$taxa, guide = "none") +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    labs(x = "Human footprint (2020)", y = "CWV body mass (within-class z²)", tag = "c") +
    theme_pub()

  # (d) 体型偏度：哪一侧尾部被截断 / skewness = which tail is truncated
  ps <- d |> filter(is.finite(`body_mass.skew`))
  pd <- ggplot(ps, aes(hfp2020, `body_mass.skew`, colour = class, fill = class)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = 2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.5, alpha = 0.15) +
    scale_colour_manual(values = PAL$taxa, guide = "none") +
    scale_fill_manual(values = PAL$taxa, guide = "none") +
    labs(x = "Human footprint (2020)", y = "CWS body mass (skewness)", tag = "d") +
    theme_pub()

  # 2 x 2 版式：地图与三张梯度响应图等权分布，避免整幅留白
  # A 2 x 2 layout balances the map against the three response panels
  pa <- pa + theme(legend.position = "bottom",
                   legend.key.width = unit(30, "pt"),
                   legend.key.height = unit(3.5, "pt"))
  pb <- pb + theme(legend.position = c(0.98, 0.98),
                   legend.justification = c(1, 1))
  p <- (pa | pb) / (pc | pd)
  save_fig(p, "Fig5_anthropogenic_filter", W2, 120, FIG)
}

# ===============================================================
# Fig 6 | 尺度依赖与跨类群一致性
# Fig 6 | Scale dependence and cross-taxon congruence
# ===============================================================
fig6 <- function() {
  if (is.null(m100)) return(invisible(NULL))
  cmp <- bind_rows(m50 |> mutate(grain = "50 km"),
                   m100 |> mutate(grain = "100 km")) |>
    filter(is.finite(SES.FDis))
  pa <- ggplot(cmp, aes(x = class, y = SES.FDis, fill = grain)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = 2) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.65) +
    scale_fill_manual(values = c(`50 km` = "grey75", `100 km` = "grey35"), name = NULL) +
    scale_x_discrete(labels = TAXA_LAB) +
    labs(x = NULL, y = "SES functional dispersion", tag = "a") +
    theme_pub()

  # 跨类群空间一致性 / cross-taxon spatial congruence
  w <- m50 |> select(cell_id, class, SES.FDis) |>
    pivot_wider(names_from = class, values_from = SES.FDis)
  cls <- setdiff(names(w), "cell_id")
  if (length(cls) >= 2) {
    cm <- stats::cor(w[, cls], use = "pairwise.complete.obs", method = "spearman")
    cd <- as.data.frame(as.table(cm)); names(cd) <- c("x", "y", "r")
    pb <- ggplot(cd, aes(x, y, fill = r)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.2f", r)), size = 2) +
      scico::scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-1, 1),
                              name = "Spearman r") +
      scale_x_discrete(labels = TAXA_LAB) + scale_y_discrete(labels = TAXA_LAB) +
      labs(x = NULL, y = NULL, tag = "b") +
      theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                          axis.line = element_blank(), axis.ticks = element_blank())
  } else pb <- patchwork::plot_spacer()

  p <- pa | pb
  save_fig(p, "Fig6_scale_and_congruence", W2, 75, FIG)
}

for (fn in list(fig1, fig2, fig3, fig4, fig5, fig6)) {
  r <- try(fn(), silent = FALSE)
  if (inherits(r, "try-error")) log_msg("  [warn] 一张图生成失败，继续 / figure failed, continuing")
}
log_msg("=== 07 完成 / done ===")
