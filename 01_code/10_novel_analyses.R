# ============================================================
# 科学问题 / Scientific questions:
#   1. 性状**变异度**的地理格局是什么？它与均值格局是否解耦？
#      What is the geography of assemblage trait VARIANCE, and is it
#      decoupled from the geography of the mean?
#   2. 在中国的不同地区，被过滤的**是哪一个性状**？
#      Which trait is filtered where — does the identity of the filter
#      itself vary geographically?
#   3. 中国的脊椎动物群落可以按性状结构划分为哪些**功能区**？
#      四个纲划出的功能区是否一致？
#      Can assemblages be regionalised by trait structure, and do the four
#      classes delineate the same functional regions?
#   4. 性状–环境耦合强度在空间上是均匀的，还是在特定区域才出现？
#      恒温/变温的差距是恒定的，还是在热胁迫最强处才拉开？
#      Is trait-environment coupling spatially uniform, and does the
#      endotherm-ectotherm gap widen where thermal stress is greatest?
#
# 分析目标 / Objective:
#   利用 05 已算好的「每格 × 每性状 × 四阶矩」的 SES，做四项此前未做的
#   探索性分析并出图。全部基于已有中间结果，不重新计算零模型。
#   Four exploratory analyses built on the per-cell, per-trait, per-moment
#   SES values already produced by script 05. No null models are recomputed.
#
# 输入数据 / Input:
#   03_data_derived/metrics_<class>_<grain>.rds, model_data.rds, grid_*.gpkg
#
# 主要流程 / Workflow:
#   Fig 7  性状变异地图 + 均值-方差解耦类型图
#   Fig 8  「主导被过滤性状」地图
#   Fig 9  基于性状结构的功能分区，及跨纲一致性
#   Fig 10 地理加权的性状-环境耦合强度
#
# 预期输出 / Expected output:
#   05_figures/Fig7 .. Fig10 (.pdf + .png)
#   04_results/tables/table_10 .. table_13 (.csv)
#
# 关键假设 / Key assumptions:
#   - SES 已剔除物种丰富度的机械影响，故可跨格、跨纲直接比较。
#   - 地理加权回归采用 k 近邻固定样本量而非固定带宽，以避免稀疏区
#     （西部）与密集区的样本量差异造成的系数不可比。
#   - 功能分区的簇数由轮廓系数选择，而非人为指定。
#
# 主要包 / Main packages: sf, dplyr, cluster, ggplot2, patchwork
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(cluster); library(parallel)
})

log_msg("=== 10 新增探索分析 / Novel exploratory analyses ===")
BM  <- load_basemap(PATH$derived)
FIG <- PATH$figures
N_CORES <- max(1L, parallel::detectCores() - 2L)

CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
          "habitat_breadth", "range_size")
TRAIT_LAB <- c(body_mass = "Body mass", body_length = "Body length",
               nocturnality = "Nocturnality", verticality = "Verticality",
               habitat_breadth = "Habitat breadth", range_size = "Range size")
PRED_LAB <- c(ax_thermal = "Thermal energy", ax_water = "Water availability",
              ax_productivity = "Productivity", ax_structure = "Habitat structure",
              ax_human = "Human pressure")

g50 <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE) |> drop_inset("50km")
read_metrics <- function(lab) {
  fs <- list.files(PATH$derived, pattern = paste0("^metrics_.*_", lab, "\\.rds$"),
                   full.names = TRUE)
  lapply(fs, readRDS) |> bind_rows() |> drop_inset(lab) |>
    left_join(TAXA[, c("class", "thermal")], by = "class")
}
m50 <- read_metrics("50km")
env <- readRDS(file.path(PATH$derived, "model_data.rds"))[["50km"]] |>
  distinct(cell_id, .keep_all = TRUE) |> drop_inset("50km")
join_grid <- function(d) g50 |> inner_join(d, by = "cell_id")

# 胡焕庸线（黑河—腾冲线，Hu 1935）：中国人口与自然地理的经典分界。
# 作为独立参照叠加，用来判断由性状结构自动划出的功能区是否与之重合。
# The Hu Huanyong line (Heihe-Tengchong, Hu 1935) is China's classic
# population and physiographic divide; overlaid as an independent reference.
HU_LINE <- st_sfc(st_linestring(rbind(c(127.5, 50.2), c(98.5, 25.0))),
                  crs = CRS_WGS84) |> st_transform(st_crs(g50))

# 两栖类性状 84% 为上游插补，其群落指标不可靠（见 06d / ED9），
# 因此本脚本的所有推断性分析只用鸟、兽、爬行三个纲；两栖仅在需要
# 展示四纲对照时出现，并在图注中标注不作解读。
# Amphibian traits are 84% imputed upstream and their assemblage metrics are
# not reliable, so inference here uses birds, mammals and reptiles only.
CLS3 <- c("Aves", "Mammalia", "Reptilia")

# ===============================================================
# Fig 7 | 性状变异的地理格局与均值-方差解耦
# Fig 7 | Geography of trait variance and mean-variance decoupling
# ===============================================================
fig7 <- function() {
  # (a-d) 体型变异度 SES 的地图，四纲共用色标
  # (a-d) SES of body-mass variance, one shared colour scale
  d <- m50 |> select(cell_id, class, v = SES.body_mass.var) |> filter(is.finite(v))
  LIM <- c(-1, 1) * max(abs(stats::quantile(d$v, c(0.02, 0.98), na.rm = TRUE)))
  maps <- lapply(TAXA$class, function(cl) {
    dd <- d |> filter(class == cl) |> select(cell_id, v)
    if (!nrow(dd)) return(NULL)
    map_china(join_grid(dd), "v", BM, "div", title = TAXA_LAB[cl],
              legend = "SES of body-mass variance", limits = LIM, inset = FALSE)
  }) |> Filter(f = Negate(is.null))
  for (i in seq_along(maps))
    maps[[i]] <- maps[[i]] + labs(tag = letters[i])
  top <- wrap_plots(maps, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom", legend.key.width = unit(34, "pt"),
          legend.key.height = unit(3.5, "pt"))

  # (e) 均值 SES 与方差 SES 的空间相关：若两者高度相关则方差不含新信息
  # (e) Spatial correlation between mean-SES and variance-SES per trait
  mv <- lapply(CLS3, function(cl) {
    x <- m50 |> filter(class == cl)
    do.call(rbind, lapply(CORE, function(tt) {
      a <- x[[paste0("SES.", tt, ".mean")]]; b <- x[[paste0("SES.", tt, ".var")]]
      ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 100) return(NULL)
      data.frame(class = cl, trait = tt,
                 r = stats::cor(a[ok], b[ok], method = "spearman"))
    }))
  }) |> bind_rows()
  mv$trait_f <- factor(TRAIT_LAB[mv$trait], levels = rev(unname(TRAIT_LAB)))
  pe <- ggplot(mv, aes(r, trait_f, fill = class)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    scale_fill_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(x = "Spearman correlation between SES of mean and SES of variance",
         y = NULL, tag = "e",
         subtitle = "Range size is the one trait whose mean and variance move in opposite directions in every class") +
    theme_pub() + theme(legend.position = "bottom",
                        legend.key.size = unit(6, "pt"),
                        legend.text = element_text(size = 5),
                        plot.subtitle = element_text(size = 4.6, colour = "grey30"))

  # (f) 过滤机制类型图：按 (均值 SES, 方差 SES) 的符号把每格分为四类
  #     这把「方向性过滤」与「过滤强度」组合成一张可读的机制地图
  # (f) Filtering-regime map: each cell classified by the signs of its
  #     mean-SES and variance-SES, combining direction with strength
  reg <- m50 |> filter(class == "Aves") |>
    transmute(cell_id,
              mu = SES.body_mass.mean, va = SES.body_mass.var) |>
    filter(is.finite(mu), is.finite(va)) |>
    mutate(regime = case_when(
      mu > 0 & va < 0 ~ "Converged on large",
      mu < 0 & va < 0 ~ "Converged on small",
      mu > 0 & va > 0 ~ "Expanded, large-shifted",
      TRUE            ~ "Expanded, small-shifted"))
  REG_COL <- c("Converged on large"      = "#8C2D04",
               "Converged on small"      = "#2C7FB8",
               "Expanded, large-shifted" = "#FDBB84",
               "Expanded, small-shifted" = "#A6D8EF")
  pf <- ggplot() +
    geom_sf(data = BM$land, fill = "grey93", colour = NA) +
    geom_sf(data = join_grid(reg), aes(fill = regime), colour = NA) +
    add_basemap_lines(BM) +
    scale_fill_manual(values = REG_COL, name = NULL) +
    coord_sf(xlim = CN_XLIM, ylim = CN_YLIM, expand = FALSE) +
    labs(tag = "f", title = "Body-mass filtering regime (birds)") +
    theme_map() + theme(legend.position = "bottom",
                        legend.text = element_text(size = 5),
                        legend.key.size = unit(6, "pt"))

  tab <- m50 |> filter(class %in% CLS3) |>
    group_by(class) |>
    summarise(across(all_of(paste0("SES.", CORE, ".var")),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  write_table(tab, "table_10_variance_ses_by_class")
  write_table(mv, "table_11_mean_variance_decoupling")

  p <- top / (pe | pf) + plot_layout(heights = c(1, 1.15))
  save_fig(p, "Fig7_trait_variance_geography", W2, 118, FIG)
}

# ===============================================================
# Fig 8 | 主导被过滤性状：过滤器的身份随地理变化
# Fig 8 | The dominant filtered trait varies geographically
# ===============================================================
fig8 <- function() {
  # 对每格取 |SES.均值| 最大的核心性状作为「主导被过滤性状」
  # For each cell, the core trait with the largest |SES of the mean|
  dom <- lapply(CLS3, function(cl) {
    x <- m50 |> filter(class == cl)
    M <- as.matrix(x[, paste0("SES.", CORE, ".mean")])
    ok <- rowSums(is.finite(M)) == length(CORE)
    if (!any(ok)) return(NULL)
    M <- M[ok, , drop = FALSE]
    j <- max.col(abs(M), ties.method = "first")
    data.frame(cell_id = x$cell_id[ok], class = cl,
               dom_trait = CORE[j],
               dom_ses = M[cbind(seq_len(nrow(M)), j)])
  }) |> bind_rows()
  # 固定因子水平，否则各纲缺少的性状会让 patchwork 无法合并图例
  # Fixed factor levels; otherwise traits absent from a class split the legend
  dom$dom_trait <- factor(dom$dom_trait, levels = CORE, labels = TRAIT_LAB[CORE])

  TCOL <- setNames(
    c("#B2182B", "#EF8A62", "#4D4D4D", "#1B7837", "#2166AC", "#9970AB"),
    CORE)
  maps <- lapply(CLS3, function(cl) {
    dd <- dom |> filter(class == cl) |> select(cell_id, dom_trait)
    ggplot() +
      geom_sf(data = BM$land, fill = "grey93", colour = NA) +
      geom_sf(data = join_grid(dd), aes(fill = dom_trait), colour = NA) +
      add_basemap_lines(BM) +
      scale_fill_manual(values = setNames(TCOL, TRAIT_LAB[CORE]), name = NULL,
                        drop = FALSE) +
      coord_sf(xlim = CN_XLIM, ylim = CN_YLIM, expand = FALSE) +
      labs(title = TAXA_LAB[cl]) + theme_map()
  })
  for (i in seq_along(maps)) maps[[i]] <- maps[[i]] + labs(tag = letters[i])
  top <- wrap_plots(maps, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom", legend.text = element_text(size = 5),
          legend.key.size = unit(6, "pt"))

  # (d) 各性状成为主导过滤性状的网格占比
  # (d) Share of cells in which each trait dominates
  sh <- dom |> count(class, dom_trait, .drop = FALSE) |> group_by(class) |>
    mutate(pct = 100 * n / sum(n)) |> ungroup()
  sh$trait_f <- factor(sh$dom_trait, levels = unname(TRAIT_LAB))
  sh$class_f <- factor(TAXA_LAB[sh$class], levels = unname(TAXA_LAB))
  pd <- ggplot(sh, aes(class_f, pct, fill = trait_f)) +
    geom_col(width = 0.62) +
    scale_fill_manual(values = setNames(TCOL, TRAIT_LAB[CORE]), name = NULL,
                      drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Cells where the trait dominates (%)", tag = "d") +
    theme_pub() + theme(legend.position = "none")

  # (e) 主导性状随热能量梯度的变化 —— 过滤器身份是否随气候切换
  # (e) Does the identity of the dominant filter switch along the thermal gradient
  de <- dom |> left_join(env[, c("cell_id", "ax_thermal")], by = "cell_id") |>
    filter(is.finite(ax_thermal)) |>
    mutate(bin = cut(ax_thermal, breaks = stats::quantile(ax_thermal,
                     seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE,
                     labels = paste0(seq(5, 95, 10), "%"))) |>
    count(class, bin, dom_trait) |> group_by(class, bin) |>
    mutate(pct = 100 * n / sum(n)) |> ungroup()
  de$trait_f <- factor(de$dom_trait, levels = unname(TRAIT_LAB))
  de$class_f <- factor(TAXA_LAB[de$class], levels = unname(TAXA_LAB))
  pe <- ggplot(de, aes(bin, pct, fill = trait_f)) +
    geom_col(width = 0.95) +
    facet_wrap(~class_f, nrow = 1) +
    scale_fill_manual(values = setNames(TCOL, TRAIT_LAB[CORE]), name = NULL,
                      drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Thermal-energy percentile (cold → warm)",
         y = "Cells (%)", tag = "e") +
    theme_pub() + theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1))

  write_table(sh, "table_12_dominant_filtered_trait")
  p <- top / ((pd | pe) + plot_layout(widths = c(1, 2.4), guides = "collect") &
                theme(legend.position = "right",
                      legend.text = element_text(size = 5),
                      legend.key.size = unit(6, "pt"))) +
    plot_layout(heights = c(1.3, 1))
  save_fig(p, "Fig8_dominant_filtered_trait", W2, 125, FIG)
}

# ===============================================================
# Fig 9 | 基于性状结构的功能分区
# Fig 9 | Functional regionalisation from trait structure
# ===============================================================
fig9 <- function() {
  regionalise <- function(cl, kmax = 8) {
    x <- m50 |> filter(class == cl)
    V <- paste0("SES.", rep(CORE, each = 2), c(".mean", ".var"))
    V <- intersect(V, names(x))
    M <- as.matrix(x[, V])
    ok <- rowSums(is.finite(M)) == ncol(M)
    M <- M[ok, , drop = FALSE]; ids <- x$cell_id[ok]
    if (nrow(M) < 200) return(NULL)
    M <- apply(M, 2, zscore)
    # 用轮廓系数选簇数，而非人为指定 / choose k by average silhouette width
    set.seed(SEED)
    sub <- sample.int(nrow(M), min(1500, nrow(M)))
    sil <- vapply(2:kmax, function(k) {
      km <- stats::kmeans(M[sub, ], centers = k, nstart = 10, iter.max = 50)
      mean(cluster::silhouette(km$cluster, dist(M[sub, ]))[, 3])
    }, numeric(1))
    k <- (2:kmax)[which.max(sil)]
    set.seed(SEED)
    km <- stats::kmeans(M, centers = k, nstart = 25, iter.max = 100)
    list(df = data.frame(cell_id = ids, region = factor(km$cluster)),
         k = k, sil = max(sil), centers = km$centers)
  }
  res <- lapply(CLS3, regionalise) |> setNames(CLS3)
  res <- res[!vapply(res, is.null, logical(1))]

  # 量化功能分区与胡焕庸线的吻合度，而非仅凭目视断言。
  # 分区是无监督的，胡线未进入任何输入，故一致率是真正的外部验证。
  # Quantify agreement with the Hu line rather than asserting it visually.
  # The partition is unsupervised and the line was never an input, so the
  # agreement is genuine external validation.
  hu_side <- st_drop_geometry(g50) |>
    transmute(cell_id,
              se = ((lon - 98.5) * (50.2 - 25.0) -
                    (lat - 25.0) * (127.5 - 98.5)) > 0)
  hu_agree <- vapply(names(res), function(cl) {
    x <- m50 |> filter(class == cl)
    V <- intersect(paste0("SES.", rep(CORE, each = 2), c(".mean", ".var")), names(x))
    M <- as.matrix(x[, V]); ok <- rowSums(is.finite(M)) == ncol(M)
    M <- apply(M[ok, , drop = FALSE], 2, zscore); set.seed(SEED)
    d <- data.frame(cell_id = x$cell_id[ok],
                    reg = stats::kmeans(M, 2, nstart = 25)$cluster) |>
      inner_join(hu_side, by = "cell_id")
    tb <- table(d$reg, d$se)
    max((tb[1, 1] + tb[2, 2]) / sum(tb), (tb[1, 2] + tb[2, 1]) / sum(tb))
  }, numeric(1))

  maps <- lapply(names(res), function(cl) {
    r <- res[[cl]]
    ggplot() +
      geom_sf(data = BM$land, fill = "grey93", colour = NA) +
      geom_sf(data = join_grid(r$df), aes(fill = region), colour = NA) +
      add_basemap_lines(BM) +
      geom_sf(data = HU_LINE, colour = "grey10", linewidth = 0.4, linetype = 2) +
      scale_fill_brewer(palette = "Set2", name = "Functional\nregion") +
      coord_sf(xlim = CN_XLIM, ylim = CN_YLIM, expand = FALSE) +
      labs(title = sprintf("%s (k = %d)", TAXA_LAB[cl], r$k),
           subtitle = sprintf("%.0f%% agreement with the Hu line",
                              100 * hu_agree[[cl]])) +
      theme_map() + theme(legend.position = "none",
                          plot.subtitle = element_text(size = 4.6,
                                                       colour = "grey30",
                                                       hjust = 0.5))
  })
  for (i in seq_along(maps)) maps[[i]] <- maps[[i]] + labs(tag = letters[i])

  # (d) 跨纲一致性：任意两纲的分区是否划出相同的边界（调整兰德指数）
  # (d) Cross-class congruence of the regionalisations (adjusted Rand index)
  ari <- function(a, b) {
    tb <- table(a, b); n <- sum(tb)
    s <- function(x) sum(choose(x, 2))
    idx <- s(tb); ea <- s(rowSums(tb)); eb <- s(colSums(tb))
    exp <- ea * eb / choose(n, 2); mx <- (ea + eb) / 2
    (idx - exp) / (mx - exp)
  }
  # 各纲的最优 k 不同（鸟 2、兽 3、爬行 2），会把「分区数不同」混进 ARI。
  # 因此额外做一次强制 k = 2 的对照，分离「边界位置不同」与「分区数不同」。
  # Optimal k differs between classes, which confounds the ARI; a forced k = 2
  # comparison separates a different boundary from a different number of regions.
  fixed_k2 <- lapply(names(res), function(cl) {
    x <- m50 |> filter(class == cl)
    V <- intersect(paste0("SES.", rep(CORE, each = 2), c(".mean", ".var")), names(x))
    M <- as.matrix(x[, V]); ok <- rowSums(is.finite(M)) == ncol(M)
    M <- apply(M[ok, , drop = FALSE], 2, zscore)
    set.seed(SEED)
    data.frame(cell_id = x$cell_id[ok],
               region = factor(stats::kmeans(M, 2, nstart = 25)$cluster))
  }) |> setNames(names(res))

  cn <- names(res); pr <- list()
  for (i in seq_along(cn)) for (j in seq_along(cn)) {
    if (i == j) { pr[[length(pr) + 1]] <- data.frame(x = cn[i], y = cn[j],
                                                     ari = NA, ari_k2 = NA); next }
    mm <- inner_join(res[[cn[i]]]$df, res[[cn[j]]]$df, by = "cell_id")
    m2 <- inner_join(fixed_k2[[cn[i]]], fixed_k2[[cn[j]]], by = "cell_id")
    pr[[length(pr) + 1]] <- data.frame(
      x = cn[i], y = cn[j],
      ari    = if (nrow(mm) > 100) ari(mm$region.x, mm$region.y) else NA,
      ari_k2 = if (nrow(m2) > 100) ari(m2$region.x, m2$region.y) else NA)
  }
  pr <- bind_rows(pr)
  pr$x_f <- factor(TAXA_LAB[pr$x], levels = unname(TAXA_LAB))
  pr$y_f <- factor(TAXA_LAB[pr$y], levels = rev(unname(TAXA_LAB)))
  pd <- ggplot(pr, aes(x_f, y_f, fill = ari)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(ari), "",
                                 sprintf("%.2f\n(%.2f)", ari, ari_k2))),
              size = 1.9, lineheight = 0.95) +
    scale_fill_scico(palette = "davos", direction = -1, limits = c(-0.1, 1),
                     na.value = "grey92", name = "Adjusted\nRand index") +
    labs(x = NULL, y = NULL, tag = "d",
         subtitle = "Do the classes delineate the same regions?\nUpper value: own optimal k. In brackets: k forced to 2") +
    theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                        plot.subtitle = element_text(size = 4.6, colour = "grey30"))

  # (e) 每个区的性状特征 —— 让分区可解释，而不只是色块
  # (e) Trait signature of each region, so the map is interpretable
  sigs <- lapply(names(res), function(cl) {
    ct <- res[[cl]]$centers
    data.frame(class = cl, region = rep(rownames(ct), ncol(ct)),
               var = rep(colnames(ct), each = nrow(ct)),
               val = as.vector(ct))
  }) |> bind_rows() |>
    filter(grepl("\\.mean$", var)) |>
    mutate(trait = sub("^SES\\.", "", sub("\\.mean$", "", var)))
  sigs$trait_f <- factor(TRAIT_LAB[sigs$trait], levels = rev(unname(TRAIT_LAB)))
  sigs$class_f <- factor(TAXA_LAB[sigs$class], levels = unname(TAXA_LAB))
  lim <- max(abs(sigs$val), na.rm = TRUE)
  pe <- ggplot(sigs, aes(region, trait_f, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_wrap(~class_f, nrow = 1, scales = "free_x") +
    scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim),
                     name = "Cluster centre\n(z of SES mean)") +
    labs(x = "Functional region", y = NULL, tag = "e") +
    theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                        axis.text.y = element_text(size = 4.6),
                        legend.position = "right",
                        legend.key.width = unit(5, "pt"))

  write_table(pr, "table_13_regionalisation_congruence")
  write_table(data.frame(class = names(hu_agree),
                         agreement_with_hu_line = round(hu_agree, 3)),
              "table_16_hu_line_agreement")
  p <- wrap_plots(maps, nrow = 1) / (pd | pe) +
    plot_layout(heights = c(1.3, 1), widths = c(1, 2.3))
  save_fig(p, "Fig9_functional_regionalisation", W2, 125, FIG)
}

# ===============================================================
# Fig 10 | 地理加权的性状-环境耦合强度
# Fig 10 | Geographically weighted trait-environment coupling
# ===============================================================
fig10 <- function() {
  K <- 200L   # 每个局部窗口的近邻数 / neighbours per local window
  RESP <- "SES.body_mass.mean"

  local_coupling <- function(cl) {
    x <- m50 |> filter(class == cl) |>
      select(cell_id, y = all_of(RESP)) |>
      inner_join(env[, c("cell_id", "x_albers", "y_albers", names(PRED_LAB))],
                 by = "cell_id") |>
      filter(is.finite(y))
    if (nrow(x) < 3 * K) return(NULL)
    XY <- as.matrix(x[, c("x_albers", "y_albers")])
    P  <- apply(as.matrix(x[, names(PRED_LAB)]), 2, zscore)
    yv <- zscore(x$y)
    n  <- nrow(x)

    # 固定近邻数（而非固定带宽），使西部稀疏区与东部密集区的估计可比
    # Fixed neighbour count keeps estimates comparable between the sparse west
    # and the dense east, which a fixed bandwidth would not
    out <- parallel::mclapply(seq_len(n), function(i) {
      d2 <- (XY[, 1] - XY[i, 1])^2 + (XY[, 2] - XY[i, 2])^2
      idx <- order(d2)[seq_len(K)]
      Xi <- cbind(1, P[idx, , drop = FALSE]); yi <- yv[idx]
      fit <- try(.lm.fit(Xi, yi), silent = TRUE)
      if (inherits(fit, "try-error")) return(rep(NA_real_, 2))
      b <- fit$coefficients[-1]
      c(mean(abs(b), na.rm = TRUE), abs(b[1]))   # 总体耦合, 热能量耦合
    }, mc.cores = N_CORES)
    o <- do.call(rbind, out)
    data.frame(cell_id = x$cell_id, class = cl,
               coup_all = o[, 1], coup_thermal = o[, 2])
  }

  gw <- lapply(CLS3, local_coupling) |> bind_rows()
  if (!nrow(gw)) { log_msg("  [skip] Fig10"); return(invisible(NULL)) }

  LIM <- c(0, stats::quantile(gw$coup_thermal, 0.98, na.rm = TRUE))
  maps <- lapply(CLS3, function(cl) {
    dd <- gw |> filter(class == cl) |> select(cell_id, coup_thermal)
    if (!nrow(dd)) return(NULL)
    map_china(join_grid(dd), "coup_thermal", BM, "seq", title = TAXA_LAB[cl],
              legend = "Local coupling to thermal energy",
              limits = LIM, inset = FALSE)
  }) |> Filter(f = Negate(is.null))
  for (i in seq_along(maps)) maps[[i]] <- maps[[i]] + labs(tag = letters[i])
  top <- wrap_plots(maps, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom", legend.key.width = unit(34, "pt"),
          legend.key.height = unit(3.5, "pt"))

  # (d) 局部耦合差值（变温 − 恒温）。用差值而非比值：比值是两个含噪估计
  #     之商，分母接近 0 时会爆炸，其均值被少数极端格拉高（1.54）而
  #     中位数仅 0.78——差值没有这个病态。
  # (d) Difference, not ratio: a ratio of two noisy estimates blows up when the
  #     denominator approaches zero (mean 1.54 versus median 0.78); the
  #     difference is well behaved.
  wide <- gw |> select(cell_id, class, coup_thermal) |>
    pivot_wider(names_from = class, values_from = coup_thermal) |>
    mutate(endo = rowMeans(cbind(Aves, Mammalia), na.rm = TRUE),
           ecto = Reptilia, diff = ecto - endo) |>
    filter(is.finite(diff))
  dl <- c(-1, 1) * max(abs(stats::quantile(wide$diff, c(0.02, 0.98), na.rm = TRUE)))
  pd <- map_china(join_grid(wide |> select(cell_id, diff)), "diff", BM, "div",
                  title = "Local coupling: reptiles minus birds and mammals",
                  legend = "Difference in local coupling",
                  limits = dl, inset = FALSE) +
    labs(tag = "d") +
    theme(legend.position = "bottom", legend.key.width = unit(30, "pt"),
          legend.key.height = unit(3.5, "pt"))

  # (e) 为什么局部对比不能直接与科级结果（图 4）对读：回归衰减诊断。
  #     标准化斜率会被响应变量中的噪声按 sqrt(信号/(信号+噪声)) 压低；
  #     每格物种数越少，SES 越噪，斜率被压得越狠。爬行类每格中位仅 13 种。
  # (e) Why the local contrast must not be read against the family-level result
  #     in Fig. 4: standardised slopes are attenuated by noise in the response,
  #     and species-poor assemblages give noisier SES. Reptiles hold a median of
  #     13 species per cell, against 129 for birds.
  rich <- lapply(list.files(PATH$derived, pattern = "^metrics_.*_50km\\.rds$",
                            full.names = TRUE), readRDS) |> bind_rows()
  att <- gw |> left_join(rich[, c("cell_id", "class", "richness")],
                         by = c("cell_id", "class")) |>
    filter(is.finite(coup_thermal), is.finite(richness))
  rr <- att |> group_by(class) |>
    summarise(r = stats::cor(coup_thermal, richness, method = "spearman"),
              med_rich = stats::median(richness), .groups = "drop")
  pe <- ggplot(att, aes(richness, coup_thermal, colour = class)) +
    geom_point(size = 0.15, alpha = 0.10) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                linewidth = 0.6, se = FALSE) +
    geom_text(data = rr, inherit.aes = FALSE,
              aes(x = Inf, y = Inf, colour = class,
                  label = sprintf("%s: rho = %+.2f (median %d spp.)",
                                  TAXA_LAB[class], r, as.integer(med_rich))),
              hjust = 1.02, vjust = seq(1.7, by = 1.4, length.out = nrow(rr)),
              size = 1.8, show.legend = FALSE) +
    scale_x_continuous(trans = "log10") +
    scale_colour_manual(values = PAL$taxa, guide = "none") +
    labs(x = expression("Species richness of the cell ("*log[10]*" scale)"),
         y = "Local coupling to thermal energy", tag = "e",
         subtitle = "Sparse reptile assemblages give noisy SES and attenuated local slopes") +
    theme_pub() + theme(plot.subtitle = element_text(size = 5, colour = "grey30"))

  write_table(rr, "table_15_local_coupling_attenuation")

  summ <- gw |> group_by(class) |>
    summarise(mean_coupling_all = mean(coup_all, na.rm = TRUE),
              mean_coupling_thermal = mean(coup_thermal, na.rm = TRUE),
              sd_coupling_thermal = sd(coup_thermal, na.rm = TRUE),
              .groups = "drop")
  write_table(summ, "table_14_geographically_weighted_coupling")
  saveRDS(gw, file.path(PATH$derived, "gw_coupling_50km.rds"))
  log_msg("  局部配对检验 / local paired test: reptiles higher in ",
          round(100 * mean(wide$diff > 0, na.rm = TRUE), 1), "% of cells")

  p <- top / (pd | pe) + plot_layout(heights = c(1.25, 1))
  save_fig(p, "Fig10_geographically_weighted_coupling", W2, 130, FIG)
}

for (fn in list(fig7, fig8, fig9, fig10)) {
  r <- try(fn(), silent = FALSE)
  if (inherits(r, "try-error")) log_msg("  [warn] 一张新图失败，继续 / one figure failed")
}
log_msg("=== 10 完成 / done ===")
