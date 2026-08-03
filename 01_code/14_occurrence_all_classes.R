# ============================================================
# 科学问题 / Scientific question:
#   §4 的核心结论「变温类群对当代人类压力的耦合是恒温类群的 2.05 倍」
#   建立在**专家范围图**之上，而 §5f 已证明范围图把人为效应压低约 5.6 倍。
#   那么用**实际观测记录**重做，这一不对称是否仍然成立？
#   The headline anthropogenic asymmetry rests on range maps, which section 5f
#   showed understate human effects roughly fivefold. Does the asymmetry
#   survive when assemblages are built from observed occurrences instead?
#
# 分析目标 / Objective:
#   把 §5f 的鸟类分析扩展到四个类群，在同一网格、同一模型框架下
#   并列比较「范围图群落」与「观测群落」的人为效应量。
#
# 输入数据 / Input:
#   02_data_raw/gbif/gbif_{Mammalia,Amphibia,Reptilia}_china.csv（GBIF API 下载）
#   04_results/cache/gbif_birds_gridded_50km.rds（鸟类，由脚本 13 缓存）
#   03_data_derived/comm_<class>_50km.rds, traits_imputed.rds, model_data.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_23_occurrence_coverage_by_class.csv
#   04_results/tables/table_24_human_effect_four_classes.csv
#   05_figures/Fig13_occurrence_four_classes.pdf/.png
#
# 关键假设 / Key assumptions:
#   - 四个类群的记录密度相差两个数量级（鸟 458 万 vs 爬行 3.5 万），因此
#     统一的完整度阈值会让变温类群几乎无格可用。本脚本对每个类群报告
#     不同阈值下的可用网格数，并以能保证 >= 150 格的最宽松阈值建模，
#     同时把完整度作为协变量。阈值的选择及其后果全部显式报告。
#   - 观测群落遗漏稀有种的方向与人为过滤同向，故必须控制完整度；
#     若控制后效应减弱，则原效应可能是采样偏倚。
#
# 主要包 / Main packages: sf, dplyr, spdep, spatialreg
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(spdep); library(spatialreg)
  library(matrixStats)
})

log_msg("=== 14 四类群观测群落 / Occurrence assemblages, four classes ===")
BM  <- load_basemap(PATH$derived)
FIG <- PATH$figures
g50 <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE) |> drop_inset("50km")
CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
          "habitat_breadth", "range_size")
AX <- c("ax_thermal", "ax_water", "ax_productivity", "ax_structure", "ax_human")

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))
env <- readRDS(file.path(PATH$derived, "model_data.rds"))[["50km"]] |>
  distinct(cell_id, .keep_all = TRUE) |> drop_inset("50km")

# ---------------------------------------------------------------
# 1. 读取并网格化各类群的观测记录 / Grid the occurrence records
# ---------------------------------------------------------------
grid_occ <- function(cl) {
  if (cl == "Aves") {
    f <- file.path(PATH$cache, "gbif_birds_gridded_50km.rds")
    if (!file.exists(f)) return(NULL)
    return(readRDS(f))
  }
  f <- file.path(PATH$raw, "gbif", sprintf("gbif_%s_china.csv", cl))
  if (!file.exists(f)) return(NULL)
  d <- data.table::fread(f)
  d <- d[!is.na(d$species) & d$species != "" &
           is.finite(d$decimalLongitude) & is.finite(d$decimalLatitude), ]
  # 与鸟类一致的过滤：坐标不确定度 <= 10 km（缺失视为通过）
  # Same filter as for birds: coordinate uncertainty at most 10 km
  u <- suppressWarnings(as.numeric(d$coordinateUncertaintyInMeters))
  d <- d[!is.finite(u) | u <= 10000, ]
  pts <- st_as_sf(d[, c("species", "decimalLongitude", "decimalLatitude")],
                  coords = c("decimalLongitude", "decimalLatitude"),
                  crs = CRS_WGS84) |> st_transform(st_crs(g50))
  st_join(pts, g50["cell_id"], join = st_within) |> st_drop_geometry() |>
    filter(!is.na(cell_id)) |> mutate(species = norm_name(species)) |>
    group_by(cell_id, species) |> summarise(n = n(), .groups = "drop")
}

occ_all <- lapply(TAXA$class, function(cl) {
  o <- cache_rds(paste0("occ_gridded_", cl), grid_occ(cl))
  if (is.null(o) || !nrow(o)) return(NULL)
  o$class <- cl; o
}) |> Filter(f = Negate(is.null)) |> bind_rows()

# ---------------------------------------------------------------
# 2. Chao 完整度与阈值诊断 / Coverage and threshold diagnostics
# ---------------------------------------------------------------
comp <- occ_all |> group_by(class, cell_id) |>
  summarise(n_rec = sum(n), n_sp = n(),
            f1 = sum(n == 1), f2 = sum(n == 2), .groups = "drop") |>
  mutate(coverage = ifelse(
    n_rec > 0 & f1 > 0,
    1 - (f1 / n_rec) * (((n_rec - 1) * f1) / ((n_rec - 1) * f1 + 2 * pmax(f2, 1))),
    ifelse(n_rec > 0, 1, NA_real_)))

# 记录密度相差两个数量级，故必须逐类群报告不同阈值下的可用网格数
# Record density spans two orders of magnitude, so report cell counts at
# several thresholds rather than imposing one and hoping for the best
thr_grid <- expand.grid(cov_min = c(0.5, 0.6, 0.7, 0.8),
                        sp_min = c(5, 10, 20), stringsAsFactors = FALSE)
diag <- lapply(TAXA$class, function(cl) {
  cc <- comp |> filter(class == cl)
  if (!nrow(cc)) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(thr_grid)), function(i) {
    data.frame(class = cl, cov_min = thr_grid$cov_min[i],
               sp_min = thr_grid$sp_min[i],
               n_cells = sum(cc$coverage >= thr_grid$cov_min[i] &
                               cc$n_sp >= thr_grid$sp_min[i], na.rm = TRUE))
  }))
}) |> bind_rows()
write_table(diag, "table_23_occurrence_coverage_by_class")
log_msg("--- 各阈值下的可用网格数 / cells passing each threshold ---")
print(diag |> pivot_wider(names_from = sp_min, values_from = n_cells,
                          names_prefix = "sp>=") |> as.data.frame())

MIN_CELLS <- 150L
pick_threshold <- function(cl) {
  d <- diag |> filter(class == cl, n_cells >= MIN_CELLS)
  if (!nrow(d)) return(NULL)
  # 在满足最小网格数的组合中取最严格的一组 / strictest that still qualifies
  d |> arrange(desc(cov_min), desc(sp_min)) |> slice(1)
}

# ---------------------------------------------------------------
# 3. 性状矩 / Trait moments
# ---------------------------------------------------------------
moments_for <- function(comm, cl) {
  tr <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  sp <- intersect(unique(comm$species), tr$species)
  if (length(sp) < 20) return(NULL)
  tr <- tr[match(sp, tr$species), ]
  use <- CORE[vapply(CORE, function(v) {
    x <- suppressWarnings(as.numeric(tr[[v]])); mean(is.finite(x)) > 0.95
  }, logical(1))]
  if (length(use) < 4) return(NULL)
  TM <- apply(as.matrix(sapply(tr[, use], as.numeric)), 2, zscore)
  rownames(TM) <- sp
  cm <- comm[comm$species %in% sp, ]
  by_cell <- split(match(cm$species, sp), cm$cell_id)
  by_cell <- by_cell[lengths(by_cell) >= 5]
  if (length(by_cell) < 50) return(NULL)
  out <- lapply(names(by_cell), function(cid) {
    M <- TM[by_cell[[cid]], , drop = FALSE]
    setNames(as.list(c(length(by_cell[[cid]]), colMeans(M), colVars(M))),
             c("richness", paste0(use, ".mean"), paste0(use, ".var")))
  })
  cbind(data.frame(cell_id = names(by_cell)), bind_rows(out))
}

fit_sar <- function(md, resp, extra = NULL) {
  if (is.null(md) || !resp %in% names(md)) return(NULL)
  d <- md |> inner_join(env[, c("cell_id", "x_albers", "y_albers", AX)],
                        by = "cell_id")
  if (!is.null(extra)) d <- d |> left_join(extra, by = "cell_id")
  vars <- c(AX, if (!is.null(extra)) setdiff(names(extra), "cell_id"))
  d$y <- d[[resp]]
  d <- d[stats::complete.cases(d[, c("y", vars, "x_albers", "y_albers")]), ]
  if (nrow(d) < 100) return(NULL)
  for (v in vars) d[[v]] <- zscore(d[[v]])
  d$y <- zscore(d$y)
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(d[, c("x_albers", "y_albers")]), k = 8))
  lw <- spdep::nb2listw(nb, style = "W")
  f <- stats::as.formula(paste("y ~", paste(vars, collapse = " + ")))
  m <- try(spatialreg::errorsarlm(f, data = d, listw = lw, zero.policy = TRUE),
           silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$Coef
  data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2],
             p = s[, 4], n = nrow(d), row.names = NULL)
}

# ---------------------------------------------------------------
# 4. 逐类群比较 / Per-class comparison
# ---------------------------------------------------------------
RESP <- c("body_mass.mean", "habitat_breadth.mean", "range_size.mean",
          "habitat_breadth.var")
rows <- list(); used_thr <- list()
for (cl in TAXA$class) {
  th <- pick_threshold(cl)
  if (is.null(th)) { log_msg("  ", cl, ": 无阈值可满足最小网格数，跳过"); next }
  used_thr[[length(used_thr) + 1]] <-
    data.frame(class = cl, cov_min = th$cov_min, sp_min = th$sp_min,
               n_cells = th$n_cells)
  good <- comp |> filter(class == cl, coverage >= th$cov_min, n_sp >= th$sp_min)
  oc <- occ_all |> filter(class == cl, cell_id %in% good$cell_id)
  m_occ <- moments_for(oc, cl)
  fmap <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  m_map <- if (file.exists(fmap)) moments_for(readRDS(fmap), cl) else NULL
  cov_cov <- good[, c("cell_id", "coverage")]
  log_msg("  ", cl, ": 阈值 coverage>=", th$cov_min, " 物种>=", th$sp_min,
          " -> 观测 ", ifelse(is.null(m_occ), 0, nrow(m_occ)), " 格")
  # 决定性对照：有记录的网格本身人类压力显著更高（+0.31 至 +0.53 vs -0.21 至 -0.34，
  # 全部 P < 2e-16）。若不控制这一点，「观测群落效应更大」可能只是因为分析被限制在
  # 人类压力更高、关系更陡的那一段。因此在**完全相同的网格子集**上重跑范围图模型。
  # Cells with records carry substantially higher human pressure than cells
  # without. Refitting the range-map models on exactly the same cell subset
  # separates a data-source effect from a cell-selection effect.
  m_map_same <- if (is.null(m_map) || is.null(m_occ)) NULL else
    m_map |> filter(cell_id %in% m_occ$cell_id)
  for (r in RESP) {
    a <- fit_sar(m_map, r);   if (!is.null(a)) { a$class <- cl; a$source <- "Range map";  a$response <- r; rows[[length(rows)+1]] <- a }
    a2 <- fit_sar(m_map_same, r)
    if (!is.null(a2)) { a2$class <- cl; a2$source <- "Range map, same cells"; a2$response <- r; rows[[length(rows)+1]] <- a2 }
    b <- fit_sar(m_occ, r);   if (!is.null(b)) { b$class <- cl; b$source <- "Occurrence"; b$response <- r; rows[[length(rows)+1]] <- b }
    cc <- fit_sar(m_occ, r, extra = cov_cov)
    if (!is.null(cc)) { cc$class <- cl; cc$source <- "Occurrence + coverage"; cc$response <- r; rows[[length(rows)+1]] <- cc }
  }
}
thr_used <- bind_rows(used_thr)
eff <- bind_rows(rows) |> filter(term %in% c(AX, "coverage")) |>
  left_join(TAXA[, c("class", "thermal")], by = "class")
eff$p_adj <- stats::p.adjust(eff$p, "BH")
write_table(eff, "table_24_human_effect_four_classes")

hum <- eff |> filter(term == "ax_human")
log_msg("--- 人为压力轴效应量 / effect of the human axis ---")
print(hum |> select(class, response, source, estimate, se, p_adj, n) |>
        mutate(across(where(is.numeric), ~ round(.x, 4))) |> as.data.frame())

# 核心检验：观测数据下恒温 vs 变温的人为耦合 / the asymmetry in occurrence data
asym <- hum |> filter(source %in% c("Range map", "Range map, same cells",
                                   "Occurrence + coverage")) |>
  group_by(source, thermal) |>
  summarise(mean_abs_effect = mean(abs(estimate), na.rm = TRUE),
            n_models = n(), .groups = "drop") |>
  pivot_wider(names_from = thermal, values_from = c(mean_abs_effect, n_models)) |>
  mutate(ratio_ecto_endo = mean_abs_effect_Ectotherm / mean_abs_effect_Endotherm)
write_table(asym, "table_25_human_asymmetry_by_source")
log_msg("--- 人为耦合的恒温/变温不对称 / asymmetry by data source ---")
print(as.data.frame(asym))

# ---------------------------------------------------------------
# 5. 制图 / Figure
# ---------------------------------------------------------------
join_grid <- function(d) g50 |> inner_join(d, by = "cell_id")
maps <- lapply(TAXA$class, function(cl) {
  d <- comp |> filter(class == cl) |> select(cell_id, coverage)
  if (!nrow(d)) return(NULL)
  th <- pick_threshold(cl)
  map_china(join_grid(d), "coverage", BM, "seq", title = TAXA_LAB[cl],
            legend = "Chao coverage", limits = c(0.3, 1), inset = FALSE) +
    labs(subtitle = if (is.null(th)) "too few usable cells" else
      sprintf("%d cells at coverage >= %.1f, >= %d spp.",
              th$n_cells, th$cov_min, th$sp_min))
}) |> Filter(f = Negate(is.null))
for (i in seq_along(maps))
  maps[[i]] <- maps[[i]] + labs(tag = letters[i]) +
    theme(plot.subtitle = element_text(size = 4.4, colour = "grey30", hjust = 0.5))
row1 <- wrap_plots(maps, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.key.width = unit(30, "pt"),
        legend.key.height = unit(3.5, "pt"))

RESP_LAB <- c(body_mass.mean = "Body mass, mean",
              habitat_breadth.mean = "Habitat breadth, mean",
              habitat_breadth.var = "Habitat breadth, variance",
              range_size.mean = "Range size, mean")
hum$response_f <- factor(RESP_LAB[hum$response], levels = rev(unname(RESP_LAB)))
hum$class_f <- factor(TAXA_LAB[hum$class], levels = unname(TAXA_LAB))
pe <- ggplot(hum, aes(estimate, response_f, colour = source)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                  position = position_dodge(width = 0.65), size = 0.18,
                  linewidth = 0.3) +
  facet_wrap(~class_f, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c(`Range map` = "grey65",
                                 `Range map, same cells` = "grey25",
                                 `Occurrence` = "#B5541F",
                                 `Occurrence + coverage` = "#2E6E8E"),
                      name = NULL) +
  labs(x = "Standardised effect of the human-pressure axis", y = NULL, tag = "e") +
  theme_pub() + theme(legend.position = "bottom",
                      legend.text = element_text(size = 4.6),
                      axis.text.y = element_text(size = 4.6),
                      strip.text = element_text(size = 6)) +
  guides(colour = guide_legend(nrow = 2))

p <- row1 / pe + plot_layout(heights = c(1.05, 1))
save_fig(p, "Fig13_occurrence_four_classes", W2, 125, FIG)

log_msg("=== 14 完成 / done ===")
