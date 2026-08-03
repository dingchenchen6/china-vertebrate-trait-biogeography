# ============================================================
# 科学问题 / Scientific questions:
#   1. 中国自然保护区网络是否捕捉到了脊椎动物**功能多样性**的地理格局？
#      Does China's nature reserve network capture the geography of
#      vertebrate functional diversity?
#   2. 保护空缺在哪里，四个类群的空缺是否重合？
#      若鸟类功能多样性不能替代爬行/两栖（Fig 6），则以鸟类为主要依据的
#      保护布局可能系统性遗漏变温类群——这是可检验的推论。
#      Where are the gaps, and do they coincide across classes? If avian
#      functional diversity is a poor surrogate (Fig. 6), a network shaped
#      around it should systematically miss ectotherms — a testable corollary.
#   3. 保护区内外的群落性状是否不同？（有效性的相关性证据，非因果）
#      Do assemblage traits differ inside and outside reserves?
#
# 分析目标 / Objective:
#   把性状框架接到保护规划上：覆盖度、代表性、空缺、连通性、内外对照。
#
# 输入数据 / Input:
#   全国自然保护区名录 + 矢量边界（1,028 个，2021 版，WGS84）
#   03_data_derived/metrics_<class>_50km.rds, model_data.rds, grid_50km.gpkg
#
# 主要流程 / Workflow:
#   1. 保护区 -> 50 km 网格的面积覆盖比例（按级别分层）
#   2. 代表性：功能多样性高的网格是否覆盖度也高（空间自回归）
#   3. 空缺：功能多样性高 + 覆盖度低的网格，逐类群制图并做跨类群重合度
#   4. 连通性：到最近保护区的距离、邻域内保护区面积
#   5. 内外对照：按气候/生产力/地形最近邻匹配，比较性状矩
#
# 预期输出 / Expected output:
#   03_data_derived/pa_coverage_50km.rds
#   04_results/tables/table_26 .. table_29
#   05_figures/Fig14_conservation_gaps.pdf/.png
#
# 关键假设 / Key assumptions:
#   - 数据仅含**自然保护区**，不含国家公园、风景名胜区、森林公园等其他
#     保护地类型，因此覆盖度是中国保护地体系的**下界**。
#   - 内外对照做的是匹配后的相关性比较，不能解释为因果效应：保护区选址
#     本身非随机（多在偏远、人类压力低处），匹配只能削弱而非消除该偏倚。
#
# 主要包 / Main packages: sf, dplyr, spdep, spatialreg
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(spdep); library(spatialreg)
})
sf_use_s2(FALSE)

log_msg("=== 15 保护规划分析 / Conservation analyses ===")
BM  <- load_basemap(PATH$derived)
FIG <- PATH$figures
PA_SHP <- file.path("/Users/dingchenchen/Downloads",
                    "List and Vector Boundaries of Nature Reserves in China",
                    "全国自然保护区名录+矢量边界", "保护区.shp")
stopifnot(file.exists(PA_SHP))

g50 <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE) |>
  drop_inset("50km")
env <- readRDS(file.path(PATH$derived, "model_data.rds"))[["50km"]] |>
  distinct(cell_id, .keep_all = TRUE) |> drop_inset("50km")
AX <- c("ax_thermal", "ax_water", "ax_productivity", "ax_structure", "ax_human")

read_metrics <- function(lab = "50km") {
  fs <- list.files(PATH$derived, pattern = paste0("^metrics_.*_", lab, "\\.rds$"),
                   full.names = TRUE)
  lapply(fs, readRDS) |> bind_rows() |> drop_inset(lab) |>
    left_join(TAXA[, c("class", "thermal")], by = "class")
}
m50 <- read_metrics()
# 两栖类性状 84% 为上游插补，其群落指标不可靠（ED9），故不参与保护推断
# Amphibian metrics rest on largely imputed traits and are excluded from inference
CLS3 <- c("Aves", "Mammalia", "Reptilia")

# ---------------------------------------------------------------
# 1. 保护区覆盖度 / Protected-area coverage per cell
# ---------------------------------------------------------------
pa_cov <- cache_rds("pa_coverage_50km", {
  pa <- st_read(PA_SHP, quiet = TRUE) |> st_make_valid()
  pa <- pa[!st_is_empty(pa), ]
  pa <- st_transform(pa, st_crs(g50)) |> st_make_valid()
  pa$lvl <- ifelse(grepl("国家", pa$`级别`), "national", "subnational")
  log_msg("  保护区 / reserves: ", nrow(pa),
          "（国家级 ", sum(pa$lvl == "national"), "）")

  inter <- function(p) {
    x <- st_intersection(g50["cell_id"], st_union(st_geometry(p))) |> st_make_valid()
    data.frame(cell_id = x$cell_id, a = as.numeric(st_area(x))) |>
      group_by(cell_id) |> summarise(a = sum(a), .groups = "drop")
  }
  all_a <- inter(pa)
  nat_a <- inter(pa[pa$lvl == "national", ])

  data.frame(cell_id = g50$cell_id,
             cell_m2 = as.numeric(st_area(g50))) |>
    left_join(all_a |> rename(pa_m2 = a), by = "cell_id") |>
    left_join(nat_a |> rename(pa_nat_m2 = a), by = "cell_id") |>
    mutate(across(c(pa_m2, pa_nat_m2), ~ ifelse(is.na(.x), 0, .x)),
           pa_frac = pmin(pa_m2 / cell_m2, 1),
           pa_nat_frac = pmin(pa_nat_m2 / cell_m2, 1))
})
saveRDS(pa_cov, file.path(PATH$derived, "pa_coverage_50km.rds"))
log_msg("全国自然保护区覆盖率 / national coverage: ",
        round(100 * sum(pa_cov$pa_m2) / sum(pa_cov$cell_m2), 2), "%（国家级 ",
        round(100 * sum(pa_cov$pa_nat_m2) / sum(pa_cov$cell_m2), 2), "%）")

# ---------------------------------------------------------------
# 2. 连通性 / Connectivity of the reserve network
# ---------------------------------------------------------------
conn <- cache_rds("pa_connectivity_50km", {
  ctr <- st_centroid(st_geometry(g50))
  prot <- pa_cov$cell_id[pa_cov$pa_frac > 0.01]
  pctr <- ctr[g50$cell_id %in% prot]
  d <- st_distance(ctr, pctr)
  # 到最近受保护网格的距离；邻域 200 km 内的受保护面积
  # Distance to the nearest protected cell, and protected area within 200 km
  nn <- apply(as.matrix(d), 1, min) / 1000
  xy <- cbind(g50$x_albers, g50$y_albers)
  nb <- lapply(seq_len(nrow(xy)), function(i) {
    dd <- sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
    which(dd <= 2e5)
  })
  neigh_pa <- vapply(nb, function(ix)
    sum(pa_cov$pa_m2[match(g50$cell_id[ix], pa_cov$cell_id)], na.rm = TRUE) / 1e6,
    numeric(1))
  data.frame(cell_id = g50$cell_id, dist_to_pa_km = nn,
             pa_km2_within_200km = neigh_pa)
})

# ---------------------------------------------------------------
# 3. 代表性：功能多样性高的地方保护得更好吗
# 3. Representation: is functional diversity better protected?
# ---------------------------------------------------------------
fit_sar <- function(d, yvar, xvars) {
  d <- d[stats::complete.cases(d[, c(yvar, xvars, "x_albers", "y_albers")]), ]
  if (nrow(d) < 200) return(NULL)
  dd <- d
  dd[[yvar]] <- zscore(dd[[yvar]])
  for (v in xvars) dd[[v]] <- zscore(dd[[v]])
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(dd[, c("x_albers", "y_albers")]), k = 8))
  lw <- spdep::nb2listw(nb, style = "W")
  f <- stats::as.formula(paste(yvar, "~", paste(xvars, collapse = " + ")))
  m <- try(spatialreg::errorsarlm(f, data = dd, listw = lw, zero.policy = TRUE),
           silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$Coef
  data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], p = s[, 4],
             n = nrow(dd), row.names = NULL)
}

rep_rows <- list()
for (cl in CLS3) {
  d <- m50 |> filter(class == cl) |>
    select(cell_id, richness, SES.FRic, SES.FDis, trait_vol) |>
    inner_join(pa_cov[, c("cell_id", "pa_frac")], by = "cell_id") |>
    inner_join(env[, c("cell_id", "x_albers", "y_albers", AX)], by = "cell_id")
  for (yv in c("richness", "SES.FRic", "trait_vol")) {
    r <- fit_sar(d, "pa_frac", c(yv, AX))
    if (is.null(r)) next
    r <- r |> filter(term == yv)
    r$class <- cl; r$predictor <- yv
    rep_rows[[length(rep_rows) + 1]] <- r
  }
}
repres <- bind_rows(rep_rows)
repres$p_adj <- stats::p.adjust(repres$p, "BH")
repres <- repres |> left_join(TAXA[, c("class", "thermal")], by = "class")
write_table(repres, "table_26_pa_representation")
log_msg("--- 代表性：保护覆盖度对功能多样性的响应 ---")
print(repres |> select(class, predictor, estimate, se, p_adj) |>
        mutate(across(where(is.numeric), ~ signif(.x, 3))) |> as.data.frame())

# ---------------------------------------------------------------
# 4. 保护空缺 / Conservation gaps
# ---------------------------------------------------------------
# 定义：功能多样性位于全国前 25% 而保护覆盖度位于后 25% 的网格
# Cells in the top quartile of functional diversity and the bottom quartile
# of protected-area coverage
gap_tbl <- list(); gap_cells <- list()
for (cl in CLS3) {
  d <- m50 |> filter(class == cl) |> select(cell_id, trait_vol, SES.FRic) |>
    inner_join(pa_cov[, c("cell_id", "pa_frac")], by = "cell_id") |>
    filter(is.finite(trait_vol))
  q_fd <- stats::quantile(d$trait_vol, 0.75, na.rm = TRUE)
  q_pa <- stats::quantile(d$pa_frac, 0.25, na.rm = TRUE)
  d$gap <- d$trait_vol >= q_fd & d$pa_frac <= q_pa
  gap_cells[[cl]] <- d |> filter(gap) |> pull(cell_id)
  gap_tbl[[length(gap_tbl) + 1]] <- data.frame(
    class = cl, n_cells = nrow(d), n_gap = sum(d$gap),
    pct_gap = round(100 * mean(d$gap), 1),
    mean_pa_frac_gap = round(mean(d$pa_frac[d$gap]), 4),
    mean_pa_frac_all = round(mean(d$pa_frac), 4))
}
gaps <- bind_rows(gap_tbl)
write_table(gaps, "table_27_conservation_gaps")
print(as.data.frame(gaps))

# 跨类群空缺的重合度 / overlap of gaps between classes
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
ov <- expand.grid(x = CLS3, y = CLS3, stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(jaccard = ifelse(x == y, NA_real_, jac(gap_cells[[x]], gap_cells[[y]]))) |>
  ungroup()
write_table(ov, "table_28_gap_overlap")
log_msg("--- 空缺网格的跨类群 Jaccard 重合度 ---")
print(ov |> filter(!is.na(jaccard)) |> as.data.frame())

# ---------------------------------------------------------------
# 5. 保护区内外的性状对照（匹配后）/ Matched inside-outside comparison
# ---------------------------------------------------------------
# 保护区选址非随机——多设在偏远、人类压力低、地形复杂处。因此直接比较
# 内外会把选址偏倚误读为保护成效。这里对每个受保护网格，在环境空间中
# 找最相似的未保护网格作为对照（最近邻匹配，不放回）。
# Reserves are not sited at random, so a raw comparison confounds siting with
# effectiveness. Each protected cell is matched to its nearest unprotected
# neighbour in environmental space, without replacement.
match_compare <- function(cl, resp) {
  d <- m50 |> filter(class == cl) |> select(cell_id, all_of(resp)) |>
    inner_join(pa_cov[, c("cell_id", "pa_frac")], by = "cell_id") |>
    inner_join(env[, c("cell_id", AX)], by = "cell_id") |>
    filter(is.finite(.data[[resp]]), stats::complete.cases(across(all_of(AX))))
  if (nrow(d) < 200) return(NULL)
  d$protected <- d$pa_frac >= 0.10
  if (sum(d$protected) < 30 || sum(!d$protected) < 30) return(NULL)
  P <- as.matrix(d[, AX]); P <- apply(P, 2, zscore)
  ip <- which(d$protected); iu <- which(!d$protected)
  used <- logical(length(iu)); pairs <- list()
  for (k in ip) {
    dd <- sqrt(colSums((t(P[iu, , drop = FALSE]) - P[k, ])^2))
    dd[used] <- Inf
    j <- which.min(dd)
    if (!is.finite(dd[j])) next
    used[j] <- TRUE
    pairs[[length(pairs) + 1]] <- c(d[[resp]][k], d[[resp]][iu[j]], dd[j])
  }
  if (length(pairs) < 30) return(NULL)
  M <- do.call(rbind, pairs)
  tt <- stats::wilcox.test(M[, 1], M[, 2], paired = TRUE)
  data.frame(class = cl, response = resp, n_pairs = nrow(M),
             mean_protected = mean(M[, 1]), mean_unprotected = mean(M[, 2]),
             diff = mean(M[, 1] - M[, 2]),
             mean_match_dist = mean(M[, 3]),
             p_wilcoxon = tt$p.value)
}
eff <- lapply(CLS3, function(cl)
  lapply(c("richness", "SES.FDis", "trait_vol"),
         function(r) match_compare(cl, r)) |> bind_rows()) |> bind_rows()
if (nrow(eff)) {
  eff$p_adj <- stats::p.adjust(eff$p_wilcoxon, "BH")
  write_table(eff, "table_29_pa_matched_comparison")
  log_msg("--- 匹配后的保护区内外对照 ---")
  print(eff |> mutate(across(where(is.numeric), ~ signif(.x, 3))) |> as.data.frame())
}

# ---------------------------------------------------------------
# 5b. 保护区的主保护对象构成 / What the reserves were designated for
# ---------------------------------------------------------------
# 代表性分析的结果（覆盖度只与哺乳类功能多样性正相关）需要一个机制解释。
# 保护区名录的「主保护对象」字段记录了设立依据，可直接检验。
# The representation result needs a mechanism; the designation records state
# what each reserve was established to protect, so it can be checked directly.
targets <- {
  pa <- st_read(PA_SHP, quiet = TRUE) |> st_drop_geometry()
  tgt <- ifelse(is.na(pa$`主保护对象`), "", pa$`主保护对象`)
  key <- list(Mammalia = "兽|哺乳|猴|熊猫|羚|鹿|虎|豹|象|猿|貂|獭",
              Aves     = "鸟|禽|鹤|雁|鹭|雉|鸨",
              Reptilia = "爬行|蛇|龟|鳄|蜥|鳖",
              Amphibia = "两栖|蛙|鲵|蟾")
  data.frame(class = names(key),
             n_reserves = vapply(key, function(k) sum(grepl(k, tgt)), numeric(1)),
             pct_reserves = round(100 * vapply(key, function(k) mean(grepl(k, tgt)),
                                               numeric(1)), 1),
             n_total = nrow(pa), row.names = NULL) |>
    left_join(TAXA[, c("class", "thermal")], by = "class")
}
write_table(targets, "table_30_reserve_designation_targets")
log_msg("--- 保护区主保护对象中提及各类群的比例 ---")
print(as.data.frame(targets))

# ---------------------------------------------------------------
# 6. 制图 / Figure
# ---------------------------------------------------------------
join_grid <- function(d) g50 |> inner_join(d, by = "cell_id")

pa_lim <- c(0, stats::quantile(pa_cov$pa_frac, 0.98, na.rm = TRUE))
p_cov <- map_china(join_grid(pa_cov[, c("cell_id", "pa_frac")]), "pa_frac", BM,
                   "seq", title = "Nature reserve coverage",
                   legend = "Fraction of cell protected", limits = pa_lim) +
  labs(tag = "a") +
  theme(legend.position = "bottom", legend.key.width = unit(28, "pt"),
        legend.key.height = unit(3.5, "pt"))

p_conn <- map_china(join_grid(conn[, c("cell_id", "dist_to_pa_km")]),
                    "dist_to_pa_km", BM, "seq",
                    title = "Isolation from the reserve network",
                    legend = "Distance to nearest reserve (km)",
                    limits = c(0, stats::quantile(conn$dist_to_pa_km, 0.98))) +
  labs(tag = "b") +
  theme(legend.position = "bottom", legend.key.width = unit(28, "pt"),
        legend.key.height = unit(3.5, "pt"))

gap_maps <- lapply(seq_along(CLS3), function(i) {
  cl <- CLS3[i]
  d <- data.frame(cell_id = g50$cell_id,
                  gap = factor(ifelse(g50$cell_id %in% gap_cells[[cl]],
                                      "Gap", "Other"), levels = c("Gap", "Other")))
  map_china_cat(join_grid(d), "gap", BM,
                values = c(Gap = "#B2182B", Other = "grey88"),
                title = TAXA_LAB[cl], legend = NULL) +
    labs(tag = letters[i + 2]) + theme(legend.position = "none")
})

pr <- repres |> filter(predictor == "trait_vol")
pr$class_f <- factor(TAXA_LAB[pr$class], levels = unname(TAXA_LAB))
p_rep <- ggplot(pr, aes(estimate, class_f, colour = class)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                  size = 0.25, linewidth = 0.4) +
  scale_colour_manual(values = PAL$taxa, guide = "none") +
  labs(x = "Effect of assemblage trait volume\non reserve coverage", y = NULL,
       tag = "f",
       subtitle = "Positive means functionally rich cells are better protected") +
  theme_pub() + theme(plot.subtitle = element_text(size = 4.6, colour = "grey30"))

ovp <- ov |> filter(!is.na(jaccard))
ovp$x_f <- factor(TAXA_LAB[ovp$x], levels = unname(TAXA_LAB))
ovp$y_f <- factor(TAXA_LAB[ovp$y], levels = rev(unname(TAXA_LAB)))
p_ov <- ggplot(ovp, aes(x_f, y_f, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.1) +
  scale_fill_scico(palette = "davos", direction = -1, limits = c(0, 1),
                   name = "Jaccard") +
  labs(x = NULL, y = NULL, tag = "g",
       subtitle = "Overlap of conservation-gap cells between classes") +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      plot.subtitle = element_text(size = 4.6, colour = "grey30"))

targets$class_f <- factor(TAXA_LAB[targets$class], levels = unname(TAXA_LAB))
p_tgt <- ggplot(targets, aes(pct_reserves, class_f, fill = class)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%  (n = %d)", pct_reserves, n_reserves)),
            hjust = -0.06, size = 1.8) +
  scale_fill_manual(values = PAL$taxa, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.55)),
                     breaks = c(0, 5, 10, 15)) +
  labs(x = "Reserves naming the class\nas a protection target (%)", y = NULL,
       tag = "h",
       subtitle = sprintf("Of %d nature reserves nationwide", targets$n_total[1])) +
  theme_pub() + theme(plot.subtitle = element_text(size = 4.6, colour = "grey30"))

p <- (p_cov | p_conn) /
     wrap_plots(gap_maps, nrow = 1) /
     (p_rep | p_ov | p_tgt) +
  plot_layout(heights = c(1.15, 1.05, 0.85))
save_fig(p, "Fig14_conservation_gaps", W2, 195, FIG)

log_msg("=== 15 完成 / done ===")
