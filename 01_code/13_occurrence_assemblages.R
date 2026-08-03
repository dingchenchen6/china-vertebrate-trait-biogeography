# ============================================================
# 科学问题 / Scientific question:
#   §4 的人为过滤效应基于**专家范围图**。范围图是「潜在分布」，对细尺度
#   土地利用不敏感——一片被彻底农田化的网格，只要落在物种的范围多边形内，
#   该物种就仍被计为存在。因此人为效应必然被系统性低估。
#   用**实际观测记录**构建的「现实群落」，能否给出更强的人为信号？
#   Range maps record potential occupancy and are insensitive to fine-scale
#   land use, so anthropogenic filtering is necessarily understated. Do
#   assemblages built from observed occurrences show a stronger signal?
#
# 分析目标 / Objective:
#   1. 由 GBIF/eBird 中国 2000-2025 的 458 万条鸟类记录构建观测群落
#   2. 逐格计算 Chao 采样完整度，并把它作为协变量显式建模
#   3. 在同一网格上对比「范围图群落」与「观测群落」的人为效应量
#
# 输入数据 / Input:
#   02_data_raw/gbif/0025648-250827131500795.csv（GBIF simple CSV，1.4 GB）
#   03_data_derived/comm_Aves_50km.rds, traits_imputed.rds, model_data.rds
#
# 主要流程 / Workflow:
#   1. duckdb 流式读取并聚合到 (物种, 0.02 度栅格) —— 1.4 GB 无需载入内存
#   2. 投影到 Albers 并归入 50 km 网格
#   3. Chao 覆盖度 + 性状矩 + 空间自回归
#   4. 与范围图结果并列比较
#
# 预期输出 / Expected output:
#   03_data_derived/comm_Aves_occurrence_50km.rds, completeness_50km.rds
#   04_results/tables/table_21_occurrence_vs_rangemap.csv
#   04_results/tables/table_22_anthropogenic_effect_comparison.csv
#   05_figures/Fig12_occurrence_vs_rangemap.pdf/.png
#
# 关键假设 / Key assumptions:
#   - 只保留 2000 年以后、坐标不确定度 <= 10 km 的记录，使观测窗口与
#     人类足迹 2020 的时间尺度大致匹配。
#   - 采样完整度用 Chao (2      ) 的样本覆盖度估计量；它只需单例与双例数，
#     不需要真实丰富度，适合机会性数据。
#   - 观测群落必然遗漏稀有种，因此其性状矩偏向常见种；这一偏倚的方向
#     （偏向广布、耐受性强的物种）恰与人为过滤同向，故比较时必须
#     控制完整度，否则会把采样偏倚误读为人为效应。
#
# 主要包 / Main packages: duckdb, sf, dplyr, spatialreg
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(duckdb); library(tidyr); library(dplyr)
  library(spdep); library(spatialreg); library(matrixStats)
})

log_msg("=== 13 观测记录群落 / Occurrence-based assemblages ===")
BM  <- load_basemap(PATH$derived)
FIG <- PATH$figures
CSV <- file.path(PATH$raw, "gbif", "0025648-250827131500795.csv")
stopifnot(file.exists(CSV))

# ---------------------------------------------------------------
# 1. duckdb 聚合 / Aggregate with duckdb
# ---------------------------------------------------------------
# 1.4 GB 的制表符分隔文件不载入内存，直接在 duckdb 里过滤并聚合到
# 0.02 度（约 2 km）的经纬网格，把 458 万条压到可在 R 中处理的规模。
# The 1.4 GB file is filtered and aggregated inside duckdb, never in memory.
agg <- cache_rds("gbif_birds_agg_002deg", {
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  # 列类型由抽样推断，部分数值列会被判为 VARCHAR，故一律 TRY_CAST；
  # TRY_CAST 对无法转换的值返回 NULL 而不报错，正是脏数据所需。
  # Types are inferred from a sample and some numeric columns come back as
  # VARCHAR, so everything numeric is TRY_CAST; it yields NULL rather than
  # failing on unparseable values.
  q <- sprintf("
    WITH raw AS (
      SELECT species,
             TRY_CAST(decimalLongitude AS DOUBLE) AS lon0,
             TRY_CAST(decimalLatitude  AS DOUBLE) AS lat0,
             TRY_CAST(year AS INTEGER) AS yr,
             TRY_CAST(coordinateUncertaintyInMeters AS DOUBLE) AS unc,
             basisOfRecord
      FROM read_csv_auto('%s', delim='\t', header=true, sample_size=-1,
                         ignore_errors=true, all_varchar=true)
    )
    SELECT species,
           round(lon0 / 0.02) * 0.02 AS lon,
           round(lat0 / 0.02) * 0.02 AS lat,
           COUNT(*) AS n
    FROM raw
    WHERE species IS NOT NULL AND species <> ''
      AND lon0 IS NOT NULL AND lat0 IS NOT NULL
      AND yr >= 2000
      AND (unc IS NULL OR unc <= 10000)
      AND basisOfRecord IN ('HUMAN_OBSERVATION','MACHINE_OBSERVATION',
                            'PRESERVED_SPECIMEN','OBSERVATION','OCCURRENCE')
    GROUP BY 1,2,3", CSV)
  dbGetQuery(con, q)
})
log_msg("聚合后 / after aggregation: ", format(nrow(agg), big.mark = ","),
        " 行, ", format(length(unique(agg$species)), big.mark = ","), " 物种, ",
        format(sum(agg$n), big.mark = ","), " 条记录")

# ---------------------------------------------------------------
# 2. 归入 50 km 网格 / Assign to the 50-km grid
# ---------------------------------------------------------------
g50 <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE) |> drop_inset("50km")
occ <- cache_rds("gbif_birds_gridded_50km", {
  pts <- st_as_sf(agg, coords = c("lon", "lat"), crs = CRS_WGS84) |>
    st_transform(st_crs(g50))
  hit <- st_join(pts, g50["cell_id"], join = st_within)
  st_drop_geometry(hit) |> filter(!is.na(cell_id)) |>
    mutate(species = norm_name(species)) |>
    group_by(cell_id, species) |> summarise(n = sum(n), .groups = "drop")
})
log_msg("网格内物种-记录对 / species-cell records: ",
        format(nrow(occ), big.mark = ","))

# ---------------------------------------------------------------
# 3. Chao 采样完整度 / Chao sample coverage
# ---------------------------------------------------------------
# Ĉ = 1 - (f1/n) * [ (n-1)f1 / ((n-1)f1 + 2f2) ]
# f1 单例种、f2 双例种、n 记录总数。只依赖样本本身，不需要真实丰富度。
comp <- occ |> group_by(cell_id) |>
  summarise(n_rec = sum(n), n_sp = n(),
            f1 = sum(n == 1), f2 = sum(n == 2), .groups = "drop") |>
  mutate(coverage = ifelse(
    n_rec > 0 & f1 > 0,
    1 - (f1 / n_rec) * (((n_rec - 1) * f1) / ((n_rec - 1) * f1 + 2 * pmax(f2, 1))),
    ifelse(n_rec > 0, 1, NA_real_)))
saveRDS(comp, file.path(PATH$derived, "completeness_50km.rds"))
log_msg("完整度 / coverage: 中位 ", round(median(comp$coverage, na.rm = TRUE), 3),
        "，>0.8 的网格占 ", round(100 * mean(comp$coverage > 0.8, na.rm = TRUE), 1), "%")

# 只保留采样足够充分的网格 / retain adequately surveyed cells
GOOD <- comp$cell_id[is.finite(comp$coverage) & comp$coverage >= 0.7 &
                       comp$n_sp >= 20]
log_msg("纳入分析的网格 / cells retained: ", length(GOOD), " / ", nrow(g50))
occ_use <- occ |> filter(cell_id %in% GOOD)
saveRDS(occ_use, file.path(PATH$derived, "comm_Aves_occurrence_50km.rds"))

# ---------------------------------------------------------------
# 4. 性状矩 / Trait moments
# ---------------------------------------------------------------
CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
          "habitat_breadth", "range_size")
traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds")) |>
  filter(class == "Aves") |> distinct(species, .keep_all = TRUE)

moments_for <- function(comm) {
  sp <- intersect(unique(comm$species), traits$species)
  tr <- traits[match(sp, traits$species), ]
  TM <- apply(as.matrix(sapply(tr[, CORE], as.numeric)), 2, zscore)
  rownames(TM) <- sp
  cm <- comm[comm$species %in% sp, ]
  by_cell <- split(match(cm$species, sp), cm$cell_id)
  by_cell <- by_cell[lengths(by_cell) >= 5]
  out <- lapply(names(by_cell), function(cid) {
    M <- TM[by_cell[[cid]], , drop = FALSE]
    setNames(as.list(c(length(by_cell[[cid]]), colMeans(M), colVars(M))),
             c("richness", paste0(CORE, ".mean"), paste0(CORE, ".var")))
  })
  cbind(data.frame(cell_id = names(by_cell)), bind_rows(out))
}
m_occ <- moments_for(occ_use)
m_map <- moments_for(readRDS(file.path(PATH$derived, "comm_Aves_50km.rds")))
log_msg("观测群落 ", nrow(m_occ), " 格；范围图群落 ", nrow(m_map), " 格")

# ---------------------------------------------------------------
# 5. 两套群落的对比 / Compare the two assemblage definitions
# ---------------------------------------------------------------
cmp <- inner_join(m_occ |> select(cell_id, richness_occ = richness,
                                  ends_with(".mean")) |>
                    rename_with(~ paste0(.x, "_occ"), ends_with(".mean")),
                  m_map |> select(cell_id, richness_map = richness,
                                  ends_with(".mean")) |>
                    rename_with(~ paste0(.x, "_map"), ends_with(".mean")),
                  by = "cell_id") |>
  left_join(comp[, c("cell_id", "coverage", "n_rec")], by = "cell_id")

detect <- cmp |>
  summarise(n_cells = n(),
            mean_richness_occ = mean(richness_occ),
            mean_richness_map = mean(richness_map),
            detection_ratio = mean(richness_occ / richness_map),
            r_richness = cor(richness_occ, richness_map, method = "spearman"),
            r_body_mass_cwm = cor(body_mass.mean_occ, body_mass.mean_map,
                                  method = "spearman"),
            r_habitat_breadth_cwm = cor(habitat_breadth.mean_occ,
                                        habitat_breadth.mean_map,
                                        method = "spearman")) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_table(detect, "table_21_occurrence_vs_rangemap")
print(as.data.frame(detect))

# ---------------------------------------------------------------
# 6. 人为效应量的对比 / Anthropogenic effect sizes side by side
# ---------------------------------------------------------------
env <- readRDS(file.path(PATH$derived, "model_data.rds"))[["50km"]] |>
  distinct(cell_id, .keep_all = TRUE) |> drop_inset("50km")
AX <- c("ax_thermal", "ax_water", "ax_productivity", "ax_structure", "ax_human")

fit_sar <- function(md, resp, extra = NULL) {
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

RESP <- c("body_mass.mean", "habitat_breadth.mean", "range_size.mean",
          "body_mass.var", "habitat_breadth.var")
cov_cov <- comp[, c("cell_id", "coverage")]
rows <- list()
for (r in RESP) {
  a <- fit_sar(m_map, r);                    if (!is.null(a)) { a$source <- "Range map";  a$response <- r; rows[[length(rows)+1]] <- a }
  b <- fit_sar(m_occ, r);                    if (!is.null(b)) { b$source <- "Occurrence"; b$response <- r; rows[[length(rows)+1]] <- b }
  cc <- fit_sar(m_occ, r, extra = cov_cov);  if (!is.null(cc)) { cc$source <- "Occurrence + coverage control"; cc$response <- r; rows[[length(rows)+1]] <- cc }
}
eff <- bind_rows(rows) |> filter(term %in% c(AX, "coverage"))
eff$p_adj <- stats::p.adjust(eff$p, "BH")
write_table(eff, "table_22_anthropogenic_effect_comparison")

hum <- eff |> filter(term == "ax_human") |>
  select(response, source, estimate, se, p_adj) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))
log_msg("--- 人为压力轴的效应量 / effect of the human axis ---")
print(as.data.frame(hum))

# ---------------------------------------------------------------
# 7. 制图 / Figure
# ---------------------------------------------------------------
join_grid <- function(d) g50 |> inner_join(d, by = "cell_id")
pa <- map_china(join_grid(comp |> select(cell_id, coverage)), "coverage", BM, "seq",
                title = "Sampling completeness", legend = "Chao coverage",
                limits = c(0.5, 1), inset = FALSE) +
  labs(tag = "a",
       subtitle = sprintf("Grey: no records or below threshold (%d of %d cells retained)",
                          length(GOOD), nrow(g50))) +
  theme(legend.position = "bottom", legend.key.width = unit(28, "pt"),
        legend.key.height = unit(3.5, "pt"),
        plot.subtitle = element_text(size = 4.6, colour = "grey30", hjust = 0.5))

pb <- ggplot(cmp, aes(richness_map, richness_occ)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25,
              colour = "grey45") +
  geom_point(size = 0.2, alpha = 0.15, colour = PAL$taxa[["Aves"]]) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
              linewidth = 0.6, colour = "grey15", se = FALSE) +
  labs(x = "Richness from range maps", y = "Richness from occurrences", tag = "b",
       subtitle = sprintf("Occurrences detect %.0f%% of mapped species",
                          100 * detect$detection_ratio)) +
  theme_pub() + theme(plot.subtitle = element_text(size = 5, colour = "grey30"))

RESP_LAB <- c(body_mass.mean = "Body mass, mean",
              body_mass.var = "Body mass, variance",
              habitat_breadth.mean = "Habitat breadth, mean",
              habitat_breadth.var = "Habitat breadth, variance",
              range_size.mean = "Range size, mean")
hum$response_f <- factor(RESP_LAB[hum$response], levels = rev(unname(RESP_LAB)))
pc <- ggplot(hum, aes(estimate, response_f, colour = source)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                  position = position_dodge(width = 0.6), size = 0.22,
                  linewidth = 0.35) +
  scale_colour_manual(values = c(`Range map` = "grey55",
                                 `Occurrence` = "#B5541F",
                                 `Occurrence + coverage control` = "#2E6E8E"),
                      name = NULL) +
  labs(x = "Standardised effect of the human-pressure axis", y = NULL, tag = "c",
       subtitle = sprintf(paste("Range-map models use %d cells, occurrence models %d;",
                                "wider intervals reflect the smaller sample, not weaker effects"),
                          nrow(m_map), nrow(m_occ))) +
  theme_pub() + theme(legend.position = "bottom",
                      legend.text = element_text(size = 5),
                      axis.text.y = element_text(size = 5),
                      plot.subtitle = element_text(size = 4.6, colour = "grey30"))

p <- (pa | pb) / pc + plot_layout(heights = c(1.2, 1))
save_fig(p, "Fig12_occurrence_vs_rangemap", W2, 125, FIG)

log_msg("=== 13 完成 / done ===")
