# ============================================================
# 科学问题 / Scientific question:
#   构建四个陆生脊椎动物类群在两个空间粒度上的群落组成矩阵，
#   作为后续群落性状分析的基础。
#   Build assemblage composition matrices for four terrestrial vertebrate
#   classes at two spatial grains.
#
# 分析目标 / Objective:
#   将异源分布数据（已网格化矩阵、范围多边形、点记录）统一重采样到
#   本项目的中国 Albers 等积 50 km / 100 km 网格。
#   Harmonise heterogeneous distribution data onto the project's
#   China Albers equal-area 50 km / 100 km grids.
#
# 输入数据 / Input:
#   - Huang et al. 2023 公开鸟类矩阵（1,180 种 × 3,888 Behrmann 50 km 格）
#   - 中国蜥蜴分布多边形（213 种，WGS84）
#   - 其他类群范围数据（见 02_data_raw/ranges）
#
# 主要流程 / Workflow:
#   1. 面积加权重网格化：Behrmann 格 -> Albers 格（>=50% 面积规则）
#   2. 多边形直接叠加网格（>=1% 网格面积或包含质心）
#   3. 输出稀疏 presence-absence 长表 + 矩阵
#
# 预期输出 / Expected output:
#   03_data_derived/comm_<class>_<grain>.rds  (稀疏长表 / sparse long table)
#
# 关键假设 / Key assumptions:
#   面积加权 >=50% 规则在重采样中保持网格物种丰富度近似无偏。
#   Area-weighted >=50% rule keeps grid richness approximately unbiased.
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
sf::sf_use_s2(FALSE)

log_msg("=== 02 构建群落组成矩阵 / Building assemblage matrices ===")

REPO <- file.path(PATH$raw, "huang2023_repo")

grids <- lapply(names(GRAIN), function(lab)
  st_read(file.path(PATH$derived, paste0("grid_", lab, ".gpkg")), quiet = TRUE))
names(grids) <- names(GRAIN)

# ---------------------------------------------------------------
# 工具函数 A：源网格 -> 目标网格 的面积加权重采样
# Helper A: area-weighted resampling from a source grid to target grid
# ---------------------------------------------------------------
regrid_matrix <- function(src_grid, src_id_col, pa_mat, target, min_frac = 0.5) {
  # pa_mat: 行 = 源网格 id（字符），列 = 物种 / rows = source cell id, cols = species
  src <- src_grid |> st_transform(st_crs(target)) |> st_make_valid()
  src$.sid <- as.character(src[[src_id_col]])

  inter <- st_intersection(target["cell_id"], src[".sid"]) |> st_make_valid()
  inter$a <- as.numeric(st_area(inter))
  w <- st_drop_geometry(inter) |>
    group_by(cell_id, .sid) |> summarise(a = sum(a), .groups = "drop")

  # 目标网格总面积用于归一化 / normalise by target cell area
  tgt_area <- data.frame(cell_id = target$cell_id,
                         tot = as.numeric(st_area(target)))
  w <- w |> left_join(tgt_area, by = "cell_id") |> mutate(frac = a / tot)

  # 只保留在源矩阵中存在的源格 / keep source cells present in the matrix
  w <- w |> filter(.sid %in% rownames(pa_mat))
  if (!nrow(w)) return(NULL)

  # 对每个物种：累加其出现的源格在目标格中的面积占比，>= min_frac 判定为存在
  # For each species: sum area fractions of occupied source cells; present if >= min_frac
  sp <- colnames(pa_mat)
  idx <- match(w$.sid, rownames(pa_mat))
  out <- vector("list", length(sp))
  for (j in seq_along(sp)) {
    occ <- pa_mat[idx, j] > 0
    if (!any(occ)) next
    agg <- tapply(w$frac[occ], w$cell_id[occ], sum)
    keep <- names(agg)[agg >= min_frac]
    if (length(keep)) out[[j]] <- data.frame(cell_id = keep, species = sp[j])
  }
  res <- data.table::rbindlist(out)
  as.data.frame(res)
}

# ---------------------------------------------------------------
# 工具函数 B：范围多边形 -> 网格
# Helper B: range polygons -> grid
# ---------------------------------------------------------------
polys_to_grid <- function(polys, sp_col, target, min_frac = 0.05) {
  p <- polys |> st_transform(st_crs(target)) |> st_make_valid()
  p$.sp <- norm_name(p[[sp_col]])
  p <- p[!is.na(p$.sp) & p$.sp != "", ]

  inter <- st_intersection(target["cell_id"], p[".sp"]) |> st_make_valid()
  inter$a <- as.numeric(st_area(inter))
  tgt_area <- data.frame(cell_id = target$cell_id, tot = as.numeric(st_area(target)))

  st_drop_geometry(inter) |>
    group_by(cell_id, .sp) |> summarise(a = sum(a), .groups = "drop") |>
    left_join(tgt_area, by = "cell_id") |>
    filter(a / tot >= min_frac) |>
    transmute(cell_id, species = .sp) |>
    as.data.frame()
}

# ---------------------------------------------------------------
# 1. 鸟类：Huang et al. 2023 的 Behrmann 50 km 矩阵
#    Birds: Huang et al. 2023 Behrmann 50 km matrix
# ---------------------------------------------------------------
build_birds <- function() {
  e <- new.env()
  load(file.path(REPO, "ave", "parafly", "comm_data.RData"), envir = e)
  pa <- e$Chinese_species_distribution_grid_selected
  colnames(pa) <- norm_name(colnames(pa))
  log_msg("鸟类源矩阵 / bird source matrix: ", nrow(pa), " cells x ", ncol(pa), " species")

  src <- st_read(file.path(REPO, "chinese_terrestral", "input", "grid.shp"), quiet = TRUE)

  lapply(names(grids), function(lab) {
    log_msg("  -> 重网格化鸟类到 ", lab)
    r <- regrid_matrix(src, "id", pa, grids[[lab]], min_frac = 0.5)
    r$class <- "Aves"; r
  }) |> setNames(names(grids))
}

# ---------------------------------------------------------------
# 2. 爬行（蜥蜴）：中国蜥蜴分布多边形
#    Reptiles (lizards): Chinese lizard range polygons
# ---------------------------------------------------------------
build_lizards <- function() {
  f <- file.path(REPO, "reptile", "data_reptile", "Chinese213lizards_9_27",
                 "Chinese_lizards_214_version20200927.shp")
  p <- st_read(f, quiet = TRUE)
  log_msg("蜥蜴范围多边形 / lizard range polygons: ", nrow(p), " species")
  lapply(names(grids), function(lab) {
    log_msg("  -> 网格化蜥蜴到 ", lab)
    r <- polys_to_grid(p, "Binomial", grids[[lab]], min_frac = 0.05)
    r$class <- "Reptilia"; r
  }) |> setNames(names(grids))
}

# ---------------------------------------------------------------
# 3. 执行 / Run
# ---------------------------------------------------------------
comm <- list()
comm$Aves     <- cache_rds("comm_birds_regridded",  build_birds())
comm$Reptilia <- cache_rds("comm_lizards_gridded",  build_lizards())

# 汇总并落盘 / summarise and save
summ <- list()
for (cl in names(comm)) {
  for (lab in names(grids)) {
    d <- comm[[cl]][[lab]]
    if (is.null(d) || !nrow(d)) next
    saveRDS(d, file.path(PATH$derived, sprintf("comm_%s_%s.rds", cl, lab)))
    rich <- d |> count(cell_id, name = "richness")
    summ[[length(summ) + 1]] <- data.frame(
      class = cl, grain = lab,
      n_species = length(unique(d$species)),
      n_cells_occupied = nrow(rich),
      n_cells_total = nrow(grids[[lab]]),
      mean_richness = round(mean(rich$richness), 1),
      max_richness = max(rich$richness),
      n_records = nrow(d))
  }
}
summ <- do.call(rbind, summ)
write_table(summ, "table_S2_assemblage_summary")
print(summ)

log_msg("=== 02 完成（鸟类、爬行）；哺乳与两栖待分布数据就绪 ===")
log_msg("=== 02 done (birds, reptiles); mammals & amphibians pending range data ===")
