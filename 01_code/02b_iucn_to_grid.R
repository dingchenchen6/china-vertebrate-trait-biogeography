# ============================================================
# 科学问题 / Scientific question:
#   以权威专家范围图定义"潜在群落"，用于跨类群群落性状比较。
#   Define "potential assemblages" from expert range maps for
#   cross-taxon comparison of assemblage trait structure.
#
# 分析目标 / Objective:
#   将 IUCN 红色名录空间数据（两栖、爬行、哺乳）按中国范围裁剪、
#   按标准准则筛选，并叠加到 50 km / 100 km 等积网格。
#
# 输入数据 / Input:
#   02_data_raw/ranges/iucn/{AMPHIBIANS,REPTILES,MAMMALS}_PART{1,2}.shp
#
# 主要流程 / Workflow:
#   1. 用中国外接矩形做空间过滤读取（避免载入全球多边形）
#   2. 筛选：presence 1-2；origin 1,2,6；seasonal 1,2；陆生
#   3. 栅格化到 5 km 亚网格，再按 50/100 km 网格统计覆盖比例
#      （比多边形求交快 1-2 个数量级，精度 5 km 对 50 km 粒度足够）
#   4. 覆盖比例 >= 1% 判定物种存在，输出长表
#
# 预期输出 / Expected output:
#   03_data_derived/comm_<class>_<grain>.rds
#
# 关键假设 / Key assumptions:
#   - presence 1-2 = 现存/可能现存；origin 1,2,6 = 本土/重引入/辅助定殖
#   - seasonal 1-2 = 留居/繁殖季，排除仅越冬与迁徙过境
#   - 5 km 亚网格离散化误差在 50 km 粒度下可忽略
#
# 主要包 / Main packages: sf, terra, dplyr, data.table
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0)

log_msg("=== 02b IUCN 范围图 -> 网格（栅格化法）/ IUCN ranges -> grid (raster method) ===")

IUCN_DIR      <- file.path(PATH$raw, "ranges", "iucn")
CHINA_WKT     <- "POLYGON((73 17,136 17,136 54,73 54,73 17))"
SUBRES        <- 5000      # 亚网格分辨率 (m) / sub-grid resolution
MIN_CELL_FRAC <- 0.01      # 网格内范围覆盖比例阈值 / min range coverage per cell

grids <- lapply(names(GRAIN), function(lab)
  st_read(file.path(PATH$derived, paste0("grid_", lab, ".gpkg")), quiet = TRUE))
names(grids) <- names(GRAIN)

# ---------------------------------------------------------------
# 0. 构建 5 km 亚网格模板，并预先建立"亚网格 -> 各粒度网格"的查找表
#    Build a 5 km sub-grid template and a lookup from sub-cells to grid cells
# ---------------------------------------------------------------
# 注意：SpatRaster 的 C++ 指针无法序列化，缓存后会失效；
# 因此模板每次由参数重建，缓存中只存纯 R 向量。
# NOTE: a SpatRaster's C++ pointer cannot be serialised, so the template is
# rebuilt from parameters each run and only plain vectors are cached.
make_template <- function() {
  bb <- st_bbox(grids[["100km"]])
  terra::rast(xmin = bb["xmin"], xmax = bb["xmax"],
              ymin = bb["ymin"], ymax = bb["ymax"],
              resolution = SUBRES, crs = CRS_ALBERS)
}
TMPL <- make_template()

build_lookup <- function() {
  lut <- list()
  for (lab in names(grids)) {
    g <- grids[[lab]]
    g$.idx <- seq_len(nrow(g))
    r <- terra::rasterize(terra::vect(g), TMPL, field = ".idx")
    lut[[lab]] <- list(idx = terra::values(r)[, 1], ids = g$cell_id,
                       n_sub = (GRAIN[[lab]] / SUBRES)^2)
    log_msg("  查找表 ", lab, ": ", sum(!is.na(lut[[lab]]$idx)), " 个 5 km 亚网格")
  }
  lut
}
LUT <- cache_rds("subgrid_lookup_5km_v2", build_lookup())

# ---------------------------------------------------------------
# 1. 读取并筛选某一类群的中国范围多边形
# ---------------------------------------------------------------
read_iucn_china <- function(prefix) {
  fs <- list.files(IUCN_DIR, pattern = paste0("^", prefix, "_PART.*shp$"), full.names = TRUE)
  stopifnot(length(fs) > 0)
  out <- lapply(fs, function(f) {
    lyr <- tools::file_path_sans_ext(basename(f))
    q <- paste0("SELECT id_no, sci_name, presence, origin, seasonal, category,",
                " marine, terrestria, freshwater, order_, family, genus FROM \"", lyr, "\"",
                " WHERE presence IN (1,2) AND origin IN (1,2,6) AND seasonal IN (1,2)")
    x <- st_read(f, query = q, wkt_filter = CHINA_WKT, quiet = TRUE)
    log_msg("   ", basename(f), ": ", nrow(x), " polygons in China bbox")
    x
  })
  x <- do.call(rbind, out)

  ter <- tolower(as.character(x$terrestria)); mar <- tolower(as.character(x$marine))
  is_ter <- ter %in% c("true", "t", "yes", "1")
  is_mar <- mar %in% c("true", "t", "yes", "1")
  x <- x[is_ter | (!is_mar & is.na(ter)), ]          # 保留陆生，剔除纯海洋 / terrestrial only

  x$species <- norm_name(x$sci_name)
  x <- x[!is.na(x$species) & nzchar(x$species), ]
  log_msg("   -> 筛选后 / after filters: ", nrow(x), " polygons, ",
          length(unique(x$species)), " species")
  x
}

# ---------------------------------------------------------------
# 2. 栅格化：单物种范围 -> 各粒度网格的覆盖比例
#    Rasterise one species' range -> coverage fraction per grid cell
# ---------------------------------------------------------------
species_to_cells <- function(v_sp) {
  r <- terra::rasterize(v_sp, TMPL, field = 1, background = NA)
  hit <- which(!is.na(terra::values(r)[, 1]))
  if (!length(hit)) return(NULL)
  res <- list()
  for (lab in names(grids)) {
    gi <- LUT[[lab]]$idx[hit]
    gi <- gi[!is.na(gi)]
    if (!length(gi)) next
    tb <- table(gi)
    frac <- as.numeric(tb) / LUT[[lab]]$n_sub
    keep <- frac >= MIN_CELL_FRAC
    if (!any(keep)) next
    res[[lab]] <- data.frame(cell_id = LUT[[lab]]$ids[as.integer(names(tb))[keep]],
                             frac = frac[keep])
  }
  res
}

# 内存安全策略 / Memory-safe strategy:
#   先把「投影 + 简化」后的范围写入磁盘 GPKG，再逐物种用 SQL 查询读取。
#   IUCN 多边形顶点极多，整体载入会耗尽内存；2 km 简化容差在 50 km
#   粒度下的离散化误差可忽略。
#   Project + simplify once, write to an on-disk GPKG, then stream one
#   species at a time. A 2 km simplification tolerance is negligible at
#   a 50 km grain but cuts memory an order of magnitude.
SIMPLIFY_TOL <- 2000

prepare_ranges_gpkg <- function(polys, cl) {
  f <- file.path(PATH$derived, paste0("ranges_albers_", cl, ".gpkg"))
  if (file.exists(f)) { log_msg("   [已存在] ", basename(f)); return(f) }
  log_msg("   投影 + 简化 (", SIMPLIFY_TOL/1000, " km 容差) 并写入磁盘 ...")
  p <- st_transform(polys[, "species"], CRS_ALBERS)
  p <- st_simplify(p, dTolerance = SIMPLIFY_TOL, preserveTopology = TRUE)
  p <- st_make_valid(p)
  p <- p[!st_is_empty(p), ]
  st_write(p, f, delete_dsn = TRUE, quiet = TRUE)
  rm(p); gc()
  f
}

ranges_to_grid_all <- function(polys, cl) {
  gpkg <- prepare_ranges_gpkg(polys, cl)
  lyr  <- tools::file_path_sans_ext(basename(gpkg))
  rm(polys); gc()

  sps <- st_read(gpkg, query = paste0("SELECT DISTINCT species FROM \"", lyr, "\""),
                 quiet = TRUE)$species
  sps <- sort(sps[!is.na(sps)])
  log_msg("   栅格化 ", length(sps), " 个物种 / rasterising ", length(sps), " species")

  acc <- setNames(vector("list", length(grids)), names(grids))
  t0 <- Sys.time()
  for (i in seq_along(sps)) {
    sp <- sps[i]
    q  <- paste0("SELECT geom FROM \"", lyr, "\" WHERE species = '",
                 gsub("'", "''", sp), "'")
    gsp <- try(st_read(gpkg, query = q, quiet = TRUE), silent = TRUE)
    if (inherits(gsp, "try-error") || !nrow(gsp)) next
    res <- species_to_cells(terra::vect(gsp))
    rm(gsp)
    if (!is.null(res)) for (lab in names(res)) {
      res[[lab]]$species <- sp
      acc[[lab]][[length(acc[[lab]]) + 1]] <- res[[lab]]
    }
    if (i %% 250 == 0) {
      gc(verbose = FALSE)
      log_msg("     ", i, "/", length(sps), "  (",
              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
    }
  }
  lapply(acc, function(x) as.data.frame(data.table::rbindlist(x)))
}

# ---------------------------------------------------------------
# 3. 执行 / Run
# ---------------------------------------------------------------
CLASS_MAP <- c(AMPHIBIANS = "Amphibia", REPTILES = "Reptilia", MAMMALS = "Mammalia")

summ <- list()
for (pref in names(CLASS_MAP)) {
  cl <- CLASS_MAP[[pref]]
  log_msg("-- ", cl, " --")

  polys <- cache_rds(paste0("iucn_china_", cl), read_iucn_china(pref))

  splist <- st_drop_geometry(polys) |>
    distinct(species, .keep_all = TRUE) |>
    transmute(class = cl, species, order = order_, family, genus,
              iucn_category = category, iucn_id = id_no)
  saveRDS(splist, file.path(PATH$derived, paste0("splist_", cl, ".rds")))

  gridded <- cache_rds(paste0("comm_", cl, "_gridded_raster"), ranges_to_grid_all(polys, cl))

  for (lab in names(grids)) {
    d <- gridded[[lab]]
    if (is.null(d) || !nrow(d)) next
    d$class <- cl
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
    log_msg("   ", lab, ": ", length(unique(d$species)), " spp, mean richness ",
            round(mean(rich$richness), 1), ", max ", max(rich$richness))
  }
}

summ <- do.call(rbind, summ)
write_table(summ, "table_S2b_iucn_assemblage_summary")
print(summ)
log_msg("=== 02b 完成 / done ===")
