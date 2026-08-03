# ============================================================
# 科学问题 / Scientific question:
#   建立可比的空间分析单元，以检验群落性状格局的尺度依赖性。
#   Establish comparable spatial analysis units to test scale-dependence
#   of assemblage trait patterns.
#
# 分析目标 / Objective:
#   基于中国官方审图号底图，构建 50 km 与 100 km 等积网格，
#   计算陆地覆盖比例、省份归属、经纬度与地理坐标，输出标准网格。
#   Build 50 km and 100 km equal-area grids over China from the official
#   approved basemap; compute land fraction, province, coordinates.
#
# 输入数据 / Input:
#   GS(2023)2767 审图号省面 shapefile（Albers 等积投影）
#
# 主要流程 / Workflow:
#   1. 读取并溶解省界 -> 中国陆域多边形
#   2. 生成规则方格网（50 km / 100 km）
#   3. 计算每格陆地面积比例，按阈值筛选
#   4. 附加省份、质心经纬度、到海岸距离
#   5. 导出 gpkg + rds
#
# 预期输出 / Expected output:
#   03_data_derived/grid_50km.gpkg, grid_100km.gpkg, china_land.gpkg
#
# 关键假设 / Key assumptions:
#   陆地覆盖 >= 50% 的网格才纳入分析，避免边缘网格的物种数被低估。
#   Cells with >=50% land are retained to avoid edge-cell richness artefacts.
#   落在南海诸岛附图框内的网格用 scs_inset 标记并在分析中排除（见正文注释）。
#   Cells inside the South China Sea inset frame are flagged and excluded.
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
sf::sf_use_s2(FALSE)

log_msg("=== 01 构建中国等积网格 / Building China equal-area grids ===")

# ---------------------------------------------------------------
# 1. 中国陆域多边形 / China land polygon
# ---------------------------------------------------------------
prov <- st_read(file.path(EXT$basemap_dir, "省面.shp"), quiet = TRUE) |>
  st_transform(CRS_ALBERS) |>
  st_make_valid()

# 统一省份字段名 / harmonise province field names
prov <- prov |>
  transmute(prov_code = as.character(DZM),
            prov_cn   = as.character(NAME),
            prov_en   = as.character(Yname),
            prov_full = as.character(`省全名`))

log_msg("省级单元 / provincial units: ", nrow(prov))

china_land <- prov |> st_union() |> st_make_valid()
land_area_km2 <- as.numeric(st_area(china_land)) / 1e6
log_msg("中国陆域面积 / China land area: ", format(round(land_area_km2), big.mark = ","), " km2")

st_write(st_sf(geom = china_land), file.path(PATH$derived, "china_land.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
st_write(prov, file.path(PATH$derived, "china_provinces.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

# 辅助制图要素（九/十段线、国界）/ auxiliary cartographic layers
for (nm in c("十段线", "国界", "省界线", "南海附图框")) {
  f <- file.path(EXT$basemap_dir, paste0(nm, ".shp"))
  if (file.exists(f)) {
    st_read(f, quiet = TRUE) |> st_transform(CRS_ALBERS) |>
      st_write(file.path(PATH$derived, paste0("basemap_", nm, ".gpkg")),
               delete_dsn = TRUE, quiet = TRUE)
  }
}
log_msg("底图图层已导出 / basemap layers exported")

# ---------------------------------------------------------------
# 2. 生成网格 / Generate grids
# ---------------------------------------------------------------
build_grid <- function(cellsize, label) {

  log_msg("-- 构建 ", label, " 网格 / building ", label, " grid --")

  # 对齐到 cellsize 的整数倍原点，保证 50 km 网格可嵌套进 100 km 网格
  # Snap origin to a multiple of cellsize so 50 km nests inside 100 km
  bb  <- st_bbox(china_land)
  off <- c(floor(bb["xmin"] / cellsize) * cellsize,
           floor(bb["ymin"] / cellsize) * cellsize)

  g <- st_make_grid(china_land, cellsize = cellsize, offset = off, what = "polygons")
  g <- st_sf(cell_id_raw = seq_along(g), geom = g, crs = st_crs(china_land))

  # 只保留与陆域相交的网格 / keep cells intersecting land
  hit <- lengths(st_intersects(g, china_land)) > 0
  g   <- g[hit, ]
  log_msg("   与陆域相交网格 / cells intersecting land: ", nrow(g))

  # 陆地面积比例 / land fraction
  inter <- st_intersection(g, china_land) |> st_make_valid()
  la <- data.frame(cell_id_raw = inter$cell_id_raw,
                   land_m2     = as.numeric(st_area(inter))) |>
    group_by(cell_id_raw) |> summarise(land_m2 = sum(land_m2), .groups = "drop")

  g <- g |>
    left_join(la, by = "cell_id_raw") |>
    mutate(land_m2    = ifelse(is.na(land_m2), 0, land_m2),
           cell_area_km2 = (cellsize^2) / 1e6,
           land_area_km2 = land_m2 / 1e6,
           land_frac  = land_m2 / (cellsize^2))

  g <- g |> filter(land_frac >= MIN_LAND_FRAC)
  log_msg("   陆地占比>=", MIN_LAND_FRAC, " 的网格 / cells retained: ", nrow(g))

  # 主导省份（按面积最大交集）/ dominant province by largest intersecting area
  pi_ <- st_intersection(g["cell_id_raw"], prov) |> st_make_valid()
  pi_$a <- as.numeric(st_area(pi_))
  dom <- st_drop_geometry(pi_) |>
    group_by(cell_id_raw) |> slice_max(a, n = 1, with_ties = FALSE) |> ungroup() |>
    select(cell_id_raw, prov_code, prov_cn, prov_en)
  g <- g |> left_join(dom, by = "cell_id_raw")

  # 质心坐标（投影坐标与经纬度）/ centroid in projected and geographic CRS
  ctr <- st_centroid(st_geometry(g))
  xy  <- st_coordinates(ctr)
  ll  <- st_coordinates(st_transform(ctr, CRS_WGS84))
  g$x_albers <- xy[, 1]; g$y_albers <- xy[, 2]
  g$lon <- ll[, 1];      g$lat <- ll[, 2]

  # 稳定的网格编号（自西南向东北）/ stable cell id ordered SW -> NE
  ord <- order(g$y_albers, g$x_albers)
  g <- g[ord, ]
  g$cell_id <- sprintf("%s_%05d", label, seq_len(nrow(g)))

  # 官方审图号底图把南海诸岛**预置在附图框内**，而非其真实地理位置。
  # 这些位置生成的网格，其反投影经纬度落在西太平洋（129-136 E, 25-27 N），
  # 于是「正确的物种组成」被配上「错误位置的气候」，必须排除。
  # 这里用标记而非直接删除：删除会让后续所有 cell_id 重新编号，
  # 使已算好的群落矩阵、环境表、指标表全部失效。
  # The approved basemap pre-places the South China Sea islands inside the
  # inset frame; cells there back-project to the western Pacific and pair
  # correct species lists with the wrong climate. They are flagged rather than
  # deleted, because deleting would renumber every downstream cell_id.
  g$scs_inset <- FALSE
  fr <- file.path(PATH$derived, "basemap_南海附图框.gpkg")
  if (file.exists(fr)) {
    box <- st_sf(geom = st_as_sfc(st_bbox(st_read(fr, quiet = TRUE))))
    g$scs_inset <- lengths(st_intersects(st_centroid(st_geometry(g)), box)) > 0
    log_msg("   标记为南海附图框内的伪网格 / cells flagged as inset artefacts: ",
            sum(g$scs_inset))
  }

  g <- g |> select(cell_id, cell_id_raw, prov_code, prov_cn, prov_en,
                   cell_area_km2, land_area_km2, land_frac, scs_inset,
                   x_albers, y_albers, lon, lat, geom)

  st_write(g, file.path(PATH$derived, paste0("grid_", label, ".gpkg")),
           delete_dsn = TRUE, quiet = TRUE)
  saveRDS(st_drop_geometry(g), file.path(PATH$derived, paste0("grid_", label, "_attr.rds")))
  log_msg("   已保存 / saved grid_", label, ".gpkg")
  g
}

grids <- lapply(names(GRAIN), function(lab) build_grid(GRAIN[[lab]], lab))
names(grids) <- names(GRAIN)

# ---------------------------------------------------------------
# 3. 网格嵌套关系（50 km -> 100 km）/ nesting lookup
# ---------------------------------------------------------------
# 用于多尺度分析中的层级随机效应 / used as a hierarchical random effect
nest <- st_join(st_centroid(grids[["50km"]]["cell_id"]),
                grids[["100km"]]["cell_id"],
                join = st_within, left = TRUE) |>
  st_drop_geometry() |>
  rename(cell_50km = cell_id.x, cell_100km = cell_id.y)
saveRDS(nest, file.path(PATH$derived, "grid_nesting_50_in_100.rds"))
log_msg("嵌套关系 / nesting: ", sum(!is.na(nest$cell_100km)), " of ",
        nrow(nest), " 50 km cells fall inside a retained 100 km cell")

# ---------------------------------------------------------------
# 4. 汇总 / Summary
# ---------------------------------------------------------------
summ <- do.call(rbind, lapply(names(grids), function(lab) {
  g <- grids[[lab]]
  data.frame(grain = lab,
             n_cells = nrow(g),
             n_analysable = sum(!g$scs_inset),
             cell_area_km2 = unique(g$cell_area_km2),
             total_land_km2 = round(sum(g$land_area_km2)),
             mean_land_frac = round(mean(g$land_frac), 3),
             n_provinces = length(unique(na.omit(g$prov_cn))))
}))
write_table(summ, "table_S1_grid_summary")
print(summ)

log_msg("=== 01 完成 / done ===")
