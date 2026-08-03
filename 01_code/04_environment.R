# ============================================================
# 科学问题 / Scientific question:
#   群落性状组成受哪一类环境过滤器支配——热能量、水分、生产力、
#   生境结构，还是当代人类压力？各类群是否不同？
#   Which environmental filter governs assemblage trait composition —
#   thermal energy, water, productivity, habitat structure, or
#   contemporary human pressure — and does this differ among classes?
#
# 分析目标 / Objective:
#   为每个网格提取五组预测变量（热、水、生产力、生境结构/地形、人为压力），
#   包含格内均值与格内异质性（标准差/极差）。
#
# 输入数据 / Input:
#   WorldClim 2.1 (2.5', 30" DEM)、MODIS NPP、ESA WorldCover、
#   EarthEnv 生境纹理、Human Footprint、CLCD 土地利用
#
# 主要流程 / Workflow:
#   1. 将网格投影到各栅格自身坐标系（比重投影大栅格快得多）
#   2. exactextractr 分区统计：均值、标准差、极差
#   3. 合并为网格 × 变量表，输出并做相关性/VIF 诊断
#
# 预期输出 / Expected output:
#   03_data_derived/env_<grain>.rds
#   04_results/tables/table_S4_env_summary.csv
#
# 关键假设 / Key assumptions:
#   - 2.5' (~4.6 km) 气候栅格在 50 km 网格内提供约 120 个像元，
#     足以刻画格内气候异质性。
#   - 人类足迹 (HFP) 综合了人口、建设、道路、耕地等多重压力，
#     作为当代人为过滤强度的总体代理。
#
# 主要包 / Main packages: terra, exactextractr, sf, dplyr
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({ library(exactextractr) })
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0)

log_msg("=== 04 环境变量提取 / Environmental predictor extraction ===")

grids <- lapply(names(GRAIN), function(lab)
  st_read(file.path(PATH$derived, paste0("grid_", lab, ".gpkg")), quiet = TRUE))
names(grids) <- names(GRAIN)

#' 分区统计辅助函数 / zonal statistics helper
#' @param r SpatRaster (单层或多层 / single or multi-layer)
#' @param g sf 网格 / grid
#' @param stats 需要的统计量 / statistics required
zonal_stats <- function(r, g, stats = "mean", prefix = "") {
  gp <- st_transform(g["cell_id"], terra::crs(r))
  v  <- exactextractr::exact_extract(r, gp, fun = stats, progress = FALSE)
  v  <- as.data.frame(v)
  if (ncol(v) == 1 && nzchar(prefix)) names(v) <- prefix
  else if (nzchar(prefix)) names(v) <- paste0(prefix, "_", names(v))
  cbind(cell_id = g$cell_id, v)
}

extract_all <- function(lab) {
  g <- grids[[lab]]
  log_msg("-- 粒度 ", lab, " (", nrow(g), " 格) --")
  out <- data.frame(cell_id = g$cell_id)

  ## ---- A. 气候：WorldClim 2.1, 2.5 arc-min --------------------------------
  wcdir <- file.path(PATH$raw, "env", "wc2.1_2.5m")
  bios  <- file.path(wcdir, paste0("wc2.1_2.5m_bio_", 1:19, ".tif"))
  bios  <- bios[file.exists(bios)]
  if (length(bios)) {
    r <- terra::rast(bios)
    # 从 "wc2.1_2.5m_bio_12.tif" 提取编号 12 / extract the index after "bio_"
    names(r) <- paste0("bio", sub("^.*_bio_(\\d+)\\.tif$", "\\1", basename(bios)))
    m <- zonal_stats(r, g, "mean")
    names(m)[-1] <- paste0(names(r), "_mean")
    out <- left_join(out, m, by = "cell_id")
    # 格内气候异质性（温度与降水）/ within-cell climatic heterogeneity
    rh <- r[[c("bio1", "bio12")]]
    s <- zonal_stats(rh, g, "stdev")
    names(s)[-1] <- c("bio1_sd", "bio12_sd")
    out <- left_join(out, s, by = "cell_id")
    log_msg("   + WorldClim bio1-19 (", length(bios), " 层)")
  }

  ## ---- B. 地形：WorldClim 30" DEM -----------------------------------------
  dem <- file.path(PATH$raw, "env", "wc2.1_30s", "wc2.1_30s_elev.tif")
  if (file.exists(dem)) {
    r <- terra::rast(dem)
    gp <- st_transform(g["cell_id"], terra::crs(r))
    v <- exactextractr::exact_extract(r, gp, fun = c("mean", "stdev", "min", "max"),
                                      progress = FALSE)
    out <- left_join(out, data.frame(cell_id = g$cell_id,
                                     elev_mean  = v$mean,
                                     elev_sd    = v$stdev,
                                     elev_range = v$max - v$min),
                     by = "cell_id")
    log_msg("   + 地形 / topography (elev mean, sd, range)")
  }

  ## ---- C. 生产力：MODIS NPP ------------------------------------------------
  if (file.exists(EXT$npp_tif)) {
    r <- terra::rast(EXT$npp_tif)
    out <- left_join(out, zonal_stats(r, g, "mean", "npp"), by = "cell_id")
    log_msg("   + NPP")
  }

  ## ---- D. 土地覆被占比：ESA WorldCover ------------------------------------
  wc <- list.files(EXT$worldcover_dir, pattern = "tif$", full.names = TRUE)
  if (length(wc)) {
    for (f in wc) {
      nm <- tolower(gsub("WorldCover_|_30s\\.tif", "", basename(f)))
      r <- terra::rast(f)
      out <- left_join(out, zonal_stats(r, g, "mean", paste0("lc_", nm)), by = "cell_id")
    }
    log_msg("   + WorldCover 占比 (", length(wc), " 类)")
  }

  ## ---- E. 生境异质性：EarthEnv 纹理 ---------------------------------------
  ee <- list.files(EXT$earthenv_dir, pattern = "tif$", full.names = TRUE)
  if (length(ee)) {
    for (f in ee) {
      nm <- tolower(sub("_.*$", "", basename(f)))
      r <- terra::rast(f)
      out <- left_join(out, zonal_stats(r, g, "mean", paste0("het_", nm)), by = "cell_id")
    }
    log_msg("   + EarthEnv 生境纹理 (", length(ee), " 指标)")
  }

  ## ---- F. 人为压力：Human Footprint ---------------------------------------
  for (yr in c(2000, 2020)) {
    f <- file.path(EXT$hfp_dir, paste0("hfp", yr, ".tif"))
    if (!file.exists(f)) next
    r <- terra::rast(f)
    out <- left_join(out, zonal_stats(r, g, "mean", paste0("hfp", yr)), by = "cell_id")
  }
  if (all(c("hfp2000", "hfp2020") %in% names(out)))
    out$hfp_change <- out$hfp2020 - out$hfp2000
  log_msg("   + 人类足迹 HFP 2000/2020 及变化")

  ## ---- G. 网格几何属性 / grid geometry ------------------------------------
  out <- left_join(out, st_drop_geometry(g)[, c("cell_id", "lon", "lat",
                                                "x_albers", "y_albers",
                                                "land_area_km2", "prov_cn")],
                   by = "cell_id")
  out
}

env <- list()
for (lab in names(grids)) {
  env[[lab]] <- cache_rds(paste0("env_", lab), extract_all(lab))
  saveRDS(env[[lab]], file.path(PATH$derived, paste0("env_", lab, ".rds")))
  log_msg("   -> env_", lab, ".rds  (", nrow(env[[lab]]), " 格 x ",
          ncol(env[[lab]]) - 1, " 变量)")
}

# ---------------------------------------------------------------
# 汇总与诊断 / Summary and diagnostics
# ---------------------------------------------------------------
summ <- lapply(names(env), function(lab) {
  d <- env[[lab]]
  num <- d[, sapply(d, is.numeric), drop = FALSE]
  data.frame(grain = lab, variable = names(num),
             n_valid = sapply(num, function(x) sum(is.finite(x))),
             mean = round(sapply(num, mean, na.rm = TRUE), 3),
             sd   = round(sapply(num, sd, na.rm = TRUE), 3),
             min  = round(sapply(num, min, na.rm = TRUE), 3),
             max  = round(sapply(num, max, na.rm = TRUE), 3),
             row.names = NULL)
}) |> bind_rows()
write_table(summ, "table_S4_env_summary")
print(head(summ, 30))

log_msg("=== 04 完成 / done ===")
