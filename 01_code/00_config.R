# ============================================================
# 项目 / Project:
#   中国陆生脊椎动物群落性状的空间格局与驱动机制
#   Spatial patterns and drivers of community-level trait structure
#   in Chinese terrestrial vertebrates
#
# 本文件 / This file:
#   全局配置：路径、坐标系、常量、通用函数
#   Global configuration: paths, CRS, constants, shared helpers
#
# 科学问题 / Scientific question:
#   恒温类群（鸟、兽）与变温类群（爬行、两栖）的群落性状组成
#   是否受到不同环境过滤过程的驱动？该差异是否具有尺度依赖性？
#   Do endothermic (birds, mammals) and ectothermic (reptiles, amphibians)
#   assemblages differ systematically in how climate, productivity and
#   human land use filter their trait composition, and is this scale-dependent?
#
# 关键假设 / Key assumptions:
#   1. 等积网格（Albers）可保证网格面积可比，避免纬度面积偏倚。
#      Equal-area grid removes latitudinal area bias.
#   2. 采用官方审图号底图（GS(2023)2767），符合中国出版制图规范。
#      Official approved basemap for Chinese publication compliance.
#
# 主要包 / Main packages: sf, terra, dplyr, data.table
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(data.table)
})

# ---------------------------------------------------------------
# 1. 路径 / Paths
# ---------------------------------------------------------------
PROJ_ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"

PATH <- list(
  root       = PROJ_ROOT,
  code       = file.path(PROJ_ROOT, "01_code"),
  raw        = file.path(PROJ_ROOT, "02_data_raw"),
  raw_trait  = file.path(PROJ_ROOT, "02_data_raw", "traits"),
  raw_range  = file.path(PROJ_ROOT, "02_data_raw", "ranges"),
  raw_phylo  = file.path(PROJ_ROOT, "02_data_raw", "phylo"),
  raw_env    = file.path(PROJ_ROOT, "02_data_raw", "env"),
  raw_bound  = file.path(PROJ_ROOT, "02_data_raw", "boundary"),
  raw_check  = file.path(PROJ_ROOT, "02_data_raw", "checklists"),
  derived    = file.path(PROJ_ROOT, "03_data_derived"),
  results    = file.path(PROJ_ROOT, "04_results"),
  tables     = file.path(PROJ_ROOT, "04_results", "tables"),
  models     = file.path(PROJ_ROOT, "04_results", "models"),
  cache      = file.path(PROJ_ROOT, "04_results", "cache"),
  figures    = file.path(PROJ_ROOT, "05_figures"),
  ms         = file.path(PROJ_ROOT, "06_manuscript"),
  logs       = file.path(PROJ_ROOT, "07_logs")
)
invisible(lapply(PATH, function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)))

# 外部既有数据（复用用户本地已下载资源）/ Reuse of user's existing local data
EXT <- list(
  # 官方审图号底图 / Official approved Chinese basemap
  basemap_dir   = "/Users/dingchenchen/Documents/SDMs/GS(2023)2767审图号",
  # WorldClim 2.1 10 arc-min 生物气候变量 / bioclimatic variables
  worldclim_dir = "/Users/dingchenchen/Documents/New project/bird_full_community_analysis/data/external/climate/wc2.1_10m",
  # 人类足迹指数 / Human footprint index
  hfp_dir       = "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/data/external/hfi",
  # 中国土地覆被数据 CLCD (Yang & Huang 2021) / China land cover dataset
  clcd_dir      = "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis/data/external/clcd",
  # ESA WorldCover 类别占比 / land-cover fractions
  worldcover_dir= "/Users/dingchenchen/Documents/New project/bird_grid_community_analysis/data/external/landcover/landuse",
  # EarthEnv 生境异质性纹理 / habitat heterogeneity texture
  earthenv_dir  = "/Users/dingchenchen/Documents/New project/bird_grid_community_analysis/data/external/earthenv_texture",
  # 净初级生产力 / net primary productivity
  npp_tif       = "/Users/dingchenchen/Documents/New project/bird_grid_community_analysis/data/external/npp/npp_mean_china.tif",
  # 国家重点保护物种分布点 / national key protected species occurrence points
  cn_points_shp = "/Users/dingchenchen/Documents/SDMs/物种分布/物种分布.shp"
)

# ---------------------------------------------------------------
# 2. 坐标系与网格 / CRS and grid definition
# ---------------------------------------------------------------
# 中国 Albers 等积圆锥投影（与官方审图号底图一致）
# China Albers Equal Area Conic, identical to the official basemap
CRS_ALBERS <- "+proj=aea +lat_0=0 +lon_0=110 +lat_1=25 +lat_2=47 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
CRS_WGS84  <- "EPSG:4326"

# 两个分析粒度 / two analysis grains
GRAIN <- c(`50km` = 50000, `100km` = 100000)

# 网格纳入分析的最小陆地覆盖比例（剔除边缘破碎网格）
# Minimum land-cover fraction for a grid cell to enter the analysis
MIN_LAND_FRAC <- 0.5

# ---------------------------------------------------------------
# 3. 类群定义 / Taxon definitions
# ---------------------------------------------------------------
TAXA <- data.frame(
  class      = c("Aves", "Mammalia", "Reptilia", "Amphibia"),
  class_cn   = c("鸟类", "哺乳类", "爬行类", "两栖类"),
  thermal    = c("Endotherm", "Endotherm", "Ectotherm", "Ectotherm"),
  thermal_cn = c("恒温", "恒温", "变温", "变温"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------
# 4. 可复现性 / Reproducibility
# ---------------------------------------------------------------
SEED <- 20260801
set.seed(SEED)
N_NULL <- 999L   # 零模型随机化次数 / null-model randomisations
N_BOOT <- 999L   # bootstrap 次数 / bootstrap replicates

# ---------------------------------------------------------------
# 5. 通用函数 / Shared helpers
# ---------------------------------------------------------------

#' 带时间戳的日志输出 / Timestamped console log
log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  utils::flush.console()
}

#' 缓存化求值：若缓存文件存在则读取，否则计算并写入
#' Memoised evaluation: read cache if present, else compute and write
cache_rds <- function(name, expr, refresh = FALSE) {
  f <- file.path(PATH$cache, paste0(name, ".rds"))
  if (file.exists(f) && !refresh) {
    log_msg("  [cache hit] ", name)
    return(readRDS(f))
  }
  log_msg("  [computing] ", name)
  val <- force(expr)
  saveRDS(val, f, compress = "xz")
  val
}

#' 学名标准化：去首尾空格、统一大小写、去亚种、下划线转空格
#' Normalise scientific names: trim, capitalise, drop subspecies/underscores
norm_name <- function(x) {
  x <- as.character(x)
  # 清理非法/非 ASCII 编码（部分性状库含 latin1 与不间断空格）
  # Sanitise invalid / non-ASCII encodings (some trait databases carry
  # latin1 bytes and non-breaking spaces)
  x <- iconv(x, from = "UTF-8", to = "ASCII", sub = " ")
  x[is.na(x)] <- NA_character_
  x <- gsub("[^A-Za-z .-]", " ", x)
  x <- gsub("_", " ", x)
  x <- gsub("[[:space:]]+", " ", trimws(x))
  x <- sub("^([A-Za-z]+) ([a-z-]+).*$", "\\1 \\2", x)   # 属 + 种加词 / genus + epithet
  out <- paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
  out[is.na(x) | !nzchar(x)] <- NA_character_
  out
}

#' 安全 z 标准化（常数向量返回 0）/ Safe z-standardisation
zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' 可用于分析的网格编号（排除南海附图框内的伪网格）
#' Cell ids usable for analysis, excluding the South China Sea inset artefacts
#'
#' 官方审图号底图把南海诸岛预置在附图框内，那里生成的网格反投影到西太平洋，
#' 会把正确的物种组成与错误位置的气候配对。建网格时已用 scs_inset 标记，
#' 此处提供统一的过滤入口，避免各脚本各自实现而遗漏。
#' The approved basemap pre-places the South China Sea islands inside the inset
#' frame; cells there pair correct species lists with the wrong climate. One
#' shared accessor keeps every script consistent.
analysable_cells <- function(grain = "50km") {
  f <- file.path(PATH$derived, paste0("grid_", grain, "_attr.rds"))
  if (!file.exists(f)) return(NULL)
  a <- readRDS(f)
  if (!"scs_inset" %in% names(a)) return(a$cell_id)
  a$cell_id[!a$scs_inset]
}

#' 按可用网格过滤任意含 cell_id 的数据框
#' Filter any data frame with a cell_id column to the analysable cells
drop_inset <- function(d, grain = "50km") {
  ok <- analysable_cells(grain)
  if (is.null(ok) || !"cell_id" %in% names(d)) return(d)
  d[d$cell_id %in% ok, , drop = FALSE]
}

#' 写出表格并记录 / Write a results table with logging
write_table <- function(df, name) {
  f <- file.path(PATH$tables, paste0(name, ".csv"))
  data.table::fwrite(df, f)
  log_msg("  -> wrote ", basename(f), " (", nrow(df), " rows)")
  invisible(f)
}

log_msg("配置已加载 / config loaded — root: ", PROJ_ROOT)
