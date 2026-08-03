# ============================================================
# 科学问题 / Scientific question:
#   要在恒温与变温类群之间比较"性状–环境耦合"，必须先有一套
#   四类群定义一致、可直接比较的核心性状。
#   Comparing trait-environment coupling between endotherms and
#   ectotherms requires a core trait set defined identically across
#   all four classes.
#
# 分析目标 / Objective:
#   以 TetrapodTraits v2.0.1 为统一骨架，整合 AVONET / COMBINE /
#   ReptTraits / AmphiBIO / Etard / Meiri / GlobTherm / AnAge，
#   构建 (a) 四类群通用核心性状集 与 (b) 类群特异性状集，
#   并用系统发育信息辅助插补缺失值。
#
# 输入数据 / Input:
#   02_data_raw/traits/*  ；03_data_derived/splist_*.rds
#
# 主要流程 / Workflow:
#   1. 读取骨架，按中国物种名录取子集
#   2. 合并补充数据库（学名标准化后左连接）
#   3. 定义核心性状（8 个，四类群齐全）与扩展性状
#   4. 报告插补前覆盖率
#   5. missForest + 系统发育特征向量插补，报告 OOB 误差
#   6. 输出性状表与覆盖率表
#
# 预期输出 / Expected output:
#   03_data_derived/traits_core.rds, traits_imputed.rds
#   04_results/tables/table_S3_trait_coverage.csv
#
# 关键假设 / Key assumptions:
#   - 通用核心性状在四类群间语义一致（如"窝仔数"= 单次繁殖后代数）
#   - 插补仅用于 <40% 缺失的性状；更高缺失率的性状仅作敏感性分析
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({ library(tidyr) })

log_msg("=== 03 性状整合 / Trait integration ===")
TR <- file.path(PATH$raw, "traits")

# ---------------------------------------------------------------
# 1. 中国物种名录（来自分布数据）/ Chinese species list from distributions
# ---------------------------------------------------------------
comm_files <- list.files(PATH$derived, pattern = "^comm_.*_50km\\.rds$", full.names = TRUE)
cn_species <- lapply(comm_files, function(f) {
  d <- readRDS(f); data.frame(class = d$class[1], species = unique(d$species))
}) |> bind_rows() |> distinct()
log_msg("中国物种（50 km 网格出现）/ Chinese species in 50 km grids: ", nrow(cn_species))
print(table(cn_species$class))

# ---------------------------------------------------------------
# 2. 骨架：TetrapodTraits v2.0.1 / Backbone
# ---------------------------------------------------------------
tt <- data.table::fread(file.path(TR, "TetrapodTraits_v2.0.1.csv"), showProgress = FALSE)
tt[, species := norm_name(Scientific.Name)]
log_msg("TetrapodTraits: ", nrow(tt), " species globally")

# --- 名录协调 / Taxonomic harmonisation -------------------------------------
# TetrapodTraits 采用系统发育名录 (Jetz/Upham/Tonini)，与 IUCN 名录常有差异；
# 幸而其自带 IUCN_Binomial 列，可作第二把钥匙；仍不匹配者再查 GBIF 骨架。
# TetrapodTraits follows phylogeny taxonomies, which often differ from IUCN
# names; its IUCN_Binomial column provides a second key, and GBIF's backbone
# resolves the remainder.
tt[, species_iucn := norm_name(IUCN_Binomial)]

harmonise <- function(cn_species, tt) {
  key1 <- setNames(tt$species, tt$species)                    # 直接匹配 / direct
  key2 <- setNames(tt$species, tt$species_iucn)               # 经 IUCN 名 / via IUCN name
  key2 <- key2[!is.na(names(key2)) & nzchar(names(key2))]

  cn_species$tt_key <- unname(key1[cn_species$species])
  need <- is.na(cn_species$tt_key)
  cn_species$tt_key[need] <- unname(key2[cn_species$species[need]])
  cn_species$match_route <- ifelse(!is.na(unname(key1[cn_species$species])), "direct",
                            ifelse(!is.na(cn_species$tt_key), "iucn_binomial", NA))

  # 仍未匹配者：查 GBIF 骨架获取被接受名 / resolve remainder via GBIF backbone
  todo <- cn_species$species[is.na(cn_species$tt_key)]
  if (length(todo)) {
    log_msg("  GBIF 骨架解析 ", length(todo), " 个未匹配学名 ...")
    resolved <- vapply(todo, function(s) {
      u <- paste0("https://api.gbif.org/v1/species/match?name=",
                  utils::URLencode(s, reserved = TRUE), "&strict=false")
      r <- try(jsonlite::fromJSON(u), silent = TRUE)
      if (inherits(r, "try-error") || is.null(r$species)) return(NA_character_)
      norm_name(r$species)
    }, character(1), USE.NAMES = FALSE)
    acc <- unname(key1[resolved]); acc2 <- unname(key2[resolved])
    fill <- ifelse(!is.na(acc), acc, acc2)
    idx <- which(is.na(cn_species$tt_key))
    cn_species$tt_key[idx]      <- fill
    cn_species$match_route[idx] <- ifelse(!is.na(fill), "gbif_backbone", NA)
  }
  cn_species
}
cn_species <- cache_rds("taxonomy_harmonisation", harmonise(cn_species, tt))
log_msg("名录协调 / harmonisation: ",
        sum(!is.na(cn_species$tt_key)), " / ", nrow(cn_species), " matched")
print(table(cn_species$match_route, useNA = "ifany"))
write_table(cn_species, "table_S3c_taxonomy_harmonisation")

# ---------------------------------------------------------------
# 3. 定义核心性状（四类群语义一致）/ Core traits, identical across classes
# ---------------------------------------------------------------
# 说明 / Notes:
#  body_mass       体重 (g, log10)             — 代谢与生态位的主轴
#  body_length     体长 (mm, log10)            — 体型的第二度量
#  nocturnality    夜行性 (0-1)                — 时间生态位
#  verticality     垂直生境位置 (0-1)          — 穴居(0) -> 空中(1)
#  diet_breadth    食性宽度 (类别数)           — 营养生态位宽度
#  diet_vert       脊椎动物食性占比            — 营养级代理
#  diet_plant      植物性食物占比              — 营养级代理
#  habitat_breadth IUCN 主要生境类型数          — 生境生态位宽度
#  litter_size     窝仔/窝卵数                 — 繁殖投入（生活史"快慢"轴）
#  max_longevity   最大寿命 (月)               — 生活史"快慢"轴
#  range_size      分布范围面积 (km2, log10)   — 地理稀有性
CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
          "diet_breadth", "diet_vert", "diet_plant", "habitat_breadth",
          "litter_size", "max_longevity", "range_size")

core <- tt[, .(
  species,
  tree_taxon      = TreeTaxon,
  class           = Class,
  order           = Order,
  family          = Family,
  body_mass       = suppressWarnings(as.numeric(BodyMass_g)),
  body_length     = suppressWarnings(as.numeric(BodyLength_mm)),
  nocturnality    = suppressWarnings(as.numeric(Nocturnality)),
  verticality     = suppressWarnings(as.numeric(Verticality)),
  diet_breadth    = suppressWarnings(as.numeric(DietBreadth)),
  diet_vert       = suppressWarnings(as.numeric(DietVend) + as.numeric(DietVect) +
                                     as.numeric(DietVfish) + as.numeric(DietVunk)),
  diet_plant      = suppressWarnings(as.numeric(DietFruit) + as.numeric(DietNect) +
                                     as.numeric(DietSeed)  + as.numeric(DietPlant)),
  habitat_breadth = suppressWarnings(as.numeric(MajorHabitatSum)),
  litter_size     = suppressWarnings(as.numeric(LitterSize)),
  max_longevity   = suppressWarnings(as.numeric(MaxLongevity)),
  range_size      = suppressWarnings(as.numeric(RangeSize_km2)),
  # 保留插补标记以便敏感性分析 / keep imputation flags for sensitivity analysis
  flag_mass = ImputedMass, flag_len = ImputedLength, flag_act = ImputedActTime,
  flag_hab = ImputedHabitat, flag_mhab = ImputedMajorHabitat
)]

# 对数变换（右偏性状）/ log10-transform right-skewed traits
for (v in c("body_mass", "body_length", "range_size", "max_longevity", "litter_size")) {
  core[[v]] <- ifelse(is.finite(core[[v]]) & core[[v]] > 0, log10(core[[v]]), NA_real_)
}

# ---------------------------------------------------------------
# 4. 补充数据库 / Supplementary databases
# ---------------------------------------------------------------
add_left <- function(base, df, keycol, cols, prefix = "") {
  df <- as.data.frame(df)
  df$species <- norm_name(df[[keycol]])
  df <- df[!duplicated(df$species), c("species", cols), drop = FALSE]
  if (nzchar(prefix)) names(df)[-1] <- paste0(prefix, names(df)[-1])
  dplyr::left_join(base, df, by = "species")
}

ext <- as.data.frame(core)

## 4.1 AVONET —— 鸟类形态与扩散能力 / bird morphology & dispersal
f <- file.path(TR, "AVONET_supp1.xlsx")
if (file.exists(f)) {
  sheets <- readxl::excel_sheets(f)
  sh <- grep("BirdLife", sheets, value = TRUE)[1]
  av <- readxl::read_excel(f, sheet = sh)
  nm <- names(av)
  pick <- c(`hwi` = grep("^Hand-Wing", nm, value = TRUE)[1],
            `tarsus` = grep("^Tarsus", nm, value = TRUE)[1],
            `wing` = grep("^Wing.Length", nm, value = TRUE)[1],
            `beak_len` = grep("^Beak.Length_Culmen", nm, value = TRUE)[1],
            `trophic_level` = grep("^Trophic.Level", nm, value = TRUE)[1],
            `trophic_niche` = grep("^Trophic.Niche", nm, value = TRUE)[1],
            `migration` = grep("^Migration", nm, value = TRUE)[1],
            `habitat_density` = grep("^Habitat.Density", nm, value = TRUE)[1])
  pick <- pick[!is.na(pick)]
  av2 <- av[, c(grep("Species1", nm, value = TRUE)[1], unname(pick))]
  names(av2) <- c("sp", names(pick))
  ext <- add_left(ext, av2, "sp", names(pick))
  log_msg("  + AVONET: ", length(pick), " bird traits")
}

## 4.2 COMBINE —— 哺乳类生活史 / mammal life history
f <- file.path(TR, "COMBINE_imputed.csv")
if (file.exists(f)) {
  cb <- data.table::fread(f, showProgress = FALSE)
  want <- intersect(c("gestation_length_d", "maturity_d", "litters_per_year_n",
                      "generation_length_d", "dispersal_km", "home_range_km2",
                      "density_n_km2", "hibernation_torpor", "fossoriality",
                      "brain_mass_g", "habitat_breadth_n"), names(cb))
  key <- if ("iucn2020_binomial" %in% names(cb)) "iucn2020_binomial" else names(cb)[1]
  ext <- add_left(ext, cb, key, want, prefix = "mam_")
  log_msg("  + COMBINE: ", length(want), " mammal traits")
}

## 4.3 AmphiBIO —— 两栖类繁殖与生态 / amphibian reproduction & ecology
f <- file.path(TR, "AmphiBIO_v1.zip")
if (file.exists(f)) {
  td <- file.path(tempdir(), "amphibio"); dir.create(td, showWarnings = FALSE)
  utils::unzip(f, exdir = td)
  cs <- list.files(td, pattern = "csv$", recursive = TRUE, full.names = TRUE)
  cs <- cs[which.max(file.size(cs))]
  ab <- data.table::fread(cs, showProgress = FALSE)
  want <- intersect(c("Age_at_maturity_min_y", "Age_at_maturity_max_y",
                      "Body_size_mm", "Size_at_maturity_min_mm",
                      "Longevity_max_y", "Litter_size_min_n", "Litter_size_max_n",
                      "Reproductive_output_y", "Offspring_size_min_mm",
                      "Dir", "Lar", "Viv", "Diu", "Noc", "Crepu",
                      "Fos", "Ter", "Aqu", "Arb"), names(ab))
  key <- if ("Species" %in% names(ab)) "Species" else names(ab)[1]
  ext <- add_left(ext, ab, key, want, prefix = "amp_")
  log_msg("  + AmphiBIO: ", length(want), " amphibian traits")
}

## 4.4 Meiri 2018 —— 蜥蜴体温与生态 / lizard body temperature & ecology
f <- file.path(TR, "Meiri2018_lizards.csv")
if (file.exists(f)) {
  mz <- data.table::fread(f, showProgress = FALSE)
  nm <- names(mz)
  want <- nm[grepl("body temperature|activity time|foraging mode|reproductive mode|clutch size|substrate|diet",
                   nm, ignore.case = TRUE)]
  want <- head(want, 12)
  key <- grep("binomial|species", nm, ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(key) && length(want)) {
    ext <- add_left(ext, mz, key, want, prefix = "liz_")
    log_msg("  + Meiri 2018: ", length(want), " lizard traits")
  }
}

## 4.5 GlobTherm —— 热耐受 / thermal tolerance
f <- file.path(TR, "GlobTherm.csv")
if (file.exists(f)) {
  gt <- data.table::fread(f, showProgress = FALSE)
  nm <- names(gt)
  gen <- grep("^Genus$", nm, value = TRUE); sp <- grep("^Species$", nm, value = TRUE)
  if (length(gen) && length(sp)) {
    gt$binom <- paste(gt[[gen]], gt[[sp]])
    want <- intersect(c("Tmax", "tmin", "max_metric", "min_metric"), nm)
    ext <- add_left(ext, gt, "binom", want, prefix = "th_")
    log_msg("  + GlobTherm: ", length(want), " thermal traits")
  }
}

## 4.6 Etard 2020 —— IUCN 生境宽度与人工生境利用 / habitat breadth & artificial habitat use
for (cl in c("birds", "mammals", "reptiles", "amphibians")) {
  f <- file.path(TR, paste0("Etard2020_", cl, ".csv"))
  if (!file.exists(f)) next
  et <- data.table::fread(f, showProgress = FALSE)
  want <- intersect(c("Habitat_breadth_IUCN", "Artificial_habitat_use",
                      "Specialisation", "Primary_diet", "Diel_activity",
                      "Trophic_level", "Generation_length_d"), names(et))
  key <- grep("Best_guess_binomial|Binomial|Species", names(et), value = TRUE)[1]
  if (is.na(key) || !length(want)) next
  tmp <- as.data.frame(et)[, c(key, want), drop = FALSE]
  tmp$species <- norm_name(tmp[[key]]); tmp[[key]] <- NULL
  tmp <- tmp[!duplicated(tmp$species), ]
  names(tmp)[names(tmp) != "species"] <- paste0("et_", names(tmp)[names(tmp) != "species"])
  # 只填补尚缺的物种 / fill only species not yet covered
  ext <- dplyr::rows_patch(ext, tmp, by = "species", unmatched = "ignore") |>
    suppressWarnings() |> tryCatch(error = function(e) dplyr::left_join(ext, tmp, by = "species"))
}
log_msg("  + Etard 2020 habitat breadth merged")

# ---------------------------------------------------------------
# 4.7 核心性状补缺 / Gap-filling of core traits from primary sources
#     骨架中 max_longevity 与两栖类的食性/窝仔数缺失严重，
#     用 AnAge、Amniote Life History DB、AmphiBIO 的原始记录补齐。
#     The backbone is sparse for longevity and for amphibian diet/litter;
#     fill from AnAge, the Amniote database and AmphiBIO.
# ---------------------------------------------------------------
fill_from <- function(base, df, key, map) {
  df <- as.data.frame(df)
  df$species <- norm_name(df[[key]])
  df <- df[!is.na(df$species) & !duplicated(df$species), ]
  for (target in names(map)) {
    src <- map[[target]]
    if (!all(src %in% names(df))) next
    val <- if (length(src) == 1) suppressWarnings(as.numeric(df[[src]]))
           else rowMeans(sapply(src, function(s) suppressWarnings(as.numeric(df[[s]]))), na.rm = TRUE)
    v <- setNames(val, df$species)
    hit <- v[base$species]
    n_before <- sum(is.na(base[[target]]))
    base[[target]] <- ifelse(is.na(base[[target]]) & is.finite(hit), hit, base[[target]])
    n_filled <- n_before - sum(is.na(base[[target]]))
    if (n_filled > 0) log_msg("     补缺 ", target, ": +", n_filled, " 种")
  }
  base
}

## AnAge：全脊椎动物最大寿命（年 -> 月）/ AnAge maximum longevity (years -> months)
f <- file.path(TR, "AnAge_dataset.zip")
if (file.exists(f)) {
  td <- file.path(tempdir(), "anage"); dir.create(td, showWarnings = FALSE)
  utils::unzip(f, exdir = td)
  ff <- list.files(td, pattern = "txt$|csv$", recursive = TRUE, full.names = TRUE)
  if (length(ff)) {
    an <- data.table::fread(ff[which.max(file.size(ff))], showProgress = FALSE)
    if (all(c("Genus", "Species") %in% names(an))) {
      an$binom <- paste(an$Genus, an$Species)
      lc <- grep("Maximum longevity", names(an), value = TRUE)[1]
      if (!is.na(lc)) {
        # TetrapodTraits 的 MaxLongevity 单位为「年」，此处保持一致
        # TetrapodTraits stores MaxLongevity in YEARS; keep the same unit
        an$long_y <- suppressWarnings(as.numeric(an[[lc]]))
        ext <- fill_from(ext, an, "binom", list(max_longevity = "long_y"))
        log_msg("  + AnAge 寿命补缺 / longevity gap-filling")
      }
    }
  }
}

## Amniote Life History Database：鸟/兽/爬行 寿命与窝仔数
f <- file.path(TR, "Amniote_LifeHistory.zip")
if (file.exists(f)) {
  td <- file.path(tempdir(), "amniote"); dir.create(td, showWarnings = FALSE)
  utils::unzip(f, exdir = td)
  ff <- list.files(td, pattern = "csv$", recursive = TRUE, full.names = TRUE)
  if (length(ff)) {
    am <- data.table::fread(ff[which.max(file.size(ff))], showProgress = FALSE)
    am[am == -999] <- NA
    if (all(c("genus", "species") %in% names(am))) {
      am$binom <- paste(am$genus, am$species)
      mp <- list()
      if ("longevity_y" %in% names(am)) { am$long_m <- suppressWarnings(as.numeric(am$longevity_y))*12
                                          mp$max_longevity <- "long_m" }
      lz <- intersect(c("litter_or_clutch_size_n"), names(am))
      if (length(lz)) mp$litter_size <- lz
      if (length(mp)) {
        ext <- fill_from(ext, am, "binom", mp)
        log_msg("  + Amniote DB 补缺 / gap-filling")
      }
    }
  }
}

## AmphiBIO：两栖类寿命、窝仔数、食性 / amphibian longevity, litter, diet
if ("amp_Longevity_max_y" %in% names(ext)) {
  ext$max_longevity <- ifelse(is.na(ext$max_longevity) & is.finite(ext$amp_Longevity_max_y),
                              ext$amp_Longevity_max_y * 12, ext$max_longevity)
}
lz <- intersect(c("amp_Litter_size_min_n", "amp_Litter_size_max_n"), names(ext))
if (length(lz)) {
  amp_lit <- rowMeans(as.data.frame(ext[, lz, drop = FALSE]), na.rm = TRUE)
  ext$litter_size <- ifelse(is.na(ext$litter_size) & is.finite(amp_lit), amp_lit, ext$litter_size)
}
# 两栖类食性：AmphiBIO 的二元食物类别 -> 食性宽度与营养组成
amp_diet <- intersect(paste0("amp_", c("Leaves","Flowers","Seeds","Fruits","Arthro","Vert")), names(ext))
if (length(amp_diet) >= 3) {
  D <- as.data.frame(ext[, amp_diet, drop = FALSE])
  D[] <- lapply(D, function(z) ifelse(is.na(suppressWarnings(as.numeric(z))), 0, 1))
  nb <- rowSums(D)
  ext$diet_breadth <- ifelse(is.na(ext$diet_breadth) & nb > 0, nb, ext$diet_breadth)
}
log_msg("  + AmphiBIO 补缺 / gap-filling")

# 补缺后再做对数变换（这些源为原始尺度）/ re-log the newly filled raw-scale values
for (v in c("max_longevity", "litter_size")) {
  big <- which(is.finite(ext[[v]]) & ext[[v]] > 10)   # 已取对数的值不会 >10
  if (length(big)) ext[[v]][big] <- log10(ext[[v]][big])
}

# ---------------------------------------------------------------
# 5. 取中国物种子集并统计覆盖率 / Subset to China and summarise coverage
# ---------------------------------------------------------------
# ext 以 TetrapodTraits 学名为键；改用协调后的 tt_key 连接，
# 并把物种名换回分布数据所用的名字，以便与群落矩阵对接。
# ext is keyed by the TetrapodTraits name; join via the harmonised key and
# rename back to the distribution-data name so it matches the assemblages.
link <- cn_species |> filter(!is.na(tt_key)) |>
  select(dist_species = species, class_dist = class, tt_key, match_route)

cn <- ext |>
  inner_join(link, by = c("species" = "tt_key"), relationship = "many-to-many") |>
  mutate(species_tt = species, species = dist_species) |>
  select(-dist_species) |>
  distinct(species, .keep_all = TRUE)
# 类群以分布数据为准（IUCN 的 Reptilia 含龟鳖/鳄，与树分类一致）
cn$class <- ifelse(is.na(cn$class), cn$class_dist, cn$class)

missing_sp <- cn_species |> filter(is.na(tt_key))
log_msg("骨架匹配 / matched to backbone: ", nrow(cn), " / ", nrow(cn_species),
        "  (未匹配 ", nrow(missing_sp), ")")
print(table(cn$class))

cov <- cn |>
  select(class, all_of(CORE)) |>
  pivot_longer(-class, names_to = "trait", values_to = "v") |>
  group_by(class, trait) |>
  summarise(n = n(), n_obs = sum(!is.na(v)),
            coverage = round(100 * mean(!is.na(v)), 1), .groups = "drop") |>
  arrange(trait, class)

# 诚实报告：TetrapodTraits 上游已插补的比例（Moura et al. 2024 提供标记）
# Honest reporting: share already imputed upstream by TetrapodTraits,
# using the ImputedX flags supplied by Moura et al. (2024)
FLAG_MAP <- c(body_mass = "flag_mass", body_length = "flag_len",
              nocturnality = "flag_act", verticality = "flag_hab",
              habitat_breadth = "flag_mhab")
up <- lapply(names(FLAG_MAP), function(tr) {
  fl <- FLAG_MAP[[tr]]
  if (!fl %in% names(cn)) return(NULL)
  cn |> group_by(class) |>
    summarise(trait = tr,
              upstream_imputed_pct = round(100 * mean(as.numeric(.data[[fl]]) %in% 1, na.rm = TRUE), 1),
              .groups = "drop")
}) |> bind_rows()
cov <- left_join(cov, up, by = c("class", "trait"))

# 标注性状集归属 / label which trait belongs to which analysis set
CORE_ALL4 <- c("body_mass", "body_length", "nocturnality", "verticality",
               "habitat_breadth", "range_size")
cov$trait_set <- ifelse(cov$trait %in% CORE_ALL4, "core (all 4 classes)",
                 ifelse(cov$trait == "max_longevity", "excluded (data-deficient)",
                        "extended (birds/mammals/reptiles)"))
write_table(cov, "table_S3_trait_coverage_before_imputation")
print(as.data.frame(cov))

saveRDS(cn, file.path(PATH$derived, "traits_china_raw.rds"))
saveRDS(missing_sp, file.path(PATH$derived, "traits_unmatched_species.rds"))
write_table(missing_sp, "table_S3b_unmatched_species")

log_msg("=== 03 完成（插补见 03b）/ done (imputation in 03b) ===")
