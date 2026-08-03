# ============================================================
# 科学问题 / Scientific question:
#   §5d 指出，中国两栖类的群落功能结构无法可靠判定，因为其核心性状
#   大多来自模型插补而非实测。若换用**中国本土编纂**的性状数据集，
#   实测覆盖能提高多少？结论会不会改变？
#   Section 5d concluded that Chinese amphibian assemblage structure could not
#   be determined because most core traits were modelled rather than measured.
#   Do China-specific trait compilations change that?
#
# 分析目标 / Objective:
#   整合王彦平团队在《生物多样性》发表的中国本土性状数据集，
#   与全球库（TetrapodTraits / AmphiBIO / Etard / Meiri / ReptTraits）取并集，
#   逐性状记录数据来源，并重算实测覆盖率。
#
# 输入数据 / Input:
#   02_data_raw/traits/china_wang/
#     China_amphibian_traits_Song2022.xlsx   Song, Chen & Wang 2022
#       生物多样性 30:22053  doi:10.17520/biods.2022053  （591 种 × 28 性状）
#     China_lizards_traits_Zhong2022.xlsx    Zhong, Chen & Wang 2022
#       生物多样性 30:22071  doi:10.17520/biods.2022071  （250 种 × 24 性状）
#     China_birds_2021201 / China_snakes_2023126（若已下载）
#       Wang, Song & Zhong 2021 doi:10.17520/biods.2021201
#       Wang J. et al. 2023    doi:10.17520/biods.2023126
#
# 预期输出 / Expected output:
#   03_data_derived/traits_observed_china.rds      逐种逐性状的实测值与来源
#   04_results/tables/table_30_china_trait_coverage.csv
#
# 关键假设 / Key assumptions:
#   - 这些数据集的缺失值写作字符串 "NA" 而非真正的 NA，必须显式处理，
#     否则会把缺失误读为有值（初次统计即因此得到 100% 的假覆盖率）。
#   - 两栖类以雄性体长为体长代表（雄性覆盖 96% 高于雌性 81%）；
#     两性体长都有时取均值。
#   - 垂直生境位置由「成体微生境」映射到与主分析同义的 0（穴居）至
#     3（树栖）刻度；映射表在代码中显式给出，便于复核。
#
# 主要包 / Main packages: readxl, dplyr, data.table
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({ library(readxl); library(dplyr) })

log_msg("=== 03c 中国本土性状数据集整合 / China-specific trait compilations ===")
WD <- file.path(PATH$raw_trait, "china_wang")

# 这些表把缺失写成字符串 "NA"；不显式处理会得到假的 100% 覆盖率
# Missing values are the literal string "NA" in these files
nz  <- function(x) { x <- as.character(x); !is.na(x) & trimws(x) != "" & x != "NA" }
#' 解析可能写成区间的数值
#'
#' 这几份数据集把测量值写成 "1010~1210" 这样的区间。若用 gsub 剥离非数字
#' 字符，"1010~1210" 会被拼成 10101210——即 10 公里长的蛇。必须先按分隔符
#' 切分再取中点。分隔符包括 ~ 、连字符、各种破折号；只有夹在数字之间的
#' 连字符才视为区间，以免误伤负号。
#' These files record measurements as ranges such as "1010~1210". Stripping
#' non-numeric characters concatenates them into 10101210 - a ten-metre snake.
#' Split on the separator and take the midpoint instead.
num <- function(x) {
  x <- as.character(x)
  x[!nz(x)] <- NA
  x <- trimws(gsub("[[:space:]\u00a0]+", "", x, perl = TRUE))
  x <- gsub("(?<=[0-9])[~\u2013\u2014\u2212-](?=[0-9])", "|", x, perl = TRUE)
  vapply(x, function(v) {
    if (is.na(v) || !nzchar(v)) return(NA_real_)
    parts <- suppressWarnings(as.numeric(strsplit(v, "|", fixed = TRUE)[[1]]))
    parts <- parts[is.finite(parts)]
    if (!length(parts)) return(NA_real_)
    mean(parts)
  }, numeric(1), USE.NAMES = FALSE)
}

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))

# ---------------------------------------------------------------
# 1. 两栖类 / Amphibians (Song, Chen & Wang 2022)
# ---------------------------------------------------------------
amph <- NULL
f_amp <- file.path(WD, "China_amphibian_traits_Song2022.xlsx")
if (file.exists(f_amp)) {
  W <- read_excel(f_amp, sheet = 1)
  W$sp <- norm_name(W[["种拉丁名"]])

  # 体长：雄性覆盖 96% 优于雌性 81%；两性皆有时取均值
  # Body length: male coverage (96%) exceeds female (81%); average when both
  ml <- num(W[["雄性体长/mm"]]); fl <- num(W[["雌性体长/mm"]])
  blen <- ifelse(is.finite(ml) & is.finite(fl), (ml + fl) / 2,
                 ifelse(is.finite(ml), ml, fl))

  # 活动模式 -> 夜行性 0/1（全日性记 0.5）
  act <- as.character(W[["活动模式"]]); act[!nz(act)] <- NA
  noct <- dplyr::case_when(grepl("夜行", act) ~ 1,
                           grepl("昼行", act) ~ 0,
                           grepl("全日", act) ~ 0.5,
                           TRUE ~ NA_real_)

  # 成体微生境 -> 垂直生境位置，刻度与主分析一致：穴居 0 / 地面-水生 1 / 树栖 3
  # Adult microhabitat mapped onto the same verticality scale as the main analysis
  mh <- as.character(W[["成体微生境"]]); mh[!nz(mh)] <- NA
  vert <- dplyr::case_when(grepl("穴|洞|地下", mh) ~ 0,
                           grepl("树|灌|攀", mh)   ~ 3,
                           grepl("水|溪|流|塘", mh) ~ 1,
                           grepl("陆|地面|草", mh) ~ 1,
                           TRUE ~ NA_real_)

  amph <- data.frame(species = W$sp, class = "Amphibia",
                     cn_body_length = blen,
                     cn_nocturnality = noct,
                     cn_verticality = vert,
                     cn_litter_size = num(W[["窝卵数"]]),
                     stringsAsFactors = FALSE) |>
    filter(!is.na(species)) |> distinct(species, .keep_all = TRUE)
  log_msg("两栖 / amphibians: ", nrow(amph), " 种（Song et al. 2022）")
}

# ---------------------------------------------------------------
# 2. 蜥蜴 / Lizards (Zhong, Chen & Wang 2022)
# ---------------------------------------------------------------
liz <- NULL
f_liz <- file.path(WD, "China_lizards_traits_Zhong2022.xlsx")
if (file.exists(f_liz)) {
  L <- read_excel(f_liz, sheet = 1)
  L$sp <- norm_name(L[["种拉丁名"]])
  at <- as.character(L[["活动时间"]]); at[!nz(at)] <- NA
  hb <- as.character(L[["栖息生境"]]); hb[!nz(hb)] <- NA
  liz <- data.frame(species = L$sp, class = "Reptilia",
                    cn_body_length = num(L[["平均体长(mm)"]]),
                    cn_body_mass   = num(L[["平均体重(g)"]]),
                    cn_nocturnality = dplyr::case_when(grepl("夜", at) ~ 1,
                                                       grepl("昼|日", at) ~ 0,
                                                       TRUE ~ NA_real_),
                    cn_verticality = dplyr::case_when(grepl("穴|洞|地下", hb) ~ 0,
                                                      grepl("树|灌|攀", hb) ~ 3,
                                                      grepl("岩|石", hb) ~ 2,
                                                      grepl("地面|草|沙|荒", hb) ~ 1,
                                                      TRUE ~ NA_real_),
                    cn_habitat_breadth = num(L[["栖息地宽度"]]),
                    cn_litter_size = num(L[["窝卵数"]]),
                    stringsAsFactors = FALSE) |>
    filter(!is.na(species)) |> distinct(species, .keep_all = TRUE)
  log_msg("蜥蜴 / lizards: ", nrow(liz), " 种（Zhong et al. 2022）")
}

# ---------------------------------------------------------------
# 3. 鸟类 / Birds (Wang, Song & Zhong 2021)
# ---------------------------------------------------------------
# 该数据集本地已有副本（随 Huang et al. 2023 的仓库分发），故不重复下载。
# A copy ships with the Huang et al. (2023) repository, so it is not re-fetched.
bird <- NULL
f_bird <- file.path(PATH$raw, "huang2023_repo", "ave", "Chinesebirdsdata",
                    "Chinesebirdsdata.csv")
if (file.exists(f_bird)) {
  B <- data.table::fread(f_bird)
  B$sp <- norm_name(B[["iucn_2021_binomial"]])
  mid <- function(a, b) { a <- num(B[[a]]); b <- num(B[[b]])
    ifelse(is.finite(a) & is.finite(b), (a + b) / 2,
           ifelse(is.finite(a), a, b)) }
  bird <- data.frame(species = B$sp, class = "Aves",
                     cn_body_mass = {
                       m <- mid("male_mass_min_g", "male_mass_max_g")
                       f <- mid("female_mass_min_g", "female_mass_max_g")
                       ifelse(is.finite(m) & is.finite(f), (m + f) / 2,
                              ifelse(is.finite(m), m, f)) },
                     cn_body_length = {
                       m <- mid("male_size_min_mm", "male_size_max_mm")
                       f <- mid("female_size_min_mm", "female_size_max_mm")
                       ifelse(is.finite(m) & is.finite(f), (m + f) / 2,
                              ifelse(is.finite(m), m, f)) },
                     cn_litter_size = mid("clutch_size_min", "clutch_size_max"),
                     stringsAsFactors = FALSE) |>
    filter(!is.na(species), species != "") |> distinct(species, .keep_all = TRUE)
  log_msg("鸟类 / birds: ", nrow(bird), " 种（Wang et al. 2021，本地副本）")
}

# 蛇类数据集（Wang J. et al. 2023, doi:10.17520/biods.2023126）未能取得：
# 期刊站点对本机 IP 限流，多次重试均连接失败。爬行类的全球库覆盖已达
# 83-100%，故缺此项不影响结论；获取方式见 DATA.md。
# The snake dataset could not be retrieved (journal site rate-limited this IP).
# Reptile coverage from global sources is already 83-100%, so nothing hinges
# on it; DATA.md records how to obtain it.

# ---------------------------------------------------------------
# 3b. 哺乳类 / Mammals (Ding et al. 2022)
# ---------------------------------------------------------------
# 该数据集**不分性别记录量度**，因此哺乳类无法进入性二态性分析（脚本 16）。
# 其形态附肢性状（耳长 6.7%、尾长 4.3%、前臂长 3.7%）在本项目物种中过于
# 稀疏，不足以支撑群落尺度分析（例如 Allen 法则检验）。
# 本项目哺乳类的核心性状本已有 92.7-98.5% 的实测率，故这份数据的价值主要
# 在**交叉验证**而非补缺——见 table_34。
# This compilation does not separate the sexes, so mammals cannot enter the
# sexual-dimorphism analysis. Its appendage traits are too sparse here for
# assemblage-level work. Its value is cross-validation, not gap-filling.
mam <- NULL
f_mam <- file.path(WD, "China_mammals_traits_Ding2022.xlsx")
if (file.exists(f_mam)) {
  M <- read_excel(f_mam, sheet = 1)
  sci <- grep("Scientific name", names(M), value = TRUE)[1]
  pick <- function(pat) { v <- grep(pat, names(M), value = TRUE)[1]
                          if (is.na(v)) rep(NA_real_, nrow(M)) else num(M[[v]]) }
  ac <- as.character(M[[grep("Activity_cycle", names(M), value = TRUE)[1]]])
  ac[!nz(ac)] <- NA
  mam <- data.frame(species = norm_name(M[[sci]]), class = "Mammalia",
                    cn_body_mass = pick("Body_mass_g"),
                    cn_body_length = pick("Body_length_mm"),
                    cn_litter_size = pick("Litter_size_n"),
                    cn_habitat_breadth = pick("Habitat_ ?breadth"),
                    cn_nocturnality = dplyr::case_when(
                      grepl("夜|Noct", ac) ~ 1,
                      grepl("昼|Diur", ac) ~ 0,
                      grepl("晨昏|Crepus", ac) ~ 0.5,
                      TRUE ~ NA_real_),
                    stringsAsFactors = FALSE) |>
    filter(!is.na(species), species != "") |> distinct(species, .keep_all = TRUE)
  log_msg("哺乳 / mammals: ", nrow(mam), " 种（Ding et al. 2022）")
}

cn <- bind_rows(amph, liz, bird, mam)
saveRDS(cn, file.path(PATH$derived, "traits_observed_china.rds"))

# ---------------------------------------------------------------
# 4. 覆盖率：加入本土数据前后 / Coverage before and after
# ---------------------------------------------------------------
FLAGS <- c(body_mass = "flag_mass", body_length = "flag_len",
           nocturnality = "flag_act", verticality = "flag_hab",
           habitat_breadth = "flag_mhab")
rows <- list()
for (cl in c("Amphibia", "Reptilia", "Aves", "Mammalia")) {
  tr <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  cc <- cn |> filter(class == cl)
  j  <- cc[match(tr$species, cc$species), ]
  for (v in c("body_length", "nocturnality", "verticality",
              "habitat_breadth", "body_mass", "litter_size")) {
    fl <- if (v %in% names(FLAGS)) FLAGS[[v]] else NULL
    # 全球库中「实测」= 存在且未被 TetrapodTraits 标记为插补
    was <- if (!is.null(fl) && fl %in% names(tr))
      is.finite(suppressWarnings(as.numeric(tr[[v]]))) &
        !(suppressWarnings(as.numeric(tr[[fl]])) %in% 1) else
      is.finite(suppressWarnings(as.numeric(tr[[v]])))
    now_col <- paste0("cn_", v)
    add <- if (now_col %in% names(j)) is.finite(j[[now_col]]) else rep(FALSE, nrow(tr))
    rows[[length(rows) + 1]] <- data.frame(
      class = cl, trait = v, n_species = nrow(tr),
      pct_global_observed = round(100 * mean(was), 1),
      pct_china_dataset   = round(100 * mean(add), 1),
      pct_union           = round(100 * mean(was | add), 1),
      gain_pp             = round(100 * (mean(was | add) - mean(was)), 1))
  }
}
cov <- bind_rows(rows)
write_table(cov, "table_30_china_trait_coverage")
print(as.data.frame(cov))

# ---------------------------------------------------------------
# 5. 交叉验证：中国本土值 vs 全球库值 / Cross-validation
# ---------------------------------------------------------------
# 这些数据集是独立编纂的，因此二者的一致程度是对**双方**的检验。
# 量表不同的性状（如生境宽度）只报告相关，不做合并。
# The compilations are independent, so agreement tests both sources. Traits on
# different scales are reported as correlations only, never merged.
val <- list()
for (cl in c("Amphibia", "Reptilia", "Aves", "Mammalia")) {
  tr <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  cc <- cn |> filter(class == cl)
  if (!nrow(cc)) next
  j <- cc[match(tr$species, cc$species), ]
  for (v in c("body_mass", "body_length", "habitat_breadth", "litter_size")) {
    a <- suppressWarnings(as.numeric(j[[paste0("cn_", v)]]))
    b <- suppressWarnings(as.numeric(tr[[v]]))
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 30) next
    # traits_imputed.rds 中体型、窝仔数、分布范围、寿命**已经是 log10**；
    # 中国本土表是原始尺度。因此只对中国值取一次对数，绝不对已取对数的
    # 全球值再取一次（那会扭曲 Pearson，Spearman 则因单调而不受影响）。
    # The stored global values are already log10 for these traits; log only the
    # China side. Double-logging distorts Pearson though not Spearman.
    LOGGED <- c("body_mass", "body_length", "litter_size", "range_size",
                "max_longevity")
    lg <- v %in% LOGGED
    x <- if (lg) log10(pmax(a[ok], 1e-9)) else a[ok]   # 中国值 -> log10
    y <- b[ok]                                          # 全球值：已是所需尺度
    val[[length(val) + 1]] <- data.frame(
      class = cl, trait = v, n = sum(ok),
      scale = if (lg) "log10 (both sides)" else "raw",
      pearson = round(stats::cor(x, y), 3),
      spearman = round(stats::cor(x, y, method = "spearman"), 3),
      mean_china = round(mean(x), 2), mean_global = round(mean(y), 2),
      mean_diff = round(mean(x) - mean(y), 3))
  }
}
if (length(val)) {
  val <- bind_rows(val)
  write_table(val, "table_34_china_global_validation")
  log_msg("--- 中国本土值与全球库值的一致性 ---")
  print(as.data.frame(val))
}

log_msg("=== 03c 完成 / done ===")
