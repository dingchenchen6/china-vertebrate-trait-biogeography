# ============================================================
# 科学问题 / Scientific question:
#   性二态性（sexual size dimorphism, SSD）是全球性状库普遍缺失的维度——
#   AVONET、COMBINE、TetrapodTraits 都只给单一体型值。王彦平团队的三份
#   中国数据集**分性别记录量度**，使得在群落尺度上提问成为可能：
#     1. 群落平均性二态性在中国有地理格局吗？
#     2. 它是否与体型的地理格局一致——即 Rensch 法则在群落尺度上成立吗？
#     3. 恒温与变温类群的 SSD 格局是否受同一批环境轴驱动？
#   Sexual size dimorphism is absent from the global trait databases, which
#   report a single body size per species. Three China-specific compilations
#   record measurements by sex, which makes assemblage-level questions possible.
#
#   哺乳类缺席：Ding et al. (2022, 生物多样性 30:21520) 的中国哺乳动物数据集
#   记录体重、脑容量、体长、尾长、耳长等，但**不分性别**，因此哺乳类无法
#   进入本分析。这不是遗漏，而是现有中文与全球数据源的共同空白。
#   Mammals are absent because no available source separates the sexes.
#
# 分析目标 / Objective:
#   1. 由分性别量度计算每个物种的 SSD 指数
#   2. 计算逐格的群落加权平均 SSD 与其方差
#   3. 制图，并用与主分析同一套空间自回归框架检验环境驱动
#   4. 在群落尺度上检验 Rensch 法则
#
# 输入数据 / Input:
#   02_data_raw/traits/china_wang/
#     China_amphibian_traits_Song2022.xlsx  雄性/雌性体长
#     China_snakes_traits_Wang2023.xls      雄性/雌性标准体长（双层表头，skip = 1）
#     02_data_raw/huang2023_repo/ave/Chinesebirdsdata/Chinesebirdsdata.csv
#                                           雌雄体重、体长、嘴峰、翅长、尾长、跗蹠
#
# 预期输出 / Expected output:
#   03_data_derived/ssd_species.rds
#   04_results/tables/table_31_ssd_coverage.csv
#   04_results/tables/table_32_ssd_drivers.csv
#   05_figures/Fig15_sexual_dimorphism.pdf/.png
#
# 关键假设 / Key assumptions:
#   - SSD 指数用 log10(雌/雄)：以 0 为对称中心，正值表示雌性更大，
#     且对「哪一性别作分母」不敏感（比值指数没有这个性质）。
#   - 鸟类以体重为首选量度（体重比线性量度更能反映能量投入），
#     缺失时退回体长。
#   - 蛇类用标准体长而非最大体长：最大值是极端观测，样本量依赖强。
#   - 群落加权在存在-缺失数据下即等权平均，与主分析一致。
#
# 主要包 / Main packages: readxl, dplyr, sf, spdep, spatialreg
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(readxl); library(tidyr); library(dplyr)
  library(spdep); library(spatialreg); library(matrixStats)
})

log_msg("=== 16 性二态性 / Sexual size dimorphism ===")
BM  <- load_basemap(PATH$derived)
FIG <- PATH$figures
WD  <- file.path(PATH$raw_trait, "china_wang")

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

#' SSD 指数：log10(雌/雄)
#' 以 0 对称，正值 = 雌性更大；换分母只改变符号而不改变量级，
#' 而 (大/小 - 1) 型指数不具备这一性质。
ssd_index <- function(female, male) {
  ok <- is.finite(female) & is.finite(male) & female > 0 & male > 0
  out <- rep(NA_real_, length(female))
  out[ok] <- log10(female[ok] / male[ok])
  out
}

# ---------------------------------------------------------------
# 1. 逐类群计算物种级 SSD / Species-level SSD per class
# ---------------------------------------------------------------
ssd <- list()

# 两栖类 —— 体长 / amphibians, body length
f <- file.path(WD, "China_amphibian_traits_Song2022.xlsx")
if (file.exists(f)) {
  A <- read_excel(f, sheet = 1)
  ssd$Amphibia <- data.frame(
    species = norm_name(A[["种拉丁名"]]), class = "Amphibia",
    measure = "body length",
    ssd = ssd_index(num(A[["雌性体长/mm"]]), num(A[["雄性体长/mm"]])),
    size = num(A[["雄性体长/mm"]]), stringsAsFactors = FALSE)
}

# 蛇类 —— 标准体长（双层表头）/ snakes, standard length; two-row header
f <- file.path(WD, "China_snakes_traits_Wang2023.xls")
if (file.exists(f)) {
  S <- read_excel(f, sheet = 1, skip = 1)
  sci <- grep("Scientific name", names(S), value = TRUE)[1]
  ml  <- grep("Standard length in mm \\(male\\)",   names(S), value = TRUE)[1]
  fl  <- grep("Standard length in mm \\(female\\)", names(S), value = TRUE)[1]
  if (!is.na(sci) && !is.na(ml) && !is.na(fl)) {
    ssd$Reptilia <- data.frame(
      species = norm_name(gsub("[*+]", "", S[[sci]])), class = "Reptilia",
      measure = "standard length",
      ssd = ssd_index(num(S[[fl]]), num(S[[ml]])),
      size = num(S[[ml]]), stringsAsFactors = FALSE)
  }
}

# 鸟类 —— 体重优先，缺失时退回体长 / birds, mass first then size
f <- file.path(PATH$raw, "huang2023_repo", "ave", "Chinesebirdsdata",
               "Chinesebirdsdata.csv")
if (file.exists(f)) {
  B <- data.table::fread(f)
  mid <- function(a, b) { a <- num(B[[a]]); b <- num(B[[b]])
    ifelse(is.finite(a) & is.finite(b), (a + b) / 2, ifelse(is.finite(a), a, b)) }
  mm <- mid("male_mass_min_g", "male_mass_max_g")
  fm <- mid("female_mass_min_g", "female_mass_max_g")
  ms <- mid("male_size_min_mm", "male_size_max_mm")
  fs <- mid("female_size_min_mm", "female_size_max_mm")
  s_mass <- ssd_index(fm, mm); s_size <- ssd_index(fs, ms)
  ssd$Aves <- data.frame(
    species = norm_name(B[["iucn_2021_binomial"]]), class = "Aves",
    measure = ifelse(is.finite(s_mass), "body mass", "body length"),
    ssd = ifelse(is.finite(s_mass), s_mass, s_size),
    size = ifelse(is.finite(mm), mm, ms), stringsAsFactors = FALSE)
}

SSD <- bind_rows(ssd) |>
  filter(!is.na(species), species != "", is.finite(ssd)) |>
  distinct(species, .keep_all = TRUE)

# 爬行类的分性别量度只来自蛇类数据集——蜥蜴数据集不分性别记录体型。
# 因此下文的「Reptilia」实为蛇类，图表标注必须如实反映。
# The reptile SSD comes entirely from the snake dataset; the lizard compilation
# does not separate the sexes. Label it as snakes, not reptiles.
LAB_SSD <- c(Aves = "Birds", Reptilia = "Snakes", Amphibia = "Amphibians")
saveRDS(SSD, file.path(PATH$derived, "ssd_species.rds"))
log_msg("有 SSD 的物种 / species with SSD: ", nrow(SSD))

# 覆盖率 / coverage against the analysed species pool
cov <- lapply(names(ssd), function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  pool <- unique(readRDS(f)$species)
  s <- SSD |> filter(class == cl)
  data.frame(class = cl, n_in_dataset = nrow(s),
             n_analysed = length(pool),
             n_matched = sum(pool %in% s$species),
             pct = round(100 * mean(pool %in% s$species), 1),
             mean_ssd = round(mean(s$ssd, na.rm = TRUE), 4),
             pct_female_larger = round(100 * mean(s$ssd > 0, na.rm = TRUE), 1))
}) |> bind_rows()
write_table(cov, "table_31_ssd_coverage")
print(as.data.frame(cov))

# ---------------------------------------------------------------
# 2. 群落加权 SSD / Assemblage-level SSD
# ---------------------------------------------------------------
g50 <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE) |>
  drop_inset("50km")
env <- readRDS(file.path(PATH$derived, "model_data.rds"))[["50km"]] |>
  distinct(cell_id, .keep_all = TRUE) |> drop_inset("50km")
AX <- c("ax_thermal", "ax_water", "ax_productivity", "ax_structure", "ax_human")

cw <- lapply(names(ssd), function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  comm <- drop_inset(readRDS(f), "50km")
  s <- SSD |> filter(class == cl) |> select(species, ssd, size)
  d <- comm |> inner_join(s, by = "species") |>
    group_by(cell_id) |>
    summarise(n_sp_ssd = n(),
              cwm_ssd  = mean(ssd, na.rm = TRUE),
              cwv_ssd  = stats::var(ssd, na.rm = TRUE),
              cwm_size = mean(log10(size), na.rm = TRUE), .groups = "drop") |>
    filter(n_sp_ssd >= 5)
  d$class <- cl; d
}) |> bind_rows()
saveRDS(cw, file.path(PATH$derived, "ssd_assemblage_50km.rds"))
log_msg("有群落 SSD 的网格 / cells with assemblage SSD: ",
        paste(sprintf("%s %d", cov$class,
                      sapply(cov$class, function(x) sum(cw$class == x))),
              collapse = "; "))

# ---------------------------------------------------------------
# 3. 环境驱动 / Environmental drivers, same SAR framework as the main analysis
# ---------------------------------------------------------------
fit_sar <- function(d, resp) {
  d <- d |> inner_join(env[, c("cell_id", "x_albers", "y_albers", AX)],
                       by = "cell_id")
  d$y <- d[[resp]]
  d <- d[stats::complete.cases(d[, c("y", AX, "x_albers", "y_albers")]), ]
  if (nrow(d) < 100) return(NULL)
  for (v in AX) d[[v]] <- zscore(d[[v]])
  d$y <- zscore(d$y)
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(d[, c("x_albers", "y_albers")]), k = 8))
  m <- try(spatialreg::errorsarlm(
    stats::as.formula(paste("y ~", paste(AX, collapse = " + "))),
    data = d, listw = spdep::nb2listw(nb, style = "W"), zero.policy = TRUE),
    silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$Coef
  data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2],
             p = s[, 4], n = nrow(d), row.names = NULL)
}

rows <- list()
for (cl in unique(cw$class)) for (r in c("cwm_ssd", "cwv_ssd")) {
  o <- fit_sar(cw |> filter(class == cl), r)
  if (!is.null(o)) { o$class <- cl; o$response <- r; rows[[length(rows) + 1]] <- o }
}
eff <- bind_rows(rows) |> filter(term %in% AX) |>
  left_join(TAXA[, c("class", "thermal")], by = "class")
eff$p_adj <- stats::p.adjust(eff$p, "BH")
write_table(eff, "table_32_ssd_drivers")
log_msg("--- SSD 的环境驱动（BH 校正后显著项）---")
print(eff |> filter(p_adj < 0.05) |>
        select(class, response, term, estimate, se, p_adj) |>
        mutate(across(where(is.numeric), ~ signif(.x, 3))) |> as.data.frame())

# ---------------------------------------------------------------
# 4. 群落尺度的 Rensch 法则 / Rensch's rule at the assemblage level
# ---------------------------------------------------------------
# Rensch 法则的常规检验在**物种**尺度。群落尺度的对应问题是群落平均 SSD
# 是否随群落平均体型变化——但这可能只是物种尺度关系的数学继承，因为群落
# 值就是物种值的平均。因此必须同时给出物种尺度的相关，二者的差距才是
# 「空间共现是否额外贡献了信息」的证据。
# The assemblage-level relationship could be inherited arithmetically from the
# species level, since assemblage values are species averages. Reporting the
# species-level correlation alongside is what makes the comparison informative.
sp_rensch <- SSD |> filter(is.finite(ssd), is.finite(size), size > 0) |>
  group_by(class) |>
  summarise(n_species = n(),
            species_slope = stats::coef(stats::lm(ssd ~ log10(size)))[2],
            species_r = stats::cor(ssd, log10(size)), .groups = "drop")

rensch <- lapply(unique(cw$class), function(cl) {
  d <- cw |> filter(class == cl) |>
    inner_join(env[, c("cell_id", "x_albers", "y_albers")], by = "cell_id")
  d <- d[stats::complete.cases(d[, c("cwm_ssd", "cwm_size", "x_albers", "y_albers")]), ]
  if (nrow(d) < 100) return(NULL)
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(d[, c("x_albers", "y_albers")]), k = 8))
  m <- try(spatialreg::errorsarlm(zscore(cwm_ssd) ~ zscore(cwm_size), data = d,
                                  listw = spdep::nb2listw(nb, style = "W"),
                                  zero.policy = TRUE), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$Coef
  data.frame(class = cl, n_cells = nrow(d),
             assemblage_slope = s[2, 1], se = s[2, 2], p = s[2, 4],
             assemblage_r = stats::cor(d$cwm_ssd, d$cwm_size, method = "spearman"))
}) |> bind_rows()
if (nrow(rensch)) {
  rensch <- rensch |> left_join(sp_rensch, by = "class") |>
    mutate(amplification = round(abs(assemblage_r) / pmax(abs(species_r), 1e-9), 2))
  rensch$p_adj <- stats::p.adjust(rensch$p, "BH")
  write_table(rensch, "table_33_rensch_assemblage")
  log_msg("--- 群落尺度的 Rensch 法则 ---")
  print(rensch |> mutate(across(where(is.numeric), ~ signif(.x, 3))) |> as.data.frame())
}

# ---------------------------------------------------------------
# 5. 制图 / Figure
# ---------------------------------------------------------------
join_grid <- function(d) g50 |> inner_join(d, by = "cell_id")
# 覆盖率：鸟 94.3%、两栖 62.7%、蛇 16.0%。蛇类的分性别量度只有 83 个物种
# 进入分析（且蜥蜴数据集不分性别），故蛇类结果的解释力弱于另外两类。
# Coverage is 94.3% for birds, 62.7% for amphibians and 16.0% for snakes;
# the snake result therefore carries less weight than the other two.
CLS <- intersect(c("Aves", "Reptilia", "Amphibia"), unique(cw$class))
LIM <- c(-1, 1) * max(abs(stats::quantile(cw$cwm_ssd, c(0.02, 0.98), na.rm = TRUE)))

maps <- lapply(CLS, function(cl) {
  d <- cw |> filter(class == cl) |> select(cell_id, cwm_ssd)
  map_china(join_grid(d), "cwm_ssd", BM, "div", title = LAB_SSD[cl],
            legend = expression("Assemblage mean SSD, "*log[10]*"(female / male)"),
            limits = LIM, inset = FALSE)
})
for (i in seq_along(maps)) maps[[i]] <- maps[[i]] + labs(tag = letters[i])
row1 <- wrap_plots(maps, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.key.width = unit(34, "pt"),
        legend.key.height = unit(3.5, "pt"))

# (d) 物种级 SSD 的分布 / species-level SSD distribution
pd <- ggplot(SSD |> filter(class %in% CLS), aes(ssd, fill = class)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, colour = "grey45") +
  geom_density(alpha = 0.45, linewidth = 0.2) +
  scale_fill_manual(values = PAL$taxa, labels = LAB_SSD, name = NULL) +
  labs(x = expression(log[10]*"(female / male)"), y = "Density", tag = "d",
       subtitle = sprintf("Right of zero: females larger.\nCoverage: birds %.0f%%, amphibians %.0f%%, snakes %.0f%%",
                          cov$pct[cov$class == "Aves"], cov$pct[cov$class == "Amphibia"],
                          cov$pct[cov$class == "Reptilia"])) +
  theme_pub() + theme(legend.position = c(0.98, 0.98),
                      legend.justification = c(1, 1),
                      plot.subtitle = element_text(size = 4.6, colour = "grey30",
                                                  lineheight = 1.1))

# (e) 环境效应 / environmental effects on assemblage SSD
PRED_LAB <- c(ax_thermal = "Thermal energy", ax_water = "Water availability",
              ax_productivity = "Productivity", ax_structure = "Habitat structure",
              ax_human = "Human pressure")
ef <- eff |> filter(response == "cwm_ssd")
ef$term_f <- factor(PRED_LAB[ef$term], levels = rev(unname(PRED_LAB)))
pe <- ggplot(ef, aes(estimate, term_f, colour = class)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                  position = position_dodge(width = 0.6), size = 0.25,
                  linewidth = 0.35) +
  scale_colour_manual(values = PAL$taxa, labels = LAB_SSD, name = NULL) +
  labs(x = "Standardised effect on assemblage mean SSD", y = NULL, tag = "e",
       caption = "Mammals absent: no source records mammal measurements by sex") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.caption = element_text(size = 4.6, colour = "grey40"))

# (f) 群落尺度 Rensch 法则 / assemblage-level Rensch
# 三个类群的体型单位不同（鸟为克，蛇与两栖为毫米），共用横轴会造成误读；
# 因此在每个类群内部对体型做 z 标准化后再叠加。
# The three classes use different size units, so a shared axis is misleading;
# size is z-standardised within class before plotting.
cwz <- cw |> filter(class %in% CLS) |> group_by(class) |>
  mutate(size_z = zscore(cwm_size)) |> ungroup()
pf <- ggplot(cwz, aes(size_z, cwm_ssd, colour = class)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, colour = "grey45") +
  geom_point(size = 0.2, alpha = 0.12) +
  geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.6, se = FALSE) +
  scale_colour_manual(values = PAL$taxa, guide = "none") +
  labs(x = "Assemblage mean body size (z, within class)",
       y = "Assemblage mean SSD", tag = "f",
       subtitle = if (nrow(rensch))
         paste0("Assemblage slope (species-level r)\n",
                paste(sprintf("%s %+.2f%s (%+.2f)", LAB_SSD[rensch$class],
                              rensch$assemblage_slope,
                              ifelse(rensch$p_adj < 0.05, "*", ""),
                              rensch$species_r),
                      collapse = "   ")) else NULL) +
  theme_pub() + theme(plot.subtitle = element_text(size = 4.6, colour = "grey30",
                                                  lineheight = 1.1))

p <- row1 / (pd | pe | pf) + plot_layout(heights = c(1.25, 1))
save_fig(p, "Fig15_sexual_dimorphism", W2, 130, FIG)

log_msg("=== 16 完成 / done ===")
