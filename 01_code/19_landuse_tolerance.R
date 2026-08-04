# ============================================================
# 科学问题 / Scientific question:
#   被系统性保护不足的性状端（慢生活史、窄生态位、大型爬行动物），
#   是否正是对土地利用最敏感的那一端？若是，保护缺口就不只是「不均衡」，
#   而是**恰好落在最需要保护的地方**。
#   Are the trait ends we found under-protected also the ones least tolerant
#   of land use? If so, the gap is not merely uneven — it falls exactly where
#   protection matters most.
#
# 为什么现在能问这个 / Why this is answerable now:
#   Etard & Newbold (2024, Conservation Biology, doi:10.1111/cobi.14208) 论证
#   生态特征与物种对土地利用的响应相关。其底层性状表即本研究已整合的
#   Etard et al. (2020, GEB)，其中有一个变量我们此前未用：
#   **人工生境利用**（artificial habitat use）——物种是否被记录使用
#   人工生境，是土地利用耐受性在性状层面最直接的度量。
#   覆盖率：鸟 80.0%、兽 84.5%、两栖 67.4%、爬行 21.4%。
#   爬行类覆盖过低，只作描述不进入推断。
#
# 三个检验 / Three tests:
#   1. 土地利用耐受性本身与三个功能轴的关系（哪种性状的物种更耐受？）
#   2. 群落尺度：耐受性的地理格局，及其与人类足迹的关系（过滤是否可见？）
#   3. 关键检验：逐物种「受保护比例」与「土地利用耐受性」的关系。
#      若不耐受的物种同时受保护更少，则保护缺口与压力敏感性**同向叠加**。
#
# 输入 / Input:
#   03_data_derived/traits_imputed.rds（含 et_Artificial_habitat_use）
#   functional_axes_species.rds, comm_*_50km.rds, pa_coverage_50km.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_49_landuse_tolerance_traits.csv
#   04_results/tables/table_50_protection_vs_tolerance.csv
#   05_figures/Fig10_landuse_tolerance.pdf/.png
#
# 关键假设 / Key assumptions:
#   - 「人工生境利用」是二元记录，缺失并非随机：研究较少的物种更可能缺失，
#     而研究少往往与分布窄相关。因此所有模型都控制分布区大小，
#     且结论只在覆盖率 > 60% 的三个纲中给出。
#   - 该变量记录的是「是否**能**使用」而非「是否偏好」，因此它是耐受性的
#     宽松上界：被记为不使用人工生境的物种，几乎确定是不耐受的。
#
# 主要包 / Main packages: dplyr, Matrix
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(Matrix); library(sf)
})

log_msg("=== 19 土地利用耐受性 / Land-use tolerance ===")

AXES <- c("SIZE", "PACE", "NICHE")
AXLAB <- c(SIZE = "Body size", PACE = "Slow–fast pace", NICHE = "Niche breadth")
MIN_COV <- 60      # 覆盖率下限（%），低于此只作描述 / inference threshold

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))
FA <- readRDS(file.path(PATH$derived, "functional_axes_species.rds"))
attr_g <- readRDS(file.path(PATH$derived, "grid_50km_attr.rds"))
pa <- readRDS(file.path(PATH$derived, "pa_coverage_50km.rds"))

# ---------------------------------------------------------------
# 1. 提取并编码人工生境利用 / extract and code artificial habitat use
# ---------------------------------------------------------------
stopifnot("et_Artificial_habitat_use" %in% names(traits))
# 该变量在本项目的性状表中已是逻辑型（TRUE/FALSE/NA）。
# 仍写成宽容的解析，是因为上游若改回 0/1 或 Yes/No 也不会静默失效。
# The variable is logical in our trait table; the tolerant parser guards
# against an upstream change back to 0/1 or Yes/No.
tol <- traits |>
  distinct(species, .keep_all = TRUE) |>
  transmute(species, class, raw = et_Artificial_habitat_use)
tol$tolerant <- if (is.logical(tol$raw)) as.numeric(tol$raw) else {
  r <- trimws(as.character(tol$raw))
  dplyr::case_when(grepl("^(1|yes|true|y)$", r, ignore.case = TRUE) ~ 1,
                   grepl("^(0|no|false|n)$", r, ignore.case = TRUE) ~ 0,
                   TRUE ~ NA_real_)
}
cov_tab <- tol |> group_by(class) |>
  summarise(n = n(), n_coded = sum(is.finite(tolerant)),
            pct_coded = round(100 * mean(is.finite(tolerant)), 1),
            pct_tolerant = round(100 * mean(tolerant, na.rm = TRUE), 1),
            .groups = "drop") |>
  mutate(used_for_inference = pct_coded >= MIN_COV)
write_table(cov_tab, "table_49b_tolerance_coverage")
log_msg("人工生境利用的覆盖率 / coverage of artificial habitat use:")
print(as.data.frame(cov_tab))
USE_CL <- cov_tab$class[cov_tab$used_for_inference]
if (!length(USE_CL)) { log_msg("无可用类群，脚本结束"); quit(save = "no") }
log_msg("进入推断的类群 / classes used for inference: ", paste(USE_CL, collapse = ", "))

SP <- FA |> select(species, class, all_of(AXES)) |>
  inner_join(tol |> select(species, tolerant), by = "species") |>
  filter(is.finite(tolerant))

# ---------------------------------------------------------------
# 2. 检验一：耐受性与三个功能轴 / tolerance vs the three axes
# ---------------------------------------------------------------
# 逻辑回归：耐受(1/0) ~ 三轴。控制分布区大小由 NICHE 部分承担，
# 但仍单列 log(分布格数) 以免把「广布」误读为「耐受」。
comm <- bind_rows(lapply(TAXA$class, function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  drop_inset(readRDS(f), "50km") |> mutate(class = cl)
}))
ncell <- comm |> count(species, name = "n_cells")
SP <- SP |> left_join(ncell, by = "species") |> filter(is.finite(n_cells), n_cells > 0)

t1 <- bind_rows(lapply(USE_CL, function(cl) {
  d <- SP |> filter(class == cl)
  if (nrow(d) < 50 || length(unique(d$tolerant)) < 2) return(NULL)
  m <- stats::glm(tolerant ~ SIZE + PACE + NICHE + log10(n_cells),
                  family = stats::binomial(), data = d)
  s <- summary(m)$coefficients
  data.frame(class = cl, term = rownames(s), estimate = s[, 1], se = s[, 2],
             p = s[, 4], n = nrow(d), row.names = NULL)
})) |> filter(term %in% AXES)
t1$p_adj <- stats::p.adjust(t1$p, "BH")
write_table(t1, "table_49_landuse_tolerance_traits")
log_msg("耐受性 ~ 功能轴 / tolerance vs axes:")
print(t1 |> mutate(across(c(estimate, se), ~ round(.x, 3)),
                   p_adj = signif(p_adj, 2)) |> as.data.frame())

# ---------------------------------------------------------------
# 3. 检验二：群落尺度的耐受性与人类足迹 / assemblage tolerance vs footprint
# ---------------------------------------------------------------
env <- readRDS(file.path(PATH$derived, "env_50km.rds")) |> drop_inset("50km")
asm_tol <- comm |> filter(class %in% USE_CL) |>
  inner_join(SP |> select(species, tolerant), by = "species") |>
  group_by(class, cell_id) |>
  summarise(n_sp = n(), frac_tolerant = mean(tolerant), .groups = "drop") |>
  filter(n_sp >= 5) |>
  inner_join(env, by = "cell_id")
t2 <- bind_rows(lapply(USE_CL, function(cl) {
  d <- asm_tol |> filter(class == cl)
  if (nrow(d) < 100) return(NULL)
  r <- stats::cor.test(d$frac_tolerant, d$hfp2020, method = "spearman",
                       exact = FALSE)
  data.frame(class = cl, n_cells = nrow(d),
             mean_frac_tolerant = round(mean(d$frac_tolerant), 3),
             rho_with_footprint = round(unname(r$estimate), 3),
             p = signif(r$p.value, 3))
}))
write_table(t2, "table_49c_assemblage_tolerance")
log_msg("群落耐受比例 与 人类足迹 / assemblage tolerance vs human footprint:")
print(as.data.frame(t2))

# ---------------------------------------------------------------
# 4. 检验三（关键）：受保护比例 与 土地利用耐受性
# Key test: protection versus land-use tolerance
# ---------------------------------------------------------------
keep_cells <- attr_g$cell_id[!attr_g$scs_inset]
area <- setNames(attr_g$land_area_km2, attr_g$cell_id)[keep_cells]
area[!is.finite(area) | area < 0] <- 0
pfrac <- setNames(pa$pa_frac, pa$cell_id)[keep_cells]
pfrac[!is.finite(pfrac)] <- 0

cc <- comm |> filter(cell_id %in% keep_cells, species %in% SP$species)
spp <- sort(unique(cc$species)); cells <- keep_cells
M <- sparseMatrix(i = match(cc$species, spp), j = match(cc$cell_id, cells), x = 1,
                  dims = c(length(spp), length(cells)))
prot_km2  <- as.numeric(M %*% (pfrac * area))
range_km2 <- as.numeric(M %*% area)
# SP 已带 n_cells（全国占据格数）；此处矩阵算出的是**可分析格内**的占据数，
# 两者含义不同，join 前先去掉 SP 的那一列，避免出现 n_cells.x / n_cells.y。
PS <- data.frame(species = spp, prot_km2 = prot_km2, range_km2 = range_km2,
                 n_cells = as.numeric(Matrix::rowSums(M))) |>
  inner_join(SP |> select(-n_cells), by = "species") |>
  filter(range_km2 > 0, n_cells > 0)

t3 <- bind_rows(lapply(USE_CL, function(cl) {
  d <- PS |> filter(class == cl)
  if (nrow(d) < 50 || length(unique(d$tolerant)) < 2) return(NULL)
  # quasi-binomial：响应为受保护/未受保护面积，容许 0，按分布区加权
  m <- stats::glm(cbind(prot_km2, range_km2 - prot_km2) ~
                    tolerant + log10(n_cells),
                  family = stats::quasibinomial(), data = d)
  s <- summary(m)$coefficients
  # 同时给出两组的原始受保护比例，便于解读系数方向
  pr <- d |> mutate(p = prot_km2 / range_km2) |> group_by(tolerant) |>
    summarise(med = stats::median(p), .groups = "drop")
  data.frame(class = cl, n = nrow(d),
             n_intolerant = sum(d$tolerant == 0), n_tolerant = sum(d$tolerant == 1),
             median_prot_intolerant = round(100 * pr$med[pr$tolerant == 0], 2),
             median_prot_tolerant   = round(100 * pr$med[pr$tolerant == 1], 2),
             estimate = s["tolerant", 1], se = s["tolerant", 2],
             p = s["tolerant", 4], row.names = NULL)
}))
t3$p_adj <- stats::p.adjust(t3$p, "BH")
write_table(t3, "table_50_protection_vs_tolerance")
log_msg("受保护比例 ~ 土地利用耐受性 / protection vs land-use tolerance:")
print(t3 |> mutate(across(c(estimate, se), ~ round(.x, 3)),
                   p_adj = signif(p_adj, 3)) |> as.data.frame())

# ---------------------------------------------------------------
# 5. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures
t1$class_f <- factor(TAXA_LAB[t1$class], levels = unname(TAXA_LAB))
t1$axis_f  <- factor(AXLAB[t1$term], levels = unname(AXLAB))
pa1 <- ggplot(t1, aes(axis_f, estimate, fill = class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.78), width = 0.68) +
  geom_errorbar(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
                position = position_dodge(width = 0.78), width = 0.16,
                linewidth = 0.25) +
  geom_text(data = t1 |> filter(p_adj < 0.05),
            aes(y = estimate + sign(estimate) * (1.96 * se + 0.06), label = "*"),
            position = position_dodge(width = 0.78), size = 2.1,
            show.legend = FALSE) +
  scale_fill_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
  labs(x = NULL, y = "Effect on using artificial habitat (log-odds)", tag = "a",
       subtitle = "Which trait ends tolerate human-modified habitat? Range size controlled. * BH-adjusted P < 0.05.") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.5, colour = "grey30"))

# (b) 群落耐受比例 vs 人类足迹
asm_tol$class_f <- factor(TAXA_LAB[asm_tol$class], levels = unname(TAXA_LAB))
pb <- ggplot(asm_tol, aes(hfp2020, frac_tolerant, colour = class)) +
  geom_point(size = 0.18, alpha = 0.25) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
              se = FALSE, linewidth = 0.5) +
  scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
  labs(x = "Human footprint (2020)", y = "Fraction of the assemblage\ntolerating artificial habitat",
       tag = "b",
       subtitle = "The anthropogenic filter, read directly off a trait rather than inferred") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.5, colour = "grey30"))

# (c) 关键面板：耐受与不耐受物种的受保护比例
PS$class_f <- factor(TAXA_LAB[PS$class], levels = unname(TAXA_LAB))
PS$tol_f <- factor(ifelse(PS$tolerant == 1, "Tolerates artificial habitat",
                          "Does not"),
                   levels = c("Does not", "Tolerates artificial habitat"))
pc <- ggplot(PS |> filter(class %in% USE_CL),
             aes(class_f, 100 * prot_km2 / range_km2, fill = tol_f)) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.25, width = 0.66,
               position = position_dodge(width = 0.78)) +
  scale_fill_manual(values = c(`Does not` = "#B5541F",
                               `Tolerates artificial habitat` = "grey72"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 30), oob = scales::squish) +
  labs(x = NULL, y = "% of range protected", tag = "c",
       subtitle = "Do the species least able to live in modified habitat also receive less protection?") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.5, colour = "grey30"))

p <- pa1 / (pb | pc) + plot_layout(heights = c(1, 1))
save_fig(p, "Fig10_landuse_tolerance", W2, 150, FIG)

log_msg("=== 19 完成 / done ===")
