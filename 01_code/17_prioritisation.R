# ============================================================
# 科学问题 / Scientific question:
#   如果保护规划的目标是物种（现行做法），功能性状空间会被漏掉多少？
#   哪一类性状被系统性地保护不足？把目标换成性状，代价是多少土地？
#   If conservation planning targets species — as it does now — how much of
#   functional trait space is left out, which traits are systematically
#   under-protected, and what would it cost to close that gap?
#
# 为什么这是新问题 / Why this is new:
#   前序结果显示：(i) 变温类群受环境过滤更强、(ii) 鸟类功能多样性不能替代
#   爬行与两栖（Fig 6 一致性分析）。若两点成立，则一个以物种数（实际上由
#   鸟类主导）为依据的保护网络，应当在功能空间上出现**可预测的、有方向性的**
#   缺口——不是随机遗漏，而是系统性偏向某一端的性状。本脚本检验这一推论，
#   并给出补齐缺口的空间方案。这是两篇参考文献未涉及的部分。
#
# 分析目标 / Objective:
#   1. 互补性优先区选择：以物种为目标 vs 以性状空间为目标，比较两条积累曲线。
#   2. 现行保护区网络落在这条效率前沿的什么位置（对照随机与富集度贪心）。
#   3. 逐物种保护比例 ~ 功能轴：识别被系统性保护不足的性状端。
#   4. 昆明–蒙特利尔 30x30 目标下的缺口与补齐方案。
#
# 输入数据 / Input:
#   03_data_derived/comm_*_50km.rds, functional_axes_species.rds,
#   pa_coverage_50km.rds, grid_50km_attr.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_41 .. table_44
#   05_figures/Fig6_prioritisation.pdf/.png
#   03_data_derived/prioritisation_50km.rds
#
# 关键假设 / Key assumptions:
#   - 贪心互补性算法给出的是近似最优解，不是精确最优。对 maximal-coverage
#     问题贪心有 (1-1/e) 的最坏情况保证，且与整数规划解通常相差 <2%，
#     用于比较**目标之间**的差异是充分的。
#   - 保护比例按网格内保护区面积比例加权，未考虑保护区内的实际适宜生境，
#     因此是保护程度的**上界**；这一偏倚对所有物种同向，不影响性状间的比较。
#   - 现行网络与贪心方案在**同等面积**下比较，但现行保护地分散在大量网格中、
#     每格只占小比例，贪心则整格选取。因此二者的差距既反映选址原则的差异，
#     也反映集中度的差异；这正是「设计」与「历史积累」的区别，但不应读作
#     现行保护地全部错置。
#     The gap partly reflects that existing reserves are thinly spread across
#     many cells whereas the greedy takes whole cells; this is the difference
#     between design and historical accretion, not evidence that reserves are
#     simply misplaced.
#   - 性状空间用三个先验功能轴的三维网格离散化；分箱数经敏感性检验（表 41b）。
#
# 主要包 / Main packages: Matrix, dplyr, sf, ggplot2
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr); library(sf)
})

log_msg("=== 17 互补性优先区与性状缺口 / Complementarity prioritisation ===")

N_BIN  <- 6L    # 每个功能轴的分箱数 / bins per functional axis
N_RAND <- 99L   # 随机基准的重复次数 / random benchmark replicates
AXES   <- c("SIZE", "PACE", "NICHE")

FA   <- readRDS(file.path(PATH$derived, "functional_axes_species.rds"))
attr_g <- readRDS(file.path(PATH$derived, "grid_50km_attr.rds"))
pa   <- readRDS(file.path(PATH$derived, "pa_coverage_50km.rds"))

# 只用可分析格（排除南海附图框内的重复格）/ analysable cells only
keep_cells <- attr_g$cell_id[!attr_g$scs_inset]
attr_g <- attr_g |> filter(cell_id %in% keep_cells)
pa <- pa |> filter(cell_id %in% keep_cells)

comm <- bind_rows(lapply(TAXA$class, function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  readRDS(f) |> filter(cell_id %in% keep_cells) |>
    select(cell_id, species) |> mutate(class = cl)
}))
comm <- comm |> semi_join(FA, by = "species")   # 只保留有功能轴分数的物种
log_msg("规划矩阵 / planning matrix: ", length(unique(comm$cell_id)), " 格 x ",
        length(unique(comm$species)), " 种")

# ---------------------------------------------------------------
# 1. 构建稀疏的 物种 x 网格 矩阵 / sparse species-by-cell matrix
# ---------------------------------------------------------------
cells <- sort(unique(comm$cell_id))
spp   <- sort(unique(comm$species))
M <- sparseMatrix(i = match(comm$species, spp), j = match(comm$cell_id, cells),
                  x = 1, dims = c(length(spp), length(cells)),
                  dimnames = list(spp, as.character(cells)))
sp_class <- comm$class[match(spp, comm$species)]

area <- attr_g$land_area_km2[match(cells, attr_g$cell_id)]
area[!is.finite(area) | area <= 0] <- 0
tot_area <- sum(area)
pa_frac <- pa$pa_frac[match(cells, pa$cell_id)]
pa_frac[!is.finite(pa_frac)] <- 0

# ---------------------------------------------------------------
# 2. 性状空间的离散化 / discretising trait space
# ---------------------------------------------------------------
# 三个功能轴各切 N_BIN 段，一个物种落入一个三维格。一个网格「覆盖」了它所含
# 物种落入的全部三维格。选址即 maximal coverage：让被覆盖的三维格最多。
# 分箱在**纲内**进行，因为轴分数已在纲内标准化；再拼上纲标签，保证
# 「大型鸟」与「大型两栖」是不同的功能格，不会互相顶替。
bin_of <- FA |> filter(species %in% spp) |>
  group_by(class) |>
  mutate(across(all_of(AXES),
                ~ cut(.x, breaks = stats::quantile(.x, seq(0, 1, length.out = N_BIN + 1),
                                                   na.rm = TRUE),
                      include.lowest = TRUE, labels = FALSE))) |>
  ungroup() |>
  mutate(bin = paste(class, SIZE, PACE, NICHE, sep = "_")) |>
  select(species, class, bin)
bins <- sort(unique(bin_of$bin))
log_msg("功能格数 / functional bins occupied: ", length(bins),
        "（每轴 ", N_BIN, " 段，四纲分开）")

# 功能格 x 网格 的覆盖矩阵 / bin-by-cell coverage
b_idx <- bin_of$bin[match(spp, bin_of$species)]
S2B <- sparseMatrix(i = match(b_idx, bins), j = seq_along(spp), x = 1,
                    dims = c(length(bins), length(spp)))
BC <- (S2B %*% M) > 0     # 逻辑：该功能格是否出现在该网格
BC <- as(as(BC, "dMatrix"), "generalMatrix")

# ---------------------------------------------------------------
# 3. 目标制贪心互补性 / target-based greedy complementarity
# ---------------------------------------------------------------
# 为什么用「目标制」而非「至少一格代表」/ Why targets, not mere representation:
#   以「每个物种至少被一个网格代表」为目标时，1,200 个 50 km 网格（约 10%
#   国土）就能覆盖 100% 的物种与 100% 的功能格——所有方案都饱和，无法区分。
#   更重要的是，「代表一次」与物种能否长期存续无关（Rodrigues & Gaston 2001）。
#   昆明–蒙特利尔目标 3 问的是**分布区的 30%**，因此这里的特征目标设为
#   「该特征总分布面积的 30%」，贪心每步最大化「单位面积削减的缺口量」。
#   Representation targets saturate at ~10% of land and say nothing about
#   persistence; targets are therefore set at 30% of each feature's total area.
TARGET <- 0.30

#' 目标制贪心 / target-based greedy: maximise shortfall reduction per unit area
#' @param A 特征 x 网格 的**面积**稀疏矩阵 / feature-by-cell amount matrix
greedy_target <- function(A, n_step, area_w = area, target = TARGET) {
  tot   <- as.numeric(Matrix::rowSums(A))
  goal  <- target * tot
  short <- goal                                   # 初始缺口 = 目标
  s <- Matrix::summary(as(A, "dgCMatrix"))        # 三元组 (i, j, x)
  alive <- rep(TRUE, ncol(A))
  sel <- integer(0)
  for (k in seq_len(n_step)) {
    if (!nrow(s)) break
    # 每格的贡献 = sum_f min(该格提供的量, 该特征剩余缺口)
    contrib <- pmin(s$x, short[s$i])
    g <- vapply(split(contrib, factor(s$j, levels = seq_len(ncol(A)))), sum,
                numeric(1))
    g[!alive] <- 0
    eff <- g / pmax(area_w, 1)
    if (max(eff) <= 0) break
    j <- which.max(eff)
    got <- s[s$j == j, , drop = FALSE]
    short[got$i] <- pmax(short[got$i] - got$x, 0)
    alive[j] <- FALSE
    # 已达标的特征退出，矩阵随选择进行而收缩 / drop satisfied features
    s <- s[short[s$i] > 0 & s$j != j, , drop = FALSE]
    sel <- c(sel, j)
  }
  sel
}

# 面积矩阵：物种 x 网格、功能格 x 网格 / amount matrices
Ma <- M  %*% Matrix::Diagonal(x = area)
Ba <- BC %*% Matrix::Diagonal(x = area)

#' 沿选址顺序累计「达标特征比例」与「代表比例」
#' Accumulate both the share of features meeting the 30% target and the share
#' merely represented at least once.
accumulate <- function(ord) {
  step_of <- function(A) {
    goal <- TARGET * as.numeric(Matrix::rowSums(A))
    cum  <- rep(0, nrow(A))
    out  <- integer(length(ord)); met <- rep(FALSE, nrow(A))
    Ao <- as(A[, ord, drop = FALSE], "dgCMatrix")
    p <- Ao@p; xi <- Ao@i + 1L; xv <- Ao@x
    for (k in seq_along(ord)) {
      idx <- seq_len(p[k + 1] - p[k]) + p[k]
      if (length(idx)) {
        ii <- xi[idx]; cum[ii] <- cum[ii] + xv[idx]
        new <- ii[!met[ii] & cum[ii] >= goal[ii]]
        met[new] <- TRUE; out[k] <- length(new)
      }
    }
    cumsum(out)
  }
  rep_of <- function(X) {
    s <- Matrix::summary(as(X[, ord, drop = FALSE], "dgCMatrix"))
    s <- s[order(s$i, s$j), ]
    cumsum(tabulate(s$j[!duplicated(s$i)], nbins = length(ord)))
  }
  data.frame(step = seq_along(ord),
             area_pct = 100 * cumsum(area[ord]) / tot_area,
             species_pct = 100 * step_of(Ma) / nrow(Ma),
             bin_pct     = 100 * step_of(Ba) / nrow(Ba),
             species_rep_pct = 100 * rep_of(M) / nrow(M),
             bin_rep_pct     = 100 * rep_of(BC) / nrow(BC))
}

N_STEP <- min(length(cells), 1800L)
log_msg("目标制贪心，选择 ", N_STEP, " 格 / target-based greedy, ", N_STEP, " cells")

set.seed(SEED)
g_sp  <- list(order = greedy_target(Ma, N_STEP))   # 目标：物种分布区的 30%
g_bin <- list(order = greedy_target(Ba, N_STEP))   # 目标：功能格分布的 30%
ord_rich <- order(Matrix::colSums(M), decreasing = TRUE)[seq_len(N_STEP)]

curves <- bind_rows(
  accumulate(g_sp$order)  |> mutate(strategy = "Species complementarity"),
  accumulate(g_bin$order) |> mutate(strategy = "Trait-space complementarity"),
  accumulate(ord_rich)    |> mutate(strategy = "Species richness (non-complementary)"))

# 随机基准 / random benchmark
set.seed(SEED)
rnd <- bind_rows(lapply(seq_len(N_RAND), function(i)
  accumulate(sample.int(length(cells), N_STEP)) |> mutate(rep = i)))
rnd_band <- rnd |> group_by(step) |>
  summarise(area_pct = mean(area_pct),
            sp_lo = stats::quantile(species_pct, 0.025),
            sp_hi = stats::quantile(species_pct, 0.975),
            bin_lo = stats::quantile(bin_pct, 0.025),
            bin_hi = stats::quantile(bin_pct, 0.975), .groups = "drop")

# ---------------------------------------------------------------
# 4. 两个目标之间的效率损失 / efficiency loss between targets
# ---------------------------------------------------------------
# 在同一土地面积下，以物种为目标的方案覆盖了多少性状空间，反之亦然。
at_area <- function(df, a) df[which.min(abs(df$area_pct - a)), ]
loss <- bind_rows(lapply(c(5, 10, 17, 30), function(a) {
  s <- at_area(curves |> filter(strategy == "Species complementarity"), a)
  b <- at_area(curves |> filter(strategy == "Trait-space complementarity"), a)
  r <- at_area(curves |> filter(strategy == "Species richness (non-complementary)"), a)
  data.frame(area_target_pct = a,
             species_by_speciesplan = round(s$species_pct, 1),
             bins_by_speciesplan    = round(s$bin_pct, 1),
             bins_by_traitplan      = round(b$bin_pct, 1),
             species_by_traitplan   = round(b$species_pct, 1),
             species_by_richness    = round(r$species_pct, 1),
             bins_by_richness       = round(r$bin_pct, 1),
             trait_shortfall_of_speciesplan = round(b$bin_pct - s$bin_pct, 1))
}))
write_table(loss, "table_41_target_efficiency")
log_msg("不同目标下的覆盖效率 / coverage efficiency by target:")
print(as.data.frame(loss))

# 分箱数敏感性 / sensitivity to the number of bins
# 若结论随分箱数改变，则「性状空间缺口」只是离散化的产物，不是真实格局。
sens <- bind_rows(lapply(c(4L, 5L, 6L, 8L), function(nb) {
  bo <- FA |> filter(species %in% spp) |> group_by(class) |>
    mutate(across(all_of(AXES), ~ cut(.x,
      breaks = stats::quantile(.x, seq(0, 1, length.out = nb + 1), na.rm = TRUE),
      include.lowest = TRUE, labels = FALSE))) |> ungroup() |>
    mutate(bin = paste(class, SIZE, PACE, NICHE, sep = "_"))
  bb <- sort(unique(bo$bin))
  s2b <- sparseMatrix(i = match(bo$bin[match(spp, bo$species)], bb),
                      j = seq_along(spp), x = 1, dims = c(length(bb), length(spp)))
  bc <- as(as((s2b %*% M) > 0, "dMatrix"), "generalMatrix")
  ba <- bc %*% Matrix::Diagonal(x = area)
  goal <- TARGET * as.numeric(Matrix::rowSums(ba))
  met_at <- function(ord) {                 # 30% 国土时达标的功能格比例
    n <- which.min(abs(100 * cumsum(area[ord]) / tot_area - 30))
    100 * mean(as.numeric(Matrix::rowSums(ba[, ord[seq_len(n)], drop = FALSE])) >= goal)
  }
  gg <- greedy_target(ba, 1200L)
  data.frame(n_bins = nb, n_functional_bins = nrow(bc),
             bins_met_by_speciesplan = round(met_at(g_sp$order), 1),
             bins_met_by_traitplan   = round(met_at(gg), 1),
             bins_met_by_richness    = round(met_at(ord_rich), 1))
}))
write_table(sens, "table_41b_bin_sensitivity")
log_msg("分箱数敏感性 / sensitivity to bin number:"); print(as.data.frame(sens))

# ---------------------------------------------------------------
# 5. 现行保护区网络的位置 / where the existing PA network sits
# ---------------------------------------------------------------
# 现行网络的达标率按**实际保护比例**计算（而非把整格算作已保护），
# 因此与贪心方案在同一货币下比较：都是「特征分布面积的 30%」。
# The network's attainment uses the actual protected fraction of each cell, so
# it is measured in the same currency as the greedy solutions.
pa_area <- sum(pa_frac * area)
pa_pct  <- 100 * pa_area / tot_area
#' 在每格保护比例为 f 时，达到 30% 目标的特征比例
#' A 的元素已是面积，故 A %*% f 即各特征的受保护面积。
met_frac <- function(A, f) {
  goal <- TARGET * as.numeric(Matrix::rowSums(A))
  100 * mean(as.numeric(A %*% f) >= goal)
}
pa_now <- data.frame(
  strategy = "Existing reserve network", area_pct = pa_pct,
  species_pct = met_frac(Ma, pa_frac),
  bin_pct     = met_frac(Ba, pa_frac))
opt_here <- bind_rows(
  at_area(curves |> filter(strategy == "Species complementarity"), pa_pct),
  at_area(curves |> filter(strategy == "Trait-space complementarity"), pa_pct))
eff_tab <- data.frame(
  area_pct = round(pa_pct, 2),
  network_species_pct = round(pa_now$species_pct, 1),
  network_bin_pct = round(pa_now$bin_pct, 1),
  optimal_species_pct = round(max(opt_here$species_pct), 1),
  optimal_bin_pct = round(max(opt_here$bin_pct), 1),
  random_bin_median = round(stats::median(
    rnd$bin_pct[abs(rnd$area_pct - pa_pct) < 0.5], na.rm = TRUE), 1))
write_table(eff_tab, "table_42_network_efficiency")
log_msg("现行网络 vs 最优 / existing network vs optimal:")
print(as.data.frame(eff_tab))

# ---------------------------------------------------------------
# 6. 哪一端的性状被系统性保护不足 / which traits are under-protected
# ---------------------------------------------------------------
# 逐物种：受保护比例 = sum(分布格内 保护区面积) / sum(分布格 陆域面积)
prot_km2  <- as.numeric(M %*% (pa_frac * area))
range_km2 <- as.numeric(M %*% area)
SPP <- data.frame(species = spp, class = sp_class,
                  prot_km2 = prot_km2, range_km2 = range_km2,
                  prot = prot_km2 / pmax(range_km2, 1),
                  n_cells = as.numeric(Matrix::rowSums(M))) |>
  left_join(FA |> select(species, all_of(AXES)), by = "species")

# 每纲一个模型：受保护面积占比 ~ 三个功能轴，控制分布区大小。
# 用 quasi-binomial（响应为「受保护/未受保护面积」的比例）而非对 logit 后的
# 比例做普通回归：不少物种的受保护比例恰为 0，logit 前必须截断，截断值的
# 选择会直接左右系数。quasi-binomial 天然容许 0，并按分布区面积加权，
# 分布区大的物种因而获得应有的权重；quasi- 前缀吸收过度离散。
# Many species have exactly zero protection, so a logit-transformed OLS needs an
# arbitrary clip that drives the coefficients. A quasi-binomial GLM on the
# protected/unprotected area split admits zeros and weights by range extent.
# 分布区大小必须控制：窄域种天然更容易被完全覆盖或完全遗漏。
und <- bind_rows(lapply(TAXA$class, function(cl) {
  d <- SPP |> filter(class == cl, is.finite(prot), range_km2 > 0)
  if (nrow(d) < 50) return(NULL)
  m <- stats::glm(cbind(prot_km2, range_km2 - prot_km2) ~
                    SIZE + PACE + NICHE + log10(n_cells),
                  family = stats::quasibinomial(), data = d)
  s <- summary(m)$coefficients
  data.frame(class = cl, term = rownames(s), estimate = s[, 1], se = s[, 2],
             p = s[, 4], n = nrow(d), row.names = NULL)
})) |> filter(term %in% AXES)
und$p_adj <- stats::p.adjust(und$p, "BH")
und <- und |> left_join(TAXA[, c("class", "thermal")], by = "class")
write_table(und, "table_43_underprotected_traits")
log_msg("受保护比例 ~ 功能轴 / protection vs functional axes:")
print(und |> mutate(across(c(estimate, se), ~ round(.x, 3)),
                    p_adj = signif(p_adj, 2)) |> as.data.frame())

# ---------------------------------------------------------------
# 7. 30x30 目标下的缺口与补齐 / the 30x30 shortfall and how to close it
# ---------------------------------------------------------------
# 昆明–蒙特利尔目标 3：到 2030 年有效保护 30% 陆域。这里问的是：
# 若把 30% 的国土按不同原则布置，多少物种能达到「分布区 30% 受保护」，
# 多少功能格能达标？现状一行用实际保护比例，情景各行把选中格记为全保护。
shortfall <- function(ord) {
  n <- which.min(abs(100 * cumsum(area[ord]) / tot_area - 100 * TARGET))
  f <- rep(0, length(cells)); f[ord[seq_len(n)]] <- 1
  p <- as.numeric(Ma %*% f) / pmax(range_km2, 1)
  data.frame(pct_species_meeting_30 = round(100 * mean(p >= TARGET), 1),
             pct_bins_meeting_30 = round(met_frac(Ba, f), 1),
             median_protection = round(100 * stats::median(p), 1))
}
gbf <- bind_rows(
  cbind(scenario = "Current reserve network",
        data.frame(pct_species_meeting_30 = round(100 * mean(SPP$prot >= TARGET), 1),
                   pct_bins_meeting_30 = round(met_frac(Ba, pa_frac), 1),
                   median_protection = round(100 * stats::median(SPP$prot), 1))),
  cbind(scenario = "30% by species complementarity",  shortfall(g_sp$order)),
  cbind(scenario = "30% by trait-space complementarity", shortfall(g_bin$order)),
  cbind(scenario = "30% by richness ranking",         shortfall(ord_rich)))
# 逐纲的现状达标率 / current attainment by class
gbf_cl <- SPP |> group_by(class) |>
  summarise(n = n(), pct_meeting_30 = round(100 * mean(prot >= TARGET), 1),
            median_protection = round(100 * stats::median(prot), 1), .groups = "drop") |>
  left_join(TAXA[, c("class", "thermal")], by = "class")
write_table(gbf, "table_44_gbf_3030")
write_table(gbf_cl, "table_44b_gbf_by_class")
log_msg("30x30 情景 / 30x30 scenarios:"); print(as.data.frame(gbf))
log_msg("逐纲现状 / current attainment by class:"); print(as.data.frame(gbf_cl))

# 优先区图层 / priority layer for mapping
n30 <- which.min(abs(100 * cumsum(area[g_bin$order]) / tot_area - 30))
prio <- data.frame(cell_id = cells,
                   rank_trait = match(seq_along(cells), g_bin$order),
                   rank_species = match(seq_along(cells), g_sp$order),
                   pa_frac = pa_frac, area = area)
prio$selected_trait <- prio$rank_trait <= n30 & !is.na(prio$rank_trait)
n30s <- which.min(abs(100 * cumsum(area[g_sp$order]) / tot_area - 30))
prio$selected_species <- prio$rank_species <= n30s & !is.na(prio$rank_species)
prio$category <- with(prio, ifelse(selected_trait & selected_species, "Both targets",
                            ifelse(selected_trait, "Trait target only",
                            ifelse(selected_species, "Species target only", "Neither"))))
saveRDS(prio, file.path(PATH$derived, "prioritisation_50km.rds"))
log_msg("两种目标的空间重合 / spatial overlap of the two targets:")
print(prio |> filter(category != "Neither") |> count(category) |> as.data.frame())

# ---------------------------------------------------------------
# 8. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures
PAL_ST <- c(`Species complementarity` = "#2E6E8E",
            `Trait-space complementarity` = "#B5541F",
            `Species richness (non-complementary)` = "grey55")

# (a) 性状空间积累曲线
pa1 <- ggplot(curves, aes(area_pct, bin_pct, colour = strategy)) +
  geom_ribbon(data = rnd_band, inherit.aes = FALSE,
              aes(area_pct, ymin = bin_lo, ymax = bin_hi),
              fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.55) +
  geom_point(data = pa_now, aes(area_pct, bin_pct), inherit.aes = FALSE,
             size = 1.6, shape = 21, fill = "#D9A441", colour = "grey15",
             stroke = 0.35) +
  annotate("text", x = pa_now$area_pct + 1.2, y = pa_now$bin_pct - 6,
           label = "Existing\nreserves", size = 1.75, hjust = 0, colour = "grey20",
           lineheight = 0.95) +
  geom_vline(xintercept = 30, linetype = "22", linewidth = 0.28, colour = "grey45") +
  annotate("text", x = 30.6, y = 8, label = "GBF 30%", size = 1.75,
           hjust = 0, colour = "grey35") +
  scale_colour_manual(values = PAL_ST, name = NULL) +
  coord_cartesian(xlim = c(0, 45)) +
  labs(x = "Land area selected (%)", y = "Functional bins meeting the 30% target (%)",
       tag = "a",
       subtitle = "At a tight budget (10–17% of land) the trait target covers 10–12 more percentage points of functional space than the species target; by 30% they converge. Grey band: 95% of random selections.") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.6, colour = "grey30")) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

# (b) 物种积累曲线
pa2 <- ggplot(curves, aes(area_pct, species_pct, colour = strategy)) +
  geom_ribbon(data = rnd_band, inherit.aes = FALSE,
              aes(area_pct, ymin = sp_lo, ymax = sp_hi),
              fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.55) +
  geom_point(data = pa_now, aes(area_pct, species_pct), inherit.aes = FALSE,
             size = 1.6, shape = 21, fill = "#D9A441", colour = "grey15",
             stroke = 0.35) +
  geom_vline(xintercept = 30, linetype = "22", linewidth = 0.28, colour = "grey45") +
  scale_colour_manual(values = PAL_ST, name = NULL, guide = "none") +
  coord_cartesian(xlim = c(0, 45)) +
  labs(x = "Land area selected (%)", y = "Species meeting the 30% target (%)", tag = "b",
       subtitle = "Complementarity beats richness ranking by ~9 points at the 30% budget; the existing network reaches 3%.") +
  theme_pub() + theme(plot.subtitle = element_text(size = 4.6, colour = "grey30"))

# (c) 被系统性保护不足的性状端
und$class_f <- factor(TAXA_LAB[und$class], levels = unname(TAXA_LAB))
und$axis_f <- factor(und$term, levels = AXES,
                     labels = c("Body size", "Slow–fast pace", "Niche breadth"))
pa3 <- ggplot(und, aes(axis_f, estimate, fill = class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_errorbar(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
                position = position_dodge(width = 0.78), width = 0.18,
                linewidth = 0.25) +
  geom_text(data = und |> filter(p_adj < 0.05),
            aes(y = estimate + sign(estimate) * (1.96 * se + 0.035), label = "*"),
            position = position_dodge(width = 0.78), size = 2.1,
            vjust = ifelse((und |> filter(p_adj < 0.05))$estimate > 0, 0.1, 0.85),
            show.legend = FALSE) +
  scale_fill_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
  labs(x = NULL, y = "Effect on protected fraction (logit scale)", tag = "c",
       subtitle = "Negative = the high end of that axis is under-protected. Quasi-binomial GLM weighted by range area, controlling range size. * BH-adjusted P < 0.05.") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.4, colour = "grey30"))

# (d) 30x30 情景
gg <- gbf |> mutate(scenario = factor(scenario, levels = rev(scenario)))
pa4 <- ggplot(gg, aes(scenario, pct_species_meeting_30)) +
  geom_col(fill = "#2E6E8E", width = 0.62) +
  geom_text(aes(label = sprintf("%.0f%%", pct_species_meeting_30)),
            hjust = -0.18, size = 1.9, colour = "grey20") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Species with ≥30% of range protected (%)", tag = "d",
       subtitle = "Kunming–Montreal Target 3 assessed per species") +
  theme_pub() + theme(axis.text.y = element_text(size = 5),
                      plot.subtitle = element_text(size = 4.6, colour = "grey30"))

# (e) 优先区地图
gsf <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE) |>
  filter(cell_id %in% cells) |>
  left_join(prio[, c("cell_id", "category")], by = "cell_id")
gsf$category <- factor(gsf$category,
                       levels = c("Both targets", "Trait target only",
                                  "Species target only", "Neither"))
BM <- load_basemap(PATH$derived)
pa5 <- map_china_cat(gsf, "category", BM,
                     values = c(`Both targets` = "#1F4E5F",
                                `Trait target only` = "#B5541F",
                                `Species target only` = "#7FA8BC",
                                Neither = "grey90"),
                     legend = NULL,
                     title = "Priority areas under the two targets (30% of land)") +
  labs(tag = "e")

p <- (pa1 | pa2) / (pa3 | pa4) / pa5 + plot_layout(heights = c(1, 0.92, 1.5))
save_fig(p, "Fig6_prioritisation", W2, 210, FIG)

log_msg("=== 17 完成 / done ===")
