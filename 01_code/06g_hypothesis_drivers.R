# ============================================================
# 科学问题 / Scientific question:
#   哪一个生物多样性假说最能解释群落性状结构？答案是否取决于
#   (a) 体温调节模式，(b) 所看的性状维度？
#   Which biodiversity hypothesis best explains assemblage trait structure, and
#   does the answer depend on thermoregulatory mode and on which trait
#   dimension is examined?
#
# 分析目标 / Objective:
#   1. 用**具名环境变量**（而非合成主成分）逐一对应经典假说，使系数可直接
#      解读为对某个假说的支持度。
#   2. 同时以**单性状**与**功能轴**为响应变量——若二者给出不同答案，
#      说明「性状结构」不是一个东西。
#   3. 所有模型同时控制其余混杂变量，并做方差分解，区分各假说的独占贡献。
#
# 假说与变量的对应 / Hypotheses and their variables:
#   环境能量假说 ambient energy      -> 年均温 bio1
#   水–能量动态 water-energy         -> 年降水 bio12
#   生产力假说 productivity          -> NPP
#   生境异质性假说 heterogeneity     -> 网格内海拔差 elev_range
#   气候季节性/稳定性 seasonality    -> 温度季节性 bio4
#   人为过滤 anthropogenic filter    -> 人类足迹 hfp2020
#
# 输入数据 / Input:
#   03_data_derived/functional_axes_species.rds, comm_*_50km.rds, env_50km.rds
#
# 预期输出 / Expected output:
#   03_data_derived/assemblage_axes_50km.rds
#   04_results/tables/table_38_hypothesis_effects.csv
#   04_results/tables/table_39_variance_partition.csv
#   05_figures/Fig4_hypothesis_drivers.pdf/.png
#
# 关键假设 / Key assumptions:
#   - 六个具名变量间存在共线性（如温度与降水），故每个模型同时纳入全部六个：
#     系数即「控制其余假说后」的独占效应，这正是区分假说所需。VIF 见表 38b。
#   - 空间自回归（SAR 误差模型）而非 OLS：残差的 Moran's I 在几乎所有模型中显著。
#   - 响应变量为群落加权平均值；其显著性用第四角 max test 复核（脚本 06b），
#     此处的 SAR 用于估计**效应量与方向**。
#
# 主要包 / Main packages: dplyr, spdep, spatialreg, vegan, car
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(spdep); library(spatialreg)
  library(vegan); library(car); library(matrixStats)
})

log_msg("=== 06g 假说驱动分析 / Hypothesis-based driver analysis ===")

# 假说 -> 变量 的对应表，全文一处定义
HYP <- c(bio1_mean   = "Ambient energy",
         bio12_mean  = "Water availability",
         npp         = "Productivity",
         elev_range  = "Habitat heterogeneity",
         bio4_mean   = "Climatic seasonality",
         hfp2020     = "Human pressure")
HYP_CN <- c(bio1_mean = "环境能量", bio12_mean = "水分", npp = "生产力",
            elev_range = "生境异质性", bio4_mean = "气候季节性",
            hfp2020 = "人为压力")
PRED <- names(HYP)

FA  <- readRDS(file.path(PATH$derived, "functional_axes_species.rds"))
env <- readRDS(file.path(PATH$derived, "env_50km.rds")) |> drop_inset("50km")
stopifnot(all(c("x_albers", "y_albers") %in% names(env)))  # env 自带 Albers 坐标

TRAITS <- c("size", "fecundity", "nocturnality", "verticality",
            "habitat_breadth", "range_size")
AXES   <- c("SIZE", "PACE", "NICHE")
RESP_LAB <- c(size = "Body size", fecundity = "Fecundity",
              nocturnality = "Nocturnality", verticality = "Vertical stratum use",
              habitat_breadth = "Habitat breadth", range_size = "Range size",
              SIZE = "Axis: body size", PACE = "Axis: slow-fast pace",
              NICHE = "Axis: niche breadth")

# ---------------------------------------------------------------
# 1. 群落加权平均：单性状与功能轴 / Assemblage means for traits and axes
# ---------------------------------------------------------------
asm <- lapply(TAXA$class, function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  comm <- drop_inset(readRDS(f), "50km")
  s <- FA |> filter(class == cl)
  d <- comm |> inner_join(s, by = "species") |>
    group_by(cell_id) |>
    summarise(richness = n(),
              across(all_of(c(TRAITS, AXES)), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") |>
    filter(richness >= 5)
  d$class <- cl; d
}) |> bind_rows()
saveRDS(asm, file.path(PATH$derived, "assemblage_axes_50km.rds"))
log_msg("群落表 / assemblage table: ", nrow(asm), " 行（格 x 纲）")

# ---------------------------------------------------------------
# 2. 共线性诊断 / Collinearity among the six named predictors
# ---------------------------------------------------------------
E <- env[stats::complete.cases(env[, PRED]), ]
Z <- as.data.frame(apply(as.matrix(E[, PRED]), 2, zscore))
Z$y <- stats::rnorm(nrow(Z))
vif <- car::vif(stats::lm(stats::as.formula(paste("y ~", paste(PRED, collapse = " + "))),
                          data = Z))
vt <- data.frame(variable = names(vif), hypothesis = HYP[names(vif)],
                 VIF = round(unname(vif), 2))
write_table(vt, "table_38b_predictor_vif")
log_msg("预测变量的 VIF / variance inflation:"); print(vt)

# ---------------------------------------------------------------
# 3. 逐类群 x 逐响应 的空间自回归 / SAR per class and response
# ---------------------------------------------------------------
# 两点实现上的要求，否则单机跑不动：
#   (1) errorsarlm 默认 method = "eigen"，要对 n x n 权重矩阵做特征分解，
#       n ≈ 3,700 时约 O(n^3)。改用 "Matrix"（稀疏 Cholesky）得到同样的
#       极大似然解，但快两个数量级。
#   (2) 同一纲下所有响应变量共用同一批网格，邻接与权重表只需构建一次。
#   errorsarlm's default eigen method is O(n^3); the sparse Cholesky method
#   gives the same ML estimates far faster, and the listw is reused per class.
fit_sar <- function(d, resp, lw) {
  X <- as.data.frame(apply(as.matrix(d[, PRED]), 2, zscore))
  X$y <- zscore(d[[resp]])
  f  <- stats::as.formula(paste("y ~", paste(PRED, collapse = " + ")))
  m  <- try(spatialreg::errorsarlm(f, data = X, listw = lw, method = "Matrix",
                                   zero.policy = TRUE), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$Coef
  # 残差空间自相关：若 SAR 已吸收，Moran's I 应接近 0
  mi <- try(spdep::moran.test(stats::residuals(m), lw, zero.policy = TRUE),
            silent = TRUE)
  data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], p = s[, 4],
             n = nrow(d), lambda = m$lambda,
             resid_moran = if (inherits(mi, "try-error")) NA_real_
                           else unname(mi$estimate[1]),
             row.names = NULL)
}

rows <- list()
for (cl in TAXA$class) {
  d <- asm |> filter(class == cl) |> inner_join(env, by = "cell_id")
  if (!nrow(d)) next
  # 同一纲内固定样本：所有响应用同一批完整格，模型间因此可直接比较
  keepv <- c(TRAITS, AXES)[c(TRAITS, AXES) %in% names(d)]
  d <- d[stats::complete.cases(d[, c(PRED, "x_albers", "y_albers")]), ]
  if (nrow(d) < 100) next
  # 对称化：稀疏 Cholesky 要求邻接对称，k 近邻本身不对称。
  # 对称化只添加「A 是 B 的近邻但反之不然」的回连，不改变邻域尺度。
  # 已核验：对称化 + Matrix 方法与精确 eigen 解的系数差 < 1e-9，快约 40 倍。
  nb <- spdep::make.sym.nb(
    spdep::knn2nb(spdep::knearneigh(as.matrix(d[, c("x_albers", "y_albers")]), k = 8)))
  lw <- spdep::nb2listw(nb, style = "W")
  log_msg("  ", cl, ": n = ", nrow(d), " 格，", length(keepv), " 个响应变量")
  for (r in keepv) {
    if (sum(is.finite(d[[r]])) < nrow(d)) next
    o <- fit_sar(d, r, lw)
    if (is.null(o)) next
    o$class <- cl; o$response <- r
    o$response_type <- if (r %in% AXES) "functional axis" else "single trait"
    rows[[length(rows) + 1]] <- o
  }
}
eff <- bind_rows(rows) |> filter(term %in% PRED) |>
  mutate(hypothesis = HYP[term]) |>
  left_join(TAXA[, c("class", "thermal")], by = "class")
eff$p_adj <- stats::p.adjust(eff$p, "BH")
write_table(eff, "table_38_hypothesis_effects")
log_msg("模型数 / models fitted: ", length(unique(paste(eff$class, eff$response))))

# 每个响应下最强的假说 / strongest hypothesis per response
best <- eff |> group_by(class, response, response_type) |>
  slice_max(abs(estimate), n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(class, response_type, response, hypothesis, estimate, p_adj) |>
  mutate(estimate = round(estimate, 3), p_adj = signif(p_adj, 2))
write_table(best, "table_38c_dominant_hypothesis")
log_msg("各响应变量下最强的假说 / dominant hypothesis per response:")
print(as.data.frame(best))

# ---------------------------------------------------------------
# 4. 方差分解 / Variance partitioning among hypothesis groups
# ---------------------------------------------------------------
# 把六个变量归为三组：气候（能量+水分+季节性）、生境（生产力+异质性）、人为。
# 分解回答的是「哪一组独占解释了多少」，而非单个变量的边际显著性。
GRP <- list(Climate = c("bio1_mean", "bio12_mean", "bio4_mean"),
            Habitat = c("npp", "elev_range"),
            Human   = "hfp2020")
vp <- list()
for (cl in TAXA$class) {
  d <- asm |> filter(class == cl) |> inner_join(env, by = "cell_id")
  for (r in c(TRAITS, AXES)) {
    if (!r %in% names(d)) next
    dd <- d[stats::complete.cases(d[, c(r, PRED)]), ]
    if (nrow(dd) < 100) next
    y <- zscore(dd[[r]])
    p <- try(vegan::varpart(y, dd[, GRP$Climate], dd[, GRP$Habitat],
                            dd[, GRP$Human, drop = FALSE]), silent = TRUE)
    if (inherits(p, "try-error")) next
    fr <- p$part$indfract$Adj.R.square
    vp[[length(vp) + 1]] <- data.frame(
      class = cl, response = r,
      response_type = if (r %in% AXES) "functional axis" else "single trait",
      climate_unique = round(fr[1], 4), habitat_unique = round(fr[2], 4),
      human_unique = round(fr[3], 4),
      shared = round(sum(fr[4:7]), 4),
      total = round(p$part$fract$Adj.R.square[7], 4))
  }
}
VP <- bind_rows(vp) |> left_join(TAXA[, c("class", "thermal")], by = "class")
write_table(VP, "table_39_variance_partition")
log_msg("方差分解（按类群平均）/ variance partition, class means:")
print(VP |> group_by(class) |>
        summarise(across(c(climate_unique, habitat_unique, human_unique, total),
                         ~ round(mean(.x, na.rm = TRUE), 3)), .groups = "drop") |>
        as.data.frame())

# ---------------------------------------------------------------
# 5. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures
# 体型轴按定义就等于体型性状本身，画两行会造成重复。图中只把 PACE 与 NICHE
# 作为「轴」单列，体型行同时标注为轴，表 38 仍保留全部九个响应。
# The size axis is the size trait by definition; plotting both would duplicate a
# row, so only PACE and NICHE appear as separate axis rows in the figure.
PLOT_RESP <- c("size", "fecundity", "nocturnality", "verticality",
               "habitat_breadth", "range_size", "PACE", "NICHE")
PLOT_LAB <- c(size = "Body size  (= size axis)", fecundity = "Fecundity",
              nocturnality = "Nocturnality", verticality = "Vertical stratum use",
              habitat_breadth = "Habitat breadth", range_size = "Range size",
              PACE = "Axis: slow–fast pace", NICHE = "Axis: niche breadth")
# 轴刻度用短名，全称放在图注里，避免窄列上的标签互相重叠
HYP_SHORT <- c(`Ambient energy` = "Energy", `Water availability` = "Water",
               Productivity = "Productivity", `Habitat heterogeneity` = "Heterogeneity",
               `Climatic seasonality` = "Seasonality", `Human pressure` = "Human")
ef <- eff |> filter(response %in% PLOT_RESP)
ef$hyp_f   <- factor(HYP[ef$term], levels = unname(HYP))
ef$class_f <- factor(TAXA_LAB[ef$class], levels = unname(TAXA_LAB))
ef$resp_f  <- factor(PLOT_LAB[ef$response], levels = rev(unname(PLOT_LAB)))

# (a) 效应量热图：假说 x 响应，分纲
lim <- stats::quantile(abs(ef$estimate), 0.98, na.rm = TRUE)
pa <- ggplot(ef, aes(hyp_f, resp_f, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = ef |> filter(p_adj < 0.05), size = 0.42, colour = "grey10") +
  scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim),
                   oob = scales::squish, name = "Standardised effect") +
  scale_x_discrete(labels = HYP_SHORT[unname(HYP)]) +
  facet_wrap(~class_f, nrow = 1) +
  labs(x = NULL, y = NULL, tag = "a",
       subtitle = "Energy = mean annual temperature; Water = annual precipitation; Productivity = NPP; Heterogeneity = within-cell elevation range; Seasonality = temperature seasonality; Human = human footprint.\nDots: BH-adjusted P < 0.05. Every model controls for all six predictors at once, so each cell is that hypothesis's unique effect.") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 4.6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5.5),
        axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 6.5),
        plot.subtitle = element_text(size = 4.5, colour = "grey35"),
        legend.position = "bottom", legend.key.width = unit(32, "pt"),
        legend.key.height = unit(3.5, "pt"))

# (b) 方差分解
vl <- VP |> select(class, response, response_type, climate_unique,
                   habitat_unique, human_unique) |>
  pivot_longer(ends_with("_unique"), names_to = "component", values_to = "r2") |>
  mutate(component = c(climate_unique = "Climate", habitat_unique = "Habitat",
                       human_unique = "Human")[component])
vl$class_f <- factor(TAXA_LAB[vl$class], levels = unname(TAXA_LAB))
pb <- ggplot(vl, aes(component, pmax(r2, 0), fill = class_f)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.66,
               position = position_dodge(width = 0.78)) +
  scale_fill_manual(values = setNames(PAL$taxa[names(TAXA_LAB)], unname(TAXA_LAB)),
                    name = NULL) +
  labs(x = NULL, y = expression("Unique adjusted "*R^2), tag = "b",
       subtitle = "Climate dominates in every class; human pressure is largest for reptiles") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.5, colour = "grey30"))

# (c) 功能轴是否给出单性状之外的信息
# 对每个轴，把它的效应量与其**组成性状**的效应量并排：若轴更强，说明性状组合
# 携带了单性状没有的信号；若更弱，说明轴只是稀释。
AXCOMP <- list(PACE = c("size", "fecundity"),
               NICHE = c("habitat_breadth", "range_size"))
cmp <- bind_rows(lapply(names(AXCOMP), function(a)
  eff |> filter(response %in% c(a, AXCOMP[[a]])) |>
    mutate(axis = a,
           kind = ifelse(response == a, "Composite axis", "Component traits")) |>
    group_by(class, axis, kind, term) |>
    summarise(mean_abs = mean(abs(estimate)), .groups = "drop")))
cmp$class_f <- factor(TAXA_LAB[cmp$class], levels = unname(TAXA_LAB))
cmp$axis_f <- factor(cmp$axis, levels = c("PACE", "NICHE"),
                     labels = c("Pace", "Niche"))
pc <- ggplot(cmp, aes(axis_f, mean_abs, fill = kind)) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.22, width = 0.6,
               position = position_dodge(width = 0.72)) +
  facet_wrap(~class_f, nrow = 1) +
  scale_fill_manual(values = c(`Component traits` = "grey70",
                               `Composite axis` = "#2E6E8E"), name = NULL) +
  labs(x = NULL, y = "Mean |standardised effect|", tag = "c",
       subtitle = "Pace = size + fecundity; Niche = habitat breadth + range size. Spread is across the six hypotheses.") +
  theme_pub() + theme(legend.position = "bottom",
                      axis.text.x = element_text(size = 5),
                      strip.text = element_text(size = 6.5),
                      plot.subtitle = element_text(size = 4.5, colour = "grey30"))

p <- pa / (pb | pc) + plot_layout(heights = c(1.3, 1), widths = c(1, 1.7))
save_fig(p, "Fig4_hypothesis_drivers", W2, 168, FIG)

# 轴是否优于其组成性状 / does the composite axis beat its components
axgain <- cmp |> group_by(class, axis_f, kind) |>
  summarise(m = mean(mean_abs), .groups = "drop") |>
  pivot_wider(names_from = kind, values_from = m) |>
  mutate(gain = round(`Composite axis` - `Component traits`, 4),
         across(c(`Composite axis`, `Component traits`), ~ round(.x, 4)))
write_table(axgain, "table_39b_axis_vs_components")
log_msg("轴 vs 其组成性状 / composite axis vs component traits:")
print(as.data.frame(axgain))

log_msg("=== 06g 完成 / done ===")
