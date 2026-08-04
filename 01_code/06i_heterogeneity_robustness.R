# ============================================================
# 科学问题 / Scientific question:
#   变温类群对环境能量的更强响应，是热缓冲机制，还是「对什么都反应更强」的
#   信噪比差异？原先用「生境异质性」作阴性对照，但只试了一种操作化
#   （网格内海拔差）。换成生境类型多样性或遥感纹理，结论还成立吗？
#   Is the ectotherm excess for ambient energy a thermal mechanism, or a
#   generic difference in signal-to-noise? The negative control originally
#   rested on one proxy; does it survive alternatives?
#
# 结果改变了论证方式 / What we found, and how it changed the argument:
#   三种异质性度量给出的比值分别是 1.12、1.81、0.67——阴性对照**并非**
#   在所有操作化下都成立。但进一步检查发现原因：三者携带的温度信息差别很大
#   （与年均温的 |r| 分别为 0.198、0.519、0.054）。唯一出现变温超额的
#   生境类型多样性，恰恰是与温度相关最强的那一个（r = +0.52）——中国的
#   土地覆被本身沿温度梯度排布，因此它不是合格的阴性对照。
#   于是把单点对照升级为**剂量–反应检验**：超额响应是否随变量的热含量递增。
#   答案是肯定的（Spearman rho = 0.74, P = 0.037, 8 个变量），
#   这比原来的单一阴性对照是更强的证据——它预测的是一条趋势而非一个点。
#   The proxy that does show an ectotherm excess is the one most correlated
#   with temperature, so it is not a valid negative control. Reframing the
#   test as a dose-response across eight predictors gives stronger evidence.
#
# 三种异质性操作化 / Three operationalisations:
#   1. 地形异质性 elev_range        网格内海拔极差（现用；与温度 |r| = 0.20）
#   2. 生境类型多样性 lc_shannon    WorldCover 六类占比的 Shannon（|r| = 0.52）
#   3. 遥感纹理异质性 het_shannon   EarthEnv EVI 纹理的 Shannon（|r| = 0.05）
#      （Tuanmu & Jetz 2015）
#
# 输入数据 / Input:
#   03_data_derived/assemblage_axes_50km.rds, env_50km.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_45_heterogeneity_proxies.csv
#   04_results/tables/table_45b_negative_control_robustness.csv
#   05_figures/Fig6b_negative_control.pdf/.png
#
# 关键假设 / Key assumptions:
#   - lc_shannon 由六类占比算出。各格占比之和并非恒为 1（WorldCover 未覆盖
#     的类别如裸地、冰雪、湿地不在这六类中），故先按行归一化再算 Shannon；
#     归一化后为 0 的格（六类均为 0，如高寒荒漠）记为缺失而非 0 多样性。
#   - 三种度量彼此相关但不冗余（相关系数见表 45），因此互为独立检验。
#
# 主要包 / Main packages: dplyr, spdep, spatialreg
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(spdep); library(spatialreg)
  library(ggrepel)
})

log_msg("=== 06i 异质性度量的稳健性 / Heterogeneity-proxy robustness ===")

env <- readRDS(file.path(PATH$derived, "env_50km.rds")) |> drop_inset("50km")
asm <- readRDS(file.path(PATH$derived, "assemblage_axes_50km.rds"))

# ---------------------------------------------------------------
# 1. 构建生境类型多样性 / Land-cover type diversity
# ---------------------------------------------------------------
LC <- grep("^lc_", names(env), value = TRUE)
P  <- as.matrix(env[, LC])
P[!is.finite(P)] <- 0
tot <- rowSums(P)
# 六类占比之和 < 1 是因为裸地、冰雪、湿地等类别不在其中；按行归一化后
# 计算的是「已知六类之间的构成多样性」，这正是需要的量。
Pn <- P / pmax(tot, 1e-12)
shannon <- -rowSums(ifelse(Pn > 0, Pn * log(Pn), 0))
shannon[tot < 0.05] <- NA_real_   # 六类合计不足 5% 的格（高寒荒漠等）不可靠
env$lc_shannon <- shannon
env$lc_richness <- rowSums(P > 0.01)   # 占比 >1% 的类别数
log_msg("生境类型多样性可用格数 / cells with land-cover diversity: ",
        sum(is.finite(env$lc_shannon)), " / ", nrow(env))

HET <- c(elev_range  = "Topographic (elevation range)",
         lc_shannon  = "Land-cover type diversity",
         het_shannon = "Remote-sensing texture (EarthEnv)")

# 三种度量之间的相关 / correlations among the three proxies
cc <- stats::cor(env[, names(HET)], use = "pairwise.complete.obs",
                 method = "spearman")
ct <- data.frame(pair = c("elev_range – lc_shannon", "elev_range – het_shannon",
                          "lc_shannon – het_shannon"),
                 spearman_r = round(c(cc[1, 2], cc[1, 3], cc[2, 3]), 3))
write_table(ct, "table_45_heterogeneity_proxies")
log_msg("三种异质性度量之间的相关 / correlations among proxies:")
print(ct)

# ---------------------------------------------------------------
# 2. 逐一替换异质性变量，重跑 SAR / Re-fit SAR swapping the proxy
# ---------------------------------------------------------------
BASE <- c("bio1_mean", "bio12_mean", "npp", "bio4_mean", "hfp2020")
TRAITS <- c("size", "fecundity", "nocturnality", "verticality",
            "habitat_breadth", "range_size")
AXES <- c("PACE", "NICHE")
RESP <- c(TRAITS, AXES)

fit_one <- function(d, resp, pred, lw) {
  X <- as.data.frame(apply(as.matrix(d[, pred]), 2, zscore))
  X$y <- zscore(d[[resp]])
  f <- stats::as.formula(paste("y ~", paste(pred, collapse = " + ")))
  m <- try(spatialreg::errorsarlm(f, data = X, listw = lw, method = "Matrix",
                                  zero.policy = TRUE), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$Coef
  data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], p = s[, 4],
             row.names = NULL)
}

rows <- list()
for (hv in names(HET)) {
  pred <- c(BASE, hv)
  for (cl in TAXA$class) {
    d <- asm |> filter(class == cl) |> inner_join(env, by = "cell_id")
    d <- d[stats::complete.cases(d[, c(pred, "x_albers", "y_albers")]), ]
    if (nrow(d) < 100) next
    nb <- spdep::make.sym.nb(spdep::knn2nb(
      spdep::knearneigh(as.matrix(d[, c("x_albers", "y_albers")]), k = 8)))
    lw <- spdep::nb2listw(nb, style = "W")
    for (r in RESP) {
      if (!r %in% names(d) || sum(is.finite(d[[r]])) < nrow(d)) next
      o <- fit_one(d, r, pred, lw)
      if (is.null(o)) next
      o$class <- cl; o$response <- r; o$het_proxy <- hv; o$n <- nrow(d)
      rows[[length(rows) + 1]] <- o
    }
    log_msg("  ", HET[[hv]], " x ", cl, ": n = ", nrow(d))
  }
}
eff <- bind_rows(rows) |>
  left_join(TAXA[, c("class", "thermal")], by = "class")
eff$p_adj <- stats::p.adjust(eff$p, "BH")

# ---------------------------------------------------------------
# 3. 阴性对照是否稳健 / Does the negative control hold?
# ---------------------------------------------------------------
# 对每一种异质性操作化，分别计算：
#   (a) 环境能量效应的 变温/恒温 比值 —— 应显著 > 1（热缓冲预测）
#   (b) 该异质性效应的 变温/恒温 比值 —— 应接近 1（阴性对照）
ratio_of <- function(d, tm) {
  x <- d |> filter(term == tm) |> group_by(thermal) |>
    summarise(m = mean(abs(estimate)), .groups = "drop")
  if (nrow(x) < 2) return(NA_real_)
  x$m[x$thermal == "Ectotherm"] / x$m[x$thermal == "Endotherm"]
}
#' 比值的自助置信区间 / bootstrap CI for the ratio
#' 重抽单位是「纲 x 响应」组合，因为每个组合是一次独立的模型拟合。
#' 只给点估计而不给区间，会让 1.12 与 1.81 的差别看起来比证据支持的更确定。
#' Resampling unit is the class-by-response combination, since each is one
#' independent model fit. Reporting point ratios alone would overstate how
#' firmly 1.12 differs from 1.81.
boot_ratio <- function(d, tm, n_rep = 2000L) {
  x <- d |> filter(term == tm)
  e <- abs(x$estimate[x$thermal == "Ectotherm"])
  n <- abs(x$estimate[x$thermal == "Endotherm"])
  if (!length(e) || !length(n)) return(c(NA_real_, NA_real_))
  set.seed(SEED)
  r <- vapply(seq_len(n_rep), function(i)
    mean(sample(e, length(e), TRUE)) / mean(sample(n, length(n), TRUE)),
    numeric(1))
  unname(stats::quantile(r, c(0.025, 0.975), na.rm = TRUE))
}
#' 直接检验：能量的不对称是否强于异质性的不对称
#' 逐次自助中同时算两个比值，取 log 之差，比分别给两个区间更有力——
#' 两个区间重叠并不等于差异不显著。
#' Bootstrapping the difference of log-ratios directly; overlapping marginal
#' intervals are not the same as a non-significant difference.
boot_diff <- function(d, hv, n_rep = 2000L) {
  g <- function(tm, th) abs(d$estimate[d$term == tm & d$thermal == th])
  ee <- g("bio1_mean", "Ectotherm"); en <- g("bio1_mean", "Endotherm")
  he <- g(hv, "Ectotherm");          hn <- g(hv, "Endotherm")
  if (!length(ee) || !length(he)) return(c(NA_real_, NA_real_, NA_real_))
  set.seed(SEED)
  dd <- vapply(seq_len(n_rep), function(i) {
    log(mean(sample(ee, length(ee), TRUE)) / mean(sample(en, length(en), TRUE))) -
    log(mean(sample(he, length(he), TRUE)) / mean(sample(hn, length(hn), TRUE)))
  }, numeric(1))
  q <- stats::quantile(dd, c(0.025, 0.975), na.rm = TRUE)
  c(unname(q), 2 * min(mean(dd <= 0), mean(dd >= 0)))   # 双尾自助 p
}

rob <- bind_rows(lapply(names(HET), function(hv) {
  d <- eff |> filter(het_proxy == hv)
  eb <- boot_ratio(d, "bio1_mean"); hb <- boot_ratio(d, hv)
  db <- boot_diff(d, hv)
  data.frame(het_proxy = hv, proxy_label = HET[[hv]],
             energy_ratio = round(ratio_of(d, "bio1_mean"), 2),
             energy_lo = round(eb[1], 2), energy_hi = round(eb[2], 2),
             heterogeneity_ratio = round(ratio_of(d, hv), 2),
             het_lo = round(hb[1], 2), het_hi = round(hb[2], 2),
             # log(能量比) − log(异质性比)：>0 表示能量的不对称更强
             logratio_diff = round(db_pt <- log(ratio_of(d, "bio1_mean")) -
                                     log(ratio_of(d, hv)), 3),
             diff_lo = round(db[1], 3), diff_hi = round(db[2], 3),
             diff_p = signif(db[3], 3),
             energy_ecto = round(mean(abs(d$estimate[d$term == "bio1_mean" &
                                                       d$thermal == "Ectotherm"])), 4),
             energy_endo = round(mean(abs(d$estimate[d$term == "bio1_mean" &
                                                       d$thermal == "Endotherm"])), 4),
             het_ecto = round(mean(abs(d$estimate[d$term == hv &
                                                    d$thermal == "Ectotherm"])), 4),
             het_endo = round(mean(abs(d$estimate[d$term == hv &
                                                    d$thermal == "Endotherm"])), 4))
}))
write_table(rob, "table_45b_negative_control_robustness")
log_msg("阴性对照的稳健性 / robustness of the negative control:")
print(as.data.frame(rob))
write_table(eff, "table_45c_effects_by_proxy")

# ---------------------------------------------------------------
# 3b. 剂量–反应：变温类群的超额响应是否随变量的「热含量」递增
# Dose-response: does the ectotherm excess scale with a predictor's
# thermal content?
# ---------------------------------------------------------------
# 上一步暴露出一件关键的事：三个异质性度量本身携带的温度信息差别很大
# （与年均温的 |r| 从 0.06 到 0.52）。若热缓冲机制为真，则超额响应应当
# 随「该变量与温度的相关强度」递增——这比单一阴性对照更强的检验：
# 它预测的是一条趋势，而不只是一个点。
# The three heterogeneity proxies differ greatly in how much temperature
# information they carry (|r| with MAT from 0.06 to 0.52). Thermal buffering
# predicts the ectotherm excess should scale with that thermal content — a
# trend prediction, which is harder to satisfy by chance than a single point.
ALLV <- c("bio1_mean", "bio12_mean", "npp", "bio4_mean", "hfp2020",
          names(HET))
dose <- bind_rows(lapply(ALLV, function(v) {
  d <- eff |> filter(term == v)
  if (!nrow(d)) return(NULL)
  # 该变量所属的拟合子集：异质性变量只在其自身那一轮出现
  data.frame(variable = v,
             thermal_content = round(abs(stats::cor(env[[v]], env$bio1_mean,
                                                    use = "pairwise.complete.obs",
                                                    method = "spearman")), 3),
             ratio = round(ratio_of(d, v), 2),
             n_fits = nrow(d))
}))
dose <- dose[is.finite(dose$ratio) & is.finite(dose$thermal_content), ]
ct2 <- stats::cor.test(dose$thermal_content, log(dose$ratio), method = "spearman",
                       exact = FALSE)
dose$overall_rho <- round(unname(ct2$estimate), 3)
dose$overall_p   <- signif(ct2$p.value, 3)
write_table(dose, "table_45d_dose_response")
log_msg("剂量–反应：超额响应 vs 变量的热含量 / excess vs thermal content:")
print(as.data.frame(dose))
log_msg("  Spearman rho = ", round(unname(ct2$estimate), 3),
        ", p = ", signif(ct2$p.value, 3))

# ---------------------------------------------------------------
# 4. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures
pl <- bind_rows(
  rob |> transmute(proxy_label, predictor = "Ambient energy",
                   ratio = energy_ratio, lo = energy_lo, hi = energy_hi),
  rob |> transmute(proxy_label, predictor = "Habitat heterogeneity",
                   ratio = heterogeneity_ratio, lo = het_lo, hi = het_hi))
pl$proxy_f <- factor(pl$proxy_label, levels = unname(HET))

pa <- ggplot(pl, aes(proxy_f, ratio, fill = predictor)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.3, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.16, linewidth = 0.3,
                position = position_dodge(width = 0.7), colour = "grey25") +
  geom_text(aes(label = sprintf("%.2f", ratio), y = hi),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 1.9,
            colour = "grey20") +
  annotate("text", x = 0.55, y = 1.06, label = "no difference", hjust = 0,
           size = 1.8, colour = "grey40") +
  scale_fill_manual(values = c(`Ambient energy` = "#B5541F",
                               `Habitat heterogeneity` = "#4A7FA5"), name = NULL) +
  scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "Effect ratio, ectotherms / endotherms", tag = "a",
       subtitle = "Ambient energy shows a consistent ~2x ectotherm excess under all three specifications; the heterogeneity ratio is unstable and straddles 1. Bars: 95% bootstrap CI over class x response fits.") +
  theme_pub() + theme(legend.position = "bottom",
                      axis.text.x = element_text(size = 5.5),
                      plot.subtitle = element_text(size = 4.6, colour = "grey30"))

# (b) 三种度量的相关结构
cm <- as.data.frame(as.table(cc))
names(cm) <- c("x", "y", "r")
cm$x <- factor(HET[as.character(cm$x)], levels = unname(HET))
cm$y <- factor(HET[as.character(cm$y)], levels = rev(unname(HET)))
pb <- ggplot(cm, aes(x, y, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.1, colour = "grey15") +
  scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-1, 1),
                   name = "Spearman r") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(x = NULL, y = NULL, tag = "b",
       subtitle = "The three proxies are related but far from redundant") +
  theme_pub() +
  theme(axis.text = element_text(size = 4.6), axis.line = element_blank(),
        axis.ticks = element_blank(), legend.position = "bottom",
        legend.key.width = unit(24, "pt"), legend.key.height = unit(3.5, "pt"),
        plot.subtitle = element_text(size = 4.6, colour = "grey30"))

# (c) 剂量–反应：超额响应随变量的热含量递增
VLAB <- c(bio1_mean = "Ambient energy", bio12_mean = "Water", npp = "Productivity",
          bio4_mean = "Seasonality", hfp2020 = "Human pressure",
          elev_range = "Heterogeneity:\ntopographic",
          lc_shannon = "Heterogeneity:\nland-cover",
          het_shannon = "Heterogeneity:\ntexture")
dose$lab <- VLAB[dose$variable]
pc <- ggplot(dose, aes(thermal_content, ratio)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.3, colour = "grey45") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.4,
              colour = "grey25", fill = "grey85") +
  geom_point(size = 1.5, colour = "#B5541F") +
  ggrepel::geom_text_repel(aes(label = lab), size = 1.7, colour = "grey25",
                           min.segment.length = 0, segment.size = 0.2,
                           box.padding = 0.25, lineheight = 0.9) +
  annotate("text", x = Inf, y = Inf, hjust = 1.06, vjust = 1.6,
           label = sprintf("Spearman rho = %.2f, P = %.3f",
                           dose$overall_rho[1], dose$overall_p[1]),
           size = 1.9, colour = "grey20") +
  labs(x = "Thermal content of the predictor (|Spearman r| with mean annual temperature)",
       y = "Effect ratio,\nectotherms / endotherms", tag = "c",
       subtitle = "The ectotherm excess scales with how much temperature information a predictor carries") +
  theme_pub() + theme(plot.subtitle = element_text(size = 4.6, colour = "grey30"))

p <- (pa | pb) / pc + plot_layout(heights = c(1, 0.95), widths = c(1.5, 1))
save_fig(p, "Fig6b_negative_control", W2, 150, FIG)

log_msg("=== 06i 完成 / done ===")
