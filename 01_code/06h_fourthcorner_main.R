# ============================================================
# 科学问题 / Scientific question:
#   在统计学上正确的检验下（同时置换样方与物种），哪些「假说变量 x 性状/
#   功能轴」的关联真正成立？朴素 CWM 回归会夸大到什么程度？
#   Under a statistically valid test that permutes both sites and species,
#   which hypothesis-variable by trait/axis links actually hold?
#
# 分析目标 / Objective:
#   把第四角检验从附录提到正文：用与 06g 完全相同的六个具名环境变量
#   （对应六个假说）与九个响应（六个单性状 + 三个先验功能轴）构成
#   6 x 9 的检验矩阵，逐纲给出 max test 的效应与显著性。
#   这样 06g 的 SAR 系数（效应量）与本脚本的 max test（显著性）
#   互为印证：SAR 说方向和强度，max test 说这条关联是否经得起置换。
#
# 为什么必须用 max test / Why the max test is required:
#   CWM ~ 环境 的普通回归中，所有样方共用同一个性状向量，残差因此不独立，
#   第一类错误率可膨胀到 0.5 以上（Peres-Neto et al. 2017; ter Braak 2018）。
#   modeltype = 6 取「置换样方 (model 2)」与「置换物种 (model 4)」两个
#   p 值的最大值，是目前公认正确的检验。
#
# 输入数据 / Input:
#   03_data_derived/comm_<class>_50km.rds, functional_axes_species.rds, env_50km.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_40_fourthcorner_hypothesis.csv
#   04_results/tables/table_40b_typeI_inflation.csv
#   05_figures/Fig5_fourth_corner.pdf/.png
#
# 关键假设 / Key assumptions:
#   - 样方做空间分层抽样至 1200 格：完整矩阵的双向置换在单机不可行，
#     分层抽样保持地理覆盖，且是置换检验的标准做法。
#   - 功能轴作为「性状」进入 Q 矩阵是合法的：轴是性状的线性组合，
#     第四角对 Q 的列不要求彼此独立。
#
# 主要包 / Main packages: ade4, dplyr
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({ library(ade4); library(dplyr); library(tidyr) })

log_msg("=== 06h 正文第四角检验 / Main-text fourth-corner ===")

NREP    <- 999L
N_CELLS <- 1200L

# 与 06g 完全一致的六个假说变量 / the same six hypothesis variables as 06g
HYP <- c(bio1_mean   = "Ambient energy",
         bio12_mean  = "Water availability",
         npp         = "Productivity",
         elev_range  = "Habitat heterogeneity",
         bio4_mean   = "Climatic seasonality",
         hfp2020     = "Human pressure")
PRED <- names(HYP)
TRAITS <- c("size", "fecundity", "nocturnality", "verticality",
            "habitat_breadth", "range_size")
AXES   <- c("SIZE", "PACE", "NICHE")
RESP_LAB <- c(size = "Body size", fecundity = "Fecundity",
              nocturnality = "Nocturnality", verticality = "Vertical stratum use",
              habitat_breadth = "Habitat breadth", range_size = "Range size",
              SIZE = "Axis: body size", PACE = "Axis: slow–fast pace",
              NICHE = "Axis: niche breadth")

FA <- readRDS(file.path(PATH$derived, "functional_axes_species.rds"))

run_fc <- function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  comm <- drop_inset(readRDS(f), "50km")
  env  <- drop_inset(readRDS(file.path(PATH$derived, "env_50km.rds")), "50km")
  tr   <- FA |> filter(class == cl) |> distinct(species, .keep_all = TRUE)

  sp <- intersect(unique(comm$species), tr$species)
  if (length(sp) < 30) return(NULL)
  tr <- tr[match(sp, tr$species), ]
  comm <- comm[comm$species %in% sp, ]

  cells <- sort(unique(comm$cell_id))
  L <- matrix(0L, length(cells), length(sp), dimnames = list(cells, sp))
  L[cbind(match(comm$cell_id, cells), match(comm$species, sp))] <- 1L

  # 空间分层抽样，保持地理覆盖 / spatially stratified subsample
  if (nrow(L) > N_CELLS) {
    xy <- env[match(rownames(L), env$cell_id), c("lon", "lat")]
    bin <- paste(cut(xy$lon, 12, labels = FALSE), cut(xy$lat, 12, labels = FALSE))
    set.seed(SEED)
    idx <- unlist(lapply(split(seq_len(nrow(L)), bin), function(ii)
      if (length(ii) <= ceiling(N_CELLS / length(unique(bin)))) ii
      else sample(ii, ceiling(N_CELLS / length(unique(bin))))))
    L <- L[sort(idx), , drop = FALSE]
  }
  L <- L[rowSums(L) >= 5, , drop = FALSE]
  L <- L[, colSums(L) >= 3, drop = FALSE]
  if (nrow(L) < 100 || ncol(L) < 20) return(NULL)

  R <- env[match(rownames(L), env$cell_id), PRED, drop = FALSE]
  R <- as.data.frame(lapply(R, function(z) {
    z[!is.finite(z)] <- stats::median(z, na.rm = TRUE); as.numeric(zscore(z)) }))
  rownames(R) <- rownames(L)

  Q <- tr[match(colnames(L), tr$species), c(TRAITS, AXES), drop = FALSE]
  Q <- as.data.frame(lapply(Q, function(z) {
    z <- as.numeric(z); z[!is.finite(z)] <- stats::median(z, na.rm = TRUE); z }))
  rownames(Q) <- colnames(L)
  Q <- Q[, vapply(Q, function(z) stats::sd(z) > 0, logical(1)), drop = FALSE]
  if (ncol(Q) < 4) return(NULL)

  log_msg("  ", cl, ": L = ", nrow(L), " x ", ncol(L), ", R = ", ncol(R),
          " env, Q = ", ncol(Q), " responses")

  set.seed(SEED)
  fc <- try(ade4::fourthcorner(tabR = R, tabL = as.data.frame(L), tabQ = Q,
                               modeltype = 6, p.adjust.method.G = "BH",
                               p.adjust.method.D = "BH", nrepet = NREP),
            silent = TRUE)
  if (inherits(fc, "try-error")) { log_msg("   [fail] ", cl); return(NULL) }

  out <- data.frame(class = cl,
                    env   = rep(colnames(R), times = ncol(Q)),
                    response = rep(colnames(Q), each = ncol(R)),
                    stat  = as.numeric(fc$tabD2$obs),
                    p_maxtest = as.numeric(fc$tabD2$adj.pvalue),
                    row.names = NULL)

  # 朴素 CWM 回归，量化第一类错误膨胀 / naive CWM regression for contrast
  naive <- list()
  for (tq in colnames(Q)) {
    cwm <- as.numeric(L %*% Q[[tq]] / rowSums(L))
    for (e1 in colnames(R)) {
      m <- stats::lm(cwm ~ R[[e1]])
      naive[[length(naive) + 1]] <- data.frame(
        env = e1, response = tq,
        r_naive = sign(stats::coef(m)[2]) * sqrt(summary(m)$r.squared),
        p_naive = summary(m)$coefficients[2, 4])
    }
  }
  naive <- bind_rows(naive)
  naive$p_naive_adj <- stats::p.adjust(naive$p_naive, "BH")
  # 第四角统计量本身无符号，用朴素相关的符号标注方向
  # The fourth-corner statistic is unsigned; direction is taken from the
  # naive correlation, which is unbiased for the point estimate.
  left_join(out, naive, by = c("env", "response"))
}

# 置换检验约需 9 分钟，缓存结果使制图可反复调整而不必重算
# The permutation test takes ~9 min; caching lets the figure be refined freely.
res <- cache_rds("fourthcorner_hypothesis",
  bind_rows(lapply(TAXA$class, function(cl) {
    r <- try(run_fc(cl), silent = TRUE)
    if (inherits(r, "try-error")) NULL else r
  })))
stopifnot(nrow(res) > 0)
res <- res |> left_join(TAXA[, c("class", "thermal")], by = "class") |>
  mutate(hypothesis = HYP[env],
         response_type = ifelse(response %in% AXES, "functional axis", "single trait"),
         signed_stat = sign(r_naive) * abs(stat))
write_table(res, "table_40_fourthcorner_hypothesis")

# 第一类错误膨胀 / Type I inflation
infl <- res |> group_by(class, thermal) |>
  summarise(n_tests = n(),
            sig_naive = sum(p_naive_adj < 0.05, na.rm = TRUE),
            sig_maxtest = sum(p_maxtest < 0.05, na.rm = TRUE),
            pct_naive = round(100 * mean(p_naive_adj < 0.05, na.rm = TRUE), 1),
            pct_maxtest = round(100 * mean(p_maxtest < 0.05, na.rm = TRUE), 1),
            inflation = round(sig_naive / pmax(sig_maxtest, 1), 2),
            .groups = "drop")
write_table(infl, "table_40b_typeI_inflation")
log_msg("第一类错误膨胀 / Type I inflation:"); print(as.data.frame(infl))

# 通过 max test 的关联，按假说汇总 / links surviving the max test, by hypothesis
surv <- res |> filter(p_maxtest < 0.05) |>
  count(hypothesis, response_type, name = "n_significant") |>
  arrange(desc(n_significant))
write_table(surv, "table_40c_surviving_links")
log_msg("通过 max test 的关联数 / links surviving, by hypothesis:")
print(as.data.frame(surv))

# ---------------------------------------------------------------
# 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures
res$hyp_f  <- factor(HYP[res$env], levels = unname(HYP))
res$resp_f <- factor(RESP_LAB[res$response], levels = rev(unname(RESP_LAB)))
res$class_f <- factor(TAXA_LAB[res$class], levels = unname(TAXA_LAB))
res$sig <- res$p_maxtest < 0.05

lim <- stats::quantile(abs(res$signed_stat), 0.98, na.rm = TRUE)
pa <- ggplot(res, aes(hyp_f, resp_f)) +
  geom_tile(aes(fill = signed_stat), colour = "white", linewidth = 0.35) +
  geom_point(data = res |> filter(sig), size = 0.5, colour = "grey5") +
  scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim),
                   oob = scales::squish,
                   name = "Fourth-corner statistic (signed)") +
  facet_wrap(~class_f, nrow = 1) +
  labs(x = NULL, y = NULL, tag = "a",
       subtitle = "Dots: significant after the max test (BH-adjusted P < 0.05)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5.5),
        axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 6.5),
        plot.subtitle = element_text(size = 4.8, colour = "grey30"),
        legend.position = "bottom", legend.key.width = unit(34, "pt"),
        legend.key.height = unit(3.5, "pt"))

# (b) 朴素检验 vs max test：第一类错误膨胀
il <- infl |> select(class, thermal, pct_naive, pct_maxtest) |>
  pivot_longer(starts_with("pct_"), names_to = "test", values_to = "pct") |>
  mutate(test = c(pct_naive = "Naive CWM regression",
                  pct_maxtest = "Max test (valid)")[test])
il$class_f <- factor(TAXA_LAB[il$class], levels = unname(TAXA_LAB))
pb <- ggplot(il, aes(class_f, pct, fill = test)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.66) +
  scale_fill_manual(values = c(`Naive CWM regression` = "#C2703D",
                               `Max test (valid)` = "#2E6E8E"), name = NULL) +
  labs(x = NULL, y = "% of tests declared significant", tag = "b",
       subtitle = "The naive test inflates significance because all sites share one trait vector") +
  theme_pub() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1),
        plot.subtitle = element_text(size = 4.6, colour = "grey30"))

# (c) 各假说下通过检验的**比例**，分响应类型。
# 用比例而非原始计数：单性状有 6 个、功能轴只有 3 个，两类的检验数不同，
# 直接比计数会把「性状数量更多」误读成「性状证据更强」。
# Proportions, not raw counts: there are six traits but only three axes, so
# counts would confound the number of tests with the strength of the evidence.
sc <- res |> group_by(hyp_f, response_type) |>
  summarise(pct_sig = 100 * mean(sig, na.rm = TRUE), n_tot = n(), .groups = "drop")
pc <- ggplot(sc, aes(hyp_f, pct_sig, fill = response_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_manual(values = c(`single trait` = "grey65",
                               `functional axis` = "#2E6E8E"), name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "% of tests significant (max test)", tag = "c",
       subtitle = "Pooled across the four classes: 24 tests per trait bar, 12 per axis bar") +
  theme_pub() +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5.5),
        plot.subtitle = element_text(size = 4.6, colour = "grey30"))

p <- pa / (pb | pc) + plot_layout(heights = c(1.4, 1))
save_fig(p, "Fig5_fourth_corner", W2, 158, FIG)

log_msg("=== 06h 完成 / done ===")
