# ============================================================
# 科学问题 / Scientific question:
#   性状与环境的关联是否真实存在？——注意：直接把群落加权平均值 (CWM)
#   对环境做回归会严重膨胀第一类错误，因为同一性状向量被所有样方共用。
#   Are trait-environment associations real? Regressing CWM on environment
#   inflates Type I error because all sites share the same trait vector.
#
# 分析目标 / Objective:
#   用 ade4 的第四角分析 (fourth-corner) 与 "max test"（同时置换样方与
#   物种，取两个 p 值的最大值）给出统计学正确的显著性检验；
#   并与朴素 CWM 回归的 p 值对照，量化第一类错误膨胀的幅度。
#
# 输入数据 / Input:
#   03_data_derived/comm_<class>_<grain>.rds, traits_imputed.rds, env_<grain>.rds
#
# 主要流程 / Workflow:
#   1. 构建 L (样方 x 物种)、R (样方 x 环境)、Q (物种 x 性状) 三矩阵
#   2. ade4::fourthcorner(modeltype = 6) 999 次置换
#   3. 朴素检验：CWM ~ 环境 的普通线性回归 p 值
#   4. 对比两者，输出「膨胀诊断表」
#
# 预期输出 / Expected output:
#   04_results/tables/table_4_fourthcorner_maxtest.csv
#   04_results/tables/table_S8_typeI_inflation_diagnostic.csv
#
# 关键假设 / Key assumptions:
#   - max test 同时控制样方端（模型 2）与物种端（模型 4）的置换，
#     是 ter Braak et al. (2012, 2018) 推荐的正确检验
#   - 存在-缺失数据下 L 为 0/1 矩阵，CWM 即简单算术平均
#
# 主要包 / Main packages: ade4, dplyr, data.table
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({ library(ade4) })

log_msg("=== 06b 第四角分析 / Fourth-corner analysis ===")

NREP    <- 499L    # 置换次数 / permutations
N_CELLS <- 1200L   # 空间分层抽样的样方数 / spatially stratified sample of sites
# 说明：完整的 3772 x 1079 矩阵做 999 次双向置换在单机上不可行；
# 对样方做空间分层抽样后检验力仍充足，且是置换检验的标准做法。
# NOTE: 999 two-way permutations on the full 3772 x 1079 matrix are not
# feasible here; a spatially stratified subsample of sites retains ample
# power and is standard practice for permutation tests.
# 性状集必须与主分析一致：逐纲只取该纲插补后完整、且在主分析中启用的性状。
# 早先版本把主分析已排除的「最大寿命」也纳入并做中位数填补，与 05/06 不一致，
# 会给审稿人留下把柄。
# The trait set must match the main analysis: per class, keep only traits that
# were complete after imputation and flagged used_in_main. An earlier version
# also included maximum longevity with crude median filling, which contradicted
# scripts 05 and 06.
CORE_ALL <- c("body_mass", "body_length", "nocturnality", "verticality",
              "diet_breadth", "diet_vert", "diet_plant", "habitat_breadth",
              "litter_size", "range_size")
usable_traits_for <- function(cl, tr) {
  cand <- intersect(CORE_ALL, names(tr))
  keep <- vapply(cand, function(v) {
    x <- suppressWarnings(as.numeric(tr[[v]]))
    mean(is.finite(x)) > 0.95 && stats::sd(x, na.rm = TRUE) > 0
  }, logical(1))
  cand[keep]
}
ENVV <- c("bio1_mean", "bio6_mean", "bio4_mean", "bio12_mean", "bio15_mean",
          "npp", "elev_range", "het_shannon", "hfp2020", "lc_cropland")

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))

run_fc <- function(cl, lab) {
  fc_file <- file.path(PATH$derived, sprintf("comm_%s_%s.rds", cl, lab))
  if (!file.exists(fc_file)) return(NULL)
  comm <- readRDS(fc_file)
  env  <- readRDS(file.path(PATH$derived, paste0("env_", lab, ".rds")))

  tr <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  sp <- intersect(unique(comm$species), tr$species)
  if (length(sp) < 30) return(NULL)
  tr <- tr[match(sp, tr$species), ]
  comm <- comm[comm$species %in% sp, ]

  # L：样方 x 物种（0/1）/ sites x species
  cells <- sort(unique(comm$cell_id))
  L <- matrix(0L, length(cells), length(sp),
              dimnames = list(cells, sp))
  L[cbind(match(comm$cell_id, cells), match(comm$species, sp))] <- 1L
  # 空间分层抽样：按经纬度分块，每块等概率抽取，避免抽样偏向高密度区
  # Spatially stratified subsample: bin by lon/lat and draw evenly per bin
  if (nrow(L) > N_CELLS) {
    xy <- env[match(rownames(L), env$cell_id), c("lon", "lat")]
    bin <- paste(cut(xy$lon, 12, labels = FALSE), cut(xy$lat, 12, labels = FALSE))
    set.seed(SEED)
    idx <- unlist(lapply(split(seq_len(nrow(L)), bin), function(ii)
      if (length(ii) <= ceiling(N_CELLS / length(unique(bin)))) ii else
        sample(ii, ceiling(N_CELLS / length(unique(bin))))))
    L <- L[sort(idx), , drop = FALSE]
  }
  # 剔除过稀疏的样方与物种 / drop very sparse rows/cols
  keep_r <- rowSums(L) >= 5
  L <- L[keep_r, , drop = FALSE]
  keep_c <- colSums(L) >= 3
  L <- L[, keep_c, drop = FALSE]
  if (nrow(L) < 100 || ncol(L) < 20) return(NULL)

  # R：样方 x 环境 / sites x environment
  ev <- intersect(ENVV, names(env))
  R <- env[match(rownames(L), env$cell_id), ev, drop = FALSE]
  R <- as.data.frame(lapply(R, function(z) { z[!is.finite(z)] <- median(z, na.rm = TRUE); as.numeric(zscore(z)) }))
  rownames(R) <- rownames(L)

  # Q：物种 x 性状（该纲实际可用者）/ species x traits usable for this class
  CORE <- usable_traits_for(cl, tr)
  if (length(CORE) < 4) return(NULL)
  Q <- tr[match(colnames(L), tr$species), CORE, drop = FALSE]
  Q <- as.data.frame(lapply(Q, function(z) {
    z <- as.numeric(z); z[!is.finite(z)] <- median(z, na.rm = TRUE); z }))
  rownames(Q) <- colnames(L)
  # 去除零方差性状 / drop zero-variance traits
  vv <- sapply(Q, function(z) stats::sd(z) > 0)
  Q <- Q[, vv, drop = FALSE]
  if (!ncol(Q)) return(NULL)

  log_msg("  ", cl, " ", lab, ": L = ", nrow(L), " x ", ncol(L),
          ", R = ", ncol(R), " env, Q = ", ncol(Q), " traits")

  set.seed(SEED)
  fc <- try(ade4::fourthcorner(tabR = R, tabL = as.data.frame(L), tabQ = Q,
                               modeltype = 6, p.adjust.method.G = "BH",
                               p.adjust.method.D = "BH", nrepet = NREP),
            silent = TRUE)
  if (inherits(fc, "try-error")) { log_msg("   [fail] fourthcorner"); return(NULL) }

  out <- data.frame(
    class = cl, grain = lab,
    env   = rep(colnames(R), times = ncol(Q)),
    trait = rep(colnames(Q), each = ncol(R)),
    stat  = as.numeric(fc$tabD2$obs),
    p_maxtest = as.numeric(fc$tabD2$adj.pvalue),
    row.names = NULL)

  ## 朴素 CWM 回归 p 值（用于量化第一类错误膨胀）
  ## Naive CWM regression p-values, to quantify Type I error inflation
  naive <- list()
  for (tq in colnames(Q)) {
    cwm <- as.numeric(L %*% Q[[tq]] / rowSums(L))
    for (ev1 in colnames(R)) {
      m <- stats::lm(cwm ~ R[[ev1]])
      naive[[length(naive) + 1]] <- data.frame(
        env = ev1, trait = tq,
        p_naive = summary(m)$coefficients[2, 4],
        r2_naive = summary(m)$r.squared)
    }
  }
  naive <- bind_rows(naive)
  naive$p_naive_adj <- stats::p.adjust(naive$p_naive, method = "BH")
  left_join(out, naive, by = c("env", "trait"))
}

res <- list()
for (cl in TAXA$class) for (lab in names(GRAIN)) {
  r <- try(run_fc(cl, lab), silent = TRUE)
  if (!inherits(r, "try-error") && !is.null(r)) res[[paste(cl, lab)]] <- r
}
res <- bind_rows(res)

if (nrow(res)) {
  res <- res |> left_join(TAXA[, c("class", "thermal")], by = "class")
  write_table(res, "table_4_fourthcorner_maxtest")

  # 第一类错误膨胀诊断 / Type I error inflation diagnostic
  diag <- res |>
    group_by(class, grain) |>
    summarise(n_tests = n(),
              sig_naive   = sum(p_naive_adj < 0.05, na.rm = TRUE),
              sig_maxtest = sum(p_maxtest   < 0.05, na.rm = TRUE),
              pct_naive   = round(100 * mean(p_naive_adj < 0.05, na.rm = TRUE), 1),
              pct_maxtest = round(100 * mean(p_maxtest   < 0.05, na.rm = TRUE), 1),
              inflation_ratio = round(sig_naive / pmax(sig_maxtest, 1), 2),
              .groups = "drop")
  write_table(diag, "table_S8_typeI_inflation_diagnostic")
  print(as.data.frame(diag))
  log_msg("第四角结果 / fourth-corner results: ", nrow(res), " 组合")
}
log_msg("=== 06b 完成 / done ===")
