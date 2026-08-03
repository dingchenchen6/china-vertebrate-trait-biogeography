# ============================================================
# 科学问题 / Scientific question:
#   群落性状分布的四阶矩（均值、方差、偏度、峰度）能把不同的
#   群落构建过程分离开：均值=方向性过滤；方差=过滤强度 vs 限制相似性；
#   偏度=哪一侧尾部被截断；峰度=生态位填充程度。
#   The first four moments of the assemblage trait distribution separate
#   distinct assembly processes: mean = directional filtering; variance =
#   filtering strength vs limiting similarity; skewness = which tail is
#   truncated; kurtosis = niche packing.
#
# 分析目标 / Objective:
#   为四个类群 × 两个粒度 × 每个核心性状，计算群落性状矩与多元功能
#   多样性指标，并用丰富度约束零模型计算标准化效应量 (SES)。
#
# 输入数据 / Input:
#   03_data_derived/comm_<class>_<grain>.rds, traits_imputed.rds
#
# 主要流程 / Workflow:
#   1. 观测值：每格每性状的 CWM / CWV / CWS / CWK
#   2. 多元指标：PCoA -> FRic(凸包体积), FDis, RaoQ, 性状体积/密度
#   3. 零模型：从类群物种库随机抽 S 种，重复 999 次
#      *关键优化*：零分布只依赖丰富度 S，故按唯一 S 值计算一次即可，
#      精确且比逐格模拟快 1-2 个数量级
#   4. SES = (obs - mean_null) / sd_null
#
# 预期输出 / Expected output:
#   03_data_derived/metrics_<class>_<grain>.rds
#   04_results/tables/table_S5_community_metrics_summary.csv
#
# 关键假设 / Key assumptions:
#   - 存在-缺失数据下"群落加权"= 等权，即对出现物种取简单矩
#   - 偏度/峰度需 S>=5；FRic 凸包需 S>=k+1（k 为 PCoA 轴数）
#   - 零模型固定丰富度后从物种库抽样。抽样权重的选择会实质改变结论：
#     * 等概率抽样把狭域特有种与广布种视为等可能，而真实群落由广布种主导；
#       若广布种功能上更相似，就会产生**虚假的功能聚集**。
#     * 按物种占据网格数加权（本项目的**主零模型**）消除这一来源。
#     实测：加权后平均 SES.FDis 由鸟 -1.67 -> -0.83、兽 -1.24 -> -0.91、
#     爬行 -1.87 -> -1.10，即等概率零模型把聚集夸大了 27-50%。
#     两套零模型的结果都输出，主分析用加权版。
#   - The choice of sampling weights materially changes the answer. Drawing
#     species equiprobably treats narrow endemics and widespread species alike,
#     whereas real assemblages are dominated by widespread species; if those are
#     functionally more alike, spurious clustering follows. Weighting by the
#     number of cells a species occupies removes that source and is used as the
#     primary null here. Both are computed and saved.
#
# 主要包 / Main packages: dplyr, data.table, ape, geometry, vegan
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({
  library(ape); library(geometry); library(cluster); library(matrixStats); library(parallel)
})

log_msg("=== 05 群落性状指标 / Assemblage trait metrics ===")

N_PCOA   <- 4L    # 用于多元指标的 PCoA 轴数 / PCoA axes for multivariate metrics
MIN_S_MOM <- 5L   # 计算矩的最小丰富度 / minimum richness for moments
N_CORES  <- max(1L, parallel::detectCores() - 2L)  # 并行核数 / cores for null model

# 核心性状集：四个类群定义一致且覆盖完整，用于跨类群正式对比检验
# Core set: identically defined and fully covered in all four classes;
# used for the formal cross-taxon contrasts (H1, H2c)
CORE_ALL4 <- c("body_mass", "body_length", "nocturnality", "verticality",
               "habitat_breadth", "range_size")
# 扩展性状集：鸟/兽/爬行可用（两栖数据不足），用于类群内分析
# Extended set: available for birds/mammals/reptiles only
EXTENDED  <- c("diet_breadth", "diet_vert", "diet_plant", "litter_size")

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))

#' 该类群实际可用的性状（插补后无缺失且有方差）
#' Traits actually usable for a class (complete after imputation, non-constant)
usable_traits <- function(tr) {
  cand <- c(CORE_ALL4, EXTENDED)
  cand <- cand[cand %in% names(tr)]
  keep <- vapply(cand, function(v) {
    x <- suppressWarnings(as.numeric(tr[[v]]))
    mean(is.finite(x)) > 0.95 && stats::sd(x, na.rm = TRUE) > 0
  }, logical(1))
  cand[keep]
}

# ---------------------------------------------------------------
# 矩函数（无偏估计）/ Moment functions (unbiased where standard)
# ---------------------------------------------------------------
m_mean <- function(x) mean(x)
m_var  <- function(x) stats::var(x)
m_skew <- function(x) {
  n <- length(x); s <- stats::sd(x)
  if (!is.finite(s) || s == 0 || n < 3) return(NA_real_)
  sum((x - mean(x))^3) / n / s^3
}
m_kurt <- function(x) {
  n <- length(x); s <- stats::sd(x)
  if (!is.finite(s) || s == 0 || n < 4) return(NA_real_)
  sum((x - mean(x))^4) / n / s^4 - 3   # 超额峰度 / excess kurtosis
}

#' 计算一个群落（物种索引向量）的全部指标
#' Compute all metrics for one assemblage (vector of species indices)
#'
#' 性能说明 / Performance note:
#'   物种两两距离矩阵 DM 在类群层面预计算一次后按索引取子集，
#'   避免在 999 x n_richness 次零模型中重复计算 O(S^2) 距离，
#'   这是整个流程最关键的加速点。
#'   The full species-by-species distance matrix DM is precomputed once per
#'   class and subset by index, avoiding recomputation of O(S^2) distances
#'   inside the 999 x n_richness null loop.
calc_metrics <- function(idx, TMAT, PCO, DM) {
  S <- length(idx)
  out <- c(richness = S)

  ## 单性状矩：矩阵化一次算完所有性状（原逐性状循环是主要瓶颈）
  ## Per-trait moments, vectorised across all traits at once
  Xs  <- TMAT[idx, , drop = FALSE]
  mu  <- colMeans(Xs)
  d   <- Xs - rep(mu, each = S)
  d2  <- d * d
  ss  <- colSums(d2)
  sdp <- sqrt(ss / S)                       # 总体标准差 / population sd
  vr  <- ss / (S - 1)                       # 样本方差 / sample variance
  ok  <- sdp > 0
  sk  <- rep(NA_real_, ncol(Xs)); ku <- sk
  if (S >= 3 && any(ok)) sk[ok] <- colSums(d2 * d)[ok] / S / sdp[ok]^3
  if (S >= 4 && any(ok)) ku[ok] <- colSums(d2 * d2)[ok] / S / sdp[ok]^4 - 3
  if (S < MIN_S_MOM) { mu[] <- NA_real_; vr[] <- NA_real_; sk[] <- NA_real_; ku[] <- NA_real_ }
  nm <- colnames(TMAT)
  out[paste0(nm, ".mean")] <- mu
  out[paste0(nm, ".var")]  <- vr
  out[paste0(nm, ".skew")] <- sk
  out[paste0(nm, ".kurt")] <- ku

  ## 多元功能指标 / multivariate functional metrics
  P   <- PCO[idx, , drop = FALSE]
  ctr <- colMeans(P)
  Pd  <- P - rep(ctr, each = S)
  out["FDis"]      <- mean(sqrt(rowSums(Pd * Pd)))          # 功能离散度
  out["trait_vol"] <- sum(colSums(Pd * Pd) / (S - 1))       # 性状体积（方差和）
  if (S >= 3) {
    dm <- DM[idx, idx]
    out["RaoQ"]       <- sum(dm) / (S * (S - 1))            # Rao 二次熵
    diag(dm) <- Inf
    out["trait_dens"] <- mean(matrixStats::rowMins(dm))     # 性状密度（最近邻距离）
  } else {
    out["trait_dens"] <- NA_real_; out["RaoQ"] <- NA_real_
  }
  ## 功能丰富度 = 凸包体积 / functional richness = convex hull volume
  if (S > ncol(P)) {
    v <- try(geometry::convhulln(P, options = "FA")$vol, silent = TRUE)
    out["FRic"] <- if (inherits(v, "try-error")) NA_real_ else v
  } else out["FRic"] <- NA_real_
  out
}

# ---------------------------------------------------------------
# 零模型：按唯一丰富度值计算（精确且高效）
# Null model: computed per unique richness value (exact and efficient)
# ---------------------------------------------------------------
#' @param w 抽样权重；NULL 表示等概率 / sampling weights, NULL means equiprobable
build_null <- function(S_vals, TMAT, PCO, DM, n_rep = N_NULL, w = NULL) {
  pool <- nrow(TMAT)
  tmpl <- calc_metrics(seq_len(min(5, pool)), TMAT, PCO, DM)
  S_vals <- S_vals[S_vals >= 2 & S_vals <= pool]

  one_S <- function(S) {
    # 每个丰富度值用确定性的独立种子，保证并行下结果可复现
    # Deterministic per-S seed keeps parallel results reproducible
    set.seed(SEED + S)
    sim <- vapply(seq_len(n_rep), function(i)
      calc_metrics(if (is.null(w)) sample.int(pool, S)
                   else sample.int(pool, S, prob = w), TMAT, PCO, DM),
      FUN.VALUE = tmpl)
    # rowSds 不保留行名，需显式补上（后续按指标名索引）
    # rowSds drops names; restore them since SES indexes by metric name
    list(mu = setNames(rowMeans(sim, na.rm = TRUE), names(tmpl)),
         sd = setNames(matrixStats::rowSds(sim, na.rm = TRUE), names(tmpl)))
  }

  res <- parallel::mclapply(S_vals, one_S, mc.cores = N_CORES)
  names(res) <- as.character(S_vals)
  # mclapply 出错时回退到串行 / fall back to serial if a worker failed
  bad <- !vapply(res, is.list, logical(1))
  if (any(bad)) for (i in which(bad)) res[[i]] <- one_S(S_vals[i])
  res
}

# ---------------------------------------------------------------
# 主流程 / Main loop
# ---------------------------------------------------------------
run_class_grain <- function(cl, lab) {
  f <- file.path(PATH$derived, sprintf("comm_%s_%s.rds", cl, lab))
  if (!file.exists(f)) { log_msg("  [skip] ", cl, " ", lab); return(NULL) }
  comm <- readRDS(f)

  tr <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  tr <- tr[!is.na(tr$species), ]
  sp_keep <- intersect(unique(comm$species), tr$species)
  if (length(sp_keep) < 20) { log_msg("  [skip] ", cl, " ", lab, " 物种不足"); return(NULL) }

  tr <- tr[match(sp_keep, tr$species), ]
  comm <- comm[comm$species %in% sp_keep, ]

  # 该类群可用性状 / traits usable for this class
  use <- usable_traits(tr)
  if (length(use) < 4) { log_msg("  [skip] ", cl, " ", lab, " 可用性状不足"); return(NULL) }

  # 性状矩阵：z 标准化以便跨性状可比 / z-standardise traits for comparability
  TMAT <- as.matrix(sapply(tr[, use, drop = FALSE], function(z) suppressWarnings(as.numeric(z))))
  rownames(TMAT) <- tr$species
  TMAT <- apply(TMAT, 2, zscore)
  # 多元指标只用四类群通用核心性状，保证跨类群可比
  # Multivariate metrics use only the common core set for comparability
  core_use <- intersect(CORE_ALL4, use)

  # PCoA（Gower 距离）/ PCoA on Gower distance of the common core traits
  gd  <- cluster::daisy(as.data.frame(TMAT[, core_use, drop = FALSE]), metric = "gower")
  pc  <- ape::pcoa(as.dist(gd))
  PCO <- pc$vectors[, seq_len(min(N_PCOA, ncol(pc$vectors))), drop = FALSE]
  rownames(PCO) <- tr$species
  # 预计算全物种功能距离矩阵（关键加速）/ precompute the full distance matrix
  DM <- as.matrix(stats::dist(PCO))
  var_expl <- round(100 * pc$values$Relative_eig[seq_len(ncol(PCO))], 1)
  log_msg("  ", cl, " ", lab, ": ", length(sp_keep), " 种, ", length(use),
          " 性状 (核心 ", length(core_use), "), PCoA 解释度 ",
          paste0(var_expl, "%", collapse = ", "))

  # 群落 -> 物种索引 / assemblages -> species indices
  sp_idx <- setNames(seq_len(nrow(TMAT)), rownames(TMAT))
  by_cell <- split(sp_idx[comm$species], comm$cell_id)
  by_cell <- by_cell[lengths(by_cell) >= 2]

  # 观测值 / observed
  obs <- t(vapply(by_cell, calc_metrics, TMAT = TMAT, PCO = PCO, DM = DM,
                  FUN.VALUE = calc_metrics(1:5, TMAT, PCO, DM)))
  obs <- as.data.frame(obs); obs$cell_id <- names(by_cell)

  # 零模型 / null model
  # 主零模型按物种占据的网格数加权。等概率抽样会把狭域特有种与广布种
  # 视为等可能，而真实群落由广布种主导；若广布种功能上更相似，就会读出
  # 虚假的功能聚集。加权版是保守且更可辩护的选择，等概率版作为敏感性保留。
  # The primary null weights species by the number of cells they occupy.
  # Equiprobable draws treat narrow endemics like widespread species and
  # manufacture apparent clustering; the weighted null is the conservative
  # choice, with the equiprobable version retained as a sensitivity.
  S_vals <- sort(unique(lengths(by_cell)))
  log_msg("    零模型：", length(S_vals), " 个唯一丰富度值 x ", N_NULL, " 次 x 2 套")
  occ   <- as.numeric(table(factor(comm$species, levels = rownames(TMAT))))
  w_occ <- occ / sum(occ)
  nul    <- build_null(S_vals, TMAT, PCO, DM, w = w_occ)   # 主 / primary
  nul_eq <- build_null(S_vals, TMAT, PCO, DM)              # 敏感性 / sensitivity

  # SES
  mets <- setdiff(names(obs), c("cell_id", "richness"))
  ses    <- matrix(NA_real_, nrow(obs), length(mets), dimnames = list(NULL, mets))
  ses_eq <- ses
  Svec <- as.character(obs$richness)
  for (i in seq_len(nrow(obs))) {
    ov <- as.numeric(obs[i, mets])
    nn <- nul[[Svec[i]]]
    if (!is.null(nn)) ses[i, ] <- (ov - nn$mu[mets]) / nn$sd[mets]
    ne <- nul_eq[[Svec[i]]]
    if (!is.null(ne)) ses_eq[i, ] <- (ov - ne$mu[mets]) / ne$sd[mets]
  }
  ses <- as.data.frame(ses);       names(ses)    <- paste0("SES.", mets)
  ses_eq <- as.data.frame(ses_eq); names(ses_eq) <- paste0("SESeq.", mets)

  out <- cbind(obs[, c("cell_id", "richness", mets)], ses, ses_eq)
  out$class <- cl; out$grain <- lab
  saveRDS(out, file.path(PATH$derived, sprintf("metrics_%s_%s.rds", cl, lab)))
  log_msg("    -> metrics_", cl, "_", lab, ".rds (", nrow(out), " 格 x ",
          ncol(out), " 列)")
  out
}

all_out <- list()
for (cl in TAXA$class) for (lab in names(GRAIN)) {
  r <- try(run_class_grain(cl, lab), silent = FALSE)
  if (!inherits(r, "try-error") && !is.null(r)) all_out[[paste(cl, lab)]] <- r
}

if (length(all_out)) {
  summ <- lapply(all_out, function(d)
    data.frame(class = d$class[1], grain = d$grain[1], n_cells = nrow(d),
               mean_richness = round(mean(d$richness), 1),
               mean_FDis = round(mean(d$FDis, na.rm = TRUE), 3),
               mean_SES_FDis = round(mean(d$SES.FDis, na.rm = TRUE), 3),
               mean_SES_FRic = round(mean(d$SES.FRic, na.rm = TRUE), 3),
               mean_SES_FDis_equiprob = round(mean(d$SESeq.FDis, na.rm = TRUE), 3),
               mean_SES_FRic_equiprob = round(mean(d$SESeq.FRic, na.rm = TRUE), 3))) |>
  bind_rows()
  write_table(summ, "table_S5_community_metrics_summary")
  print(summ)
}
log_msg("=== 05 完成 / done ===")
