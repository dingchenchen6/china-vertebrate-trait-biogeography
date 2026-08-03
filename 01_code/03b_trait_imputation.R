# ============================================================
# 科学问题 / Scientific question:
#   性状缺失若与类群或生态类型非随机相关，会系统性扭曲群落性状格局。
#   Non-random trait missingness would systematically distort assemblage
#   trait patterns, so it must be imputed and the imputation validated.
#
# 分析目标 / Objective:
#   用随机森林 (missForest) 结合系统发育信息（系统发育特征向量）或
#   分类阶元信息插补缺失性状；报告袋外误差与插补前后覆盖率。
#
# 输入数据 / Input:
#   03_data_derived/traits_china_raw.rds；02_data_raw/phylo/*.tre
#
# 主要流程 / Workflow:
#   1. 逐类群提取系统发育特征向量 (PVR) 前 10 轴；无树则用目/科哑变量
#   2. missForest 插补（性状 + PVR 一起进入）
#   3. 记录 OOB 归一化均方根误差 (NRMSE) 与分类误差 (PFC)
#   4. 保留原始值与插补标记，供敏感性分析
#
# 预期输出 / Expected output:
#   03_data_derived/traits_imputed.rds
#   04_results/tables/table_S6_imputation_diagnostics.csv
#
# 关键假设 / Key assumptions:
#   - 性状具系统发育信号，故近缘种可互相提供信息
#   - 缺失率 > 40% 的性状不参与主分析（仅作敏感性）
#
# 主要包 / Main packages: missForest, ape, PVR/vegan, dplyr
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({ library(missForest); library(ape) })

log_msg("=== 03b 性状插补 / Trait imputation ===")

CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
          "diet_breadth", "diet_vert", "diet_plant", "habitat_breadth",
          "litter_size", "max_longevity", "range_size")
MAX_MISS <- 0.40   # 参与主分析的最大缺失率 / max missingness for main analysis
N_PVR    <- 10L    # 系统发育特征向量数 / number of phylogenetic eigenvectors

tr <- readRDS(file.path(PATH$derived, "traits_china_raw.rds"))
log_msg("输入性状表 / input: ", nrow(tr), " 种 x ", ncol(tr), " 列")

# ---------------------------------------------------------------
# 1. 系统发育树登记（按类群）/ Phylogeny registry by class
# ---------------------------------------------------------------
# 四类群各一棵定年树 / one dated tree per class
#   两栖 Jetz & Pyron 2018；哺乳 Upham et al. 2019；
#   鸟类 Jetz et al. 2012 (Hackett backbone)；有鳞目 Tonini et al. 2016
TREE_FILES <- list(
  Amphibia = file.path(PATH$raw_phylo, "amph_shl_new_Consensus_7238.tre"),
  Mammalia = c(file.path(PATH$raw_phylo, "MamPhy_Completed_5911sp_MCC_v2_target_NEWICK.tre"),
               file.path(PATH$raw_phylo, "MamPhy_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")),
  Aves     = c(file.path(PATH$raw_phylo, "bird_Jetz2012_AllBirdsHackett1_tree1_vertlife.tre"),
               file.path(PATH$raw_phylo, "bird_Jetz2012_Hackett_tree1.tre")),
  Reptilia = c(file.path(PATH$raw_phylo, "squamate_Tonini_tree1of100.tre"))
)

load_tree <- function(paths) {
  for (p in paths) {
    if (!file.exists(p)) next
    t <- try(ape::read.tree(p), silent = TRUE)
    if (inherits(t, "try-error") || is.null(t)) t <- try(ape::read.nexus(p), silent = TRUE)
    if (inherits(t, "try-error") || is.null(t)) next
    if (inherits(t, "multiPhylo")) t <- t[[1]]
    return(t)
  }
  NULL
}

#' 标准化树叶标签为 "Genus_species"
#' Normalise tip labels to "Genus_species"
#' Upham et al. 2019 的标签形如 Zaglossus_attenboroughi_TACHYGLOSSIDAE_MONOTREMATA
clean_tips <- function(tree) {
  lab <- sub("^_+", "", tree$tip.label)
  lab <- sub("^([A-Za-z]+_[A-Za-z-]+).*$", "\\1", lab)
  tree$tip.label <- lab
  tree <- ape::drop.tip(tree, which(duplicated(lab)))
  tree
}

#' 系统发育特征向量 / Phylogenetic eigenvectors (PVR)
pvr_axes <- function(tree, spp, k = N_PVR) {
  tree <- clean_tips(tree)
  tips <- gsub(" ", "_", spp)
  keep <- intersect(tips, tree$tip.label)
  if (length(keep) < 20) return(NULL)
  st <- ape::keep.tip(tree, keep)
  D  <- ape::cophenetic.phylo(st)
  # 主坐标分析 -> 前 k 轴 / PCoA on phylogenetic distances
  pc <- try(stats::cmdscale(as.dist(D), k = min(k, nrow(D) - 1)), silent = TRUE)
  if (inherits(pc, "try-error")) return(NULL)
  colnames(pc) <- paste0("pvr", seq_len(ncol(pc)))
  out <- as.data.frame(pc)
  out$species <- gsub("_", " ", rownames(pc))
  out
}

# ---------------------------------------------------------------
# 2. 逐类群插补 / Impute per class
# ---------------------------------------------------------------
diag_list <- list(); out_list <- list()

for (cl in unique(tr$class)) {
  d <- tr[tr$class == cl, ]
  if (nrow(d) < 30) { out_list[[cl]] <- d; next }
  log_msg("-- ", cl, " (", nrow(d), " 种) --")

  # 缺失率筛选 / drop traits that are too incomplete
  miss <- sapply(d[, CORE, drop = FALSE], function(x) mean(is.na(x)))
  use  <- names(miss)[miss <= MAX_MISS]
  drop <- names(miss)[miss > MAX_MISS]
  if (length(drop))
    log_msg("   缺失率过高，排除 / dropped (>", MAX_MISS*100, "%): ",
            paste0(drop, " (", round(100*miss[drop]), "%)", collapse = ", "))
  if (!length(use)) { out_list[[cl]] <- d; next }

  # 系统发育信息 / phylogenetic information
  tree <- load_tree(TREE_FILES[[cl]])
  pv <- if (!is.null(tree)) pvr_axes(tree, d$species) else NULL
  if (!is.null(pv)) {
    d2 <- dplyr::left_join(d, pv, by = "species")
    pvcols <- grep("^pvr", names(d2), value = TRUE)
    cov_pv <- mean(!is.na(d2[[pvcols[1]]]))
    log_msg("   系统发育特征向量 / PVR: ", length(pvcols), " 轴, 覆盖 ",
            round(100 * cov_pv), "% 物种")
    # PVR 自身缺失用 0 填补（等价于置于系统发育中心）
    for (p in pvcols) d2[[p]][is.na(d2[[p]])] <- 0
  } else {
    d2 <- d; pvcols <- character(0)
    log_msg("   无可用系统发育树，改用分类阶元信息 / no tree; using taxonomy")
  }

  # 分类阶元作为辅助预测变量 / taxonomy as auxiliary predictors
  d2$order_f  <- factor(ifelse(is.na(d2$order),  "unk", d2$order))
  d2$family_f <- factor(ifelse(is.na(d2$family), "unk", d2$family))
  # 水平过多的因子会拖慢随机森林 / collapse rare levels for speed
  collapse_rare <- function(f, min_n = 5) {
    tb <- table(f); keep <- names(tb)[tb >= min_n]
    factor(ifelse(as.character(f) %in% keep, as.character(f), "other"))
  }
  d2$order_f  <- collapse_rare(d2$order_f)
  d2$family_f <- collapse_rare(d2$family_f)

  X <- d2[, c(use, pvcols, "order_f", "family_f"), drop = FALSE]
  X <- X[, sapply(X, function(z) length(unique(z[!is.na(z)])) > 1), drop = FALSE]

  set.seed(SEED)
  fit <- try(missForest::missForest(as.data.frame(X), maxiter = 8, ntree = 200,
                                    variablewise = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    log_msg("   [warn] missForest 失败，保留原值 / imputation failed, keeping raw")
    out_list[[cl]] <- d; next
  }

  imp <- fit$ximp
  err <- fit$OOBerror
  log_msg("   OOB 误差 / OOB error: ",
          paste(names(err), round(err, 4), collapse = "; "))

  # 写回：记录哪些值被插补 / write back with imputation flags
  res <- d
  for (v in use) {
    flag <- is.na(res[[v]])
    res[[v]] <- imp[[v]]
    res[[paste0("imp_", v)]] <- flag
  }
  out_list[[cl]] <- res

  diag_list[[cl]] <- data.frame(
    class = cl, n_species = nrow(d),
    trait = names(miss),
    missing_pct_before = round(100 * miss, 1),
    used_in_main = names(miss) %in% use,
    imputed_pct = round(100 * sapply(names(miss), function(v)
      if (v %in% use) mean(is.na(d[[v]])) else NA_real_), 1),
    OOB_NRMSE = round(unname(err["NRMSE"]), 4),
    row.names = NULL)
}

traits_imp <- dplyr::bind_rows(out_list)
saveRDS(traits_imp, file.path(PATH$derived, "traits_imputed.rds"))
log_msg("插补后性状表 / imputed traits: ", nrow(traits_imp), " 种")

diag_df <- dplyr::bind_rows(diag_list)
write_table(diag_df, "table_S6_imputation_diagnostics")
print(as.data.frame(diag_df))

# 插补后覆盖率 / coverage after imputation
cov_after <- traits_imp |>
  select(class, any_of(CORE)) |>
  tidyr::pivot_longer(-class, names_to = "trait", values_to = "v") |>
  group_by(class, trait) |>
  summarise(coverage_after = round(100 * mean(!is.na(v)), 1), .groups = "drop")
write_table(cov_after, "table_S6b_trait_coverage_after_imputation")

log_msg("=== 03b 完成 / done ===")
