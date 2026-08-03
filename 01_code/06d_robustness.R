# ============================================================
# 科学问题 / Scientific question:
#   核心结论（变温类群性状–环境耦合更强）是否依赖于某一个纲？
#   两栖类的功能超离散是否只是上游性状插补的产物？
#   Does the headline result depend on any single class, and is the amphibian
#   over-dispersion an artefact of upstream trait imputation?
#
# 分析目标 / Objective:
#   1. 留一纲检验（leave-one-class-out）：逐次剔除一个纲后重做 H1 检验
#   2. 上游插补敏感性：仅用实测（非 TetrapodTraits 插补）性状值的物种
#      重算功能离散度 SES
#   3. 输出可直接进正文/附录的稳健性表
#
# 输入数据 / Input:
#   04_results/tables/table_7_family_coupling_strength.csv
#   03_data_derived/traits_imputed.rds, comm_*_50km.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_8_leave_one_class_out.csv
#   04_results/tables/table_9_imputation_sensitivity.csv
#
# 关键假设 / Key assumptions:
#   - 留一纲后每组仅剩 7-8 个科，检验力必然大幅下降；因此关注的是
#     **效应方向与比值是否保持**，而非 P 值是否仍然显著。
#   - 实测子集本身有偏（研究充分、体型较大的物种更可能有实测值），
#     故该检验只能判断"结论是否可能被插补驱动"，不能替代主分析。
#
# 主要包 / Main packages: lme4, lmerTest, cluster, ape, matrixStats
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(cluster); library(ape)
  library(matrixStats); library(parallel)
})

log_msg("=== 06d 稳健性检验 / Robustness checks ===")

N_CORES <- max(1L, parallel::detectCores() - 2L)
CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
          "habitat_breadth", "range_size")

# ---------------------------------------------------------------
# 1. 留一纲检验 / Leave-one-class-out
# ---------------------------------------------------------------
fam <- data.table::fread(file.path(PATH$tables, "table_7_family_coupling_strength.csv"))
fam$thermal <- factor(fam$thermal, levels = c("Endotherm", "Ectotherm"))

one_test <- function(d, metric) {
  dd <- d[is.finite(d[[metric]]), ]; dd$y <- dd[[metric]]
  if (length(unique(dd$thermal)) < 2) return(NULL)
  w <- stats::wilcox.test(y ~ thermal, data = dd, exact = FALSE)
  m <- try(lmerTest::lmer(y ~ thermal + (1 | class/order), data = dd), silent = TRUE)
  pm <- if (inherits(m, "try-error")) NA_real_ else {
    s <- summary(m)$coefficients; pc <- grep("^Pr\\(", colnames(s))
    if (length(pc)) s["thermalEctotherm", pc[1]] else NA_real_
  }
  e <- mean(dd$y[dd$thermal == "Endotherm"]); k <- mean(dd$y[dd$thermal == "Ectotherm"])
  data.frame(metric = metric,
             n_endo = sum(dd$thermal == "Endotherm"),
             n_ecto = sum(dd$thermal == "Ectotherm"),
             endo_mean = round(e, 4), ecto_mean = round(k, 4),
             ratio = round(k / e, 2),
             p_wilcoxon = signif(w$p.value, 3),
             p_mixed = if (is.na(pm)) NA else signif(pm, 3),
             row.names = NULL)
}

loo <- list()
scenarios <- c(list(`Full data set (33 families)` = NULL),
               setNames(as.list(TAXA$class),
                        paste0("Excluding ", TAXA$class)))
for (nm in names(scenarios)) {
  drop_cl <- scenarios[[nm]]
  d <- if (is.null(drop_cl)) fam else fam[fam$class != drop_cl, ]
  for (mt in c("coupling", "coupling_thermal", "coupling_human")) {
    r <- one_test(d, mt)
    if (is.null(r)) next
    r$scenario <- nm; r$class_dropped <- if (is.null(drop_cl)) NA else drop_cl
    loo[[length(loo) + 1]] <- r
  }
}
loo <- bind_rows(loo) |>
  select(scenario, class_dropped, metric, n_endo, n_ecto,
         endo_mean, ecto_mean, ratio, p_wilcoxon, p_mixed)
write_table(loo, "table_8_leave_one_class_out")
log_msg("--- 留一纲检验 / leave-one-class-out ---")
print(as.data.frame(loo))
log_msg("方向保持（比值 > 1）的比例: ",
        round(100 * mean(loo$ratio > 1, na.rm = TRUE), 1), "%")

# ---------------------------------------------------------------
# 2. 上游插补敏感性 / Upstream-imputation sensitivity
# ---------------------------------------------------------------
# TetrapodTraits 对部分性状提供了 ImputedX 标记（1 = 该值为模型插补）。
# 两栖类的体重与夜行性中分别有 75% 与 71% 是插补值，因此必须检验
# "功能超离散"是否只在插补值上成立。
# TetrapodTraits flags which values it imputed. For amphibians, 75% of body
# masses and 71% of nocturnality values are imputed, so the over-dispersion
# result must be tested against the observed-only subset.
traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))
FLAGS <- c(body_mass = "flag_mass", body_length = "flag_len",
           nocturnality = "flag_act", verticality = "flag_hab",
           habitat_breadth = "flag_mhab")

fdis_ses <- function(tr, comm, n_rep = 499L) {
  sp <- intersect(unique(comm$species), tr$species)
  if (length(sp) < 50) return(NULL)
  tr <- tr[match(sp, tr$species), ]
  TM <- apply(as.matrix(sapply(tr[, CORE], as.numeric)), 2, zscore)
  rownames(TM) <- sp
  if (any(!is.finite(TM))) return(NULL)
  PCO <- ape::pcoa(as.dist(cluster::daisy(as.data.frame(TM), metric = "gower")))$vectors[, 1:4]
  cm <- comm[comm$species %in% sp, ]
  by_cell <- split(match(cm$species, sp), cm$cell_id)
  by_cell <- by_cell[lengths(by_cell) >= 5]
  if (length(by_cell) < 100) return(NULL)
  fdis <- function(idx) {
    P <- PCO[idx, , drop = FALSE]; ctr <- colMeans(P)
    mean(sqrt(rowSums((P - rep(ctr, each = length(idx)))^2)))
  }
  obs <- vapply(by_cell, fdis, numeric(1)); S <- lengths(by_cell)
  uS <- sort(unique(S))
  nul <- parallel::mclapply(uS, function(k) {
    set.seed(SEED + k)
    v <- replicate(n_rep, fdis(sample.int(nrow(TM), k))); c(mean(v), stats::sd(v))
  }, mc.cores = N_CORES)
  names(nul) <- as.character(uS)
  ses <- vapply(seq_along(obs), function(i) {
    n <- nul[[as.character(S[i])]]
    if (!is.numeric(n) || length(n) < 2) return(NA_real_)
    (obs[i] - n[1]) / n[2]
  }, numeric(1))
  list(n_species = nrow(TM), n_cells = length(obs),
       mean_ses = mean(ses, na.rm = TRUE),
       pct_positive = 100 * mean(ses > 0, na.rm = TRUE))
}

sens <- list()
for (cl in TAXA$class) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) next
  comm <- readRDS(f)
  tr_all <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)

  # 仅保留核心性状全部为实测值的物种 / keep species with no upstream imputation
  keep <- rep(TRUE, nrow(tr_all))
  for (fl in FLAGS) if (fl %in% names(tr_all))
    keep <- keep & !(suppressWarnings(as.numeric(tr_all[[fl]])) %in% 1)
  tr_obs <- tr_all[keep, ]

  a <- fdis_ses(tr_all, comm); b <- fdis_ses(tr_obs, comm)
  sens[[length(sens) + 1]] <- data.frame(
    class = cl,
    n_species_all = nrow(tr_all), n_species_observed_only = nrow(tr_obs),
    pct_species_retained = round(100 * nrow(tr_obs) / nrow(tr_all), 1),
    mean_SES_all = if (is.null(a)) NA else round(a$mean_ses, 2),
    mean_SES_observed_only = if (is.null(b)) NA else round(b$mean_ses, 2),
    pct_positive_all = if (is.null(a)) NA else round(a$pct_positive, 1),
    pct_positive_observed_only = if (is.null(b)) NA else round(b$pct_positive, 1),
    n_cells_observed_only = if (is.null(b)) NA else b$n_cells,
    row.names = NULL)
  log_msg("  ", cl, ": ", nrow(tr_obs), "/", nrow(tr_all), " 种为实测")
}
sens <- bind_rows(sens)
write_table(sens, "table_9_imputation_sensitivity")
log_msg("--- 上游插补敏感性 / upstream-imputation sensitivity ---")
print(as.data.frame(sens))

log_msg("=== 06d 完成 / done ===")
