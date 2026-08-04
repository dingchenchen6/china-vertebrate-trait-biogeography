# ============================================================
# 科学问题 / Scientific question:
#   H1「恒温 vs 变温」的对比若以「纲」为重复单元只有 4 个重复；
#   若以「目」为单元，中国变温类群仅 ANURA 与 SQUAMATA 两个目达标，
#   设计严重失衡。因此本脚本以「科」为重复单元——两侧数量接近，
#   并以 纲/目 嵌套随机效应吸收系统发育非独立性。
#   CLASS as the replicate gives n = 4; ORDER is badly unbalanced because
#   Chinese ectotherms are concentrated in just two orders. FAMILY gives a
#   far more balanced design, with class/order nested random effects
#   absorbing phylogenetic non-independence.
#
# 分析目标 / Objective:
#   1. 为每个物种数充足的目，独立构建群落矩阵并计算性状矩与 FD
#   2. 逐目拟合空间自回归模型，提取标准化环境效应
#   3. 以「目」为重复，检验 效应强度 ~ 体温调节模式（纲为随机效应，
#      控制目在纲内的系统发育非独立性）
#
# 输入数据 / Input:
#   03_data_derived/comm_<class>_50km.rds, traits_imputed.rds, env_50km.rds
#   03_data_derived/splist_<class>.rds（目/科信息）
#
# 预期输出 / Expected output:
#   04_results/tables/table_5_clade_level_effects.csv
#   04_results/tables/table_6_thermal_mode_contrast.csv
#
# 关键假设 / Key assumptions:
#   - 目内物种库是该目群落构建的合理零模型参照系
#   - 目需 >= MIN_SP_CLADE 个中国物种才纳入，保证矩与 FD 可估
#   - 纲作为随机效应部分吸收目间的系统发育非独立性
#
# 主要包 / Main packages: spdep, spatialreg, lme4, matrixStats
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({
  library(spdep); library(spatialreg); library(lme4); library(lmerTest)
  library(matrixStats); library(geometry); library(cluster); library(ape)
  library(parallel)
})

log_msg("=== 06c 目级重复分析 / Clade(order)-level replication ===")

CLADE_RANK   <- "family"  # 重复单元 / replicate unit
# 阈值可由环境变量覆盖，用于阈值敏感性检验（06j）。默认值不变。
# Thresholds can be overridden by environment variables so that 06j can
# re-run the whole pipeline at other cut-offs; defaults are unchanged.
MIN_SP_CLADE <- as.integer(Sys.getenv("CLADE_MIN_SP", "20"))
MIN_CELLS    <- as.integer(Sys.getenv("CLADE_MIN_CELLS", "100"))
TAG <- Sys.getenv("CLADE_TAG", "")   # 非空时给输出表加后缀，避免覆盖主结果
GRAIN_LAB    <- "50km"
N_NULL_C     <- 499L   # 目级零模型次数（目数量多，适度降低）/ nulls per clade
N_CORES      <- max(1L, parallel::detectCores() - 2L)
CORE_ALL4    <- c("body_mass", "body_length", "nocturnality", "verticality",
                  "habitat_breadth", "range_size")
PRED <- c("ax_thermal", "ax_water", "ax_productivity", "ax_structure", "ax_human")

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))
mdat   <- readRDS(file.path(PATH$derived, "model_data.rds"))[[GRAIN_LAB]]
env    <- mdat |> distinct(cell_id, .keep_all = TRUE) |>
  select(cell_id, all_of(PRED), x_albers, y_albers)

# ---------------------------------------------------------------
# 1. 物种 -> 目 的对应表 / species-to-order lookup
# ---------------------------------------------------------------
build_clade_lookup <- function() {
  out <- list()
  # IUCN 名录带目/科信息 / IUCN species lists carry order and family
  for (cl in c("Amphibia", "Reptilia", "Mammalia")) {
    f <- file.path(PATH$derived, paste0("splist_", cl, ".rds"))
    if (file.exists(f)) {
      s <- readRDS(f)
      out[[cl]] <- data.frame(class = cl, species = s$species,
                              order  = toupper(trimws(as.character(s$order))),
                              family = toupper(trimws(as.character(s$family))))
    }
  }
  # 鸟类取自 TetrapodTraits 骨架 / birds from the trait backbone
  tt <- data.table::fread(file.path(PATH$raw_trait, "TetrapodTraits_v2.0.1.csv"),
                          select = c("Scientific.Name", "Class", "Order", "Family"),
                          showProgress = FALSE)
  tt$species <- norm_name(tt$Scientific.Name)
  out$Aves <- data.frame(class = "Aves",
                         species = tt$species[tt$Class == "Aves"],
                         order   = toupper(trimws(tt$Order[tt$Class == "Aves"])),
                         family  = toupper(trimws(tt$Family[tt$Class == "Aves"])))
  bind_rows(out) |>
    filter(!is.na(family), nzchar(family), !is.na(order), nzchar(order)) |>
    distinct(class, species, .keep_all = TRUE)
}
ORD <- cache_rds("species_clade_lookup", build_clade_lookup())
ORD$clade <- ORD[[CLADE_RANK]]
log_msg("物种-科 对应: ", nrow(ORD), " 条, ", length(unique(ORD$clade)), " 个科")

# ---------------------------------------------------------------
# 2. 单目的群落指标（与 05 同法，简化为核心指标）
#    Metrics for one order (same method as script 05, core metrics only)
# ---------------------------------------------------------------
clade_metrics <- function(sp_pool, comm_sub) {
  tr <- traits[match(sp_pool, traits$species), CORE_ALL4, drop = FALSE]
  TM <- as.matrix(sapply(tr, function(z) suppressWarnings(as.numeric(z))))
  rownames(TM) <- sp_pool
  keep <- colSums(is.finite(TM)) == nrow(TM) & matrixStats::colSds(TM) > 0
  TM <- TM[, keep, drop = FALSE]
  if (ncol(TM) < 3) return(NULL)
  TM <- apply(TM, 2, zscore)

  gd  <- cluster::daisy(as.data.frame(TM), metric = "gower")
  pc  <- ape::pcoa(as.dist(gd))
  PCO <- pc$vectors[, seq_len(min(3, ncol(pc$vectors))), drop = FALSE]
  DM  <- as.matrix(stats::dist(PCO))

  one <- function(idx) {
    S <- length(idx)
    Xs <- TM[idx, , drop = FALSE]
    mu <- colMeans(Xs); d <- Xs - rep(mu, each = S); ss <- colSums(d * d)
    P <- PCO[idx, , drop = FALSE]; ctr <- colMeans(P); Pd <- P - rep(ctr, each = S)
    c(richness = S,
      setNames(mu, paste0(colnames(TM), ".mean")),
      setNames(ss / (S - 1), paste0(colnames(TM), ".var")),
      FDis = mean(sqrt(rowSums(Pd * Pd))),
      RaoQ = sum(DM[idx, idx]) / (S * (S - 1)))
  }

  sp_idx <- setNames(seq_along(sp_pool), sp_pool)
  by_cell <- split(sp_idx[comm_sub$species], comm_sub$cell_id)
  by_cell <- by_cell[lengths(by_cell) >= 5]
  if (length(by_cell) < MIN_CELLS) return(NULL)

  obs <- t(vapply(by_cell, one, FUN.VALUE = one(1:5)))
  obs <- as.data.frame(obs); obs$cell_id <- names(by_cell)

  # 零模型（按唯一丰富度）/ null model per unique richness
  # 抽样按物种占据的网格数加权，与脚本 05 的主零模型一致。等概率抽样把狭域
  # 特有种与广布种视为等可能，会读出虚假的功能聚集（实测：等概率下平均
  # SES.FDis 鸟 -1.68、兽 -1.26、爬行 -1.88；加权后为 -0.84、-0.90、-1.13）。
  # Weight draws by how many cells each species occupies, matching the primary
  # null in script 05; equiprobable draws manufacture apparent clustering.
  pool <- length(sp_pool); tmpl <- one(1:5)
  occ_w <- as.numeric(table(factor(comm_sub$species, levels = sp_pool)))
  occ_w <- occ_w / sum(occ_w)
  S_vals <- sort(unique(lengths(by_cell)))
  S_vals <- S_vals[S_vals >= 5 & S_vals <= pool]
  nul <- parallel::mclapply(S_vals, function(S) {
    set.seed(SEED + S)
    sim <- vapply(seq_len(N_NULL_C),
                  function(i) one(sample.int(pool, S, prob = occ_w)),
                  FUN.VALUE = tmpl)
    list(mu = setNames(rowMeans(sim, na.rm = TRUE), names(tmpl)),
         sd = setNames(matrixStats::rowSds(sim, na.rm = TRUE), names(tmpl)))
  }, mc.cores = N_CORES)
  names(nul) <- as.character(S_vals)

  mets <- setdiff(names(obs), c("cell_id", "richness"))
  ses <- matrix(NA_real_, nrow(obs), length(mets), dimnames = list(NULL, mets))
  Sv <- as.character(obs$richness)
  for (i in seq_len(nrow(obs))) {
    nn <- nul[[Sv[i]]]
    if (is.null(nn) || !is.list(nn)) next
    ses[i, ] <- (as.numeric(obs[i, mets]) - nn$mu[mets]) / nn$sd[mets]
  }
  ses <- as.data.frame(ses); names(ses) <- paste0("SES.", mets)
  cbind(obs[, c("cell_id", "richness")], ses)
}

# ---------------------------------------------------------------
# 3. 逐目：构建群落 -> 指标 -> SAR 效应
# ---------------------------------------------------------------
fit_sar_simple <- function(dd, resp) {
  ok <- is.finite(dd[[resp]])
  if (sum(ok) < MIN_CELLS) return(NULL)
  x <- dd[ok, ]
  f <- stats::as.formula(paste0("`", resp, "` ~ ", paste(PRED, collapse = " + ")))
  nb <- spdep::make.sym.nb(spdep::knn2nb(
    spdep::knearneigh(as.matrix(x[, c("x_albers", "y_albers")]), k = 8)))
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  sar <- try(spatialreg::errorsarlm(f, data = x, listw = lw, zero.policy = TRUE,
                                    method = "Matrix"), silent = TRUE)
  if (inherits(sar, "try-error")) return(NULL)
  co <- as.data.frame(summary(sar)$Coef)
  co$term <- rownames(co); rownames(co) <- NULL
  names(co)[1:4] <- c("estimate", "se", "z", "p")
  co$n <- nrow(x); co$lambda <- sar$lambda
  co[co$term %in% PRED, ]
}

results <- list()
for (cl in TAXA$class) {
  f <- file.path(PATH$derived, sprintf("comm_%s_%s.rds", cl, GRAIN_LAB))
  if (!file.exists(f)) next
  comm <- readRDS(f)
  lk <- ORD |> filter(class == cl)
  comm <- comm |> left_join(lk[, c("species", "order", "clade")], by = "species")
  comm <- comm[!is.na(comm$clade), ]
  comm <- comm[comm$species %in% traits$species, ]

  tab <- comm |> group_by(clade, order) |>
    summarise(n_sp = n_distinct(species), n_cells = n_distinct(cell_id), .groups = "drop") |>
    filter(n_sp >= MIN_SP_CLADE, n_cells >= MIN_CELLS)
  log_msg("-- ", cl, ": ", nrow(tab), " 个科达标 (>=", MIN_SP_CLADE, " 种, >=",
          MIN_CELLS, " 格)")

  for (i in seq_len(nrow(tab))) {
    od <- tab$clade[i]
    sub <- comm[comm$clade == od, ]
    sp_pool <- sort(unique(sub$species))
    m <- try(clade_metrics(sp_pool, sub), silent = TRUE)
    if (inherits(m, "try-error") || is.null(m)) { log_msg("   [skip] ", od); next }
    dd <- m |> left_join(env, by = "cell_id")
    resp_set <- intersect(c("SES.FDis", "SES.RaoQ", "SES.body_mass.mean",
                            "SES.body_mass.var", "SES.habitat_breadth.mean",
                            "SES.habitat_breadth.var"), names(dd))
    for (r in resp_set) {
      co <- try(fit_sar_simple(dd, r), silent = TRUE)
      if (inherits(co, "try-error") || is.null(co)) next
      co$class <- cl; co$clade <- od; co$order <- tab$order[i]; co$response <- r
      co$n_species <- tab$n_sp[i]
      results[[length(results) + 1]] <- co
    }
    log_msg("   ", od, " (", tab$n_sp[i], " 种, ", nrow(m), " 格) 完成")
  }
}

eff <- bind_rows(results)
if (!nrow(eff)) { log_msg("[warn] 无目级结果"); quit(save = "no") }
eff <- eff |> left_join(TAXA[, c("class", "thermal")], by = "class")
eff$p_adj <- stats::p.adjust(eff$p, method = "BH")
write_table(eff, paste0("table_5_clade_level_effects", TAG))
log_msg("科级效应: ", nrow(eff), " 行, ", length(unique(eff$clade)), " 个科")

# ---------------------------------------------------------------
# 4. H1 正式检验：以「目」为重复，纲为随机效应
#    Formal H1 test with ORDERS as replicates and class as a random effect
# ---------------------------------------------------------------
contrasts_out <- list()
for (r in unique(eff$response)) {
  for (tm in PRED) {
    dd <- eff |> filter(response == r, term == tm)
    if (nrow(dd) < 6 || length(unique(dd$thermal)) < 2) next
    dd$thermal <- factor(dd$thermal, levels = c("Endotherm", "Ectotherm"))
    # 以效应绝对值衡量「耦合强度」，以效应本身衡量「方向」
    # |effect| measures coupling STRENGTH; effect measures DIRECTION
    for (what in c("strength", "direction")) {
      dd$y <- if (what == "strength") abs(dd$estimate) else dd$estimate
      # 纲/目 嵌套随机效应吸收科间的系统发育非独立性
      # class/order nesting absorbs phylogenetic non-independence among families
      m <- try(lme4::lmer(y ~ thermal + (1 | class/order), data = dd, REML = TRUE),
               silent = TRUE)
      if (inherits(m, "try-error")) {
        m2 <- stats::lm(y ~ thermal, data = dd)
        s <- summary(m2)$coefficients
        src <- "lm (nested RE failed)"
      } else {
        s <- summary(m)$coefficients; src <- "lmer (1|class/order)"
      }
      if (!"thermalEctotherm" %in% rownames(s)) next
      # p 值必须按列名取：lme4 的系数表没有 p 列，按位置取会误取 t 值
      # Extract p by NAME: lme4's coefficient table has no p column, so
      # positional extraction silently returns the t statistic instead
      pcol <- grep("^Pr\\(", colnames(s))
      est <- s["thermalEctotherm", 1]; sev <- s["thermalEctotherm", 2]
      if (length(pcol)) {
        pval <- s["thermalEctotherm", pcol[1]]
      } else {
        dfree <- max(1, nrow(dd) - length(unique(dd$class)) - 1)
        pval <- 2 * stats::pt(abs(est / sev), df = dfree, lower.tail = FALSE)
      }
      contrasts_out[[length(contrasts_out) + 1]] <- data.frame(
        response = r, predictor = tm, quantity = what,
        n_clades = nrow(dd),
        n_endo = sum(dd$thermal == "Endotherm"),
        n_ecto = sum(dd$thermal == "Ectotherm"),
        endo_mean = round(mean(dd$y[dd$thermal == "Endotherm"]), 4),
        ecto_mean = round(mean(dd$y[dd$thermal == "Ectotherm"]), 4),
        contrast  = round(est, 4),
        se        = round(sev, 4),
        t         = round(est / sev, 3),
        p         = pval,
        model     = src, row.names = NULL)
    }
  }
}
ct <- bind_rows(contrasts_out)
if (nrow(ct)) {
  ct$p_adj <- stats::p.adjust(ct$p, method = "BH")
  write_table(ct, paste0("table_6_thermal_mode_contrast", TAG))
  log_msg("体温调节模式对比: ", nrow(ct), " 项, 显著 (BH<0.05): ",
          sum(ct$p_adj < 0.05, na.rm = TRUE))
  print(as.data.frame(ct |> filter(p_adj < 0.10) |> arrange(p_adj)))
}

# ---------------------------------------------------------------
# 5. 总体检验（首要推断）/ Headline test
#    上面 30 个逐项对比彼此相关且各自检验力有限；正确的整体推断是
#    先把每个科的耦合强度汇总为一个值，再以「科」为唯一重复做单次检验。
#    The 30 per-combination contrasts are correlated and individually
#    underpowered. The valid headline inference collapses each family to a
#    single coupling-strength value and runs ONE test with family as the
#    sole replicate.
# ---------------------------------------------------------------
fam <- eff |>
  group_by(class, order, clade, thermal, n_species) |>
  summarise(coupling = mean(abs(estimate), na.rm = TRUE),
            coupling_thermal = mean(abs(estimate[term == "ax_thermal"]), na.rm = TRUE),
            coupling_human   = mean(abs(estimate[term == "ax_human"]),   na.rm = TRUE),
            .groups = "drop") |>
  mutate(thermal = factor(thermal, levels = c("Endotherm", "Ectotherm")))
write_table(fam, paste0("table_7_family_coupling_strength", TAG))

headline <- list()
for (v in c("coupling", "coupling_thermal", "coupling_human")) {
  dd <- fam[is.finite(fam[[v]]), ]; dd$y <- dd[[v]]
  m <- try(lmerTest::lmer(y ~ thermal + (1 | class/order), data = dd, REML = TRUE),
           silent = TRUE)
  if (inherits(m, "try-error")) {
    m <- stats::lm(y ~ thermal, data = dd); s <- summary(m)$coefficients
    src <- "lm"; pv <- s["thermalEctotherm", 4]
  } else {
    s <- summary(m)$coefficients; src <- "lmer (1|class/order)"
    pc <- grep("^Pr\\(", colnames(s))
    pv <- if (length(pc)) s["thermalEctotherm", pc[1]] else
      2 * stats::pnorm(abs(s["thermalEctotherm", 1] / s["thermalEctotherm", 2]),
                       lower.tail = FALSE)
  }
  # 非参数备份：Wilcoxon 秩和检验 / non-parametric back-up
  wt <- stats::wilcox.test(y ~ thermal, data = dd, exact = FALSE)
  headline[[length(headline) + 1]] <- data.frame(
    metric = v, n_families = nrow(dd),
    n_endo = sum(dd$thermal == "Endotherm"), n_ecto = sum(dd$thermal == "Ectotherm"),
    endo_mean = round(mean(dd$y[dd$thermal == "Endotherm"]), 4),
    ecto_mean = round(mean(dd$y[dd$thermal == "Ectotherm"]), 4),
    ratio_ecto_endo = round(mean(dd$y[dd$thermal == "Ectotherm"]) /
                            mean(dd$y[dd$thermal == "Endotherm"]), 2),
    contrast = round(s["thermalEctotherm", 1], 4),
    se = round(s["thermalEctotherm", 2], 4),
    p_mixed = signif(pv, 3), p_wilcoxon = signif(wt$p.value, 3),
    model = src, row.names = NULL)
}
headline <- bind_rows(headline)
write_table(headline, paste0("table_6b_headline_thermal_test", TAG))
log_msg("--- 总体检验 / headline test ---")
print(as.data.frame(headline))

# 符号一致性检验：30 个对比中变温 > 恒温 的比例
# Sign-consistency test across the 30 contrasts
st <- ct |> filter(quantity == "strength")
if (nrow(st)) {
  k <- sum(st$contrast > 0); n <- nrow(st)
  bt <- stats::binom.test(k, n, p = 0.5)
  log_msg("符号一致性 / sign consistency: ", k, "/", n,
          " 的对比中变温耦合更强, binomial p = ", signif(bt$p.value, 3))
  log_msg("  注：30 个对比彼此相关，此检验仅作方向一致性的描述性佐证；",
          "首要推断以 table_6b 的科级单次检验为准。")
  write_table(data.frame(n_contrasts = n, n_ectotherm_stronger = k,
                         binomial_p = signif(bt$p.value, 4),
                         note = "correlated tests; descriptive support only"),
              "table_6c_sign_consistency")
}
log_msg("=== 06c 完成 / done ===")
