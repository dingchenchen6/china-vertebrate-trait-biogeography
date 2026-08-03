# ============================================================
# 科学问题 / Scientific question:
#   H1 热约束不对称：变温类群的群落性状与"热能量轴"耦合更强，
#        恒温类群与"生产力/生境结构轴"耦合更强？
#   H2 人为过滤器替代：随人类足迹升高，气候对性状组成的解释力下降、
#        人为因子解释力上升，且恒温类群更明显？
#   H2.5 均值–方差解耦：CWM 与 CWV 受不同环境轴支配？
#
# 分析目标 / Objective:
#   (a) 逐类群空间自回归模型，得到标准化效应量；
#   (b) 跨类群联合混合模型，正式检验 环境 × 体温调节模式 交互；
#   (c) 方差分解：气候 | 生产力 | 生境结构 | 人为 | 空间；
#   (d) 气候耦合强度随人类足迹的衰减（H2 的交互检验）；
#   (e) 随机森林非线性重要性作为稳健性。
#
# 输入数据 / Input:
#   03_data_derived/metrics_<class>_<grain>.rds, env_<grain>.rds
#
# 主要流程 / Workflow:
#   1. 预测变量分组、标准化、共线性筛选 (VIF<5)
#   2. Moran 特征向量 (MEM) 提取空间滤波器
#   3. SAR 误差模型（逐类群）+ 残差 Moran's I 检验
#   4. 联合混合模型（class 与 cell 随机效应 + MEM）检验交互
#   5. vegan::varpart 方差分解
#   6. 空间分块交叉验证
#
# 预期输出 / Expected output:
#   04_results/tables/table_2_*.csv, table_3_*.csv ...
#
# 关键假设 / Key assumptions:
#   - 响应变量使用零模型 SES，已剔除丰富度的机械效应
#   - 预测变量组内先做 PCA 降维以缓解共线性，保留可解释的第一轴
#
# 主要包 / Main packages: spdep, spatialreg, adespatial, lme4, vegan, ranger
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({
  library(spdep); library(spatialreg); library(vegan)
  library(lme4); library(lmerTest); library(ranger); library(car)
})

log_msg("=== 06 统计建模 / Statistical modelling ===")

# ---------------------------------------------------------------
# 1. 预测变量分组 / Predictor groups (五个生态学假设轴)
# ---------------------------------------------------------------
# 前缀 ax_ 用于区分「环境轴」与类群属性 thermal（体温调节模式）
# The ax_ prefix disambiguates environmental axes from the taxon attribute
# `thermal` (thermoregulatory mode)
PRED_GROUPS <- list(
  ax_thermal      = c("bio1_mean", "bio6_mean", "bio4_mean", "bio3_mean"), # 热能量与热变异
  ax_water        = c("bio12_mean", "bio15_mean", "bio17_mean"),           # 水分供给与季节性
  ax_productivity = c("npp", "lc_trees"),                                  # 生产力
  ax_structure    = c("elev_range", "het_shannon", "bio1_sd"),             # 生境结构与异质性
  ax_human        = c("hfp2020", "lc_cropland", "lc_built")                # 当代人为压力
)
ALL_PRED <- unlist(PRED_GROUPS, use.names = FALSE)

#' 组内 PCA 第一轴：缓解组内共线性并保留生态学含义
#' First PCA axis within each group: reduces collinearity, keeps interpretation
group_axes <- function(env) {
  out <- data.frame(cell_id = env$cell_id)
  load_tab <- list()
  for (g in names(PRED_GROUPS)) {
    v <- intersect(PRED_GROUPS[[g]], names(env))
    if (!length(v)) next
    X <- as.matrix(env[, v, drop = FALSE])
    X <- apply(X, 2, function(z) { z[!is.finite(z)] <- stats::median(z, na.rm = TRUE); zscore(z) })
    p <- stats::prcomp(X, center = TRUE, scale. = FALSE)
    ax <- p$x[, 1]
    # 统一方向：使第一个变量的载荷为正，便于跨类群解释
    if (p$rotation[1, 1] < 0) ax <- -ax
    out[[g]] <- as.numeric(zscore(ax))
    load_tab[[g]] <- data.frame(group = g, variable = v,
                                loading = round(p$rotation[, 1] * ifelse(p$rotation[1,1] < 0, -1, 1), 3),
                                var_explained = round(100 * p$sdev[1]^2 / sum(p$sdev^2), 1))
  }
  attr(out, "loadings") <- bind_rows(load_tab)
  out
}

# ---------------------------------------------------------------
# 2. 数据装配 / Assemble data
# ---------------------------------------------------------------
assemble <- function(lab) {
  env <- readRDS(file.path(PATH$derived, paste0("env_", lab, ".rds")))
  ax  <- group_axes(env)
  write_table(attr(ax, "loadings"), paste0("table_S7_predictor_axis_loadings_", lab))

  fs <- list.files(PATH$derived, pattern = paste0("^metrics_.*_", lab, "\\.rds$"),
                   full.names = TRUE)
  if (!length(fs)) return(NULL)
  m <- lapply(fs, readRDS) |> bind_rows()
  d <- m |>
    left_join(ax, by = "cell_id") |>
    left_join(env[, c("cell_id", "x_albers", "y_albers", "lon", "lat",
                      "hfp2020", "land_area_km2")], by = "cell_id") |>
    left_join(TAXA[, c("class", "thermal")], by = "class")
  d$thermal <- factor(d$thermal, levels = c("Endotherm", "Ectotherm"))
  d
}

# ---------------------------------------------------------------
# 3. 空间滤波器（Moran 特征向量）/ Spatial filters (MEM)
# ---------------------------------------------------------------
make_mem <- function(coords, k = 8, n_keep = 15) {
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  me <- try(adespatial::mem(lw, MEM.autocor = "positive"), silent = TRUE)
  if (inherits(me, "try-error")) return(list(lw = lw, mem = NULL))
  list(lw = lw, mem = as.data.frame(me)[, seq_len(min(n_keep, ncol(me))), drop = FALSE])
}

# ---------------------------------------------------------------
# 4. 逐类群 SAR 模型 / Per-class SAR error models
# ---------------------------------------------------------------
fit_sar <- function(d, resp, preds, lw) {
  ok <- is.finite(d[[resp]]) & stats::complete.cases(d[, preds, drop = FALSE])
  if (sum(ok) < 100) return(NULL)
  dd <- d[ok, ]
  f <- stats::as.formula(paste0("`", resp, "` ~ ", paste(preds, collapse = " + ")))
  ols <- stats::lm(f, data = dd)

  # 仅在有效行上重建空间权重 / rebuild weights on the retained rows
  # 注意：kNN 邻接本身不对称，而 SAR 的 Matrix 求解要求对称权重，
  # 故用 make.sym.nb 对称化（标准做法）。
  # NOTE: kNN neighbours are asymmetric but the Matrix solver for SAR needs
  # symmetric weights, so symmetrise with make.sym.nb (standard practice).
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(dd[, c("x_albers", "y_albers")]), k = 8))
  nb <- spdep::make.sym.nb(nb)
  lw2 <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  mi <- try(spdep::lm.morantest(ols, lw2, zero.policy = TRUE), silent = TRUE)

  sar <- try(spatialreg::errorsarlm(f, data = dd, listw = lw2, zero.policy = TRUE,
                                    method = "Matrix"), silent = TRUE)
  if (inherits(sar, "try-error"))
    sar <- try(spatialreg::errorsarlm(f, data = dd, listw = lw2, zero.policy = TRUE,
                                      method = "LU"), silent = TRUE)
  if (inherits(sar, "try-error")) return(NULL)
  s <- summary(sar)
  co <- as.data.frame(s$Coef)
  co$term <- rownames(co); rownames(co) <- NULL
  names(co)[1:4] <- c("estimate", "se", "z", "p")
  # 标量安全赋值：summary.sarlm 的某些字段可能缺失，
  # 直接赋 numeric(0) 会让整次拟合静默失败
  # Scalar-safe assignment: some summary.sarlm fields can be absent, and
  # assigning numeric(0) would silently abort the whole fit
  sc1 <- function(x) if (length(x) == 1 && is.finite(suppressWarnings(as.numeric(x))))
    as.numeric(x) else NA_real_
  co$n <- nrow(dd)
  co$lambda      <- sc1(sar$lambda)
  co$moran_ols_p <- if (inherits(mi, "try-error")) NA_real_ else sc1(mi$p.value)
  co$pseudo_r2   <- sc1(s$NK)                       # Nagelkerke 伪 R2
  co$aic_sar     <- sc1(stats::AIC(sar))
  co$aic_ols     <- sc1(stats::AIC(ols))
  co[co$term != "(Intercept)", ]
}

# ---------------------------------------------------------------
# 5. 主循环 / Main loop
# ---------------------------------------------------------------
RESP_SETS <- function(d) {
  ns <- names(d)
  list(
    # H2.5 核心：单性状矩的 SES / per-trait moment SES
    moments = grep("^SES\\.(body_mass|body_length|nocturnality|verticality|diet_breadth|diet_vert|diet_plant|habitat_breadth|litter_size|max_longevity|range_size)\\.(mean|var|skew|kurt)$", ns, value = TRUE),
    # 多元功能多样性 / multivariate FD
    fd = intersect(c("SES.FDis", "SES.FRic", "SES.RaoQ", "SES.trait_vol", "SES.trait_dens"), ns)
  )
}

coef_all <- list(); dat_all <- list()
for (lab in names(GRAIN)) {
  d <- assemble(lab)
  if (is.null(d)) { log_msg("  [skip] ", lab); next }
  dat_all[[lab]] <- d
  log_msg("-- ", lab, ": ", nrow(d), " 行, ", length(unique(d$class)), " 类群 --")

  preds <- intersect(names(PRED_GROUPS), names(d))
  rs <- RESP_SETS(d)
  resps <- c(rs$moments, rs$fd)
  log_msg("   响应变量 / responses: ", length(resps))

  for (cl in unique(d$class)) {
    dc <- d[d$class == cl, ]
    n_ok <- 0L; n_fail <- 0L; first_err <- NULL
    for (r in resps) {
      co <- try(fit_sar(dc, r, preds, NULL), silent = TRUE)
      if (inherits(co, "try-error")) {
        n_fail <- n_fail + 1L
        if (is.null(first_err)) first_err <- as.character(co)
        next
      }
      if (is.null(co)) { n_fail <- n_fail + 1L; next }
      co$class <- cl; co$grain <- lab; co$response <- r
      coef_all[[length(coef_all) + 1]] <- co
      n_ok <- n_ok + 1L
    }
    log_msg("   ", cl, ": 成功 ", n_ok, " / 失败 ", n_fail)
    if (n_ok == 0L && !is.null(first_err))
      log_msg("     [首个错误] ", substr(gsub("\n", " ", first_err), 1, 160))
  }
}

coefs <- bind_rows(coef_all)
if (nrow(coefs)) {
  coefs <- coefs |>
    left_join(TAXA[, c("class", "thermal")], by = "class") |>
    mutate(moment = sub("^SES\\.[^.]+\\.", "", response),
           trait  = sub("^SES\\.([^.]+)\\..*$", "\\1", response),
           moment = ifelse(grepl("^SES\\.(FDis|FRic|RaoQ|trait_vol|trait_dens)$", response),
                           "multivariate", moment),
           p_adj  = stats::p.adjust(p, method = "BH"))
  write_table(coefs, "table_2_sar_coefficients_all")
  log_msg("SAR 系数表 / SAR coefficient table: ", nrow(coefs), " 行")
}

# ---------------------------------------------------------------
# 6. H1 检验：环境 × 体温调节模式 交互（联合混合模型）
#    H1 test: environment x thermoregulatory mode interaction
# ---------------------------------------------------------------
test_interaction <- function(d, resp_prefix, lab) {
  ns <- grep(paste0("^SES\\.", resp_prefix), names(d), value = TRUE)
  out <- list()
  preds <- intersect(names(PRED_GROUPS), names(d))
  for (r in ns) {
    dd <- d[is.finite(d[[r]]), ]
    if (nrow(dd) < 200 || length(unique(dd$thermal)) < 2) next
    dd$y <- as.numeric(scale(dd[[r]]))
    f <- stats::as.formula(paste0("y ~ (", paste(preds, collapse = " + "),
                                  ") * thermal + (1|class) + (1|cell_id)"))
    m <- try(lme4::lmer(f, data = dd, REML = FALSE,
                        control = lme4::lmerControl(calc.derivs = FALSE)), silent = TRUE)
    if (inherits(m, "try-error")) next
    s <- summary(m)$coefficients
    ix <- grep(":thermalEctotherm$", rownames(s))
    if (!length(ix)) next
    out[[length(out) + 1]] <- data.frame(
      grain = lab, response = r,
      term = sub(":thermalEctotherm$", "", rownames(s)[ix]),
      interaction_estimate = round(s[ix, "Estimate"], 4),
      se = round(s[ix, "Std. Error"], 4),
      t = round(s[ix, "t value"], 3),
      p = if ("Pr(>|t|)" %in% colnames(s)) s[ix, "Pr(>|t|)"] else NA_real_,
      row.names = NULL)
  }
  bind_rows(out)
}

inter <- list()
for (lab in names(dat_all)) {
  d <- dat_all[[lab]]
  inter[[lab]] <- bind_rows(
    test_interaction(d, "(body_mass|litter_size|max_longevity|habitat_breadth|nocturnality)", lab),
    test_interaction(d, "(FDis|FRic|RaoQ|trait_vol|trait_dens)$", lab))
}
inter <- bind_rows(inter)
if (nrow(inter)) {
  inter$p_adj <- stats::p.adjust(inter$p, method = "BH")
  write_table(inter, "table_3_thermal_mode_interactions")
  log_msg("交互检验 / interaction tests: ", nrow(inter), " 行, 显著 (BH<0.05): ",
          sum(inter$p_adj < 0.05, na.rm = TRUE))
}

saveRDS(dat_all, file.path(PATH$derived, "model_data.rds"))
log_msg("=== 06 完成 / done ===")
