# ============================================================
# Scientific question / 科学问题:
#   主模型若不做空间校正，残差的空间自相关有多强？
#   How strong is residual spatial autocorrelation in the main driver models
#   before the spatial correction?
#
# Objective / 分析目标:
#   对 06g 的全部主模型（四纲 x 九个响应，六个具名预测变量）拟合普通最小
#   二乘，计算其残差的 Moran's I 统计量并存表。06g 只保存了 SAR 残差的
#   Moran's I（table_38 的 resid_moran），06_models 只保存了 OLS 检验的
#   p 值，因此「校正前的 Moran's I」此前没有任何可复现的出处。
#   Fit OLS for every main model in 06g and save the Moran's I statistic of its
#   residuals. Previously only the SAR residual I (table_38) and the OLS p-value
#   (table_2) were saved, so the pre-correction statistic had no traceable source.
#
# Input data / 输入数据:
#   03_data_derived/assemblage_axes_50km.rds（06g 生成）; env_50km.rds
# Main workflow / 主要流程:
#   1. 与 06g 完全相同的样本、标准化与 k = 8 对称化邻接
#   2. OLS 拟合 -> spdep::moran.test(残差)
#   3. 与 table_38 的 SAR 残差 Moran's I 并表
# Expected output / 预期输出:
#   04_results/tables/table_38d_ols_residual_moran.csv
# Key assumptions / 关键假设:
#   样本与邻接定义与 06g 逐行一致，两组 Moran's I 因此可直接对比。
# Main packages / 主要包: spdep, dplyr
# Output directory / 输出路径: 04_results/tables/
# ============================================================
source("01_code/00_config.R")
suppressPackageStartupMessages({ library(dplyr) })

PRED   <- c("bio1_mean", "bio12_mean", "npp", "elev_range", "bio4_mean", "hfp2020")
TRAITS <- c("size", "fecundity", "nocturnality", "verticality",
            "habitat_breadth", "range_size")
AXES   <- c("SIZE", "PACE", "NICHE")

asm <- readRDS(file.path(PATH$derived, "assemblage_axes_50km.rds"))
env <- readRDS(file.path(PATH$derived, "env_50km.rds")) |> drop_inset("50km")

out <- list()
for (cl in TAXA$class) {
  # 样本与 06g 一致 / identical sample to 06g
  d <- asm |> filter(class == cl) |> inner_join(env, by = "cell_id")
  d <- d[stats::complete.cases(d[, c(PRED, "x_albers", "y_albers")]), ]
  if (nrow(d) < 100) next
  nb <- spdep::make.sym.nb(
    spdep::knn2nb(spdep::knearneigh(as.matrix(d[, c("x_albers", "y_albers")]), k = 8)))
  lw <- spdep::nb2listw(nb, style = "W")
  X  <- as.data.frame(apply(as.matrix(d[, PRED]), 2, zscore))
  f  <- stats::as.formula(paste("y ~", paste(PRED, collapse = " + ")))
  for (r in intersect(c(TRAITS, AXES), names(d))) {
    if (sum(is.finite(d[[r]])) < nrow(d)) next
    X$y <- zscore(d[[r]])
    m  <- stats::lm(f, data = X)
    mi <- spdep::moran.test(stats::residuals(m), lw, zero.policy = TRUE)
    out[[length(out) + 1]] <- data.frame(
      class = cl, response = r, n = nrow(d),
      ols_resid_moran = unname(mi$estimate[1]), ols_moran_p = mi$p.value)
  }
}
res <- bind_rows(out)

# 与 SAR 残差 Moran's I 并表 / join the SAR residual statistic from table_38
sar <- utils::read.csv(file.path(PATH$tables, "table_38_hypothesis_effects.csv")) |>
  distinct(class, response, sar_resid_moran = resid_moran)
res <- res |> left_join(sar, by = c("class", "response"))
write_table(res, "table_38d_ols_residual_moran")

log_msg(sprintf("OLS 残差 Moran's I: %.3f – %.3f（%d 个模型，全部 P < %.0e）",
                min(res$ols_resid_moran), max(res$ols_resid_moran), nrow(res),
                max(res$ols_moran_p)))
log_msg(sprintf("SAR 残差 Moran's I: %.3f – %.3f",
                min(res$sar_resid_moran, na.rm = TRUE),
                max(res$sar_resid_moran, na.rm = TRUE)))
