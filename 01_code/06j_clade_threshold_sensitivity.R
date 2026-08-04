# ============================================================
# 科学问题 / Scientific question:
#   科级检验用了 33 个科，但中国四纲共有 166 个科。被排除的 133 个科
#   是怎么排除的？这个阈值是否驱动了结论？
#   The family-level test uses 33 of 166 families. What excluded the rest,
#   and does the threshold drive the result?
#
# 筛选标准（本就写在 06c 中，此处显式化并检验）/ The criterion, made explicit:
#   一个科要进入分析，必须同时满足
#     (a) 在中国有 >= 20 个物种 —— 否则该科的群落性状均值在每个网格上
#         只由极少数物种决定，噪声压过信号；
#     (b) 至少占据 100 个 50 km 网格 —— 否则「性状 ~ 环境」这条回归的
#         样本量不足以估计斜率。
#   两条都是**估计可行性**要求，与性状值或体温调节模式无关。
#
# 本脚本做三件事 / Three things this script does:
#   1. 画出纳入/排除的依据（科的规模分布），让「为什么是 33 个」可见。
#   2. 在 06c **自身的定义**下改阈值重跑（≥10/50、≥20/100、≥30/150），
#      比较结论。这些重跑由 06c 通过环境变量完成，结果表带 _th 后缀。
#   3. 换一种耦合强度的操作化（具名热变量上的 R²，而非合成热轴上的
#      |标准化系数|），检验结论是否依赖于定义。
#
# 结论摘要（写在这里，因为它限定了正文能说什么）/ What we found:
#   方向在三个阈值下一致（比值 1.33 / 1.79 / 1.71，均 > 1），但
#   **放宽到 ≥10 种时不再显著**（p = 0.23）。放宽纳入的是小科，其耦合
#   估计噪声更大，比值因回归稀释而收缩——方向未反转，符合预期，
#   但正文不能声称该结果对纳入阈值不敏感。
#   换用 R² 操作化后比值降至约 1.05，说明科级证据比网格级 SAR 证据弱。
#   因此正文把网格级结果作为主要证据，科级结果作为**方向性佐证**。
#
# 预期输出 / Expected output:
#   04_results/tables/table_46_family_inclusion.csv
#   04_results/tables/table_46b_threshold_sensitivity.csv
#   04_results/tables/table_46c_alternative_definition.csv
#   05_figures/Fig6c_family_threshold.pdf/.png
#
# 主要包 / Main packages: dplyr, ggplot2
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({ library(tidyr); library(dplyr) })

log_msg("=== 06j 科级阈值敏感性 / Family-threshold sensitivity ===")

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))
env <- readRDS(file.path(PATH$derived, "env_50km.rds")) |> drop_inset("50km")
FA  <- readRDS(file.path(PATH$derived, "functional_axes_species.rds"))

comm <- bind_rows(lapply(TAXA$class, function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  drop_inset(readRDS(f), "50km") |> mutate(class = cl)
}))
fam <- traits |> select(species, class, order, family) |>
  filter(!is.na(family), nzchar(family)) |> distinct(species, .keep_all = TRUE)
cf <- comm |> inner_join(fam, by = c("species", "class")) |>
  group_by(class, order, family) |>
  summarise(n_species = n_distinct(species), n_cells = n_distinct(cell_id),
            .groups = "drop") |>
  left_join(TAXA[, c("class", "thermal")], by = "class")
cf$kept20 <- cf$n_species >= 20 & cf$n_cells >= 100
write_table(cf, "table_46_family_inclusion")
log_msg("科总数 / families total: ", nrow(cf), "，主分析纳入 / retained: ",
        sum(cf$kept20))

drop_share <- cf |> group_by(thermal) |>
  summarise(families_total = n(), families_kept = sum(kept20),
            species_total = sum(n_species), species_kept = sum(n_species[kept20]),
            pct_species_kept = round(100 * sum(n_species[kept20]) /
                                       sum(n_species), 1), .groups = "drop")
write_table(drop_share, "table_46d_exclusion_share")
log_msg("排除的是小科，物种覆盖仍高 / exclusion drops small families:")
print(as.data.frame(drop_share))

# ---------------------------------------------------------------
# 2. 三个阈值下的科级结论 / the family-level result at three thresholds
# ---------------------------------------------------------------
# 这些表由 06c 在不同 CLADE_MIN_SP / CLADE_MIN_CELLS 下生成。
read_thr <- function(tag, lab) {
  f <- file.path(PATH$tables, paste0("table_6b_headline_thermal_test", tag, ".csv"))
  if (!file.exists(f)) { log_msg("  [缺] ", basename(f)); return(NULL) }
  d <- utils::read.csv(f); d$threshold <- lab; d
}
sens <- bind_rows(read_thr("_th10", "≥10 spp, ≥50 cells"),
                  read_thr("",      "≥20 spp, ≥100 cells (main)"),
                  read_thr("_th30", "≥30 spp, ≥150 cells"))
stopifnot(nrow(sens) > 0)
sens$threshold <- factor(sens$threshold,
                         levels = c("≥10 spp, ≥50 cells",
                                    "≥20 spp, ≥100 cells (main)",
                                    "≥30 spp, ≥150 cells"))
write_table(sens, "table_46b_threshold_sensitivity")
log_msg("阈值敏感性（06c 自身定义）/ threshold sensitivity, 06c's own definition:")
print(sens |> filter(metric == "coupling_thermal") |>
        select(threshold, n_families, n_endo, n_ecto, endo_mean, ecto_mean,
               ratio_ecto_endo, p_mixed, p_wilcoxon) |> as.data.frame())

# ---------------------------------------------------------------
# 3. 换一种耦合强度的定义 / an alternative operationalisation
# ---------------------------------------------------------------
# 06c 的耦合强度 = 合成热轴上 |标准化 SAR 系数| 的均值。
# 这里改用：三个功能轴对两个**具名**热变量（年均温、温度季节性）的
# 多元 R² 均值。若结论只在一种定义下成立，正文必须说明。
PRED_T <- c("bio1_mean", "bio4_mean")
AX <- c("SIZE", "PACE", "NICHE")
coupling_r2 <- function(fm) {
  sp <- fam$species[fam$family == fm]
  d <- comm[comm$species %in% sp, ] |>
    inner_join(FA[FA$species %in% sp, ], by = "species") |>
    group_by(cell_id) |>
    summarise(rich = n(), across(all_of(AX), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") |>
    filter(rich >= 3) |> inner_join(env, by = "cell_id")
  d <- d[stats::complete.cases(d[, c(AX, PRED_T)]), ]
  if (nrow(d) < 30) return(NULL)
  v <- vapply(AX, function(a)
    summary(stats::lm(stats::reformulate(PRED_T, a), data = d))$r.squared,
    numeric(1))
  data.frame(family = fm, n_cells_used = nrow(d), coupling_r2 = mean(v))
}
cs <- cache_rds("family_coupling_r2", bind_rows(lapply(unique(cf$family),
  function(fm) { o <- try(coupling_r2(fm), silent = TRUE)
                 if (inherits(o, "try-error")) NULL else o })))
cs <- cs |> inner_join(cf, by = "family") |> filter(kept20)
alt <- data.frame(
  definition = "Mean R2 of the three axes on named thermal variables",
  n_families = nrow(cs),
  endo_mean = round(mean(cs$coupling_r2[cs$thermal == "Endotherm"]), 4),
  ecto_mean = round(mean(cs$coupling_r2[cs$thermal == "Ectotherm"]), 4),
  ratio = round(mean(cs$coupling_r2[cs$thermal == "Ectotherm"]) /
                  mean(cs$coupling_r2[cs$thermal == "Endotherm"]), 2),
  p_wilcoxon = signif(stats::wilcox.test(coupling_r2 ~ thermal,
                                         data = cs, exact = FALSE)$p.value, 3))
write_table(alt, "table_46c_alternative_definition")
log_msg("换一种耦合定义后的结果 / under an alternative definition:")
print(alt)

# ---------------------------------------------------------------
# 4. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures
QL <- c(coupling = "All predictors", coupling_thermal = "Thermal predictors",
        coupling_human = "Human pressure")
sens$q_f <- factor(QL[sens$metric], levels = unname(QL))

pa <- ggplot(sens, aes(threshold, ratio_ecto_endo, group = q_f, colour = q_f)) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.3, colour = "grey45") +
  geom_line(linewidth = 0.45) +
  geom_point(aes(shape = p_wilcoxon < 0.05), size = 1.6) +
  geom_text(aes(label = sprintf("%.2f", ratio_ecto_endo)), vjust = -1.2,
            size = 1.8, show.legend = FALSE) +
  scale_colour_manual(values = c(`All predictors` = "#4A7FA5",
                                 `Thermal predictors` = "#B5541F",
                                 `Human pressure` = "#5B9E6E"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "P < 0.05", `FALSE` = "n.s."),
                     name = NULL) +
  scale_x_discrete(labels = function(x) gsub(", ", ",\n", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.16))) +
  labs(x = NULL, y = "Coupling ratio, ectotherms / endotherms", tag = "a",
       subtitle = "Direction holds at every threshold, but significance is lost when small families are admitted") +
  theme_pub() + theme(legend.position = "bottom", legend.box = "vertical",
                      legend.spacing.y = unit(1, "pt"),
                      axis.text.x = element_text(size = 5, lineheight = 0.9),
                      plot.subtitle = element_text(size = 4.4, colour = "grey30"))

pb <- ggplot(cf, aes(n_species, n_cells, colour = kept20, shape = thermal)) +
  geom_vline(xintercept = 20, linetype = "22", linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = 100, linetype = "22", linewidth = 0.3, colour = "grey55") +
  geom_point(size = 1, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = "#B5541F", `FALSE` = "grey72"),
                      labels = c(`TRUE` = sprintf("retained (%d)", sum(cf$kept20)),
                                 `FALSE` = sprintf("excluded (%d)", sum(!cf$kept20))),
                      name = NULL) +
  scale_shape_manual(values = c(Endotherm = 16, Ectotherm = 17), name = NULL) +
  labs(x = "Species in China (log scale)", y = "Occupied 50-km cells (log scale)",
       tag = "b",
       subtitle = "Exclusion is by estimability alone, and removes small families of both thermal modes") +
  theme_pub() + theme(legend.position = "bottom", legend.box = "vertical",
                      legend.spacing.y = unit(1, "pt"),
                      plot.subtitle = element_text(size = 4.2, colour = "grey30"))

p <- pa | pb
save_fig(p, "Fig6c_family_threshold", W2, 100, FIG)

log_msg("=== 06j 完成 / done ===")
