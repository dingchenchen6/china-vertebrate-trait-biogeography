# ============================================================
# 科学问题 / Scientific question:
#   「群落性状结构」是一个东西还是几个独立维度？若体型、生活史节奏与
#   生态位宽度各自响应不同的环境梯度，把它们压进一个多元指数会掩盖机制。
#   Is assemblage trait structure one thing or several independent dimensions?
#
# 分析目标 / Objective:
#   构建**先验定义**的功能轴，而非让 PCA 自由旋转。
#
# 为什么不用自由 PCA / Why not free PCA:
#   纲内 PCA 的 PC1 在四个纲中载荷结构不同——鸟与两栖的 PC1 由生态位宽度
#   主导，爬行的 PC1 才是体型+繁殖力，哺乳的 PC1 是垂直生境。直接跨纲比较
#   PC1 会把不同的东西当成同一个轴来比。先验轴以固定的性状组合定义，
#   在四个纲中含义一致，因而可比。数据驱动的 PCA 仍并列输出作为对照。
#   Free PCA gives axes whose meaning differs between classes, so comparing
#   "PC1" across classes compares different things. A priori axes fix the trait
#   composition so the axis means the same thing everywhere; the data-driven
#   PCA is retained alongside as a check.
#
# 三个先验功能轴 / Three a priori axes:
#   1. 体型 SIZE          = 体型（每纲一个量度）
#   2. 生活史节奏 PACE    = 体型(+) 与 繁殖力(−) 的对比 —— 慢速端 = 大而少子
#   3. 生态位宽度 NICHE   = 生境宽度(+) 与 分布范围(+)
#   （活动节律与垂直生境保留为单性状，不合成轴——二者在生物学上不同质）
#
# 输入数据 / Input:
#   03_data_derived/traits_core.rds, traits_lifehistory.rds
#
# 预期输出 / Expected output:
#   03_data_derived/functional_axes_species.rds
#   04_results/tables/table_37_axis_definition.csv
#   05_figures/Fig2_functional_axes.pdf/.png
#
# 关键假设 / Key assumptions:
#   - 繁殖力四纲实测覆盖 60–100%，寿命仅 8–55%，故快慢轴只用「体型 vs 繁殖力」；
#     世代长度（鸟 80.8%、兽 92.2%）只进入恒温类群的扩展分析。
#   - 性状在**纲内** z 标准化后再合成轴，故轴分数在纲内可比、跨纲同义。
#
# 主要包 / Main packages: dplyr, ggplot2, ggrepel
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({ library(tidyr); library(dplyr); library(ggrepel) })

log_msg("=== 03f 先验功能轴 / A priori functional axes ===")

core <- readRDS(file.path(PATH$derived, "traits_core.rds"))
lh   <- readRDS(file.path(PATH$derived, "traits_lifehistory.rds"))

# 用整合后的繁殖力替换核心表中的窝仔数（覆盖更高、来源可追溯）
core <- core |>
  left_join(lh |> select(species, fecundity, any_of("imp_fecundity")),
            by = "species") |>
  mutate(fecundity = log10(pmax(fecundity, 0.5)))

CORE <- c("size", "fecundity", "nocturnality", "verticality",
          "habitat_breadth", "range_size")
LAB <- c(size = "Body size", fecundity = "Fecundity",
         nocturnality = "Nocturnality", verticality = "Vertical stratum use",
         habitat_breadth = "Habitat breadth", range_size = "Range size")

# ---------------------------------------------------------------
# 1. 先验轴的定义 / A priori axis definitions
# ---------------------------------------------------------------
# 每个轴 = 指定性状的 z 分数按给定符号取平均。符号即生物学方向。
# PACE 用「繁殖力对体型的残差」而非「体型 − 繁殖力」。
# 后者把 size 直接写进 PACE，导致两轴相关 0.38–0.73，无法分开
# 「体型效应」与「生活史节奏效应」。取残差后二者构造上正交：
# PACE 表示「相对于其体型，该物种少生还是多生」。
# Defining PACE as size minus fecundity embeds size in it, correlating the two
# axes at 0.38-0.73. The residual of fecundity on size is orthogonal to size by
# construction and means "fewer or more offspring than its body size predicts".
AXDEF <- list(
  SIZE  = c(size = 1),
  NICHE = c(habitat_breadth = 1, range_size = 1)
)
AXLAB_EXTRA <- c(PACE = "Slow-fast pace of life")
AXLAB <- c(SIZE = "Body size", PACE = "Slow-fast pace of life",
           NICHE = "Niche breadth")
AXNOTE <- c(SIZE = "body size",
            PACE = "negative residual of fecundity on size; high = slow",
            NICHE = "habitat breadth (+) and range size (+)")

def_tab <- rbind(
  do.call(rbind, lapply(names(AXDEF), function(a)
    data.frame(axis = a, axis_label = AXLAB[[a]],
               trait = names(AXDEF[[a]]), sign = unname(AXDEF[[a]]),
               interpretation = AXNOTE[[a]]))),
  data.frame(axis = "PACE", axis_label = AXLAB[["PACE"]],
             trait = c("size", "fecundity"), sign = c(0, -1),
             interpretation = AXNOTE[["PACE"]]))
write_table(def_tab, "table_37_axis_definition")
log_msg("先验轴的定义 / axis definitions:"); print(def_tab)

# ---------------------------------------------------------------
# 2. 计算轴分数（纲内 z 标准化）/ Axis scores, standardised within class
# ---------------------------------------------------------------
build <- function(cl) {
  d <- core |> filter(class == cl)
  d <- d[stats::complete.cases(d[, CORE]), ]
  if (nrow(d) < 50) return(NULL)
  Z <- as.data.frame(apply(as.matrix(d[, CORE]), 2, zscore))
  out <- data.frame(species = d$species, class = cl,
                    order = d$order, family = d$family, Z)
  for (a in names(AXDEF)) {
    w <- AXDEF[[a]]
    out[[a]] <- zscore(rowMeans(sweep(as.matrix(Z[, names(w), drop = FALSE]),
                                      2, w, `*`)))
  }
  # PACE = −resid(繁殖力 ~ 体型)：与 SIZE 构造上正交
  out$PACE <- zscore(-stats::residuals(stats::lm(Z$fecundity ~ Z$size)))
  # 数据驱动 PCA 作为对照 / data-driven PCA retained as a check
  p <- stats::prcomp(as.matrix(Z), center = TRUE, scale. = FALSE)
  sgn <- sign(p$rotation["size", 1]); if (sgn == 0) sgn <- 1
  out$PC1 <- p$x[, 1] * sgn; out$PC2 <- p$x[, 2]
  L <- p$rotation[, 1:2, drop = FALSE]; L[, 1] <- L[, 1] * sgn
  attr(out, "loadings") <- L
  attr(out, "ve") <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)[1:2]
  out
}
AX <- lapply(TAXA$class, build) |> setNames(TAXA$class)
AX <- AX[!vapply(AX, is.null, logical(1))]
AXES3 <- c("SIZE", "PACE", "NICHE")
FA <- bind_rows(lapply(AX, function(x) x[, c("species", "class", "order", "family",
                                             CORE, AXES3, "PC1", "PC2")]))
saveRDS(FA, file.path(PATH$derived, "functional_axes_species.rds"))
log_msg("有功能轴分数的物种 / species with axis scores: ", nrow(FA))

# ---------------------------------------------------------------
# 3. 轴之间是否独立 / Are the axes independent?
# ---------------------------------------------------------------
# 若三个轴彼此高度相关，则「多个维度」的说法不成立。
cors <- lapply(names(AX), function(cl) {
  d <- AX[[cl]]
  m <- stats::cor(d[, AXES3], method = "spearman")
  data.frame(class = cl,
             SIZE_PACE  = round(m["SIZE", "PACE"], 3),
             SIZE_NICHE = round(m["SIZE", "NICHE"], 3),
             PACE_NICHE = round(m["PACE", "NICHE"], 3),
             # 先验轴与自由 PC1 的一致性 / agreement with the free PC1
             PACE_vs_PC1 = round(stats::cor(d$PACE, d$PC1, method = "spearman"), 3))
}) |> bind_rows()
write_table(cors, "table_37b_axis_independence")
log_msg("轴间相关（Spearman）/ correlations among axes:")
print(as.data.frame(cors))

# ---------------------------------------------------------------
# 4. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures

# (a) 轴定义示意 / what each axis is made of
dd <- def_tab |> mutate(
  trait_f = factor(LAB[trait], levels = rev(unname(LAB))),
  axis_f  = factor(AXLAB[axis], levels = unname(AXLAB)))
pa <- ggplot(dd, aes(axis_f, trait_f, fill = factor(sign))) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(sign > 0, "+", "−")), size = 4,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = c(`1` = "#2E6E8E", `-1` = "#B5541F", `0` = "grey75"),
                    labels = c(`1` = "positive", `-1` = "negative",
                               `0` = "used as covariate"), name = NULL) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(x = NULL, y = NULL, tag = "a",
       subtitle = "Axes are defined a priori so they mean the same thing in every class") +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.8, colour = "grey30"))

# (b) 轴间独立性 / independence
ci <- cors |> select(class, SIZE_PACE, SIZE_NICHE, PACE_NICHE) |>
  pivot_longer(-class, names_to = "pair", values_to = "r")
ci$pair_f <- factor(ci$pair, levels = c("SIZE_PACE", "SIZE_NICHE", "PACE_NICHE"),
                    labels = c("Size –\nPace", "Size –\nNiche", "Pace –\nNiche"))
ci$class_f <- factor(TAXA_LAB[ci$class], levels = unname(TAXA_LAB))
pb <- ggplot(ci, aes(pair_f, r, fill = class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  scale_fill_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = NULL, y = "Spearman correlation", tag = "b",
       subtitle = "Weak correlations mean the axes carry separate information") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 4.8, colour = "grey30"))

# (c-f) 每纲的体型–繁殖力权衡（快慢轴的核心）
tp <- lapply(TAXA$class, function(cl) {
  if (is.null(AX[[cl]])) return(NULL)
  d <- AX[[cl]]
  r <- stats::cor(d$size, d$fecundity, method = "spearman")
  ggplot(d, aes(size, fecundity)) +
    geom_point(size = 0.25, alpha = 0.2, colour = PAL$taxa[[cl]]) +
    geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.5,
                colour = "grey15", se = FALSE) +
    annotate("text", x = Inf, y = Inf, hjust = 1.08, vjust = 1.6,
             label = sprintf("rho = %+.2f", r), size = 1.9, colour = "grey20") +
    labs(title = TAXA_LAB[cl], x = "Body size (z)", y = "Fecundity (z)") +
    theme_pub()
}) |> Filter(f = Negate(is.null))
for (i in seq_along(tp)) tp[[i]] <- tp[[i]] + labs(tag = letters[i + 2])

p <- (pa | pb) / wrap_plots(tp, nrow = 1) + plot_layout(heights = c(1, 1))
save_fig(p, "Fig2_functional_axes", W2, 128, FIG)

log_msg("=== 03f 完成 / done ===")
