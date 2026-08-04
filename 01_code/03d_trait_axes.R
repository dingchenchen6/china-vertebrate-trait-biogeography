# ============================================================
# 科学问题 / Scientific question:
#   「群落性状结构」是一个东西，还是几个彼此独立的维度？
#   若体型、生活史节奏、生态位宽度各自响应不同的环境梯度，
#   那么把它们混在一个多元功能多样性指数里会掩盖机制。
#   Is assemblage trait structure one thing or several? If body size, pace of
#   life and niche breadth respond to different gradients, collapsing them into
#   a single multivariate diversity index hides the mechanism.
#
# 分析目标 / Objective:
#   1. 重建核心性状集：每个纲**只保留一个体型量度**，避免体重与体长共线；
#      两栖类用体长（其体重实测率仅 25%），其余三纲用体重。
#   2. 加入生活史性状：窝仔数（四纲）、世代长度（仅恒温类群可得）。
#   3. 对每个纲做性状 PCA，提取可解释的功能轴，并绘图。
#
# 输入数据 / Input:
#   03_data_derived/traits_imputed.rds, traits_observed_china.rds
#
# 预期输出 / Expected output:
#   03_data_derived/traits_core.rds        重建后的核心性状表
#   03_data_derived/trait_axes_species.rds 每个物种在功能轴上的坐标
#   04_results/tables/table_35_trait_axis_loadings.csv
#   05_figures/Fig2_trait_axes.pdf/.png
#
# 关键假设 / Key assumptions:
#   - **每纲只用一个体型量度**。体重与体长在四个纲中的 log-log 相关为
#     0.86-0.97，同时纳入会让「体型」在 PCA 中获得双倍权重，并使功能轴的
#     解释失真。两栖类选体长而非体重，因为其体重 75% 为模型插补。
#   - 世代长度仅鸟（80.8%）与兽（92.2%）可得，爬行与两栖为 0，
#     故它只进入恒温类群的扩展分析，不进入四纲通用的核心集。
#   - PCA 在**纲内**进行：不同纲的性状量纲与生物学含义不可直接合并。
#
# 主要包 / Main packages: dplyr, ggplot2, ggrepel
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({ library(tidyr); library(dplyr); library(ggrepel) })

log_msg("=== 03d 核心性状集与功能轴 / Core traits and functional axes ===")

traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))
f_cn <- file.path(PATH$derived, "traits_observed_china.rds")
cn <- if (file.exists(f_cn)) readRDS(f_cn) else NULL

# ---------------------------------------------------------------
# 1. 每纲选一个体型量度 / One size measure per class
# ---------------------------------------------------------------
# 先量化体重与体长的冗余程度，作为「只选一个」的依据
# Quantify the redundancy that motivates keeping only one
red <- lapply(unique(traits$class), function(cl) {
  d <- traits |> filter(class == cl)
  a <- suppressWarnings(as.numeric(d$body_mass))
  b <- suppressWarnings(as.numeric(d$body_length))
  ok <- is.finite(a) & is.finite(b)
  data.frame(class = cl, n = sum(ok),
             r_mass_length = round(stats::cor(a[ok], b[ok]), 3))
}) |> bind_rows()
log_msg("体重与体长的相关（均为 log10 尺度）/ mass-length correlation:")
print(as.data.frame(red))

# 两栖类的体重 75% 为上游插补，体长仅 1%；故两栖用体长，其余用体重。
# Amphibian body mass is 75% imputed upstream against 1% for body length.
SIZE_TRAIT <- c(Aves = "body_mass", Mammalia = "body_mass",
                Reptilia = "body_mass", Amphibia = "body_length")
log_msg("体型量度的选择 / size measure per class: ",
        paste(sprintf("%s=%s", names(SIZE_TRAIT), SIZE_TRAIT), collapse = ", "))

# ---------------------------------------------------------------
# 2. 组装核心性状表 / Assemble the core trait table
# ---------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))

core <- lapply(unique(traits$class), function(cl) {
  d <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  sz <- num(d[[SIZE_TRAIT[[cl]]]])

  # 窝仔数：两栖类的全球库覆盖不足，用中国本土数据补
  # Clutch size: top up amphibians from the China compilation
  lit <- num(d$litter_size)
  if (!is.null(cn)) {
    j <- cn[match(d$species, cn$species), ]
    if ("cn_litter_size" %in% names(j)) {
      add <- log10(pmax(num(j$cn_litter_size), 0.5))
      lit <- ifelse(is.finite(lit), lit, add)
    }
  }

  data.frame(species = d$species, class = cl, order = d$order, family = d$family,
             size = sz,
             litter_size = lit,
             generation_length = num(d$et_Generation_length_d),
             nocturnality = num(d$nocturnality),
             verticality = num(d$verticality),
             habitat_breadth = num(d$habitat_breadth),
             range_size = num(d$range_size),
             stringsAsFactors = FALSE)
}) |> bind_rows()

cov <- core |> group_by(class) |>
  summarise(across(c(size, litter_size, generation_length, nocturnality,
                     verticality, habitat_breadth, range_size),
                   ~ round(100 * mean(is.finite(.x)), 1)), .groups = "drop")
write_table(cov, "table_35a_core_trait_coverage")
log_msg("核心性状集的完整率 / completeness of the core set:")
print(as.data.frame(cov))

saveRDS(core, file.path(PATH$derived, "traits_core.rds"))

# ---------------------------------------------------------------
# 3. 功能轴：纲内 PCA / Functional axes from within-class PCA
# ---------------------------------------------------------------
# 四纲通用的核心集（世代长度因爬行与两栖为 0 而排除）
CORE6 <- c("size", "litter_size", "nocturnality", "verticality",
           "habitat_breadth", "range_size")
LAB <- c(size = "Body size", litter_size = "Clutch / litter size",
         generation_length = "Generation length", nocturnality = "Nocturnality",
         verticality = "Vertical stratum use", habitat_breadth = "Habitat breadth",
         range_size = "Range size")

axes <- list(); loads <- list(); vexp <- list()
for (cl in TAXA$class) {
  d <- core |> filter(class == cl)
  d <- d[stats::complete.cases(d[, CORE6]), ]
  if (nrow(d) < 50) next
  M <- apply(as.matrix(d[, CORE6]), 2, zscore)
  p <- stats::prcomp(M, center = TRUE, scale. = FALSE)
  ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)

  # 符号对齐：让 PC1 的体型载荷为正，使各纲的轴方向可比
  # Sign-align so body size loads positively on PC1 in every class
  sgn <- sign(p$rotation["size", 1]); if (sgn == 0) sgn <- 1
  p$rotation[, 1] <- p$rotation[, 1] * sgn
  p$x[, 1] <- p$x[, 1] * sgn

  axes[[cl]] <- data.frame(species = d$species, class = cl,
                           PC1 = p$x[, 1], PC2 = p$x[, 2], PC3 = p$x[, 3])
  loads[[cl]] <- as.data.frame(p$rotation[, 1:3]) |>
    tibble::rownames_to_column("trait") |> mutate(class = cl)
  vexp[[cl]] <- data.frame(class = cl, axis = paste0("PC", 1:3),
                           var_explained = ve[1:3], n_species = nrow(d))
}
AX <- bind_rows(axes); LD <- bind_rows(loads); VE <- bind_rows(vexp)
saveRDS(AX, file.path(PATH$derived, "trait_axes_species.rds"))
write_table(LD |> mutate(across(where(is.numeric), ~ round(.x, 3))),
            "table_35_trait_axis_loadings")
write_table(VE, "table_35b_trait_axis_variance")
log_msg("功能轴的方差解释 / variance explained:")
print(VE |> pivot_wider(names_from = axis, values_from = var_explained) |>
        as.data.frame())

# 轴的生物学解读：由载荷自动生成，避免人为附会
# Interpret each axis from its loadings rather than by assertion
interp <- LD |> pivot_longer(PC1:PC3, names_to = "axis", values_to = "load") |>
  group_by(class, axis) |>
  slice_max(abs(load), n = 2, with_ties = FALSE) |>
  summarise(top_traits = paste(sprintf("%s (%+.2f)", LAB[trait], load),
                               collapse = "; "), .groups = "drop")
write_table(interp, "table_35c_axis_interpretation")
log_msg("各轴载荷最高的两个性状 / two strongest loadings per axis:")
print(as.data.frame(interp))

# ---------------------------------------------------------------
# 4. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures

# (a-d) 每纲的 PC1-PC2 双标图 / per-class biplot
biplots <- lapply(TAXA$class, function(cl) {
  a <- AX |> filter(class == cl); l <- LD |> filter(class == cl)
  if (!nrow(a)) return(NULL)
  v <- VE |> filter(class == cl)
  sc <- 0.85 * max(abs(c(a$PC1, a$PC2)))
  ggplot(a, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_point(size = 0.25, alpha = 0.22, colour = PAL$taxa[[cl]]) +
    geom_segment(data = l, inherit.aes = FALSE,
                 aes(x = 0, y = 0, xend = PC1 * sc, yend = PC2 * sc),
                 arrow = arrow(length = unit(1.7, "pt")),
                 linewidth = 0.28, colour = "grey20") +
    ggrepel::geom_text_repel(data = l, inherit.aes = FALSE,
                             aes(PC1 * sc, PC2 * sc, label = LAB[trait]),
                             size = 1.7, colour = "grey10", segment.size = 0.15,
                             min.segment.length = 0, max.overlaps = 20) +
    labs(title = sprintf("%s (n = %d)", TAXA_LAB[cl], nrow(a)),
         x = sprintf("PC1 (%.1f%%)", v$var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", v$var_explained[2])) +
    theme_pub()
}) |> Filter(f = Negate(is.null))
for (i in seq_along(biplots)) biplots[[i]] <- biplots[[i]] + labs(tag = letters[i])

# (e) 载荷热图 —— 让「轴是什么」一目了然
LD2 <- LD |> pivot_longer(PC1:PC3, names_to = "axis", values_to = "load")
LD2$trait_f <- factor(LAB[LD2$trait], levels = rev(LAB[CORE6]))
LD2$class_f <- factor(TAXA_LAB[LD2$class], levels = unname(TAXA_LAB))
pe <- ggplot(LD2, aes(axis, trait_f, fill = load)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f", load)), size = 1.8) +
  scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-1, 1),
                   name = "Loading") +
  facet_wrap(~class_f, nrow = 1) +
  labs(x = NULL, y = NULL, tag = "e") +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      legend.position = "right",
                      legend.key.width = unit(5, "pt"))

# (f) 方差解释
VE$class_f <- factor(TAXA_LAB[VE$class], levels = unname(TAXA_LAB))
pf <- ggplot(VE, aes(class_f, var_explained, fill = axis)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f", var_explained)),
            position = position_dodge(width = 0.78), vjust = -0.35, size = 1.8) +
  scale_fill_manual(values = c(PC1 = "#2E6E8E", PC2 = "#8FBF9F", PC3 = "#E8C07D"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Variance explained (%)", tag = "f") +
  theme_pub() + theme(legend.position = "bottom")

p <- wrap_plots(biplots, nrow = 1) / (pe | pf) +
  plot_layout(heights = c(1.15, 1), widths = c(2.4, 1))
save_fig(p, "Fig2_trait_axes", W2, 135, FIG)

log_msg("=== 03d 完成 / done ===")
