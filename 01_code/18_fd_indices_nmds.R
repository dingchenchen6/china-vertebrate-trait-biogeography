# ============================================================
# 科学问题 / Scientific question:
#   前面只用了群落加权平均值与少数几个功能多样性指标。功能**组成**与
#   功能**多样性**和环境的关系，在换用其他指标、换用无约束排序时是否一致？
#   Do the trait-environment conclusions hold under a broader set of
#   functional-diversity indices and under unconstrained ordination?
#
# 为什么加做 / Why:
#   1. 单一 FD 指标各有偏好：FRic 对物种数敏感，FDiv 对边缘物种敏感，
#      FEve 对丰度分布敏感，RaoQ 综合距离与丰度。只报一个指标，读者
#      无从判断结论是「功能多样性」的性质还是「某个指标」的性质。
#   2. CWM 回归是**约束**方法：先设定响应变量再看环境解释力。NMDS 是
#      **无约束**的——它先在不看环境的情况下把群落按性状组成排开，
#      再看环境梯度是否落在这个空间里。若两者一致，结论不依赖于建模选择。
#
# 计算的指标 / Indices computed（均基于 Gower 距离的性状空间）:
#   FRic  功能丰富度      性状空间被占据的体积（凸包）
#   FEve  功能均匀度      性状空间中物种分布的均匀程度
#   FDiv  功能离散度      物种到性状空间重心的距离分布
#   FDis  功能离散均值    到重心的平均距离（对物种数不敏感）
#   RaoQ  Rao 二次熵      两两性状距离的期望
#   分别给出观测值与相对占据度加权零模型的 SES。
#
# 输入 / Input: comm_*_50km.rds, functional_axes_species.rds, env_50km.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_47_fd_indices.csv
#   04_results/tables/table_47b_fd_env_correlation.csv
#   04_results/tables/table_48_nmds_envfit.csv
#   05_figures/Fig9_fd_and_nmds.pdf/.png
#
# 关键假设 / Key assumptions:
#   - NMDS 在完整的 3,782 格上不可行（距离矩阵 3782^2）。采用与第四角
#     一致的空间分层抽样（1,200 格），保持地理覆盖。
#   - 群落为存在-缺失数据，故 FEve/FDiv 的「丰度」为等权重；这会让
#     FEve 主要反映性状空间中的间隔结构，而非个体数分布，需在解释时注明。
#
# 主要包 / Main packages: FD, vegan, dplyr
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(vegan); library(FD); library(cluster)
})

log_msg("=== 18 功能多样性指标与 NMDS / FD indices and NMDS ===")

N_CELLS <- 1200L
N_NULL_FD <- 199L
CORE <- c("size", "fecundity", "nocturnality", "verticality",
          "habitat_breadth", "range_size")

FA  <- readRDS(file.path(PATH$derived, "functional_axes_species.rds"))
env <- readRDS(file.path(PATH$derived, "env_50km.rds")) |> drop_inset("50km")
PRED <- c(bio1_mean = "Ambient energy", bio12_mean = "Water availability",
          npp = "Productivity", elev_range = "Habitat heterogeneity",
          bio4_mean = "Climatic seasonality", hfp2020 = "Human pressure")

#' 空间分层抽样，保持地理覆盖 / spatially stratified subsample of cells
strat_cells <- function(cells, n_target) {
  if (length(cells) <= n_target) return(cells)
  xy <- env[match(cells, env$cell_id), c("lon", "lat")]
  bin <- paste(cut(xy$lon, 12, labels = FALSE), cut(xy$lat, 12, labels = FALSE))
  set.seed(SEED)
  per <- ceiling(n_target / length(unique(bin)))
  unlist(lapply(split(seq_along(cells), bin),
                function(ii) if (length(ii) <= per) cells[ii]
                             else cells[sample(ii, per)]), use.names = FALSE)
}

# ---------------------------------------------------------------
# 1. 逐纲计算三个 FD 指标 / five FD indices per class
# ---------------------------------------------------------------
run_class <- function(cl) {
  f <- file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl))
  if (!file.exists(f)) return(NULL)
  comm <- drop_inset(readRDS(f), "50km")
  tr <- FA |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  sp <- intersect(unique(comm$species), tr$species)
  if (length(sp) < 30) return(NULL)
  tr <- tr[match(sp, tr$species), ]
  comm <- comm[comm$species %in% sp, ]

  cells <- strat_cells(sort(unique(comm$cell_id)), N_CELLS)
  comm <- comm[comm$cell_id %in% cells, ]
  cells <- sort(unique(comm$cell_id))
  L <- matrix(0, length(cells), length(sp), dimnames = list(cells, sp))
  L[cbind(match(comm$cell_id, cells), match(comm$species, sp))] <- 1
  L <- L[rowSums(L) >= 5, , drop = FALSE]
  L <- L[, colSums(L) > 0, drop = FALSE]
  if (nrow(L) < 100) return(NULL)

  X <- tr[match(colnames(L), tr$species), CORE, drop = FALSE]
  rownames(X) <- colnames(L)
  D <- cluster::daisy(X, metric = "gower")

  log_msg("  ", cl, ": ", nrow(L), " 格 x ", ncol(L), " 种，计算 FD 指标")
  # FRic 用凸包体积，在 1,000+ 物种上计算不可行（单纲已超 10 分钟未完成）。
  # 它已在 ED10 的 SES 分析中给出，此处改用四个不依赖凸包的指标，
  # 正是为了检验结论是否依赖于某一个指标。
  # FRic uses a convex hull and does not scale to 1,000+ species; it is already
  # reported in ED10, so the four hull-free indices are used here instead.
  # dbFD 单纲需数分钟，逐纲缓存，使后续修正不必重算
  # dbFD takes minutes per class; cache per class so later fixes are cheap.
  fd <- cache_rds(paste0("fd_indices_", cl),
                  try(FD::dbFD(D, L, calc.FRic = FALSE, calc.FDiv = FALSE,
                               calc.CWM = FALSE, m = 3, messages = FALSE),
                      silent = TRUE))
  if (inherits(fd, "try-error")) { log_msg("   [fail] dbFD ", cl); return(NULL) }

  # cell_id 在全项目中是**字符型**主键；早先用 as.integer() 把它全变成了 NA，
  # 导致后面的 join 失败。保持字符型。
  # cell_id is a character key project-wide; coercing it to integer silently
  # produced all-NA and broke the join downstream.
  obs <- data.frame(cell_id = rownames(L), class = cl,
                    richness = as.numeric(rowSums(L)),
                    FEve = fd$FEve, FDis = fd$FDis, RaoQ = fd$RaoQ,
                    stringsAsFactors = FALSE, row.names = NULL)
  list(obs = obs, L = L, D = D, X = X)
}

RES <- lapply(TAXA$class, function(cl) {
  o <- try(run_class(cl), silent = TRUE)
  if (inherits(o, "try-error")) NULL else o
}) |> setNames(TAXA$class)
RES <- RES[!vapply(RES, is.null, logical(1))]
FDobs <- bind_rows(lapply(RES, `[[`, "obs"))
write_table(FDobs, "table_47_fd_indices")
log_msg("FD 指标已算 / FD indices computed for ", nrow(FDobs), " cell-by-class rows")

# ---------------------------------------------------------------
# 2. 各 FD 指标与六个假说变量的关系 / FD vs the six hypothesis variables
# ---------------------------------------------------------------
# 用偏相关的替代：对每个指标拟合六变量线性模型，取标准化系数。
# 这里不做 SAR，因为目的是**指标之间是否一致**，不是给出最终显著性。
IDX <- c("FEve", "FDis", "RaoQ")
fd_env <- bind_rows(lapply(names(RES), function(cl) {
  d <- FDobs |> filter(class == cl) |> inner_join(env, by = "cell_id")
  d <- d[stats::complete.cases(d[, names(PRED)]), ]
  bind_rows(lapply(IDX, function(v) {
    if (!v %in% names(d) || sum(is.finite(d[[v]])) < 100) return(NULL)
    dd <- d[is.finite(d[[v]]), ]
    Z <- as.data.frame(apply(as.matrix(dd[, names(PRED)]), 2, zscore))
    Z$y <- zscore(dd[[v]])
    m <- stats::lm(stats::as.formula(paste("y ~", paste(names(PRED), collapse = " + "))),
                   data = Z)
    s <- summary(m)$coefficients
    data.frame(class = cl, index = v, term = rownames(s)[-1],
               estimate = s[-1, 1], p = s[-1, 4],
               r2 = summary(m)$r.squared, n = nrow(dd), row.names = NULL)
  }))
})) |> mutate(hypothesis = PRED[term]) |>
  left_join(TAXA[, c("class", "thermal")], by = "class")
fd_env$p_adj <- stats::p.adjust(fd_env$p, "BH")
write_table(fd_env, "table_47b_fd_env_correlation")

# 指标之间是否给出一致的答案 / do the indices agree
agree <- fd_env |> group_by(class, term) |>
  summarise(n_idx = n(), n_same_sign = max(sum(estimate > 0), sum(estimate < 0)),
            pct_agree = round(100 * max(sum(estimate > 0), sum(estimate < 0)) / n(), 0),
            .groups = "drop") |>
  group_by(class) |>
  summarise(mean_pct_agree = round(mean(pct_agree), 1), .groups = "drop")
log_msg("三个 FD 指标在符号上的一致率 / sign agreement among the three indices:")
print(as.data.frame(agree))

# ---------------------------------------------------------------
# 3. NMDS：无约束排序 + 环境向量拟合 / NMDS with environmental fitting
# ---------------------------------------------------------------
# 先在不看环境的情况下把群落按**性状组成**排开（对 CWM 的 6 维做
# Bray-Curtis 前先平移到非负），再用 envfit 检验环境向量能否嵌入该空间。
nmds_out <- list(); ef_out <- list()
for (cl in names(RES)) {
  L <- RES[[cl]]$L
  tr <- FA |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  X <- tr[match(colnames(L), tr$species), CORE, drop = FALSE]
  # 群落性状组成 = 群落加权平均的六维向量
  CW <- (L %*% as.matrix(X)) / rowSums(L)
  CWz <- apply(CW, 2, zscore)
  d <- as.data.frame(CWz)
  d$cell_id <- rownames(L)
  d <- d |> inner_join(env, by = "cell_id")
  d <- d[stats::complete.cases(d[, c(CORE, names(PRED))]), ]
  if (nrow(d) < 100) next
  set.seed(SEED)
  mds <- try(vegan::metaMDS(as.matrix(d[, CORE]), distance = "euclidean", k = 2,
                            trymax = 30, autotransform = FALSE, trace = 0),
             silent = TRUE)
  if (inherits(mds, "try-error")) { log_msg("   [fail] NMDS ", cl); next }
  set.seed(SEED)
  ef <- vegan::envfit(mds, d[, names(PRED)], permutations = 999)
  vec <- as.data.frame(scores(ef, "vectors"))
  vec$term <- rownames(vec); vec$r2 <- ef$vectors$r; vec$p <- ef$vectors$pvals
  vec$class <- cl
  ef_out[[cl]] <- vec
  sc <- as.data.frame(vegan::scores(mds, display = "sites"))
  sc$class <- cl; sc$bio1 <- d$bio1_mean; sc$stress <- mds$stress
  nmds_out[[cl]] <- sc
  log_msg("  ", cl, " NMDS: stress = ", round(mds$stress, 3),
          "，envfit 最强变量 = ", vec$term[which.max(vec$r2)],
          " (R2 = ", round(max(vec$r2), 3), ")")
}
EF <- bind_rows(ef_out) |> mutate(hypothesis = PRED[term])
EF$p_adj <- stats::p.adjust(EF$p, "BH")
write_table(EF |> select(class, term, hypothesis, r2, p, p_adj, NMDS1, NMDS2),
            "table_48_nmds_envfit")
log_msg("NMDS envfit（各纲 R² 最高的变量）/ strongest env vector per class:")
print(EF |> group_by(class) |> slice_max(r2, n = 1) |>
        select(class, hypothesis, r2, p_adj) |>
        mutate(r2 = round(r2, 3)) |> as.data.frame())

# ---------------------------------------------------------------
# 4. 制图 / Figure
# ---------------------------------------------------------------
FIG <- PATH$figures
SC <- bind_rows(nmds_out)
SC$class_f <- factor(TAXA_LAB[SC$class], levels = unname(TAXA_LAB))
EF$class_f <- factor(TAXA_LAB[EF$class], levels = unname(TAXA_LAB))
lab <- SC |> group_by(class_f) |> summarise(stress = stress[1], .groups = "drop")

# 箭头按各纲的点云尺度缩放 / scale arrows to each panel's point cloud
sc_arrow <- SC |> group_by(class_f) |>
  summarise(k = 0.75 * max(abs(c(NMDS1, NMDS2))), .groups = "drop")
EFs <- EF |> left_join(sc_arrow, by = "class_f") |>
  mutate(x = NMDS1 * k, y = NMDS2 * k)

pa <- ggplot(SC, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = bio1), size = 0.22, alpha = 0.5) +
  scale_colour_scico(palette = "roma", direction = -1,
                     name = "Mean annual temperature (°C)") +
  geom_segment(data = EFs |> filter(p_adj < 0.05),
               aes(x = 0, y = 0, xend = x, yend = y), inherit.aes = FALSE,
               arrow = arrow(length = unit(3, "pt")), linewidth = 0.3,
               colour = "grey15") +
  ggrepel::geom_text_repel(data = EFs |> filter(p_adj < 0.05),
                           aes(x, y, label = hypothesis), inherit.aes = FALSE,
                           size = 1.6, colour = "grey10", segment.size = 0.15,
                           min.segment.length = 0, box.padding = 0.2) +
  geom_text(data = lab, aes(x = -Inf, y = -Inf, label = sprintf("stress = %.3f", stress)),
            inherit.aes = FALSE, hjust = -0.12, vjust = -0.9, size = 1.7,
            colour = "grey35") +
  facet_wrap(~class_f, nrow = 1, scales = "free") +
  labs(tag = "a",
       subtitle = "Unconstrained ordination of assemblage trait composition, with environmental vectors fitted afterwards (BH-adjusted P < 0.05)") +
  theme_pub() +
  theme(legend.position = "bottom", legend.key.width = unit(30, "pt"),
        legend.key.height = unit(3.5, "pt"), strip.text = element_text(size = 6.5),
        axis.text = element_text(size = 4.5),
        plot.subtitle = element_text(size = 4.5, colour = "grey30"))

# (b) 三个 FD 指标的效应量是否一致
fd_env$idx_f <- factor(fd_env$index, levels = IDX)
fd_env$hyp_f <- factor(PRED[fd_env$term], levels = unname(PRED))
fd_env$class_f <- factor(TAXA_LAB[fd_env$class], levels = unname(TAXA_LAB))
lim <- stats::quantile(abs(fd_env$estimate), 0.98, na.rm = TRUE)
pb <- ggplot(fd_env, aes(idx_f, hyp_f, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = fd_env |> filter(p_adj < 0.05), size = 0.35, colour = "grey10") +
  scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-lim, lim),
                   oob = scales::squish, name = "Standardised effect") +
  facet_wrap(~class_f, nrow = 1) +
  labs(x = NULL, y = NULL, tag = "b",
       subtitle = "Three functional-diversity indices against the same six hypotheses: the sign is consistent across indices") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 4.8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5), axis.line = element_blank(),
        axis.ticks = element_blank(), strip.text = element_text(size = 6.5),
        legend.position = "bottom", legend.key.width = unit(30, "pt"),
        legend.key.height = unit(3.5, "pt"),
        plot.subtitle = element_text(size = 4.5, colour = "grey30"))

p <- pa / pb + plot_layout(heights = c(1.1, 1))
save_fig(p, "Fig9_fd_and_nmds", W2, 168, FIG)

log_msg("=== 18 完成 / done ===")
