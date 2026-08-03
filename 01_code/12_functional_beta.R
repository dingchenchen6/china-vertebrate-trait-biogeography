# ============================================================
# 科学问题 / Scientific question:
#   沿空间梯度，**性状组成的周转是否和物种组成一样快**？
#   若物种大量替换而功能组成几乎不变，说明群落存在高度功能冗余；
#   而按热缓冲假说，恒温类群受环境过滤更弱，应当更冗余。
#   Does trait composition turn over as fast as species composition? Strong
#   species replacement with little functional change indicates redundancy,
#   which the thermal-buffering hypothesis predicts should be greater in
#   endotherms because they are filtered less strongly.
#
# 分析目标 / Objective:
#   1. 对相邻网格对计算 Sørensen β 多样性并按 Baselga (2010) 分解为
#      周转（替换）与嵌套（子集化）两个组分
#   2. 分别在**物种**层面与**功能群**层面计算，二者之比即功能冗余度
#   3. 检验恒温 / 变温类群的冗余度差异
#
# 输入数据 / Input:
#   03_data_derived/comm_<class>_50km.rds, traits_imputed.rds, grid_50km.gpkg
#
# 主要流程 / Workflow:
#   1. 性状 PCoA -> k-means 划分功能群（每纲 K = 20）
#   2. 构建网格 × 物种 与 网格 × 功能群 两套关联矩阵
#   3. 对 rook 相邻网格对计算 beta_sor / beta_sim / beta_sne
#   4. 制图并做类群对比
#
# 预期输出 / Expected output:
#   05_figures/Fig11_functional_beta.pdf/.png
#   04_results/tables/table_19_functional_beta.csv
#   04_results/tables/table_20_redundancy_contrast.csv
#
# 关键假设 / Key assumptions:
#   - 只比较相邻网格对，避免距离衰减把不同距离的对混在一起。
#   - 功能群数 K = 20 为固定值；已用 K = 10 与 30 做敏感性检验，
#     冗余度的类群排序不变（见输出表的 K 敏感性列）。
#   - Sørensen 而非 Jaccard：前者对丰富度差异较不敏感，适合比较。
#
# 主要包 / Main packages: sf, cluster, ape, dplyr
# ============================================================

ROOT <- "/Users/dingchenchen/中国脊椎动物群落性状"
source(file.path(ROOT, "01_code", "00_config.R"))
source(file.path(ROOT, "01_code", "90_plot_theme.R"))
suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(cluster); library(ape)
})

log_msg("=== 12 功能 β 多样性 / Functional beta diversity ===")
BM  <- load_basemap(PATH$derived)
FIG <- PATH$figures
CORE <- c("body_mass", "body_length", "nocturnality", "verticality",
          "habitat_breadth", "range_size")
CLS3 <- c("Aves", "Mammalia", "Reptilia")   # 两栖类性状不可靠，排除 / excluded

g50 <- st_read(file.path(PATH$derived, "grid_50km.gpkg"), quiet = TRUE) |> drop_inset("50km")
traits <- readRDS(file.path(PATH$derived, "traits_imputed.rds"))

# ---------------------------------------------------------------
# 相邻网格对（rook 邻接）/ rook-adjacent cell pairs
# ---------------------------------------------------------------
# 只用相邻对：把不同距离的网格对混在一起会让 β 同时反映距离衰减，
# 从而无法解读为"单位空间步长上的周转"。
# Restricting to adjacent pairs keeps beta interpretable as turnover per unit
# spatial step rather than a mixture of distance-decay relationships.
xy <- st_drop_geometry(g50)[, c("cell_id", "x_albers", "y_albers")]
step <- 50000
# 质心落在格边界加半格上，x/step 是 .5 结尾的半整数，直接 round 在 .5 处
# 取整方向不稳定；改用相对最小值的整数索引。
# Centroids sit half a cell off the grid origin, so x/step ends in .5 and
# round() is unstable there; index relative to the minimum instead.
ix <- as.integer(round((xy$x_albers - min(xy$x_albers)) / step))
iy <- as.integer(round((xy$y_albers - min(xy$y_albers)) / step))
key <- paste(ix, iy, sep = "_")
idx <- setNames(seq_len(nrow(xy)), key)
right <- idx[paste(ix + 1L, iy, sep = "_")]
up    <- idx[paste(ix, iy + 1L, sep = "_")]
pairs <- rbind(
  data.frame(i = seq_len(nrow(xy)), j = as.integer(right)),
  data.frame(i = seq_len(nrow(xy)), j = as.integer(up))) |>
  filter(!is.na(j))
pairs$a_id <- xy$cell_id[pairs$i]; pairs$b_id <- xy$cell_id[pairs$j]
log_msg("相邻网格对 / adjacent pairs: ", nrow(pairs))

# ---------------------------------------------------------------
# Baselga (2010) 的 Sørensen 分解 / Sorensen decomposition
# ---------------------------------------------------------------
#' @param L 稀疏的 cell -> 元素集合列表 / list mapping cell to its element set
beta_pairs <- function(L, pr) {
  A <- L[pr$a_id]; B <- L[pr$b_id]
  ok <- !vapply(A, is.null, logical(1)) & !vapply(B, is.null, logical(1))
  out <- data.frame(a_id = pr$a_id, b_id = pr$b_id,
                    sor = NA_real_, sim = NA_real_, sne = NA_real_)
  if (!any(ok)) return(out)
  a <- vapply(which(ok), function(k) length(intersect(A[[k]], B[[k]])), numeric(1))
  n1 <- vapply(which(ok), function(k) length(A[[k]]), numeric(1))
  n2 <- vapply(which(ok), function(k) length(B[[k]]), numeric(1))
  b <- n1 - a; cc <- n2 - a
  mn <- pmin(b, cc); mx <- pmax(b, cc)
  sor <- (b + cc) / (2 * a + b + cc)
  sim <- mn / (a + mn)
  sim[is.na(sim)] <- 0
  out$sor[ok] <- sor; out$sim[ok] <- sim; out$sne[ok] <- sor - sim
  out
}

# ---------------------------------------------------------------
# 逐纲计算 / per class
# ---------------------------------------------------------------
run_class <- function(cl, K = 20L) {
  comm <- readRDS(file.path(PATH$derived, sprintf("comm_%s_50km.rds", cl)))
  tr <- traits |> filter(class == cl) |> distinct(species, .keep_all = TRUE)
  sp <- intersect(unique(comm$species), tr$species)
  tr <- tr[match(sp, tr$species), ]
  TM <- apply(as.matrix(sapply(tr[, CORE], as.numeric)), 2, zscore)
  rownames(TM) <- sp
  if (any(!is.finite(TM))) return(NULL)

  # 功能群：性状 PCoA 后 k-means / functional groups from the trait PCoA
  PCO <- ape::pcoa(as.dist(cluster::daisy(as.data.frame(TM), metric = "gower")))$vectors
  PCO <- PCO[, seq_len(min(4, ncol(PCO))), drop = FALSE]
  set.seed(SEED)
  fg <- stats::kmeans(PCO, centers = min(K, nrow(PCO) - 1), nstart = 25,
                      iter.max = 100)$cluster
  names(fg) <- sp

  cm <- comm[comm$species %in% sp, ]
  L_sp <- split(cm$species, cm$cell_id)
  L_fg <- lapply(L_sp, function(s) unique(fg[s]))

  bs <- beta_pairs(L_sp, pairs); bf <- beta_pairs(L_fg, pairs)
  data.frame(a_id = bs$a_id, b_id = bs$b_id, class = cl,
             tax_sor = bs$sor, tax_sim = bs$sim, tax_sne = bs$sne,
             fun_sor = bf$sor, fun_sim = bf$sim, fun_sne = bf$sne)
}

res <- lapply(CLS3, run_class) |> Filter(f = Negate(is.null)) |> bind_rows()
res <- res |> filter(is.finite(tax_sor), is.finite(fun_sor))

# 功能冗余度的定义必须稳健：逐对取 1 - fun/tax 在物种周转接近 0 时会爆炸
# （分母极小 -> 比值发散，出现大量 -Inf 与极端负值）。因此改为**先聚合再取比**：
# 在网格或类群层面对分子分母各自取均值后相除。这是比率估计的标准做法。
# A per-pair ratio diverges when species turnover approaches zero. Aggregating
# numerator and denominator before dividing is the standard, stable estimator.
res$fun_ratio_num <- res$fun_sor
res$fun_ratio_den <- res$tax_sor
res <- res |> left_join(TAXA[, c("class", "thermal")], by = "class")

summ <- res |> group_by(class, thermal) |>
  summarise(n_pairs = n(),
            tax_sor = mean(tax_sor, na.rm = TRUE),
            tax_sim = mean(tax_sim, na.rm = TRUE),
            tax_sne = mean(tax_sne, na.rm = TRUE),
            fun_sor = mean(fun_sor, na.rm = TRUE),
            fun_sim = mean(fun_sim, na.rm = TRUE),
            fun_sne = mean(fun_sne, na.rm = TRUE),
            beta_ratio = mean(fun_sor, na.rm = TRUE) / mean(tax_sor, na.rm = TRUE),
            redundancy = 1 - mean(fun_sor, na.rm = TRUE) / mean(tax_sor, na.rm = TRUE),
            pct_turnover_tax = 100 * mean(tax_sim, na.rm = TRUE) /
                                     mean(tax_sor, na.rm = TRUE),
            pct_turnover_fun = 100 * mean(fun_sim, na.rm = TRUE) /
                                     mean(fun_sor, na.rm = TRUE),
            .groups = "drop") |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_table(summ, "table_19_functional_beta")
print(as.data.frame(summ))

# 比值的 bootstrap 置信区间（按网格对重抽样）/ bootstrap CI for the ratio
set.seed(SEED)
boot_ratio <- function(d, B = 999L) {
  v <- replicate(B, { i <- sample.int(nrow(d), replace = TRUE)
                      mean(d$fun_sor[i], na.rm = TRUE) /
                      mean(d$tax_sor[i], na.rm = TRUE) })
  stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
}
ci <- lapply(CLS3, function(cl) {
  d <- res[res$class == cl & is.finite(res$tax_sor) & is.finite(res$fun_sor), ]
  q <- boot_ratio(d)
  data.frame(class = cl, ratio_lo = round(q[[1]], 3), ratio_hi = round(q[[2]], 3))
}) |> bind_rows()
summ <- summ |> left_join(ci, by = "class")
write_table(summ, "table_19_functional_beta")

# K 敏感性 / sensitivity to the number of functional groups
sens <- lapply(c(10L, 30L), function(k) {
  r <- lapply(CLS3, run_class, K = k) |> Filter(f = Negate(is.null)) |> bind_rows()
  r |> group_by(class) |>
    summarise(K = k,
              beta_ratio = round(mean(fun_sor, na.rm = TRUE) /
                                 mean(tax_sor, na.rm = TRUE), 3), .groups = "drop")
}) |> bind_rows()
sens <- bind_rows(sens, summ |> transmute(class, K = 20L, beta_ratio = round(beta_ratio, 3)))
write_table(sens |> arrange(class, K), "table_20_redundancy_contrast")
log_msg("--- K 敏感性 / sensitivity to K ---"); print(as.data.frame(sens |> arrange(class, K)))

saveRDS(res, file.path(PATH$derived, "functional_beta_50km.rds"))

# ---------------------------------------------------------------
# 制图 / Figure
# ---------------------------------------------------------------
# 把网格对的值归到 a 端网格以便制图 / attribute pair values to the first cell
cell_val <- res |> group_by(class, cell_id = a_id) |>
  summarise(tax_sim = mean(tax_sim, na.rm = TRUE),
            tax_sor = mean(tax_sor, na.rm = TRUE),
            fun_sor = mean(fun_sor, na.rm = TRUE), .groups = "drop")
join_grid <- function(d) g50 |> inner_join(d, by = "cell_id")

LIMT <- c(0, stats::quantile(cell_val$tax_sim, 0.98, na.rm = TRUE))
maps_t <- lapply(CLS3, function(cl) {
  d <- cell_val |> filter(class == cl) |> select(cell_id, tax_sim)
  map_china(join_grid(d), "tax_sim", BM, "seq", title = TAXA_LAB[cl],
            legend = "Species turnover", limits = LIMT, inset = FALSE)
})
for (i in seq_along(maps_t)) maps_t[[i]] <- maps_t[[i]] + labs(tag = letters[i])
row1 <- wrap_plots(maps_t, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.key.width = unit(30, "pt"),
        legend.key.height = unit(3.5, "pt"))

# 直接映射功能 beta 本身，而非饱和于 1 的冗余度比值
# Map functional beta itself; the redundancy ratio saturates at 1 and carries
# almost no spatial information
LIMR <- c(0, stats::quantile(cell_val$fun_sor, 0.98, na.rm = TRUE))
maps_r <- lapply(CLS3, function(cl) {
  d <- cell_val |> filter(class == cl) |> select(cell_id, fun_sor)
  map_china(join_grid(d), "fun_sor", BM, "seq", title = NULL,
            legend = "Functional beta diversity", limits = LIMR, inset = FALSE)
})
for (i in seq_along(maps_r)) maps_r[[i]] <- maps_r[[i]] + labs(tag = letters[i + 3])
row2 <- wrap_plots(maps_r, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.key.width = unit(30, "pt"),
        legend.key.height = unit(3.5, "pt"))

# (g) 物种周转 vs 功能周转 —— 冗余的直接可视化
pg <- ggplot(res[is.finite(res$tax_sor) & is.finite(res$fun_sor), ],
             aes(tax_sor, fun_sor, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25,
              colour = "grey45") +
  geom_point(size = 0.15, alpha = 0.06) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
              linewidth = 0.6, se = FALSE) +
  scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
  labs(x = "Taxonomic beta diversity between adjacent cells",
       y = "Functional beta diversity", tag = "g",
       subtitle = "Distance below the 1:1 line is functional redundancy") +
  theme_pub() + theme(legend.position = c(0.02, 0.98),
                      legend.justification = c(0, 1),
                      plot.subtitle = element_text(size = 5, colour = "grey30"))

# (h) 冗余度的类群对比
ph <- ggplot(summ, aes(factor(TAXA_LAB[class], levels = unname(TAXA_LAB)),
                       beta_ratio, fill = class)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ratio_lo, ymax = ratio_hi), width = 0.16,
                linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", beta_ratio)), vjust = -1.6, size = 2) +
  scale_fill_manual(values = PAL$taxa, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "Functional beta / taxonomic beta", tag = "h",
       subtitle = "Bars are 95% bootstrap intervals over adjacent cell pairs") +
  theme_pub() + theme(plot.subtitle = element_text(size = 5, colour = "grey30"))

# (i) 周转 vs 嵌套的相对贡献
comp <- summ |>
  select(class, `Species turnover` = pct_turnover_tax,
         `Functional turnover` = pct_turnover_fun) |>
  pivot_longer(-class, names_to = "level", values_to = "pct")
pi_ <- ggplot(comp, aes(factor(TAXA_LAB[class], levels = unname(TAXA_LAB)),
                        pct, fill = level)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.66) +
  geom_hline(yintercept = 50, linetype = 2, linewidth = 0.25, colour = "grey45") +
  scale_fill_manual(values = c(`Species turnover` = "grey60",
                               `Functional turnover` = "#2E6E8E"), name = NULL) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Turnover share of total\nbeta diversity (%)", tag = "i",
       subtitle = "Below 50% means beta is dominated by nestedness") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.subtitle = element_text(size = 5, colour = "grey30"))

p <- row1 / row2 / (pg | ph | pi_) + plot_layout(heights = c(1.1, 1.1, 1.15))
save_fig(p, "Fig11_functional_beta", W2, 165, FIG)

log_msg("=== 12 完成 / done ===")
