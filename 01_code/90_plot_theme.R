# ============================================================
# 本文件 / This file:
#   出版级绘图主题与中国制图辅助函数
#   Publication-quality plotting theme and China mapping helpers
#
# 设计原则 / Design principles:
#   - Nature 风格：无冗余装饰、细线条、小字号、留白充足
#   - 感知均匀配色 (viridis / scico)，色盲友好
#   - 中国地图合规：使用审图号 GS(2023)2767 底图，
#     必须包含南海诸岛附图与完整国界线
#
# 主要包 / Main packages: ggplot2, sf, patchwork, scico, ggspatial
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(patchwork)
  library(scico); library(viridis); library(grid); library(scales); library(cowplot)
})

# ---------------------------------------------------------------
# 1. 字体与基础主题 / Fonts and base theme
# ---------------------------------------------------------------
BASE_SIZE <- 7          # Nature 正文图字号 / Nature figure font size
BASE_FAM  <- if (any(grepl("Helvetica", systemfonts::system_fonts()$family))) "Helvetica" else "sans"

theme_pub <- function(base_size = BASE_SIZE, base_family = BASE_FAM) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line        = element_line(linewidth = 0.3, colour = "grey20"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "grey20"),
      axis.ticks.length= unit(1.5, "pt"),
      axis.text        = element_text(colour = "grey20", size = base_size - 0.5),
      axis.title       = element_text(colour = "black", size = base_size),
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size - 0.5, colour = "grey35"),
      plot.tag         = element_text(size = base_size + 2, face = "bold"),
      legend.key.size  = unit(7, "pt"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 0.5),
      legend.background= element_blank(),
      legend.margin    = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size, face = "bold", hjust = 0),
      panel.spacing    = unit(4, "pt"),
      plot.margin      = margin(3, 3, 3, 3)
    )
}

theme_map <- function(base_size = BASE_SIZE, base_family = BASE_FAM) {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      legend.key.height = unit(4, "pt"),
      legend.key.width  = unit(24, "pt"),
      legend.position   = "bottom",
      legend.title      = element_text(size = base_size - 0.5, vjust = 1),
      legend.text       = element_text(size = base_size - 1.5),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.tag          = element_text(size = base_size + 2, face = "bold"),
      plot.margin       = margin(2, 2, 2, 2)
    )
}

# ---------------------------------------------------------------
# 2. 配色 / Palettes
# ---------------------------------------------------------------
PAL <- list(
  # 连续量（丰富度、性状均值）/ sequential
  seq      = function(...) scico::scale_fill_scico(palette = "batlow", ...),
  # 发散量（SES、效应量）/ diverging, centred on zero
  div      = function(...) scico::scale_fill_scico(palette = "vik", midpoint = 0, ...),
  div_col  = function(...) scico::scale_colour_scico(palette = "vik", midpoint = 0, ...),
  # 类群配色：恒温=暖色，变温=冷色 / taxa: endotherms warm, ectotherms cool
  taxa     = c(Aves = "#C1442E", Mammalia = "#E8A33D",
               Reptilia = "#3B7EA1", Amphibia = "#4C9F70"),
  thermal  = c(Endotherm = "#C1442E", Ectotherm = "#3B7EA1")
)

TAXA_LAB <- c(Aves = "Birds", Mammalia = "Mammals",
              Reptilia = "Reptiles", Amphibia = "Amphibians")
TAXA_LAB_CN <- c(Aves = "鸟类", Mammalia = "哺乳类",
                 Reptilia = "爬行类", Amphibia = "两栖类")

# ---------------------------------------------------------------
# 3. 中国制图底图要素 / China basemap layers
# ---------------------------------------------------------------
#' 载入中国制图底图（南海诸岛按真实地理位置绘制）
#' Load the China basemap with the South China Sea islands at true positions
#'
#' 底图选择 / Choice of basemap:
#'   GS(2023)2767 审图号数据把南海诸岛**预置在右侧附图框内**，出图后右下角
#'   会留出一个大空框，且诸岛与大陆的相对位置不真实。本项目改用真实位置的
#'   省界数据（范围至 3.83 N）配合独立的十段线图层，得到国内期刊常见的
#'   紧凑版式：诸岛以小点形式位于台湾以南，十段线完整绘出。
#'   The approved GS(2023)2767 dataset pre-places the South China Sea islands
#'   in a framed inset, which leaves a large empty box and misrepresents their
#'   position relative to the mainland. We instead use provincial boundaries at
#'   true positions together with a separate ten-dash-line layer, giving the
#'   compact layout usual in Chinese journals.
#'
#'   分析网格仍建立在 GS(2023)2767 的大陆边界上；两套边界在大陆部分一致，
#'   差异仅在于南海诸岛的摆放方式。
#'   The analysis grid still rests on the approved mainland boundary; the two
#'   datasets agree there and differ only in how the islands are placed.
load_basemap <- function(derived_dir) {
  get1 <- function(f) {
    p <- file.path(derived_dir, f)
    if (file.exists(p)) sf::st_read(p, quiet = TRUE) else NULL
  }
  list(
    land     = get1("cn_land_true.gpkg"),
    prov     = get1("cn_prov_true.gpkg"),
    border   = get1("cn_natborder.gpkg"),
    dashline = get1("cn_dashline.gpkg")
  )
}

# 制图范围必须向南延伸到十段线末端（y 约 3.36e5），否则南海诸岛被切掉。
# The window must reach the southern end of the ten-dash line, or the South
# China Sea islands are cut off.
CN_XLIM <- c(-3.06e6, 1.90e6)
CN_YLIM <- c(3.10e5, 5.90e6)

#' 叠加省界、国界与十段线
#' Overlay provincial boundaries, the national boundary and the dash line
add_basemap_lines <- function(bm, lwd = 0.22) {
  list(
    if (!is.null(bm$prov))
      ggplot2::geom_sf(data = bm$prov, fill = NA, colour = "grey72",
                       linewidth = lwd * 0.5),
    if (!is.null(bm$land))
      ggplot2::geom_sf(data = bm$land, fill = NA, colour = "grey15",
                       linewidth = lwd),
    if (!is.null(bm$dashline))
      ggplot2::geom_sf(data = bm$dashline, fill = NA, colour = "grey15",
                       linewidth = lwd * 0.9)
  )
}

map_china <- function(g, var, bm, type = c("seq", "div"),
                      title = NULL, legend = var, limits = NULL,
                      inset = TRUE, show_legend = TRUE) {
  type <- match.arg(type)
  g <- g[is.finite(g[[var]]), ]
  # 伪网格不进入任何图件 / inset artefacts never appear on a map
  if ("scs_inset" %in% names(g)) g <- g[!g$scs_inset, ]

  p <- ggplot() +
    { if (!is.null(bm$land))
        geom_sf(data = bm$land, fill = "grey95", colour = NA) } +
    geom_sf(data = g, aes(fill = .data[[var]]), colour = NA) +
    add_basemap_lines(bm) +
    { if (type == "seq")
        scico::scale_fill_scico(palette = "batlow", name = legend, limits = limits,
                                na.value = "grey88", oob = scales::squish,
                                guide = guide_colourbar(title.position = "top"))
      else
        scico::scale_fill_scico(palette = "vik", name = legend, midpoint = 0,
                                limits = limits, na.value = "grey88",
                                oob = scales::squish,
                                guide = guide_colourbar(title.position = "top")) } +
    coord_sf(xlim = CN_XLIM, ylim = CN_YLIM, expand = FALSE) +
    labs(title = title) +
    theme_map()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

#' 离散取值的中国地图（分区、主导性状等）
#' China map for a categorical variable, sharing the same compliant frame
map_china_cat <- function(g, var, bm, values, title = NULL, legend = NULL,
                          drop = FALSE) {
  if ("scs_inset" %in% names(g)) g <- g[!g$scs_inset, ]
  ggplot() +
    { if (!is.null(bm$land))
        geom_sf(data = bm$land, fill = "grey95", colour = NA) } +
    geom_sf(data = g, aes(fill = .data[[var]]), colour = NA) +
    add_basemap_lines(bm) +
    scale_fill_manual(values = values, name = legend, drop = drop,
                      na.value = "grey88") +
    coord_sf(xlim = CN_XLIM, ylim = CN_YLIM, expand = FALSE) +
    labs(title = title) +
    theme_map()
}

#' 多面板共享一个图例（避免每格重复图例占位）
#' Collect panels under a single shared legend
panels_shared_legend <- function(plots, ncol = length(plots), legend_from = 1) {
  leg <- cowplot::get_plot_component(
    plots[[legend_from]] + theme(legend.position = "bottom"),
    "guide-box-bottom", return_all = TRUE)
  body <- patchwork::wrap_plots(
    lapply(plots, function(p) p + theme(legend.position = "none")), ncol = ncol)
  patchwork::wrap_plots(body, patchwork::wrap_elements(leg),
                        ncol = 1, heights = c(1, 0.10))
}

#' 森林图：标准化效应量按类群对比
#' Forest plot of standardised effect sizes by taxon
forest_effects <- function(df, x = "estimate", se = "se",
                           y = "term", colour = "class",
                           facet = NULL, title = NULL, xlab = "Standardised effect") {
  df$.lo <- df[[x]] - 1.96 * df[[se]]
  df$.hi <- df[[x]] + 1.96 * df[[se]]
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour]])) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60", linetype = 2) +
    geom_linerange(aes(xmin = .lo, xmax = .hi),
                   position = position_dodge(width = 0.6), linewidth = 0.35) +
    geom_point(position = position_dodge(width = 0.6), size = 0.9) +
    scale_colour_manual(values = PAL$taxa, labels = TAXA_LAB, name = NULL) +
    labs(x = xlab, y = NULL, title = title) +
    theme_pub()
  if (!is.null(facet)) p <- p + facet_wrap(stats::as.formula(paste("~", facet)), scales = "free_x")
  p
}

#' 统一保存：同时输出 PDF（矢量）与 PNG（预览）
#' Save both a vector PDF and a raster PNG preview
save_fig <- function(p, name, width, height, dir, dpi = 600) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = width, height = height,
         units = "mm", device = grDevices::cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = width, height = height,
         units = "mm", dpi = dpi)
  message("  -> ", name, ".pdf / .png  (", width, " x ", height, " mm)")
}

# Nature 单栏 89 mm、双栏 183 mm / Nature single 89 mm, double 183 mm
W1 <- 89; W2 <- 183
