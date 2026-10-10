# ============================================================
# Scientific question / 科学问题:
#   （工具脚本，不涉及科学问题）
# Objective / 分析目标:
#   按 table_00_figure_manifest.csv 把 05_figures/ 的分析图重命名并装配进投稿包。
#   Rename and assemble the analysis figures into the submission package,
#   following the figure manifest.
# Input data / 输入数据:
#   04_results/tables/table_00_figure_manifest.csv; 05_figures/*.{pdf,png}
# Expected output / 预期输出:
#   09_submission/03_main_figures/Fig1–Fig8.{pdf,png}
#   09_submission/04_extended_data/ED1–ED25.{pdf,png}
# Key assumption / 关键假设:
#   清单中的 src 列与 05_figures/ 下的文件名一一对应。
# Output directory / 输出路径: 09_submission/
# ============================================================
man <- utils::read.csv("04_results/tables/table_00_figure_manifest.csv",
                       stringsAsFactors = FALSE)
dirs <- c(Main = "09_submission/03_main_figures",
          `Extended Data` = "09_submission/04_extended_data")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# 硬链接防护：若目标文件恰是源文件的硬链接，就地覆盖或清理目标会连带改写分析原图。
# 2026-10 曾发生过：投稿包图件与 05_figures 原图共享 inode，对投稿包运行 wm-clean
# 时连带改写了原图的元数据。因此先删除目标再复制，并校验目标的硬链接数为 1。
# Hard-link guard: writing to a target that shares an inode with its source would
# rewrite the source. Unlink first, copy, then assert the target has one link.
nlink <- function(p) {
  for (args in list(c("-f", "%l"), c("-c", "%h"))) {
    out <- suppressWarnings(tryCatch(
      system2("stat", c(args, shQuote(p)), stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)))
    if (length(out) && grepl("^[0-9]+$", out[1])) return(as.integer(out[1]))
  }
  NA_integer_
}

n <- 0L; missing <- character(0)
for (i in seq_len(nrow(man))) {
  dest <- dirs[[man$panel[i]]]
  for (ext in c("pdf", "png")) {
    src <- file.path("05_figures", paste0(man$src[i], ".", ext))
    if (file.exists(src)) {
      dest_f <- file.path(dest, paste0(man$num[i], ".", ext))
      if (file.exists(dest_f)) unlink(dest_f)          # 先删除，断开可能存在的硬链接
      file.copy(src, dest_f)
      k <- nlink(dest_f)
      if (!is.na(k) && k > 1L) stop("目标文件与他处共享 inode / shared inode: ", dest_f)
      n <- n + 1L
    } else missing <- c(missing, src)
  }
}
cat(sprintf("已装配 %d 个文件（%d 张图 x 2 种格式）\n", n, nrow(man)))
if (length(missing)) { cat("缺失 / missing:\n"); cat(paste0("  ", missing, collapse = "\n"), "\n") }
