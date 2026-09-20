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

n <- 0L; missing <- character(0)
for (i in seq_len(nrow(man))) {
  dest <- dirs[[man$panel[i]]]
  for (ext in c("pdf", "png")) {
    src <- file.path("05_figures", paste0(man$src[i], ".", ext))
    if (file.exists(src)) {
      file.copy(src, file.path(dest, paste0(man$num[i], ".", ext)), overwrite = TRUE)
      n <- n + 1L
    } else missing <- c(missing, src)
  }
}
cat(sprintf("已装配 %d 个文件（%d 张图 x 2 种格式）\n", n, nrow(man)))
if (length(missing)) { cat("缺失 / missing:\n"); cat(paste0("  ", missing, collapse = "\n"), "\n") }
