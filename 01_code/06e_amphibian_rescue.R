# ============================================================
# 科学问题 / Scientific question:
#   中国两栖类群落的功能结构，能否在**不依赖模型插补**的前提下判定？
#   Can the functional structure of Chinese amphibian assemblages be
#   determined without relying on modelled trait values?
#
# 背景 / Background:
#   主分析使用的六个核心性状中，两栖类的体重有 75%、活动节律有 71% 是
#   TetrapodTraits 的上游插补值（ED1c）。限定为五项全实测的 69 个物种时，
#   平均 SES.FDis 由 +1.85 反转为 −0.35（表 9），故原"功能超离散"结论被撤回。
#   本脚本检验一条替代路径：**改用一组本身就有较高实测率的性状**，
#   而不是削减物种数。
#   Rather than shrinking the species set, this script swaps the trait set for
#   one that is well measured in amphibians to begin with.
#
# 分析目标 / Objective:
#   1. 从 AmphiBIO v1 与 Etard et al. 2020 的原始库中提取两栖类的实测性状
#   2. 构建"高实测率性状集"：体长、垂直生境位置、生境宽度、分布范围
#      （剔除体重与活动节律——这两项在中国两栖类中确实未被测量）
#   3. 用该性状集重算功能离散度 SES，与主分析及实测子集三者对照
#
# 输入数据 / Input:
#   02_data_raw/traits/amphibio/AmphiBIO_v1.csv
#   02_data_raw/traits/Etard2020_amphibians.csv
#   03_data_derived/traits_imputed.rds, comm_Amphibia_50km.rds
#
# 预期输出 / Expected output:
#   04_results/tables/table_17_amphibian_trait_rescue.csv
#   04_results/tables/table_18_amphibian_ses_comparison.csv
#   03_data_derived/amphibian_observed_traits.rds
#
# 关键假设 / Key assumptions:
#   - AmphiBIO 的 Fos/Ter/Aqu/Arb 为二元生境利用标记，可加权合成与主分析
#     同义的 verticality（穴居 0 → 空中 4）。
#   - 分布范围由 IUCN 多边形直接计算，属实测量而非插补。
#   - 该检验仍受"有实测值的物种本身偏向研究充分者"的影响，故只能判断
#     结论是否稳健，不能证明真实格局。
#
# 主要包 / Main packages: data.table, cluster, ape, parallel
# ============================================================

source(file.path("/Users/dingchenchen/中国脊椎动物群落性状", "01_code", "00_config.R"))
suppressPackageStartupMessages({
  library(cluster); library(ape); library(parallel)
})

log_msg("=== 06e 两栖类性状抢救 / Amphibian trait rescue ===")
N_CORES <- max(1L, parallel::detectCores() - 2L)

# ---------------------------------------------------------------
# 1. 读取原始性状库 / Read the primary trait databases
# ---------------------------------------------------------------
tr <- readRDS(file.path(PATH$derived, "traits_imputed.rds")) |>
  filter(class == "Amphibia") |> distinct(species, .keep_all = TRUE)
comm <- readRDS(file.path(PATH$derived, "comm_Amphibia_50km.rds"))

ab <- data.table::fread(file.path(PATH$raw_trait, "amphibio", "AmphiBIO_v1.csv"))
# AmphiBIO 的 Species 列已是完整学名，不可再与 Genus 拼接
# The Species column already holds the full binomial
ab$sp <- norm_name(ab$Species)
et <- data.table::fread(file.path(PATH$raw_trait, "Etard2020_amphibians.csv"))
et$sp <- norm_name(et$Best_guess_binomial)

A <- ab[match(tr$species, ab$sp), ]
E <- et[match(tr$species, et$sp), ]

num <- function(x) suppressWarnings(as.numeric(x))
blank_na <- function(x) { x <- as.character(x); x[x == ""] <- NA; x }

# ---------------------------------------------------------------
# 2. 合成实测性状 / Assemble observed traits
# ---------------------------------------------------------------
# 体长：AmphiBIO 优先，Etard 补充 / body length, AmphiBIO first
body_length <- ifelse(is.finite(num(A$Body_size_mm)), num(A$Body_size_mm),
                      num(E$Body_length_mm))

# 垂直生境位置：由 AmphiBIO 的四个二元生境标记合成，权重与主分析一致
# Verticality from AmphiBIO's four binary habitat flags, weighted as in the
# main analysis: fossorial 0, terrestrial/aquatic 1, arboreal 3
vw <- cbind(fos = num(A$Fos), ter = num(A$Ter), aqu = num(A$Aqu), arb = num(A$Arb))
vw[!is.finite(vw)] <- NA
has_vert <- rowSums(!is.na(vw)) > 0
wsum <- rowSums(vw * rep(c(0, 1, 1, 3), each = nrow(vw)), na.rm = TRUE)
nsum <- rowSums(vw, na.rm = TRUE)
verticality <- ifelse(has_vert & nsum > 0, wsum / nsum, NA_real_)

habitat_breadth <- num(E$Habitat_breadth_IUCN)
range_size <- num(tr$range_size)   # 由 IUCN 多边形计算，属实测 / measured

obs <- data.frame(species = tr$species, body_length, verticality,
                  habitat_breadth, range_size,
                  body_mass = ifelse(is.finite(num(A$Body_mass_g)),
                                     num(A$Body_mass_g), num(E$Body_mass_g)),
                  nocturnality = ifelse(is.finite(num(A$Noc)), num(A$Noc),
                                        ifelse(blank_na(E$Diel_activity) == "Nocturnal", 1,
                                        ifelse(is.na(blank_na(E$Diel_activity)), NA, 0))))

cov_tab <- data.frame(
  trait = c("body_length", "verticality", "habitat_breadth", "range_size",
            "body_mass", "nocturnality"),
  n_observed = c(sum(is.finite(obs$body_length)), sum(is.finite(obs$verticality)),
                 sum(is.finite(obs$habitat_breadth)), sum(is.finite(obs$range_size)),
                 sum(is.finite(obs$body_mass)), sum(is.finite(obs$nocturnality))),
  n_total = nrow(obs))
cov_tab$pct <- round(100 * cov_tab$n_observed / cov_tab$n_total, 1)
cov_tab$retained <- cov_tab$trait %in% c("body_length", "verticality",
                                         "habitat_breadth", "range_size")
write_table(cov_tab, "table_17_amphibian_trait_rescue")
print(cov_tab)

RESCUE <- c("body_length", "verticality", "habitat_breadth", "range_size")
complete <- stats::complete.cases(obs[, RESCUE])
log_msg("四项高实测率性状全有的物种: ", sum(complete), " / ", nrow(obs),
        " (", round(100 * mean(complete), 1), "%)")
saveRDS(obs, file.path(PATH$derived, "amphibian_observed_traits.rds"))

# ---------------------------------------------------------------
# 3. 用该性状集重算功能离散度 SES / Recompute SES with this trait set
# ---------------------------------------------------------------
ses_for <- function(TM, comm, n_rep = 999L, min_sp = 5L) {
  sp <- rownames(TM)
  cm <- comm[comm$species %in% sp, ]
  if (!nrow(cm)) return(NULL)
  PCO <- ape::pcoa(as.dist(cluster::daisy(as.data.frame(TM), metric = "gower")))$vectors
  PCO <- PCO[, seq_len(min(4, ncol(PCO))), drop = FALSE]
  by_cell <- split(match(cm$species, sp), cm$cell_id)
  by_cell <- by_cell[lengths(by_cell) >= min_sp]
  if (length(by_cell) < 50) return(NULL)
  fdis <- function(idx) {
    P <- PCO[idx, , drop = FALSE]; ctr <- colMeans(P)
    mean(sqrt(rowSums((P - rep(ctr, each = length(idx)))^2)))
  }
  o <- vapply(by_cell, fdis, numeric(1)); S <- lengths(by_cell)
  uS <- sort(unique(S))
  nul <- parallel::mclapply(uS, function(k) {
    set.seed(SEED + k)
    v <- replicate(n_rep, fdis(sample.int(nrow(TM), k))); c(mean(v), stats::sd(v))
  }, mc.cores = N_CORES)
  names(nul) <- as.character(uS)
  ses <- vapply(seq_along(o), function(i) {
    n <- nul[[as.character(S[i])]]
    if (!is.numeric(n) || length(n) < 2 || !is.finite(n[2]) || n[2] == 0) return(NA_real_)
    (o[i] - n[1]) / n[2]
  }, numeric(1))
  list(n_species = nrow(TM), n_cells = length(o),
       mean_ses = mean(ses, na.rm = TRUE),
       pct_positive = 100 * mean(ses > 0, na.rm = TRUE),
       ses = data.frame(cell_id = names(by_cell), SES.FDis = ses))
}

build_TM <- function(df, vars, ids) {
  M <- apply(as.matrix(sapply(df[, vars], as.numeric)), 2, zscore)
  rownames(M) <- ids; M
}

scen <- list()

# 情景 A：主分析（六个核心性状，含插补）/ main analysis, imputed
CORE6 <- c("body_mass", "body_length", "nocturnality", "verticality",
           "habitat_breadth", "range_size")
tA <- tr[stats::complete.cases(tr[, CORE6]), ]
rA <- ses_for(build_TM(tA, CORE6, tA$species), comm)
scen[[1]] <- data.frame(scenario = "A. Main analysis (6 traits, imputed)",
                        n_traits = 6, n_species = rA$n_species, n_cells = rA$n_cells,
                        mean_SES = round(rA$mean_ses, 2),
                        pct_positive = round(rA$pct_positive, 1))

# 情景 B：仅五项全实测的物种（表 9 的做法）/ observed-only species subset
flags <- c("flag_mass", "flag_len", "flag_act", "flag_hab", "flag_mhab")
keep <- rep(TRUE, nrow(tr))
for (f in flags) if (f %in% names(tr))
  keep <- keep & !(suppressWarnings(as.numeric(tr[[f]])) %in% 1)
tB <- tr[keep & stats::complete.cases(tr[, CORE6]), ]
rB <- ses_for(build_TM(tB, CORE6, tB$species), comm)
scen[[2]] <- if (is.null(rB)) NULL else
  data.frame(scenario = "B. Observed-only species (6 traits)",
             n_traits = 6, n_species = rB$n_species, n_cells = rB$n_cells,
             mean_SES = round(rB$mean_ses, 2),
             pct_positive = round(rB$pct_positive, 1))

# 情景 C：高实测率性状集（本脚本的方案）/ well-measured trait set
oC <- obs[complete, ]
rC <- ses_for(build_TM(oC, RESCUE, oC$species), comm)
scen[[3]] <- if (is.null(rC)) NULL else
  data.frame(scenario = "C. Well-measured trait set (4 traits)",
             n_traits = 4, n_species = rC$n_species, n_cells = rC$n_cells,
             mean_SES = round(rC$mean_ses, 2),
             pct_positive = round(rC$pct_positive, 1))

# 情景 D：对照——同样 4 个性状，但用主分析的（含插补）数值
# Scenario D controls for the trait set itself by using the same four traits
# taken from the imputed table, isolating the effect of measurement
tD <- tr[stats::complete.cases(tr[, RESCUE]), ]
rD <- ses_for(build_TM(tD, RESCUE, tD$species), comm)
scen[[4]] <- if (is.null(rD)) NULL else
  data.frame(scenario = "D. Same 4 traits but imputed values",
             n_traits = 4, n_species = rD$n_species, n_cells = rD$n_cells,
             mean_SES = round(rD$mean_ses, 2),
             pct_positive = round(rD$pct_positive, 1))

out <- do.call(rbind, Filter(Negate(is.null), scen))
write_table(out, "table_18_amphibian_ses_comparison")
print(out)

if (!is.null(rC))
  saveRDS(rC$ses, file.path(PATH$derived, "amphibian_ses_rescued_50km.rds"))

log_msg("=== 06e 完成 / done ===")
