# Submission statements

Ready-to-paste text for the submission system. Square-bracketed fields must be filled
before submission; nothing in them has been guessed.

---

## Data availability

The aggregated data underlying every figure and table in this paper — per-cell assemblage
metrics, community-weighted trait moments, functional-axis scores, model coefficients and
all result tables — are available at
https://github.com/dingchenchen6/china-vertebrate-trait-biogeography and archived at
[Zenodo DOI to be minted from the GitHub release].

Species range polygons were obtained from the IUCN Red List of Threatened Species
(https://www.iucnredlist.org/resources/spatial-data-download) and BirdLife International
(http://datazone.birdlife.org/species/requestdis). **Their licences prohibit
redistribution**, so the species-by-cell occurrence matrices and the species-level trait
tables derived from them are not deposited. These are not available "on request" either:
they are obtainable by anyone, free of charge, directly from the two providers under their
standard terms, and the analysis scripts rebuild both files from those downloads without
manual intervention. `DATA.md` in the repository lists every source file, its version, its
access date and the script that consumes it.

Third-party data reused without modification: WorldClim 2.1 climate surfaces
(https://worldclim.org); MODIS MOD17A3 net primary productivity
(https://lpdaac.usgs.gov); ESA WorldCover land cover (https://esa-worldcover.org);
EarthEnv habitat heterogeneity (https://www.earthenv.org/texture); the Human Footprint
2000–2020; the World Database on Protected Areas (https://www.protectedplanet.net); and
the trait compilations AVONET, COMBINE, AmphiBIO, ReptTraits and TetrapodTraits, together
with the Chinese trait datasets cited in Methods §2.3. Each is cited in the reference list
and in `DATA.md`.

Source data for each main-text figure are provided in the repository under
`04_results/tables/`, indexed by `table_00_figure_manifest.csv`; the Supplementary
Tables workbook reproduces, in order of citation, every table the manuscript cites.

## Code availability

All analysis code is available at
https://github.com/dingchenchen6/china-vertebrate-trait-biogeography and archived at
[Zenodo DOI]. The pipeline is written in R 4.5.1 and runs end to end from the raw public
sources to every figure and table. The random seed is fixed at 20260801; each script
carries a bilingual header stating its scientific question, inputs, workflow, outputs and
key assumptions. Principal packages: sf, terra, spdep, spatialreg, vegan, FD, ade4, lme4,
missForest, Matrix and cluster.

## Author contributions

C.D. designed the study, assembled the data, performed all analyses, produced the figures
and wrote the manuscript.

## Competing interests

The author declares no competing interests.

## Funding

[To be completed — funder and grant numbers, or an explicit statement that the work
received no specific funding.]

## Ethics

This study analysed previously published distribution and trait data and generated no new
records from live animals. No ethical approval was required.

---

## 中文核对 — 提交前必须确认的字段

| 字段 | 状态 | 说明 |
|---|---|---|
| Zenodo DOI | **待补** | 从 GitHub release 铸造后，填入 Data availability、Code availability 与投稿信三处 |
| 资助信息 | **待补** | 若无专项资助，需明确写"未获专项资助"，不能留空 |
| ORCID | **待补** | 投稿系统必填 |
| R 版本号 | 已填 4.5.1 | 若投稿前更新过环境，请核对 |
| IUCN / BirdLife 数据访问日期 | 见 `DATA.md` | 审稿人可能要求在正文注明具体版本年份 |

**为什么不用"可向作者索取"**：IUCN 与 BirdLife 的数据任何人都能免费直接获取，
不存在只有作者才能提供的情形。Nature 系列把"available on request"视为弱声明，
仅在存在法律、伦理、商业或第三方限制且说明了审批路径时才接受。此处属于
**第三方再分发限制**，正确的写法是指明原始获取途径并保证脚本可自动重建，
这正是上文的处理方式。
