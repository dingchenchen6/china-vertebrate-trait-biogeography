# 数据获取与再分发说明 / Data availability and redistribution

本项目的**全部代码、文档、图件与聚合结果表**均在本仓库中公开。
部分原始数据受许可限制或体积过大，不随仓库分发；本文件给出每一项的确切来源与获取方式，
使整条流程可完全复现。

All code, documentation, figures and aggregated result tables are in this
repository. Some source data are licence-restricted or too large to ship; every
item is listed below with its exact origin so the pipeline can be reproduced in full.

---

## 1. 不随仓库分发的数据 / Not redistributed

### 1.1 物种分布 / Species distributions

| 数据 | 来源 | 为何不分发 | 如何获取 |
|---|---|---|---|
| 两栖、爬行、哺乳类范围多边形 | IUCN Red List Spatial Data v2025-1 | **IUCN 使用条款禁止再分发**。物种×网格存在矩阵是范围图的重新编排而非聚合摘要，同样受限 | 在 [iucnredlist.org/resources/spatial-data-download](https://www.iucnredlist.org/resources/spatial-data-download) 免费注册后下载 |
| 鸟类网格存在矩阵（1,180 种 × 3,888 格） | Huang et al. 2023, *The Innovation* 4:100379，其公开仓库 [Huangmp1996/Paleoclimate-and-regionalization-of-China](https://github.com/Huangmp1996/Paleoclimate-and-regionalization-of-China) | 源自 BirdLife/IUCN 范围图，受同一条款约束 | 见上述仓库 `ave/parafly/comm_data.RData` |
| 中国蜥蜴分布多边形（213 种） | 同上仓库 `reptile/data_reptile/Chinese213lizards_9_27/` | 第三方数据 | 见上述仓库 |

> **本仓库包含的是逐格聚合指标**（性状矩、功能多样性、零模型 SES 等），
> 其中不含任何物种身份信息，因此不构成对范围图的再分发。

### 1.2 性状数据库 / Trait databases

全部为开放获取，脚本 `01_code/03_traits.R` 中给出了直接下载地址：

| 数据库 | 引用 | 获取 |
|---|---|---|
| **TetrapodTraits v2.0.1**（主干） | Moura et al. 2024 *PLoS Biology* | Zenodo `10.5281/zenodo.18926700` |
| AVONET | Tobias et al. 2022 *Ecology Letters* | figshare `10.6084/m9.figshare.16586228` |
| COMBINE | Soria et al. 2021 *Ecology* | figshare `10.6084/m9.figshare.13028255` |
| AmphiBIO v1 | Oliveira et al. 2017 *Scientific Data* | figshare `10.6084/m9.figshare.4644424` |
| ReptTraits | Oskyrko et al. 2024 *Scientific Data* | figshare `10.6084/m9.figshare.24572683` |
| Etard et al. 2020 | Etard, Morrill & Newbold 2020 *GEB* | figshare `10.6084/m9.figshare.10075421` |
| Meiri 2018 蜥蜴性状 | Meiri 2018 *GEB* | Zenodo `10.5281/zenodo.4946186` |
| EltonTraits 1.0 | Wilman et al. 2014 *Ecology* | figshare `3559887` |
| Amniote Life History DB | Myhrvold et al. 2015 *Ecology* | figshare `3563457` |
| AnAge / HAGR | Tacutu et al. 2018 *NAR* | genomics.senescence.info |
| GlobTherm | Bennett et al. 2018 *Scientific Data* | Zenodo `10.5281/zenodo.4976423` |

### 1.3 环境数据 / Environmental data

| 数据 | 来源 |
|---|---|
| WorldClim 2.1 生物气候变量（2.5′）与高程（30″） | worldclim.org |
| MODIS 净初级生产力 | 由用户既有项目提供，源自 MOD17A3 |
| ESA WorldCover 土地覆被占比 | esa-worldcover.org |
| EarthEnv 生境异质性纹理 | earthenv.org |
| 人类足迹指数 2020 | Venter et al. / Mu et al. 更新版 |
| 中国土地覆被数据集 CLCD | Yang & Huang 2021 *ESSD* |

### 1.4 系统发育树 / Phylogenies

鸟类 Jetz et al. 2012；哺乳类 Upham et al. 2019；有鳞目 Tonini et al. 2016；
两栖类 Jetz & Pyron 2018。经 vertlife.org 与相应公开仓库获取。

### 1.5 保护地 / Protected areas

全国自然保护区名录与矢量边界（1,028 个，2021 版）。
**注意**：该数据仅含自然保护区，不含国家公园、风景名胜区、森林公园等其他保护地类型，
因此本研究报告的覆盖率（10.31%）是中国保护地体系的**下界**。

### 1.6 观测记录 / Occurrence records

GBIF。鸟类使用 2000–2025 年中国下载集（458 万条）；
哺乳（70,450）、两栖（102,386）、爬行（35,446）由 `01_code/11a_gbif_download.py` 经 GBIF API 获取，
该脚本可直接重跑。

---

## 2. 随仓库分发的数据 / Included in this repository

| 路径 | 内容 | 体积 |
|---|---|---|
| `01_code/` | 全部分析代码（23 个脚本，中英双语注释） | 364 KB |
| `00_docs/` | 研究设计、结果报告、本文件 | 60 KB |
| `06_manuscript/` | 论文草稿、图注、表注 | 52 KB |
| `04_results/tables/` | 45 张分析输出表 | 1.0 MB |
| `04_results/tables_formatted/` | 21 张可直接排版的表 + Excel 汇编 | 152 KB |
| `05_figures/` | 正文图 14 幅 + 扩展数据图 9 幅（矢量 PDF + 600 dpi PNG） | 62 MB |
| `03_data_derived/grid_*.gpkg` | 50 km 与 100 km 等积网格（含 `scs_inset` 标记） | 1.4 MB |
| `03_data_derived/metrics_*.rds` | 逐格性状矩与功能多样性指标（**无物种身份**） | 7 MB |
| `03_data_derived/model_data.rds` | 建模用表：环境 + 指标 | 8 MB |
| `03_data_derived/pa_coverage_50km.rds` | 逐格保护区覆盖度 | < 1 MB |
| `03_data_derived/cn_*.gpkg` | 制图底图（省界、国界、十段线） | 5 MB |

---

## 3. 复现路径 / Reproduction

1. 按 §1 获取受限数据，放入 `02_data_raw/` 下对应子目录（路径见 `01_code/00_config.R` 的 `EXT` 列表）。
2. 按 `README.md` 中的顺序执行脚本。所有随机化均由 `SEED = 20260801` 固定。
3. 若只需复现**统计与图件**（不重建群落矩阵），可直接从 `03_data_derived/model_data.rds`
   与 `metrics_*.rds` 出发运行 `06_models.R` 之后的脚本。

---

## 4. 制图合规 / Cartographic compliance

分析网格建立在自然资源部标准地图 **GS(2023)2767** 的大陆边界上。
制图底图采用南海诸岛按真实地理位置绘制的省级矢量边界配合独立十段线图层。
两套边界在大陆部分一致，差异仅在于南海诸岛的摆放方式（附图框 vs 真实位置）。

> GS(2023)2767 的省面图层把南海诸岛预置在附图框内，据此建网格会在该处生成
> 32 个（50 km）伪网格，其反投影经纬度落在西太平洋。这些网格已用 `scs_inset`
> 字段标记并在全部分析与制图中排除，详见 `01_code/01_build_grid.R`。
