# 数据来源与可获取性 / Data sources and availability

本仓库公开**全部代码、文档、图件与聚合结果**。
原始数据与物种级中间产物因体积或许可限制不入库，本文件给出逐项获取方式，使分析可完整复现。

This repository contains all code, documentation, figures and aggregated results.
Raw data and species-level intermediates are excluded for size or licensing
reasons; every item is documented here so the analysis can be reproduced in full.

---

## 1. 不可再分发的数据 / Data that may not be redistributed

| 数据 | 来源 | 限制 | 获取方式 |
|---|---|---|---|
| 两栖、爬行、哺乳类分布范围多边形 | IUCN Red List Spatial Data v2025-1 | IUCN 使用条款禁止再分发 | 在 [iucnredlist.org/resources/spatial-data-download](https://www.iucnredlist.org/resources/spatial-data-download) 注册后免费下载 |
| 鸟类网格存在矩阵 | Huang et al. 2023, *The Innovation* 4:100379 配套仓库 | 源自 BirdLife/IUCN 范围图，同受上述条款约束 | [github.com/Huangmp1996/Paleoclimate-and-regionalization-of-China](https://github.com/Huangmp1996/Paleoclimate-and-regionalization-of-China) |
| 中国自然保护区矢量边界 | 全国自然保护区名录 + 矢量边界数据集 | 需向数据提供方申请 | 见 `01_code/15_conservation.R` 头部注释 |

> 由这些数据派生的**物种 × 网格存在矩阵**（`comm_*.rds`）与投影后的范围多边形
> （`ranges_albers_*.gpkg`）同样属于再分发限制范围，已在 `.gitignore` 中排除。
> 仓库中保留的是**逐格聚合指标**（性状矩、功能多样性、SES），不含物种身份，
> 属于常规研究产出。

---

## 2. 开放数据：可由代码自动获取 / Open data fetched by the code

| 数据 | 版本 | 用途 | 获取 |
|---|---|---|---|
| **TetrapodTraits** v2.0.1 | Moura et al. 2024 *PLoS Biol* | 四类群统一性状骨架（33,281 种） | Zenodo 10.5281/zenodo.18926700 |
| **AVONET** | Tobias et al. 2022 *Ecol Lett* | 鸟类形态与生态性状 | figshare 10.6084/m9.figshare.16586228 |
| **COMBINE** | Soria et al. 2021 *Ecology* | 哺乳类生活史 | figshare 10.6084/m9.figshare.13028255 |
| **AmphiBIO** v1 | Oliveira et al. 2017 *Sci Data* | 两栖类生态性状（本项目性状抢救的主要来源） | figshare 10.6084/m9.figshare.4644424 |
| **ReptTraits** | Oskyrko et al. 2024 *Sci Data* | 爬行类性状 | figshare 10.6084/m9.figshare.24572683 |
| **Etard et al. 2020** | *Glob Ecol Biogeogr* | IUCN 生境宽度（四类群，CC0） | figshare 10.6084/m9.figshare.10075421 |
| **Meiri 2018** | *Glob Ecol Biogeogr* | 蜥蜴体温与体型 | Zenodo record 4946186 |
| **GlobTherm** | Bennett et al. 2018 *Sci Data* | 热耐受 | Zenodo record 4976423 |
| **Amniote LHDB** | Myhrvold et al. 2015 *Ecology* | 生活史补充 | figshare 3563457 |
| **AnAge / HAGR** | Tacutu et al. 2018 *NAR* | 最大寿命 | genomics.senescence.info |
| **WorldClim 2.1** | Fick & Hijmans 2017 | 19 个生物气候变量（2.5′）+ 高程（30″） | geodata.ucdavis.edu |
| **ESA WorldCover** | v100 (2020) | 土地覆被占比 | esa-worldcover.org |
| **EarthEnv** | Tuanmu & Jetz 2015 | 生境异质性纹理 | earthenv.org |
| **MODIS NPP** | MOD17A3 | 净初级生产力 | LP DAAC |
| **Human Footprint** | Mu et al. 2022 | 人类压力（2000–2020） | figshare |
| **GBIF 出现记录** | 2000–2025，中国 | 观测群落（鸟 458 万、两栖 10.2 万、兽 7.0 万、爬行 3.5 万） | 由 `01_code/11a_gbif_download.py` 经 GBIF API 自动下载 |

系统发育树：鸟 Jetz et al. 2012；兽 Upham et al. 2019；有鳞目 Tonini et al. 2016；
两栖 Jetz & Pyron 2018。均可从 [vertlife.org](https://vertlife.org) 获取。

---

## 3. 制图底图 / Cartographic basemap

所有中国地图基于自然资源部标准地图服务的 **GS(2023)2767** 号标准地图，
国界、省界与南海诸岛按官方底图原样呈现，未作任何修改。
底图本身不随仓库分发，可从 [bzdt.ch.mnr.gov.cn](http://bzdt.ch.mnr.gov.cn) 下载。

> **一处必须注意的技术细节**：该底图把南海诸岛**预置在附图框内**而非其真实
> 地理位置。若直接在其上生成规则网格，附图框内会产生 32 个（50 km）伪网格，
> 其反投影经纬度落在西太平洋，从而把正确的物种组成与错误位置的气候配对。
> `01_code/01_build_grid.R` 用 `scs_inset` 标记并在所有分析与制图中排除这些网格。

---

## 4. 仓库中包含的衍生数据 / Derived data included here

| 文件 | 内容 | 可否公开 |
|---|---|---|
| `03_data_derived/grid_{50,100}km.gpkg` | 等积分析网格及属性 | ✓ 本项目构建 |
| `03_data_derived/metrics_<class>_<grain>.rds` | 逐格性状矩、功能多样性、双零模型 SES | ✓ 聚合量，不含物种身份 |
| `03_data_derived/model_data.rds` | 建模用表（环境 + 指标） | ✓ |
| `03_data_derived/env_{50,100}km.rds` | 47 个环境变量的逐格汇总 | ✓ 源自开放数据 |
| `04_results/tables/` | 全部分析输出表 | ✓ |
| `04_results/tables_formatted/` + `.xlsx` | 可直接排版的正文表与附表 | ✓ |
| `05_figures/` | 全部图件（PDF 矢量 + 600 dpi PNG） | ✓ |

---

## 5. 复现顺序 / Reproduction order

```bash
Rscript 01_code/01_build_grid.R           # 等积网格
Rscript 01_code/02_species_pools.R        # 鸟类矩阵重网格化
Rscript 01_code/02b_iucn_to_grid.R        # IUCN 范围图 -> 网格（需自备 IUCN 数据）
Rscript 01_code/03_traits.R               # 性状整合与名录协调
Rscript 01_code/03b_trait_imputation.R    # 系统发育辅助插补
Rscript 01_code/04_environment.R          # 环境变量提取
Rscript 01_code/05_community_metrics.R    # 性状矩 + 功能多样性 + 双零模型
Rscript 01_code/06_models.R               # 逐类群空间自回归
Rscript 01_code/06b_fourthcorner.R        # 第四角 max test
Rscript 01_code/06c_clade_level.R         # 科级重复 + 核心检验
Rscript 01_code/06d_robustness.R          # 留一纲 + 插补敏感性
Rscript 01_code/06e_amphibian_rescue.R    # 两栖类实测性状重做
python3 01_code/11a_gbif_download.py      # GBIF 记录下载
Rscript 01_code/13_occurrence_assemblages.R
Rscript 01_code/14_occurrence_all_classes.R
Rscript 01_code/10_novel_analyses.R       # 探索分析
Rscript 01_code/12_functional_beta.R      # 功能 β 多样性
Rscript 01_code/15_conservation.R         # 保护空缺与有效性
Rscript 01_code/07_figures.R              # 正文图
Rscript 01_code/08_extended_figures.R     # 扩展数据图
Rscript 01_code/09_tables.R               # 规范化表格
```

环境：R 4.5.1；随机种子 `SEED = 20260801`（见 `01_code/00_config.R`）。
中间结果自动缓存于 `04_results/cache/`，可分步重跑。
