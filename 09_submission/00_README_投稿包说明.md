# 投稿包 / Submission package

**Thermoregulatory mode sets the strength of environmental filtering on vertebrate
assemblage traits, and leaves a directional conservation gap**

Chenchen Ding · School of Life Sciences, Peking University

目标期刊为 *Nature Ecology & Evolution*。与该刊格式的差距（正文字数、展示项数、
文献格式、章节顺序）及改投 *Nature Communications* 的比较，见
[`07_submission_checklist.md`](07_submission_checklist.md) 第三节。

---

## 目录

| 路径 | 内容 | 投稿系统里传哪一项 |
|---|---|---|
| `01_cover_letter.md` | 投稿信 | Cover letter |
| `02_manuscript/manuscript.docx` | 正文（题名、198 词摘要、引言、方法、结果、讨论、声明、图注、参考文献） | Main manuscript |
| `02_manuscript/manuscript.md` | 同上，Markdown 源文件 | — |
| `02_manuscript/figure_legends.md` | 8 张正文图 + 25 张扩展数据图的图注 | Figure legends（多数期刊要求单独上传） |
| `02_manuscript/abstract_long_version.md` | 原约 360 词的长摘要，留作参考 | — |
| `03_main_figures/Fig1–Fig8` | 正文图，PDF（矢量）+ PNG（600 dpi） | Figures |
| `04_extended_data/ED1–ED25` | 扩展数据图，同上两种格式 | Extended Data |
| `05_supplementary/Supplementary_Tables.xlsx` | 补充表 S1–S28 + 目录页 | Supplementary Information |
| `06_statements/statements.md` | 数据/代码可得性、作者贡献、利益冲突、资助、伦理 | 各声明栏 |
| `07_submission_checklist.md` | 提交前逐项核对 | — |

## 编号规则与复现

图、扩展数据图、补充表**全部按正文首次引用顺序编号**（Nature 系列要求）。编号由脚本
生成，改动正文后重跑即可，不必手工调整：

```bash
python3 09_submission/renumber_apply.py                 # 按引用顺序重编号手稿与图注
Rscript --vanilla 01_code/95_assemble_figures.R          # 按新编号生成清单与图件
python3 09_submission/build_supplementary_tables.py      # 按新编号构建补充表册
Rscript --vanilla 09_submission/assemble_figures.R       # 装配进投稿包
```

`display_numbering.json` 记录了分析内部表号/图号到投稿编号的映射。
投稿包内的图件是 `05_figures/` 的重命名副本，已 gitignore，由最后一步生成。

## 图件规格

- 33 张图全部 600 dpi、183 mm 宽（Nature 双栏），PDF 为嵌入字体的矢量文件
- 字体 Helvetica；**色相编码类群**（鸟红、兽橙、爬行蓝、两栖绿），
  **明度编码体温调节模式**——同一颜色不承担两个含义
- 连续量用 `batlow`、发散量用 `vik`（Crameri 感知均匀、色觉障碍友好）
- 地图为中国 Albers 等积圆锥投影，南海诸岛按真实地理位置绘制，十段线为独立图层

## 提交前必须补的四项

1. **Zenodo DOI** — 从 GitHub release 铸造，填入投稿信与两处可得性声明
2. **资助信息** — 若无专项资助也要明确写出
3. **ORCID**
4. **推荐/回避审稿人** — 见投稿信末尾的说明

---

## 本轮检查发现并修正的问题

**交叉引用（上一版投稿包的硬伤）**
- 手稿引用的 18 个补充表号中 17 个在表册里不存在——手稿用的是分析内部表号，表册却是另一套编号
- 上一版表册取自 8 月 3 日生成的格式化表，早于零模型修正，**内容过时**；现改为直接从当前结果表构建
- 25 张扩展数据图中 16 张正文从未引用；Fig. 1 从未引用；Fig. 6 只被引用过一次且引错了位置
- 手稿引用了不存在的"Supplementary Table S0"（预先声明的假设表）——已补建为 S1
- 上一轮我新增的体型稳健性表沿用了另一套编号，与其余引用冲突——已统一

**事实与口径**
- OLS 残差 Moran's *I*：五处写的"高达 0.69"无任何保存的出处；补算得 0.65–0.89，新增表 S6
- 体型量度：手稿 §2.3、申请书与三问文档写成"变温类群用体长"，而代码是**只有两栖用体长**
- 图注 5 处错误（其中 Fig7 两栖类缺席的理由是我上一轮编造的，已换成可核实的实测覆盖率）
- 投稿信"三条结论被撤回且论文如实说明"言过其实，改为点名手稿中实际写明的两处更正
- 投稿信称"改投只需换前置信息"不实，已改为如实列出版式差距
- 手稿 §3.4 与 §4.2 提议的后续设计是"匹配生境异质性"，与申请书中论证的设计相矛盾（异质性沿海拔本应变化，只能作协变量），已统一

**格式**
- 摘要由约 360 词压缩为 198 词（NEE 与 *Nature Communications* 上限均为 200 词），长版另存
- 页眉字数由"约 8,000"改为实测值

**未改动**：25 张扩展数据图全部保留（作者决定）；正文字数、图数、文献格式与章节顺序
等 NEE 版式转换，待确定目标期刊后一次完成。
