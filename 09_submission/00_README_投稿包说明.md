# 投稿包 / Submission package

**Thermoregulatory mode sets the strength of environmental filtering on vertebrate
assemblage traits, and leaves a directional conservation gap**

Chenchen Ding · School of Life Sciences, Peking University

按 Nature 系列格式组织（Extended Data + Supplementary Information）。
若改投 *Nature Communications*、*Ecology Letters* 或 *Global Ecology and Biogeography*，
正文与图件不必改动，只需替换投稿信与前置信息。

---

## 目录

| 路径 | 内容 | 投稿系统里传哪一项 |
|---|---|---|
| `01_cover_letter.md` | 投稿信 | Cover letter |
| `02_manuscript/manuscript.docx` | 正文（题名、摘要、引言、方法、结果、讨论、参考文献） | Main manuscript |
| `02_manuscript/manuscript.md` | 同上，Markdown 源文件 | — |
| `02_manuscript/figure_legends.md` | 8 张正文图 + 25 张扩展数据图的完整图注 | Figure legends（多数期刊要求单独上传） |
| `03_main_figures/Fig1–Fig8` | 正文图，PDF（矢量）+ PNG（600 dpi） | Figures |
| `04_extended_data/ED1–ED25` | 扩展数据图，同上两种格式 | Extended Data |
| `05_supplementary/Supplementary_Tables.xlsx` | 22 张补充表 + 目录页 | Supplementary Information |
| `06_statements/statements.md` | 数据/代码可得性、作者贡献、利益冲突、资助、伦理 | 各声明栏 |
| `07_submission_checklist.md` | 提交前逐项核对 | — |

---

## 图件规格

- 33 张图，**全部 600 dpi、183 mm 宽**（Nature 双栏），PDF 为带嵌入字体的矢量文件
- 字体 Helvetica；**色相编码类群**（鸟红、兽橙、爬行蓝、两栖绿），
  **明度编码体温调节模式**——同一颜色不承担两个含义
- 连续量用 `batlow`、发散量用 `vik`（Crameri 感知均匀、色觉障碍友好）
- 地图为中国 Albers 等积圆锥投影，南海诸岛按真实地理位置绘制，十段线为独立图层
- 编号对应关系见 `04_results/tables/table_00_figure_manifest.csv`

## 可复现性

全部代码在 https://github.com/dingchenchen6/china-vertebrate-trait-biogeography ，
随机种子固定为 20260801。IUCN 与 BirdLife 的分布数据禁止再分发，
因此物种×网格矩阵与物种级性状表不入库；脚本可从两家机构的公开下载重建这两份文件。
详见 `06_statements/statements.md`。

---

## 提交前必须补的四项

1. **Zenodo DOI** — 从 GitHub release 铸造，填入投稿信与两处可得性声明
2. **资助信息** — 若无专项资助也要明确写出
3. **ORCID**
4. **推荐/回避审稿人** — 见投稿信末尾的说明

---

## 本版相对前一版改了什么

- 科级热学对比新增**体型量度稳健性检验**：两栖类改用体长后比值 1.72 → 1.67，
  结论不变（TableS19）
- Fig3 体型行改为**每纲一个量度**（两栖用体长），行标题由 "CWM body mass" 改为
  "CWM body size"
- Fig6 面板 b 的体温调节模式配色改为明度对比，解除"红=鸟类/红=恒温"的一色两义
- Fig8、ED22、ED24 的面板副标题过长导致**跨面板重叠**，已换行修复
- 数据可得性声明按 Nature 政策重写：不再使用"可向作者索取"，改为指明
  第三方再分发限制与原始获取路径
- 图注全部重写：旧的 `06_manuscript/figure_captions.md` 是零模型修正之前的版本，
  编号与数值均已过时，**不可用于投稿**
