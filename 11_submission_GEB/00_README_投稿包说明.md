# 投稿包 — *Global Ecology and Biogeography*

**Thermoregulatory mode sets the strength of environmental filtering on vertebrate assemblage
traits, and leaves a directional conservation gap**

## 为什么推荐这本刊

- **题材正中该刊核心**：性状生物地理、把经典多样性假说用于性状组成、群落加权均值的推断问题
  （第四角检验、零模型）——都是 *GEB* 读者日常面对的问题。
- **章节结构与主稿一致**（方法在结果之前、作者–年份制文献），三套编号无需改动，已核验。
- **首投支持自由格式**，但结构化摘要是硬性要求，已按要求改写。
- **风险**：篇幅超过该刊参考值（见下表）。投稿信已如实说明理由，并表示可按编辑意见把部分内容移入支撑信息。

## 与该刊参考要求的对照

来源：*GEB* 作者指南（Wiley 页面被 Cloudflare 拦截无法直接抓取，以下要点取自检索引擎收录的官方页面原文）。

| 项目 | *GEB* 要求 | 当前 | 状态 |
|---|---|---|---|
| 结构化摘要 | ≤ 300 词；Aim / Location / Time period / Major taxa studied / Methods / Results / Main conclusions | 297 词，七个小标题齐全 | ✅ |
| 正文 | "roughly 5000 words in the main body" | 6,477 词 | ⬜ 超出，为参考值 |
| 参考文献 | "50 literature citations" | 100 | ⬜ 超出，为参考值 |
| 图表 | "six to eight display pieces" | 8 | ✅ |
| 弹性条款 | "Papers that are shorter in one of these respects may be longer in another" | — | 投稿信已据此说明 |

## 相对主稿的改动（全部由脚本生成）

- 摘要改为 297 词的结构化摘要
- 去掉小节编号；图表改用 Wiley 惯例：Figure 1–8、Figure S1–S25、Table S1–S28
- 支撑信息分为 Appendix S1（图）与 Appendix S2（表）
- 增加 Data availability statement、Conflict of interest statement 与 "Use of AI technology" 占位节
- 补充表 S1（假设表）的"证据"列改为 Wiley 称谓

## 目录

| 路径 | 内容 |
|---|---|
| `01_cover_letter.md` / `.docx` | 投稿信 |
| `02_manuscript/manuscript.docx` | 正文（含结构化摘要与图注） |
| `03_figures/Figure1–Figure8` | 正文图，PDF 矢量 + PNG 600 dpi |
| `04_supporting_information/Appendix_S1_Figures_S1-S25.docx` | 支撑图 S1–S25（嵌入约 300 dpi 预览图）及图注 |
| `04_supporting_information/figures/` | 支撑图原图，PDF + PNG 600 dpi |
| `04_supporting_information/Appendix_S2_Tables_S1-S28.xlsx` | 支撑表 S1–S28 + 目录页 |
| `05_statements/AI_use_disclosure.md` | AI 使用披露：Wiley 政策原文、使用记录、参考草稿 |

## 提交前你必须完成

1. **AI 使用披露**：在手稿 "Use of AI technology" 一节写入定稿（Wiley 要求写明用途、是否影响核心论点、你本人如何复核），
   并核实所用 AI 服务不会获得稿件内容的权利（含训练）
2. Zenodo DOI、资助信息、ORCID、致谢、推荐审稿人
3. 在投稿系统中核对本包未能核实的项目：biosketch、关键词数量、running title 长度

## 重新生成

```bash
python3 09_submission/build_journal_packages.py
```
