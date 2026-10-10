# 投稿包 — *Nature Communications*

**Thermoregulation sets how strongly climate filters vertebrate assemblages, and reserves leave
a directional functional gap**

## 为什么推荐这本刊

- **篇幅与现稿吻合，不用删减内容**：正文 4,288 词（官网"理想上 ≤ 5,000"）、摘要 198 词（≤ 200）、
  8 个展示项（≤ 10）。三项都已对照官网核实（2026-09/10）。
- **首投格式自由**：官网原文 "flexible with regard to the format of initial submissions … We also do not
  require a particular structure or format at first submission"，文献编号制等排版细节可留到修回阶段。
- **读者面与论文的三条贡献对得上**：方法学结论（合成轴稀释信号、朴素检验膨胀显著性）面向所有用性状数据的人；
  四纲对比面向宏观生态学；对昆明–蒙特利尔目标 3 的保护地审计面向保护政策。
- **风险**：综合性期刊对"普遍意义"的门槛高；H2 只是"部分支持"，编辑可能认为机制证据不够硬。

## 相对主稿的改动（全部由脚本生成）

- 题名改为 15 词（官网上限 15 词）；原题名 18 词
- 章节改为 Introduction → Results → Discussion → **Methods（移到文末）**，去掉小节编号
- "扩展数据图"改称"补充图"（该刊没有扩展数据图这一类别）
- 方法移到文末改变了引用顺序，**25 张补充图中 24 张、28 张补充表中 26 张重新编号**；三套编号均已核验为按首次引用升序
- 补充表 1（假设表）的"证据"列随新编号改写
- Methods 末尾加入 "Use of AI tools" 小节草稿（Nature 政策要求在 Methods 中披露）
- 拆分出独立的 Data availability 与 Code availability 两节

## 目录

| 路径 | 内容 |
|---|---|
| `01_cover_letter.md` / `.docx` | 投稿信 |
| `02_manuscript/manuscript.docx` | 正文（含正文图注） |
| `03_figures/Fig1–Fig8` | 正文图，PDF 矢量 + PNG 600 dpi |
| `04_supplementary/Supplementary_Information.docx` | 补充信息：补充图 1–25（嵌入约 300 dpi 预览图）及图注 |
| `04_supplementary/figures/` | 补充图原图，PDF + PNG 600 dpi（编辑部要求高清时用） |
| `04_supplementary/Supplementary_Tables_1-28.xlsx` | 补充表 1–28 + 目录页 |
| `05_statements/AI_use_disclosure.md` | AI 使用披露：政策、使用记录、待确认事项 |
| `numbering_from_master.json` | 主稿编号 → 本刊编号的映射 |

## 提交前你必须完成

1. **AI 使用披露**：逐项确认 `05_statements/AI_use_disclosure.md`，改定 Methods 中方括号内容并删去提示
2. **题名**：确认 15 词版本，或另拟（不超过 15 词）
3. Zenodo DOI、资助信息、ORCID、致谢、推荐审稿人
4. 投稿系统中的 Reporting Summary 等表单（本包未核实其是否在首投阶段要求）

## 重新生成

```bash
python3 09_submission/build_journal_packages.py
```

主稿（`06_manuscript/manuscript_v4_full.md`）改动后重跑即可，编号会自动重算。
