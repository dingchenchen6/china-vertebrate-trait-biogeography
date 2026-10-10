# AI 使用披露 — *Global Ecology and Biogeography*（Wiley）

**状态：须作者撰写定稿后方可提交。**

按 ARS `disclosure` 模式，Wiley / *GEB* 不在其已收录的期刊政策库中，状态为
`UNKNOWN / HALTED / UNCURATED_POLICY`——即**不替你生成"已按期刊对齐"的披露文本**。
下面给出：①本轮直接抓取的 Wiley 官方政策原文要点；②预填的使用记录；③一份供你参考改写的草稿。
手稿 "Use of AI technology" 一节目前是占位符，需要你写入最终版本。

## 一、Wiley 官方政策要点

来源：Wiley *Best Practice Guidelines on Publishing Ethics*
（https://authorservices.wiley.com/ethics-guidelines/index.html ，2026-10-10 抓取）。逐条摘录：

- "Authors should document all AI Technology used, including its purpose, whether it influenced
  key arguments or conclusions, and how they personally reviewed and verified any AI-generated
  content. Authors must also disclose their use of AI Technologies when submitting a manuscript
  to a Wiley journal."
- "If an author has used AI Technology to substantially edit, develop, or translate any part of
  a manuscript, its use must be described transparently and integrated into the manuscript"
- "AI Tools used solely for spelling, grammar, and general editing are not included in these
  disclosure requirements."
- "AI Technology may be used for creating data visualizations and illustrations, however,
  important restrictions apply to photograph enhancement"
- "Authors must ensure the AI Technology and the provider of that AI Technology does not gain any
  rights over the author's underlying content, including the right to 'train' their AI
  Technology on the content"
- AI 不能署名（依据 COPE 关于作者身份与 AI 工具的立场声明）

**Wiley 比 Nature 多出的一项：** 你需要确认所用 AI 服务的条款与账户设置**不允许服务商用你的稿件训练模型**。
这取决于你的账户类型与数据设置，我无法替你判断，请在 Claude 的隐私/数据设置中核实。

## 二、使用记录（与 *Nature Communications* 包相同，逐项确认）

工具：Claude（Anthropic），经 Claude Code 使用，包括 Opus 5 与 Opus 5.5 模型；时间 2026 年 8–10 月。

| 类别 | 状态 | AI 实际做了什么 | 影响部分 | 是否影响核心论点或结论 |
|---|---|---|---|---|
| 分析 | 已使用（属于研究本身） | 编写并运行 R 分析流程 | 方法、结果、全部图表 | ⬜ 请你判断并写明 |
| 数据可视化 | 已使用 | 编写绘图 R 代码；无生成式图像 | 全部图件 | — |
| 撰写 | 已使用 | 起草手稿、图注、支撑信息、投稿信 | 全文 | ⬜ |
| 修订与一致性核查 | 已使用 | 数字溯源、编号、更正 | 全文 | — |
| 文献核查与检索 | 已使用 | Crossref 核实文献；检索数据集与期刊要求 | 参考文献、方法 | — |
| 文字编辑 | 已使用 | 润色、格式转换 | 全文 | — |
| 本项目之外的其他 AI 工具 | **未知** | — | — | ⬜ **必须由你声明** |

Wiley 特别要求写明**你本人如何复核与验证 AI 生成的内容**。从本项目记录看，可以如实写入的复核手段包括：
所有数字回溯到结果表的逐项核对、参考文献逐条对照 Crossref、关键分析的稳健性重跑（零模型、阈值、体型量度）——
但"由谁复核"必须是你本人做过的，请据实写。

## 三、供参考的草稿（未按期刊对齐，须你改写）

> **Use of AI technology.** An AI assistant (Claude; Anthropic), including the Opus 5 and
> Opus 5.5 models accessed through Claude Code between August and October 2026, was used to
> write and run the R analysis code and produce the figures, to draft and revise the manuscript
> text and supporting information, and to check references and reported values. [State whether
> this use influenced key arguments or conclusions.] [State how you reviewed and verified the
> AI-generated code, results and text.] [Any other AI tools used, or "None".] The author takes
> full responsibility for the content of this paper.

## 四、提交前你需要确认的四件事

1. 使用记录逐行属实且完整
2. 写明 AI 使用是否影响了核心论点或结论
3. 写明你本人的复核方式
4. 核实所用 AI 服务不会获得稿件内容的权利（含训练）
