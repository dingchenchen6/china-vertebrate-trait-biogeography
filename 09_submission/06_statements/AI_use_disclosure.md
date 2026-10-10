# AI 使用披露 — *Nature Ecology & Evolution*

**状态：须作者确认后方可提交。** 本文件按 ARS `disclosure` 模式整理：使用记录由本项目的
工作记录预填，每一项都要你本人核对；你在本项目之外是否用过其他 AI 工具，只有你能声明。

## 一、期刊政策（Nature Portfolio，适用于 *Nature Ecology & Evolution*）

依据 `shared/policy_data/nature_policy.md`（Nature Portfolio AI 编辑政策页的存档摘录）：

> "Use of an LLM should be properly documented in the Methods section … of the manuscript."
>
> "In all cases, there must be human accountability for the final version of the text and
> agreement from the authors that the edits reflect their original work."

- 披露位置：**Methods**（已在手稿方法部分 2.14 节放入 "Use of AI tools" 草稿）
- LLM 不能署名
- 仅用于"AI 辅助文字编辑"（对人写文本的可读性与文风修改）可免于申报——**本项目的使用远超此范围，必须申报**
- 生成式 AI 图像一般不允许发表；本文所有图均为 R 代码绘制的数据图，不涉及生成式图像

## 二、使用记录（由本项目工作记录预填，逐项确认）

工具：Claude（Anthropic），经 Claude Code 使用，包括 Opus 5 与 Opus 5.5 模型；时间 2026 年 8–10 月。

| 类别 | 状态 | AI 实际做了什么 | 影响部分 | 你的确认 |
|---|---|---|---|---|
| 分析 | 已使用（属于研究本身） | 编写并运行 R 分析流程：数据整合、缺失值插补、零模型、空间自回归模型、第四角检验、保护优先区算法 | 方法、结果、全部图表 | ⬜ |
| 数据可视化 | 已使用 | 编写绘制全部图件的 R 代码；未使用生成式图像模型 | 正文图 1–8、扩展数据图 1–25 | ⬜ |
| 撰写 | 已使用 | 起草手稿各部分、图注、补充信息、投稿信与声明 | 全文 | ⬜ |
| 修订与一致性核查 | 已使用 | 把每个数字回溯到结果表、按引用顺序编号、更正错误 | 全文 | ⬜ |
| 文献核查 | 已使用 | 对照 Crossref 核实 100 条参考文献，核对正文引用与文献表 | 参考文献 | ⬜ |
| 文献与资料检索 | 已使用 | 检索性状数据集、期刊投稿要求与相关文献 | 引言、方法 | ⬜ |
| 文字编辑与排版 | 已使用 | 语言润色、格式转换 | 全文 | ⬜ |
| 模拟审稿 | 不确定 | 以审稿人视角做过投稿前审计 | — | ⬜ 请确认是否算作此类 |
| 本项目之外的其他 AI 工具 | **未知** | — | — | ⬜ **必须由你声明：无，或列出** |

**从工作记录看，你本人做了这些事**（请核对、补充，用于写"作者角色"那句）：确定研究主题与目标；
提出关键要求（如每纲只用一个体型量度、纳入生活史性状、把驱动变量对应到经典假说、面向保护规划）；
提供或指定数据来源（王彦平团队性状数据集、你发表的哺乳类性状数据集、标准底图、保护地数据）；
提出审阅问题（如生境异质性度量、33 个科的依据、西北空白的成因）；做出范围决定（如保留全部 25 张补充图）；
审阅各轮输出。

## 三、手稿 Methods 中的草稿（已放入，括号内须你改定）

> **Use of AI tools.** An AI assistant (Claude; Anthropic), including the Opus 5 and Opus 5.5
> models accessed through Claude Code between August and October 2026, was used under the
> author's direction to write and run the R analysis code and produce the figures from that
> code; to draft and revise the manuscript text, figure legends and Supplementary Information;
> and to check references against Crossref records and reported values against the result
> tables. No generative image model was used: every figure is a plot of data rendered by R
> code. [The author's own role.] [Any AI tools used outside this workflow, or "None".]

## 四、提交前你需要确认的三件事

1. 上表每一行是否属实、是否完整（尤其"本项目之外的其他 AI 工具"一行）
2. 你能否作出 Nature 要求的承诺：对最终文本负责，并同意这些修改反映的是你本人的研究工作
3. 删去手稿中所有方括号提示，再提交
