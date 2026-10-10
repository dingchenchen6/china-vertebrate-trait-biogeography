#!/usr/bin/env python3
# ============================================================
# Scientific question / 科学问题:
#   （工具脚本）同一篇手稿按不同期刊的硬性要求生成投稿包。
#   Tooling: derive journal-specific submission packages from one master.
#
# Objective / 分析目标:
#   以 06_manuscript/manuscript_v4_full.md 为唯一来源，生成
#   Nature Communications 与 Global Ecology and Biogeography 两个投稿包，
#   各自按首次引用顺序重排正文图、补充图与补充表的编号，并装配图件与表册。
#   Build NatCommun and GEB packages from the single master manuscript,
#   renumbering every display series in order of first citation per journal.
#
# Input data / 输入数据:
#   06_manuscript/manuscript_v4_full.md（主稿）; table_00_figure_manifest.csv;
#   09_submission/display_numbering.json; 05_figures/; journal_text/*.md
# Main workflow / 主要流程:
#   1. 拆分主稿 → 2. 按期刊重组章节与图表称谓 → 3. 按首次引用重编号
#   4. 装配图件、补充信息文档与表册 → 5. pandoc 导出 docx
# Expected output / 预期输出:
#   10_submission_NatCommun/, 11_submission_GEB/
# Key assumptions / 关键假设:
#   主稿中每张图、每张表都已被正文引用（上一步 renumber_apply.py 已保证）。
#   页眉（标题、字数行）不参与引用解析。
# Main packages / 主要包: re, openpyxl（经 build_supplementary_tables）, pandoc, sips
# Output directory / 输出路径: 10_submission_NatCommun/, 11_submission_GEB/
# ============================================================
import csv, io, json, os, re, shutil, subprocess, sys

sys.path.insert(0, '09_submission')
from build_supplementary_tables import build as build_tables  # noqa: E402
import build_supplementary_tables as bst  # noqa: E402

MASTER = '06_manuscript/manuscript_v4_full.md'
MANIFEST = '04_results/tables/table_00_figure_manifest.csv'
NUMBERING = '09_submission/display_numbering.json'
TEXT = '09_submission/journal_text'


def read(p):
    return io.open(p, encoding='utf-8').read()


def write(p, s):
    os.makedirs(os.path.dirname(p) or '.', exist_ok=True)
    io.open(p, 'w', encoding='utf-8').write(s)


def words(t):
    """与前几轮相同的计词口径 / same word-count rule as previous rounds"""
    t = re.sub(r'^#+.*$', '', t, flags=re.M)
    t = re.sub(r'[*_`>|]', ' ', t)
    return len([w for w in t.split() if re.search(r'[A-Za-z0-9]', w)])


# ------------------------------------------------------------
# 1. 拆分主稿 / split the master manuscript
# ------------------------------------------------------------
def split_master(s):
    def sec(a, b):
        i = s.index(a) + len(a)
        j = s.index(b, i) if b else len(s)
        return re.sub(r'\n-{3,}\s*$', '', s[i:j].strip()).strip()
    kw = re.search(r'\*\*Keywords\*\*\s+(.+?)\n\s*\n', s, re.S).group(1)
    return dict(
        title=re.match(r'# (.+)', s).group(1).strip(),
        abstract=sec('## Abstract\n', '**Keywords**'),
        keywords=re.sub(r'\s+', ' ', kw).strip(),
        intro=sec('## 1. Introduction\n', '## 2. Materials and Methods'),
        methods=sec('## 2. Materials and Methods\n', '## 3. Results'),
        results=sec('## 3. Results\n', '## 4. Discussion'),
        discussion=sec('## 4. Discussion\n', '## Data and code availability'),
        contrib=sec('## Author contributions\n', '## Competing interests'),
        competing=sec('## Competing interests\n', '## Figure legends'),
        legends_main=sec('## Figure legends\n', '## Extended Data figure legends'),
        legends_ed=sec('## Extended Data figure legends\n', '## References'),
        refs=sec('## References\n', None),
    )


def denumber(t):
    """去掉小节编号与 § 引用 / drop section numbers and § pointers"""
    t = re.sub(r'^### \d+\.\d+ ', '### ', t, flags=re.M)
    t = t.replace('(§ Data and code availability)', '(see Data availability)')
    t = re.sub(r'\s*\(§\d+(?:\.\d+)?\)', ' (see Results)', t)
    return t


def statements():
    s = read('09_submission/06_statements/statements.md')
    def sec(name):
        m = re.search(r'## ' + name + r'\n(.*?)(?=\n## |\n---)', s, re.S)
        return m.group(1).strip()
    data = sec('Data availability').replace('Methods §2.3', 'Methods')
    code = sec('Code availability')
    return data, code


# ------------------------------------------------------------
# 2. 引用解析与重编号 / citation parsing and renumbering
# ------------------------------------------------------------
SEP = r'(?:\s*,\s*|\s+and\s+|\s*[–-]\s*)'
PANEL = r'(?:[a-z](?:[,–-][a-z])*)'


def nums_from(head, tail):
    """'2' + '–4' → [2,3,4]；'2' + ', 5 and 7' → [2,5,7]"""
    first = int(head)
    r = re.match(r'\s*[–-]\s*(\d+)\s*$', tail or '')
    if r:
        return list(range(first, int(r.group(1)) + 1))
    return [first] + [int(x) for x in re.findall(r'\d+', tail or '')]


def fmt(nums, single, plural, prefix=''):
    nums = sorted(set(nums))
    p = lambda n: f'{prefix}{n}'
    if len(nums) == 1:
        return f'{single}{p(nums[0])}'
    if nums == list(range(nums[0], nums[-1] + 1)):
        return f'{plural}{p(nums[0])}–{p(nums[-1])}'
    return f'{plural}{", ".join(p(n) for n in nums[:-1])} and {p(nums[-1])}'


class Kind:
    """一类引用：正文图 / 补充图 / 补充表"""
    def __init__(self, rx, single, plural, prefix='', guard=None):
        self.rx = re.compile(rx)
        self.single, self.plural, self.prefix, self.guard = single, plural, prefix, guard

    def hits(self, text):
        for m in self.rx.finditer(text):
            if self.guard and self.guard(text, m.start()):
                continue
            yield m


def not_main(text, pos):
    """正文图规则不能吃掉「Supplementary Fig.」或「Extended Data Fig.」"""
    return re.search(r'(Supplementary|Extended\s+Data)\s*$', text[max(0, pos - 22):pos]) is not None


def kinds(supp_word):
    return dict(
        main=Kind(r'\b(Figs?\.)\s?(\d+)(' + PANEL + r')?((?:' + SEP + r'\d+)*)',
                  'Fig. ', 'Figs ', guard=not_main),
        supf=Kind(supp_word + r'\s+Figs?\.?\s?(\d+)(' + PANEL + r')?((?:' + SEP + r'\d+)*)',
                  supp_word + ' Fig. ', supp_word + ' Figs '),
        supt=Kind(r'Supplementary\s+Tables?\s+(\d+)()((?:' + SEP + r'\d+)*)',
                  'Supplementary Table ', 'Supplementary Tables '),
    )


def first_order(text, kind):
    seen = []
    for m in kind.hits(text):
        g = m.groups()
        num_head, tail = (g[1], g[3]) if len(g) == 4 else (g[0], g[2])
        for n in nums_from(num_head, tail):
            if n not in seen:
                seen.append(n)
    return seen


def remap(text, kind, mapping, out_single=None, out_plural=None, out_prefix=''):
    """按映射替换编号；占位符避免互换时二次替换"""
    out_single = out_single or kind.single
    out_plural = out_plural or kind.plural

    def sub(m):
        if kind.guard and kind.guard(text, m.start()):
            return m.group(0)
        g = m.groups()
        if len(g) == 4:
            head, panel, tail = g[1], g[2] or '', g[3]
        else:
            head, panel, tail = g[0], g[1] or '', g[2]
        olds = nums_from(head, tail)
        news = [mapping[o] for o in olds]
        if len(news) == 1:
            return f'\u27e6{out_single}{out_prefix}{news[0]}{panel}\u27e7'
        return '\u27e6' + fmt(news, out_single, out_plural, out_prefix) + '\u27e7'
    text = kind.rx.sub(sub, text)
    return text


def unwrap(text):
    return re.sub('[\u27e6\u27e7]', '', text)


def reorder(block, header_rx):
    """按新编号给图注块排序 / sort legend blocks by their new number"""
    parts = re.split(r'(?=' + header_rx + r')', block, flags=re.M)
    head = [p for p in parts if not re.match(header_rx, p)]
    items = [p for p in parts if re.match(header_rx, p)]
    key = lambda p: int(re.search(r'(\d+)', re.match(header_rx, p).group(0)).group(1))
    return ''.join(head) + ''.join(p.rstrip() + '\n\n' for p in sorted(items, key=key))


# ------------------------------------------------------------
# 3. 图件与表册 / figures and tables
# ------------------------------------------------------------
def manifest():
    m = {}
    for r in csv.DictReader(open(MANIFEST)):
        m[r['num']] = r['src']
    return m


def master_table_ids():
    st = json.load(open(NUMBERING))['st']
    return {v: k for k, v in st.items()}          # master number -> internal id


def copy_fig(src, dest_stem):
    for ext in ('pdf', 'png'):
        shutil.copy2(f'05_figures/{src}.{ext}', f'{dest_stem}.{ext}')


def preview_png(src, dest, width=2200):
    """补充信息文档中嵌入的预览图：约 300 dpi，控制文件体积"""
    shutil.copy2(f'05_figures/{src}.png', dest)
    subprocess.run(['sips', '-Z', str(width), dest], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def pandoc(md_path, out_path, resource_dir):
    subprocess.run(['pandoc', os.path.basename(md_path), '-o', os.path.basename(out_path),
                    f'--resource-path={resource_dir}'], check=True, cwd=os.path.dirname(md_path))


def journal_s0(transform, out_csv):
    """S1 假设表的「证据」列含图表编号，须随期刊映射改写
    The hypothesis table cites display items, so its evidence column is remapped per journal."""
    rows = list(csv.DictReader(io.open('04_results/tables/table_S0_hypotheses.csv', encoding='utf-8')))
    for r in rows:
        r['evidence'] = transform(r['evidence'])
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with io.open(out_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    bst.FILE['S0'] = out_csv


# ------------------------------------------------------------
# 4. Nature Communications
# ------------------------------------------------------------
def build_natcommun(P, out='10_submission_NatCommun'):
    K = kinds('Supplementary')
    data, code = statements()
    title = read(f'{TEXT}/natcommun_title.txt').strip()
    ai = read(f'{TEXT}/natcommun_ai_methods.md').strip()

    rename = lambda t: re.sub(r'Extended(\s+)Data(\s+)Fig', r'Supplementary\2Fig', t)
    intro, results, discussion, methods = (rename(denumber(P[k])) for k in
                                           ('intro', 'results', 'discussion', 'methods'))
    data, code = rename(data), rename(code)
    methods = methods + '\n\n### Use of AI tools\n\n' + ai

    body = '\n\n'.join([
        '## Introduction\n\n' + intro, '## Results\n\n' + results,
        '## Discussion\n\n' + discussion, '## Methods\n\n' + methods,
        '## Data availability\n\n' + data, '## Code availability\n\n' + code])
    legends_main = rename(P['legends_main'])          # 正文图注里也有扩展图引用
    legends_sup = re.sub(r'^\*\*ED (\d+) \|', r'**Supplementary Fig. \1 |', rename(P['legends_ed']), flags=re.M)

    # 首次引用顺序：正文（不含页眉）→ 必须覆盖全部图；表可补自图注
    mo, so = first_order(body, K['main']), first_order(body, K['supf'])
    to = first_order(body, K['supt'])
    for n in first_order(legends_main + legends_sup, K['supt']):
        if n not in to:
            to.append(n)
    assert sorted(mo) == list(range(1, 9)), f'正文图引用不全: {mo}'
    assert sorted(so) == list(range(1, 26)), f'补充图引用不全: {sorted(set(range(1, 26)) - set(so))}'
    assert sorted(to) == list(range(1, 29)), f'补充表引用不全: {sorted(set(range(1, 29)) - set(to))}'
    mm = {o: i + 1 for i, o in enumerate(mo)}
    sm = {o: i + 1 for i, o in enumerate(so)}
    tm = {o: i + 1 for i, o in enumerate(to)}

    def apply(t):
        t = remap(t, K['supf'], sm)
        t = remap(t, K['main'], mm)
        t = remap(t, K['supt'], tm)
        return unwrap(t)
    body = apply(body)
    legends_main = apply(re.sub(r'^### Fig\. (\d+) \|',
                                lambda m: f'### Fig. \u27e6{mm[int(m.group(1))]}\u27e7 |',
                                legends_main, flags=re.M))
    legends_sup = apply(re.sub(r'^\*\*Supplementary Fig\. (\d+) \|',
                               lambda m: f'**Supplementary Fig. \u27e6{sm[int(m.group(1))]}\u27e7 |',
                               legends_sup, flags=re.M))
    legends_main = reorder(unwrap(legends_main), r'### Fig\. \d+ \|')
    legends_sup = reorder(unwrap(legends_sup), r'\*\*Supplementary Fig\. \d+ \|')

    # 字数（与前几轮同口径）
    def part(t, a, b):
        i = t.index(a); j = t.index(b, i + len(a)) if b else len(t)
        return t[i:j]
    w = dict(abstract=words(P['abstract']),
             intro=words(part(body, '## Introduction', '## Results')),
             results=words(part(body, '## Results', '## Discussion')),
             discussion=words(part(body, '## Discussion', '## Methods')),
             methods=words(part(body, '## Methods', '## Data availability')))
    w['main'] = w['intro'] + w['results'] + w['discussion']
    nref = len([q for q in re.split(r'\n\s*\n', P['refs']) if re.search(r'\(\d{4}[a-z]?\)', q)])
    header = (f'# {title}\n\n**Chenchen Ding**\n\nSchool of Life Sciences, Peking University, '
              f'Beijing 100871, China  \nCorrespondence: jialeding1220@gmail.com\n\n'
              f'**Word count** Title {len(title.split())} words | Abstract {w["abstract"]} | '
              f'Main text {w["main"]:,} (Introduction {w["intro"]:,}, Results {w["results"]:,}, '
              f'Discussion {w["discussion"]:,}) | Methods {w["methods"]:,} | Display items 8 | '
              f'Supplementary Figures 25 | Supplementary Tables 28 | References {nref}\n\n')
    ms = (header + '## Abstract\n\n' + P['abstract'] + '\n\n' + body +
          '\n\n## References\n\n' + P['refs'] +
          '\n\n## Acknowledgements\n\n[To be completed by the author: funding, data providers, '
          'and any non-author contributions.]\n\n## Author contributions\n\n' + P['contrib'] +
          '\n\n## Competing interests\n\n' + P['competing'] +
          '\n\n## Figure legends\n\n' + legends_main.strip() + '\n')
    write(f'{out}/02_manuscript/manuscript.md', ms)
    pandoc(f'{out}/02_manuscript/manuscript.md', f'{out}/02_manuscript/manuscript.docx', '.')

    # 图件：正文图与补充图
    man, inv_m, inv_s = manifest(), {v: k for k, v in mm.items()}, {v: k for k, v in sm.items()}
    for d in ('03_figures', '04_supplementary/figures', '04_supplementary/preview'):
        shutil.rmtree(f'{out}/{d}', ignore_errors=True); os.makedirs(f'{out}/{d}')
    for n in range(1, 9):
        copy_fig(man[f'Fig{inv_m[n]}'], f'{out}/03_figures/Fig{n}')
    si = ['# Supplementary Information\n',
          f'**{title}**\n', 'Chenchen Ding\n',
          'Supplementary Figures 1–25 follow, each with its legend. Supplementary Tables 1–28 are '
          'provided as a separate Excel workbook (`Supplementary_Tables_1-28.xlsx`), numbered in order '
          'of first citation in the main text.\n']
    leg = {int(m.group(1)): m.group(0) for m in
           re.finditer(r'^\*\*Supplementary Fig\. (\d+) \|.*?(?=^\*\*Supplementary Fig\. \d+ \||\Z)',
                       legends_sup, flags=re.M | re.S)}
    for n in range(1, 26):
        src = man[f'ED{inv_s[n]}']
        copy_fig(src, f'{out}/04_supplementary/figures/SupplementaryFig{n}')
        preview_png(src, f'{out}/04_supplementary/preview/SupplementaryFig{n}.png')
        si.append('```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n')
        si.append(f'![](preview/SupplementaryFig{n}.png){{width=16cm}}\n')
        si.append(leg[n].strip() + '\n')
    write(f'{out}/04_supplementary/Supplementary_Information.md', '\n'.join(si))
    pandoc(f'{out}/04_supplementary/Supplementary_Information.md',
           f'{out}/04_supplementary/Supplementary_Information.docx', '.')
    journal_s0(lambda t: apply(rename(t)), f'{out}/.build/table_S0_hypotheses.csv')
    ids = master_table_ids()
    inv_t = {v: k for k, v in tm.items()}
    build_tables([ids[inv_t[n]] for n in range(1, 29)],
                 f'{out}/04_supplementary/Supplementary_Tables_1-28.xlsx')
    json.dump({'title': title, 'words': w, 'main': mm, 'supplementary_figures': sm,
               'supplementary_tables': tm},
              open(f'{out}/numbering_from_master.json', 'w'), indent=1)
    return w, mm, sm, tm


# ------------------------------------------------------------
# 5. Global Ecology and Biogeography
# ------------------------------------------------------------
def build_geb(P, out='11_submission_GEB'):
    K = kinds('Extended Data')
    data, code = statements()
    data = re.sub(r'the Supplementary\s+Tables workbook', 'Appendix S2', data)   # Wiley 称谓
    abstract = read(f'{TEXT}/geb_structured_abstract.md').strip()
    ai = read(f'{TEXT}/geb_ai_placeholder.md').strip()
    intro, methods, results, discussion = (denumber(P[k]) for k in
                                           ('intro', 'methods', 'results', 'discussion'))
    body = '\n\n'.join(['## Introduction\n\n' + intro, '## Methods\n\n' + methods,
                        '## Results\n\n' + results, '## Discussion\n\n' + discussion])
    # GEB 的章节顺序与主稿一致，编号应当不变——核验，而不是假定
    mo, so = first_order(body, K['main']), first_order(body, K['supf'])
    to = first_order(body, K['supt'])
    for n in first_order(P['legends_main'] + P['legends_ed'], K['supt']):
        if n not in to:
            to.append(n)
    assert mo == list(range(1, 9)), f'正文图顺序与主稿不同: {mo}'
    assert so == list(range(1, 26)), f'扩展图顺序与主稿不同: {so}'
    assert to == list(range(1, 29)), f'补充表顺序与主稿不同: {to}'
    ident = {n: n for n in range(1, 29)}

    def wiley(t):
        """改为 Wiley 惯例：Figure N / Figure SN / Table SN"""
        t = remap(t, K['supf'], ident, 'Figure ', 'Figures ', 'S')
        t = remap(t, K['main'], ident, 'Figure ', 'Figures ')
        t = remap(t, K['supt'], ident, 'Table ', 'Tables ', 'S')
        return unwrap(t)
    body, data, code = wiley(body), wiley(data), wiley(code)
    legends_main = wiley(re.sub(r'^### Fig\. (\d+) \|', r'### Figure \1 |', P['legends_main'], flags=re.M))
    legends_sup = wiley(re.sub(r'^\*\*ED (\d+) \|', r'**Figure S\1 |', P['legends_ed'], flags=re.M))

    def part(t, a, b):
        i = t.index(a); j = t.index(b, i + len(a)) if b else len(t)
        return t[i:j]
    w = dict(abstract=words(abstract),
             intro=words(part(body, '## Introduction', '## Methods')),
             methods=words(part(body, '## Methods', '## Results')),
             results=words(part(body, '## Results', '## Discussion')),
             discussion=words(body[body.index('## Discussion'):]))
    w['main_body'] = w['intro'] + w['methods'] + w['results'] + w['discussion']
    nref = len([q for q in re.split(r'\n\s*\n', P['refs']) if re.search(r'\(\d{4}[a-z]?\)', q)])
    header = (f'# {P["title"]}\n\n**Running title** Thermal mode and vertebrate trait biogeography\n\n'
              f'**Chenchen Ding**\n\nSchool of Life Sciences, Peking University, Beijing 100871, China  \n'
              f'Correspondence: jialeding1220@gmail.com\n\n'
              f'**Word count** Abstract {w["abstract"]} | Main body {w["main_body"]:,} '
              f'(Introduction {w["intro"]:,}, Methods {w["methods"]:,}, Results {w["results"]:,}, '
              f'Discussion {w["discussion"]:,}) | Figures 8 | Supporting figures 25 | '
              f'Supporting tables 28 | References {nref}\n\n')
    ms = (header + '## Abstract\n\n' + abstract + '\n\n**Keywords** ' + P['keywords'] + '\n\n' +
          body + '\n\n## Acknowledgements\n\n[To be completed by the author: funding, data providers, '
          'and any non-author contributions.]\n\n## Conflict of interest statement\n\n' + P['competing'] +
          '\n\n## Author contributions\n\n' + P['contrib'] +
          '\n\n## Data availability statement\n\n' + data + '\n\n' + code +
          '\n\n## Use of AI technology\n\n' + ai +
          '\n\n## References\n\n' + P['refs'] +
          '\n\n## Figure legends\n\n' + legends_main.strip() +
          '\n\n## Supporting Information\n\nAdditional supporting information may be found online in '
          'the Supporting Information section at the end of the article.\n\n'
          '**Appendix S1** Supporting figures S1–S25 with legends.\n\n'
          '**Appendix S2** Supporting tables S1–S28 (Excel workbook).\n')
    write(f'{out}/02_manuscript/manuscript.md', ms)
    pandoc(f'{out}/02_manuscript/manuscript.md', f'{out}/02_manuscript/manuscript.docx', '.')

    man = manifest()
    for d in ('03_figures', '04_supporting_information/figures', '04_supporting_information/preview'):
        shutil.rmtree(f'{out}/{d}', ignore_errors=True); os.makedirs(f'{out}/{d}')
    for n in range(1, 9):
        copy_fig(man[f'Fig{n}'], f'{out}/03_figures/Figure{n}')
    leg = {int(m.group(1)): m.group(0) for m in
           re.finditer(r'^\*\*Figure S(\d+) \|.*?(?=^\*\*Figure S\d+ \||\Z)', legends_sup, flags=re.M | re.S)}
    si = ['# Appendix S1. Supporting figures\n', f'**{P["title"]}**\n', 'Chenchen Ding\n',
          'Figures S1–S25, numbered in order of first citation in the main text. Tables S1–S28 are in '
          'Appendix S2 (`Appendix_S2_Tables_S1-S28.xlsx`).\n']
    for n in range(1, 26):
        src = man[f'ED{n}']
        copy_fig(src, f'{out}/04_supporting_information/figures/FigureS{n}')
        preview_png(src, f'{out}/04_supporting_information/preview/FigureS{n}.png')
        si += ['```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n', f'![](preview/FigureS{n}.png){{width=16cm}}\n', leg[n].strip() + '\n']
    write(f'{out}/04_supporting_information/Appendix_S1_Figures_S1-S25.md', '\n'.join(si))
    pandoc(f'{out}/04_supporting_information/Appendix_S1_Figures_S1-S25.md',
           f'{out}/04_supporting_information/Appendix_S1_Figures_S1-S25.docx', '.')
    journal_s0(wiley, f'{out}/.build/table_S0_hypotheses.csv')
    ids = master_table_ids()
    build_tables([ids[n] for n in range(1, 29)],
                 f'{out}/04_supporting_information/Appendix_S2_Tables_S1-S28.xlsx',
                 label='Table S{n}', sheet='Table S{n}', heading='Appendix S2. Supporting tables',
                 note='Tables S1–S28, numbered in order of first citation in the main text. Column names '
                      'are those written by the analysis code, so every value can be traced to its script.')
    json.dump({'words': w}, open(f'{out}/numbering_from_master.json', 'w'), indent=1)
    return w


if __name__ == '__main__':
    P = split_master(read(MASTER))
    w1, mm, sm, tm = build_natcommun(P)
    print('Nature Communications 正文图映射:', {k: v for k, v in mm.items() if k != v} or '不变')
    print('  补充图（主稿扩展图号→NC 补充图号）变化数:', sum(k != v for k, v in sm.items()))
    print('  补充表变化数:', sum(k != v for k, v in tm.items()))
    print('  字数:', w1)
    w2 = build_geb(P)
    print('GEB 字数:', w2)
