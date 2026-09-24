#!/usr/bin/env python3
# ============================================================
# 目的 / Purpose:
#   按正文首次引用顺序重排正文图、扩展数据图与补充表的编号（Nature 系列
#   要求按首次引用顺序编号）。先以占位符替换旧编号，再映射为新编号，
#   避免 5→6、6→5 这类互换时的二次替换错误。
#   Renumber main figures, Extended Data figures and Supplementary Tables in
#   order of first citation. Old numbers are first replaced by placeholder
#   tokens, then mapped, so swaps (5->6, 6->5) cannot double-apply.
# 用法 / Usage:
#   python3 renumber_displays.py --dry-run   # 只打印映射
#   python3 renumber_displays.py             # 应用到手稿与图注
# ============================================================
import re, io, sys, json

MS  = '06_manuscript/manuscript_v4_full.md'
LEG = '09_submission/02_manuscript/figure_legends.md'
CUT = '## Figure legends'

RE_MAIN = re.compile(r'(?<![A-Za-z]\s)(?<!Data\s)(?<!Data\n)\b(Figs?\.)(\s?)(\d+)')
RE_ED   = re.compile(r'(Extended\s+Data\s+Figs?\.?)(\s?)(\d+)((?:\s*[–-]\s*\d+)?)')
RE_ST   = re.compile(r'(Supplementary\s+Tables?\s+)(S?[0-9]+[a-z]*(?:pc)?)((?:\s*[–-]\s*[a-z]\b)?)')

def first_order(text, rx, grp):
    seen=[]
    for m in rx.finditer(text):
        n=m.group(grp)
        if n not in seen: seen.append(n)
    return seen

def ed_first_order(text):
    seen=[]
    for m in RE_ED.finditer(text):
        a=int(m.group(3)); rng=m.group(4)
        r=re.search(r'(\d+)',rng) if rng else None
        for n in (range(a,int(r.group(1))+1) if r else [a]):
            if n not in seen: seen.append(n)
    return seen

def main_first_order(text):
    seen=[]
    for m in re.finditer(r'(?<!Data)(?<!Data\s)\bFigs?\.\s?(\d+)',text):
        pre=text[max(0,m.start()-14):m.start()]
        if re.search(r'Data\s*$',pre): continue
        n=int(m.group(1))
        if n not in seen: seen.append(n)
    return seen

def st_first_order(text):
    seen=[]
    for m in RE_ST.finditer(text):
        tid=m.group(2)
        tail=m.group(3)
        ids=[tid]
        if tail:   # 'S3b–c' → S3b, S3c
            base=re.match(r'(S?\d+)',tid).group(1); ids=[tid, base+tail.strip(' –-')]
        for t in ids:
            if t not in seen: seen.append(t)
    return seen

if __name__=='__main__':
    dry='--dry-run' in sys.argv
    ms=io.open(MS,encoding='utf-8').read()
    body, legends = ms.split(CUT,1)
    mo=main_first_order(body)
    eo=ed_first_order(body)
    so=st_first_order(body)+[t for t in st_first_order(legends) if t not in st_first_order(body)]
    main_map={old:i+1 for i,old in enumerate(mo)}
    ed_map={old:i+1 for i,old in enumerate(eo)}
    st_map={old:i+1 for i,old in enumerate(so)}
    print("正文图 首次引用顺序:", mo, "→ 映射", {k:v for k,v in main_map.items() if k!=v} or "无变化")
    print("扩展图 首次引用顺序:", eo)
    print("   缺失:", [n for n in range(1,26) if n not in eo])
    print("   变化:", {k:v for k,v in ed_map.items() if k!=v})
    print("补充表 首次引用顺序（内部表号）:", so)
    json.dump({'main':main_map,'ed':ed_map,'st':st_map}, open('/tmp/renumber_maps.json','w'))
