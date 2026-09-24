#!/usr/bin/env python3
# ============================================================
# 按首次引用顺序重编号：正文图、扩展数据图、补充表（Nature 系列要求）。
# 旧编号先换成占位符，再映射为新编号，互换（5<->6）不会被二次替换。
# Renumber figures, Extended Data figures and Supplementary Tables in order of
# first citation. Old numbers become placeholders first, so swaps cannot
# double-apply. Writes /tmp/renumber_maps.json for downstream steps.
# ============================================================
import re, io, json

MS='06_manuscript/manuscript_v4_full.md'; LEG='09_submission/02_manuscript/figure_legends.md'
CUT='## Figure legends'
ED_RX = re.compile(r'(Extended\s+Data\s+Figs?\.?\s?)(\d+)([a-z](?:[,–-]?[a-z])*)?((?:\s*[–-]\s*\d+)?)')
MAIN_RX = re.compile(r'\b(Figs?\.\s?)(\d+)([a-z](?:[,–-]?[a-z])*)?')
ST_RX = re.compile(r'(Supplementary\s+Tables?\s+)((?:S?\d+[a-z]*)(?:(?:\s*,\s*|\s+and\s+|\s*[–-]\s*)(?:S?\d+[a-z]*|[a-z]\b))*)')
ID_RX = re.compile(r'S?\d+[a-z]*')

def is_ed_context(text, pos):
    return re.search(r'Extended\s+Data\s*$', text[max(0,pos-20):pos]) is not None

def expand_st(group):
    """'47b, 48' -> ['47b','48'];  'S3b–c' -> ['S3b','S3c']"""
    ids=[]; prev=None
    for tok in re.split(r'(\s*,\s*|\s+and\s+|\s*[–-]\s*)', group):
        t=tok.strip()
        if not t or re.fullmatch(r',|and|[–-]',t): continue
        if re.fullmatch(r'[a-z]',t) and prev:   # S3b–c 的 c
            t=re.match(r'S?\d+',prev).group(0)+t
        ids.append(t); prev=t
    return ids

def orders(body, legends):
    mo,eo,so=[],[],[]
    for txt in (body, legends):
        for m in MAIN_RX.finditer(txt):
            if is_ed_context(txt,m.start()): continue
            n=int(m.group(2))
            if txt is body and n not in mo: mo.append(n)
        for m in ED_RX.finditer(txt):
            a=int(m.group(2)); r=re.search(r'(\d+)',m.group(4) or '')
            for n in (range(a,int(r.group(1))+1) if r else [a]):
                if txt is body and n not in eo: eo.append(n)
        for m in ST_RX.finditer(txt):
            for t in expand_st(m.group(2)):
                if t not in so: so.append(t)
    return mo,eo,so

def fmt_list(nums, prefix_single, prefix_plural):
    nums=sorted(nums)
    if len(nums)==1: return f"{prefix_single}{nums[0]}"
    if nums==list(range(nums[0],nums[-1]+1)) and len(nums)>2: return f"{prefix_plural}{nums[0]}–{nums[-1]}"
    if len(nums)==2: return f"{prefix_plural}{nums[0]} and {nums[1]}" if nums[1]!=nums[0]+1 else f"{prefix_plural}{nums[0]}–{nums[1]}"
    return f"{prefix_plural}{', '.join(map(str,nums[:-1]))} and {nums[-1]}"

def apply(text, main_map, ed_map, st_map):
    # 扩展图（先处理，避免被正文图规则误伤）
    def ed_sub(m):
        a=int(m.group(2)); r=re.search(r'(\d+)',m.group(4) or ''); panel=m.group(3) or ''
        olds=list(range(a,int(r.group(1))+1)) if r else [a]
        news=[ed_map[o] for o in olds]
        if len(news)==1: return f"⟦ED{news[0]}{panel}⟧"
        return "⟦EDL"+",".join(map(str,news))+"⟧"
    text=ED_RX.sub(ed_sub, text)
    def main_sub(m):
        n=int(m.group(2)); panel=m.group(3) or ''
        return f"{m.group(1)}⟦F{main_map.get(n,n)}⟧{panel}"
    text=MAIN_RX.sub(main_sub, text)
    def st_sub(m):
        ids=expand_st(m.group(2)); news=sorted({st_map[i] for i in ids})
        return fmt_list(news,"Supplementary Table ","Supplementary Tables ")
    text=ST_RX.sub(st_sub, text)
    # 占位符还原
    text=re.sub(r'⟦ED(\d+)([a-z,–-]*)⟧', lambda m: f"Extended Data Fig. {m.group(1)}{m.group(2)}", text)
    text=re.sub(r'⟦EDL([\d,]+)⟧', lambda m: fmt_list([int(x) for x in m.group(1).split(',')],"Extended Data Fig. ","Extended Data Figs "), text)
    text=re.sub(r'⟦F(\d+)⟧', r'\1', text)
    return text

def reorder_legends(block, kind):
    """把图注块按新编号排序。main: '### Fig. N |'；ed: '**ED N |'"""
    if kind=='main':
        parts=re.split(r'(?=^### Fig\. \d+ \|)', block, flags=re.M)
        key=lambda p: int(re.match(r'### Fig\. (\d+)',p).group(1)) if p.startswith('### Fig.') else -1
    else:
        parts=re.split(r'(?=^\*\*ED \d+ \|)', block, flags=re.M)
        key=lambda p: int(re.match(r'\*\*ED (\d+)',p).group(1)) if p.startswith('**ED') else -1
    head=[p for p in parts if key(p)==-1]; items=sorted([p for p in parts if key(p)!=-1], key=key)
    items=[p.rstrip()+"\n\n" for p in items]
    return "".join(head)+"".join(items)

if __name__=='__main__':
    ms=io.open(MS,encoding='utf-8').read(); body,legends=ms.split(CUT,1)
    mo,eo,so=orders(body,legends)
    assert sorted(mo)==list(range(1,9)), f"正文图引用不全: {mo}"
    assert sorted(eo)==list(range(1,26)), f"扩展图引用不全: {sorted(set(range(1,26))-set(eo))}"
    main_map={o:i+1 for i,o in enumerate(mo)}; ed_map={o:i+1 for i,o in enumerate(eo)}
    st_map={o:i+1 for i,o in enumerate(so)}
    # 图注标题中的编号：'### Fig. N |' 与 '**ED N |' 也要映射
    def head_sub(t):
        t=re.sub(r'^### Fig\. (\d+) \|', lambda m: f"### Fig. ⟦F{main_map[int(m.group(1))]}⟧ |", t, flags=re.M)
        t=re.sub(r'^\*\*ED (\d+) \|', lambda m: f"**ED ⟦H{ed_map[int(m.group(1))]}⟧ |", t, flags=re.M)
        return t
    def finish(t):
        t=re.sub(r'⟦F(\d+)⟧', r'\1', t); t=re.sub(r'⟦H(\d+)⟧', r'\1', t); return t
    # 手稿
    new_body=apply(body,main_map,ed_map,st_map)
    L=finish(apply(head_sub(legends),main_map,ed_map,st_map))
    mhead, rest = L.split('## Extended Data figure legends',1)
    edblock, tail = rest.split('\n## References',1)
    L=reorder_legends(mhead,'main')+'## Extended Data figure legends'+reorder_legends(edblock,'ed')+'\n## References'+tail
    io.open(MS,'w',encoding='utf-8').write(new_body+CUT+L)
    # 独立图注文件
    lg=io.open(LEG,encoding='utf-8').read()
    pre, rest = lg.split('## Main figures',1)
    mainblk, rest2 = rest.split('## Extended Data figures',1)
    edblk, sup = rest2.split('## Superseded',1)
    mainblk=finish(apply(head_sub(mainblk),main_map,ed_map,st_map))
    edblk=finish(apply(head_sub(edblk),main_map,ed_map,st_map))
    lg=pre+'## Main figures'+reorder_legends(mainblk,'main')+'---\n\n## Extended Data figures'+reorder_legends(edblk,'ed')+'---\n\n## Superseded'+sup
    io.open(LEG,'w',encoding='utf-8').write(lg)
    json.dump({'main':main_map,'ed':ed_map,'st':st_map}, open('/tmp/renumber_maps.json','w'))
    print("正文图映射:", {k:v for k,v in main_map.items() if k!=v})
    print("扩展图映射（旧→新）:", ed_map)
    print("补充表映射（内部表号→S 编号）:", st_map)
