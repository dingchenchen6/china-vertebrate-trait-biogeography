#!/usr/bin/env python3
# ============================================================
# 目的 / Purpose:
#   按手稿首次引用顺序（renumber_apply.py 生成的映射）从 04_results/tables/
#   的当前结果表构建补充表册。不从 tables_formatted/ 取数：那里的表生成于
#   8 月 3 日，早于零模型修正，已过时。
#   Build the Supplementary Tables workbook from the CURRENT result tables, in
#   order of first citation. tables_formatted/ predates the null-model
#   correction and is not used.
# 用法 / Usage: python3 09_submission/build_supplementary_tables.py
# ============================================================
import csv, io, os, re, json, glob, math, sys
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

T = '04_results/tables'

# 内部表号 -> 文件；有歧义或非标准命名的显式指定
FILE = {'S0': f'{T}/table_S0_hypotheses.csv', 'S3b': f'{T}/table_S3b_unmatched_species.csv',
        'S3c': f'{T}/table_S3c_taxonomy_harmonisation.csv',
        '30': f'{T}/table_30_reserve_designation_targets.csv',
        '6bpc': '04_results/tables_formatted/TableS19_size_measure_robustness.csv'}

def resolve(k):
    """内部表号 -> 文件路径 / internal id -> csv path"""
    if k in FILE: return FILE[k]
    hit = [f for f in glob.glob(f'{T}/table_{k}_*.csv')]
    if len(hit) != 1: sys.exit(f"表号 {k} 无法唯一解析: {hit}")
    return hit[0]

TITLE = {
 'S0': "Questions, hypotheses and refutation conditions stated in advance, with the outcome of each test",
 'S3b': "Species names that could not be matched between trait and range sources",
 'S3c': "Taxonomic harmonisation between trait and range sources",
 '36': "Observed coverage of life-history traits by class, before imputation",
 '38b': "Variance inflation factors of the six named predictors",
 '38d': "Residual spatial autocorrelation (Moran's I) of the 36 main models under OLS and under SAR",
 '6bpc': "Family-level thermal contrast under two body-size definitions (body mass for all classes versus one measure per class)",
 '41b': "Sensitivity of the trait-space prioritisation to the number of trait bins",
 '37b': "Within-class correlations among the three a priori functional axes",
 '39b': "Mean absolute effect of component traits against the composite axis they form",
 '39': "Unique and shared variance of climate, habitat and human predictor groups, by class and response",
 '38': "Standardised effects of the six named predictors on every trait and axis, by class (SAR models)",
 '45b': "Ectotherm/endotherm effect ratios under three habitat-heterogeneity proxies",
 '46b': "Family-level thermal contrast at three inclusion thresholds",
 '40b': "Type I inflation of naive CWM regression relative to the fourth-corner max test",
 '40c': "Trait-environment links surviving the max test, by predictor",
 '47b': "Correlations of three functional-diversity indices with the six predictors",
 '48': "Environmental vectors fitted to the NMDS ordination of assemblage trait composition",
 '13': "Congruence (adjusted Rand index) of class-specific functional regionalisations",
 '20': "Functional redundancy: functional relative to taxonomic beta diversity, by class",
 '28': "Jaccard overlap of conservation-gap cells between classes",
 '30': "Classes named as protection targets in reserve designation records (1,028 reserves)",
 '42': "Efficiency of the existing reserve network against an optimal layout of the same area",
 '41': "Coverage of the 30%-of-range target under species and trait-space targets across budgets",
 '43': "Effect of the functional axes on the protected fraction of species ranges (quasi-binomial GLM)",
 '49': "Land-use tolerance of species against the functional axes",
 '50': "Protected fraction of range for land-use-tolerant and intolerant species",
 '17': "Observed coverage of amphibian traits and which are retained",
}
def val(v):
    v = v.strip()
    if re.fullmatch(r'-?\d+', v): return int(v)
    try: x = float(v)
    except ValueError: return v
    if not math.isfinite(x): return v
    return float(f"{x:.3g}") if (x != 0 and abs(x) < 1e-3) else round(x, 4)


def build(order, out, label="Supplementary Table {n}", sheet="Table S{n}",
          heading="Supplementary Tables", note=None):
    """按给定顺序构建表册。order: 内部表号列表（第 i 个即第 i+1 号表）
    Build the workbook in the given order; order[i] becomes table i+1."""
    missing = [k for k in order if k not in TITLE]
    if missing: sys.exit(f"缺标题: {missing}")
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Contents"
    H = Font(bold=True, color="FFFFFF"); F = PatternFill("solid", fgColor="1F4E79")
    ws.append([heading]); ws["A1"].font = Font(bold=True, size=14)
    ws.append(["Thermoregulatory mode sets the strength of environmental filtering on vertebrate assemblage traits"])
    ws.append([note or ("Numbered in order of first citation in the manuscript. Column names are those written by "
               "the analysis code (github.com/dingchenchen6/china-vertebrate-trait-biogeography), so every "
               "value can be traced to the script that produced it.")])
    ws.append([]); ws.append(["Table", "Title", "Rows", "Source file"])
    for c in ws[5]: c.font = H; c.fill = F
    for i, k in enumerate(order, 1):
        path = resolve(k)
        rows = list(csv.reader(io.open(path, encoding='utf-8')))
        sh = wb.create_sheet(sheet.format(n=i))
        sh.append([f"{label.format(n=i)} | {TITLE[k]}"]); sh["A1"].font = Font(bold=True)
        sh.append([]); sh.append(rows[0])
        for c in sh[3]: c.font = H; c.fill = F
        for r in rows[1:]: sh.append([val(v) for v in r])
        sh.freeze_panes = "A4"
        for j, col in enumerate(sh.iter_cols(min_row=3, max_row=min(sh.max_row, 80)), 1):
            w = max((len(str(c.value)) for c in col if c.value is not None), default=8)
            sh.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 9), 48)
        ws.append([label.format(n=i), TITLE[k], len(rows) - 1, os.path.basename(path)])
    for col, w in zip("ABCD", (20, 92, 7, 44)): ws.column_dimensions[col].width = w
    for r in ws.iter_rows(min_row=6): r[1].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells("A3:D3"); ws["A3"].alignment = Alignment(wrap_text=True); ws.row_dimensions[3].height = 32
    wb.save(out)
    return out


if __name__ == '__main__':
    MAPS = sys.argv[1] if len(sys.argv) > 1 else '09_submission/display_numbering.json'
    st = json.load(open(MAPS))['st']
    order = [k for k, _ in sorted(st.items(), key=lambda kv: kv[1])]
    out = build(order, '09_submission/05_supplementary/Supplementary_Tables.xlsx')
    print(f"✓ {out}: {len(order)} tables + contents")
