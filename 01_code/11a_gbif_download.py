#!/usr/bin/env python3
"""
下载中国境内哺乳、两栖、爬行三类群的 GBIF 出现记录
Download GBIF occurrence records for Chinese mammals, amphibians and reptiles

科学问题 / Scientific question:
    基于专家范围图的群落对细尺度土地利用不敏感，因而系统性低估人为过滤效应。
    实际观测记录能否给出更强的人为信号？
    Expert range maps are insensitive to fine-scale land use and therefore
    understate anthropogenic filtering. Do observed occurrences show more?

说明 / Notes:
    - 鸟类使用本地已下载的 GBIF/eBird 中国 2000-2025 数据集（458 万条），不重复下载。
      Birds use the locally held GBIF/eBird China download; not re-fetched.
    - occurrence/search 的 offset 上限为 100,000，两栖类超过该上限，故按纬度带分段。
      The search endpoint caps offset at 100,000, so amphibians are split by latitude.
    - 每个类群写出一个 CSV，字段仅保留下游需要的列。
"""

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

OUT_DIR = "/Users/dingchenchen/中国脊椎动物群落性状/02_data_raw/gbif"
BASE = "https://api.gbif.org/v1/occurrence/search"
PAGE = 300
MAX_OFFSET = 100_000          # GBIF 硬上限 / GBIF hard limit
FIELDS = ["species", "decimalLongitude", "decimalLatitude", "year",
          "basisOfRecord", "coordinateUncertaintyInMeters", "class", "order"]

# 现行 GBIF 主干把 Squamata / Testudines / Crocodylia 处理为纲，
# 而非置于并系的 Reptilia 之下；因此按这三个键取爬行类。
# The current GBIF backbone treats these as classes rather than nesting them
# in the paraphyletic Reptilia, so reptiles are fetched under all three keys.
GROUPS = {
    "Mammalia":  [359],
    "Amphibia":  [131],
    "Reptilia":  [11592253, 11418114, 11493978],
}


def fetch(params, retries=5):
    """带重试的单页抓取 / one page with retries"""
    url = BASE + "?" + urllib.parse.urlencode(params)
    for a in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if a == retries - 1:
                print(f"    FAIL {params.get('offset')}: {e}", flush=True)
                return None
            time.sleep(2 * (a + 1))
    return None


def count_for(params):
    d = fetch({**params, "limit": 0})
    return d["count"] if d else 0


def harvest(base_params, label):
    """把一个查询的全部记录抓下来 / harvest every record of one query"""
    n = count_for(base_params)
    if n == 0:
        return []
    if n > MAX_OFFSET:
        # 按纬度带切分，而不是按年份。GBIF 没有"年份为空"的过滤条件，
        # 而中国两栖记录中有 63% 没有年份，按年切分会整体遗漏这些记录。
        # 纬度对每条有坐标的记录都存在，因此切分是完备的。
        # Split by latitude, not year: GBIF offers no "year is null" filter and
        # 63% of Chinese amphibian records carry no year, so a year split
        # silently drops them. Every georeferenced record has a latitude.
        print(f"    {label}: {n:,} exceeds the offset cap, splitting by latitude",
              flush=True)
        rows = []
        edges = [round(v, 2) for v in
                 [-90 + i * (180 / 72) for i in range(73)]]
        for lo, hi in zip(edges[:-1], edges[1:]):
            sub = {**base_params, "decimalLatitude": f"{lo},{hi}"}
            c = count_for(sub)
            if c == 0:
                continue
            if c > MAX_OFFSET:
                # 再对半切一次 / halve once more
                mid = round((lo + hi) / 2, 3)
                for a, b in ((lo, mid), (mid, hi)):
                    rows += harvest({**base_params,
                                     "decimalLatitude": f"{a},{b}"},
                                    f"{label} lat {a}-{b}")
            else:
                rows += harvest(sub, f"{label} lat {lo}-{hi}")
        return rows

    offsets = list(range(0, min(n, MAX_OFFSET), PAGE))
    rows = []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(fetch, {**base_params, "limit": PAGE, "offset": o}): o
                for o in offsets}
        for f in as_completed(futs):
            d = f.result()
            if not d:
                continue
            for r in d.get("results", []):
                if r.get("species") and r.get("decimalLongitude") is not None:
                    rows.append([r.get(k) for k in FIELDS])
    print(f"    {label}: {len(rows):,} of {n:,} retrieved", flush=True)
    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for cls, keys in GROUPS.items():
        out = os.path.join(OUT_DIR, f"gbif_{cls}_china.csv")
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            print(f"{cls}: already present, skipping", flush=True)
            continue
        print(f"{cls}: starting", flush=True)
        allrows = []
        for k in keys:
            p = {"country": "CN", "hasCoordinate": "true",
                 "hasGeospatialIssue": "false", "taxonKey": str(k)}
            allrows += harvest(p, f"{cls}/{k}")
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(FIELDS)
            w.writerows(allrows)
        print(f"{cls}: wrote {len(allrows):,} rows -> {os.path.basename(out)}",
              flush=True)


if __name__ == "__main__":
    main()
