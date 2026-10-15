#!/usr/bin/env python3
"""
从 GTF 生成 EasyHap regions.tsv（基因本体 + 启动子区间）

对每个基因输出两行区间：
  - body     : 基因本体（全基因跨度，含外显子+内含子）
  - promoter : 启动子（ATG 起始密码子上游 N bp，默认 2000 bp）

链方向约定：
  - + 链：ATG = 最小 CDS start；promoter = [ATG - N, ATG - 1]
  - - 链：ATG = 最大 CDS end；  promoter = [ATG + 1, ATG + N]

用法:
  python generate_regions.py --gtf gene.a4.gtf --genes genes.txt --upstream 2000 \
      --out regions.tsv --map regions.gene_map.tsv

  genes.txt: 每行一个基因 ID（gene_name，如 Glyma_13G107900），# 开头为注释
  输出 regions.tsv: chr<TAB>start<TAB>end（无表头，3 列，EasyHap 直接可用）
  输出 regions.gene_map.tsv: 区间 <-> 基因/类型的对照表（便于人工核对）
"""

import argparse
import re
import sys
from collections import defaultdict


def parse_gtf(gtf_path):
    """按 gene_name 索引 GTF 记录。返回 {gene_name: [rows]}"""
    genes = defaultdict(list)
    with open(gtf_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            m = re.search(r'gene_name "([^"]+)"', cols[8])
            if not m:
                continue
            genes[m.group(1)].append(
                (cols[0], cols[2], int(cols[3]), int(cols[4]), cols[6])
            )
    return genes


def gene_regions(rows, upstream):
    """返回 (body_start, body_end, prom_start, prom_end, strand, atg)"""
    chroms = {r[0] for r in rows}
    if len(chroms) != 1:
        raise ValueError(f"基因跨越多个染色体: {chroms}")
    chrom = next(iter(chroms))

    strands = {r[4] for r in rows}
    if len(strands) != 1:
        raise ValueError(f"基因链方向不一致: {strands}")
    strand = next(iter(strands))

    gene_rows = [r for r in rows if r[1] == "gene"]
    if gene_rows:
        body_start = min(r[2] for r in gene_rows)
        body_end = max(r[3] for r in gene_rows)
    else:
        body_start = min(r[2] for r in rows)
        body_end = max(r[3] for r in rows)

    cds = [r for r in rows if r[1] == "CDS"]
    if not cds:
        raise ValueError("未找到 CDS，无法确定 ATG 位置")

    if strand == "-":
        atg = max(r[3] for r in cds)          # - 链：起始密码子在最大坐标端
        prom_start, prom_end = atg + 1, atg + upstream
    else:
        atg = min(r[2] for r in cds)          # + 链：起始密码子在最小坐标端
        prom_start, prom_end = atg - upstream, atg - 1

    return chrom, body_start, body_end, prom_start, prom_end, strand, atg


def main():
    ap = argparse.ArgumentParser(description="从 GTF 生成 EasyHap regions.tsv")
    ap.add_argument("--gtf", required=True, help="GTF/GFF3 注释文件")
    ap.add_argument("--genes", required=True, help="基因 ID 列表文件（每行一个 gene_name）")
    ap.add_argument("--upstream", type=int, default=2000, help="启动子 = ATG 上游 N bp（默认 2000）")
    ap.add_argument("--out", default="regions.tsv", help="输出 regions.tsv")
    ap.add_argument("--map", default="regions.gene_map.tsv", help="输出区间-基因对照表")
    args = ap.parse_args()

    genes = parse_gtf(args.gtf)

    with open(args.genes, "r", encoding="utf-8") as f:
        targets = [l.strip() for l in f if l.strip() and not l.strip().startswith("#")]

    out_lines = []
    map_lines = ["region\tgene\ttype\tstrand\tatg\tlength_bp"]

    for t in targets:
        key = None
        for k in genes:
            if k.upper() == t.upper():
                key = k
                break
        if key is None:
            print(f"[警告] GTF 中未找到基因: {t}", file=sys.stderr)
            continue

        chrom, bs, be, ps, pe, strand, atg = gene_regions(genes[key], args.upstream)

        body_region = f"{chrom}:{bs}-{be}"
        prom_region = f"{chrom}:{ps}-{pe}"

        out_lines.append(f"{chrom}\t{bs}\t{be}")
        out_lines.append(f"{chrom}\t{ps}\t{pe}")

        map_lines.append(f"{body_region}\t{key}\tbody\t{strand}\t{atg}\t{be - bs + 1}")
        map_lines.append(f"{prom_region}\t{key}\tpromoter({args.upstream}bp)\t{strand}\t{atg}\t{pe - ps + 1}")
        print(f"[OK] {key}: body {body_region} ({be - bs + 1}bp) | promoter {prom_region}")

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + "\n")
    with open(args.map, "w", encoding="utf-8") as f:
        f.write("\n".join(map_lines) + "\n")

    print(f"\n已生成: {args.out}（{len(out_lines)} 个区间）、{args.map}")


if __name__ == "__main__":
    main()
