#!/bin/bash
# ===========================================================================
# EasyHap 单倍型分析 — 主配置文件
# ===========================================================================

# 基础路径（根据实际部署修改）
BASE_DIR="/share/org/YZWL/yzwl_zhangtao/01.project/haplo_analysis"
RAW_DIR="${BASE_DIR}/01.rawdata"      # 原始输入
INPUT_DIR="${BASE_DIR}/02.inputs"     # 分析输入表
SCRIPT_DIR="${BASE_DIR}/03.scripts"   # 脚本
OUT_DIR="${BASE_DIR}/04.results"      # 输出

# 环境（conda 环境建在 GPFS 上，计算节点可见）
export PATH="$HOME/easyhap_env/bin:$PATH"

# ---------- 输入文件 ----------
VCF="${RAW_DIR}/cohort.filtered.vcf.gz"   # 已 BGZF 压缩 + tabix 索引
GTF="${RAW_DIR}/gene.a4.gtf"              # Phytozome Wm82.a4.v1 注释（染色体命名 Gm01…）
TRAITS_RAW="${RAW_DIR}/traits.txt"        # 原始性状表（可能含重复样本 ID）
TRAITS="${RAW_DIR}/traits.txt.clean"      # 去重合并后的性状表（脚本自动生成）

# ---------- 分析区间文件 ----------
# 格式：chr<TAB>start<TAB>end（无表头，1-based 含两端）
# 单区间也可用 --region "Gm01:1428163-1430162" 直接指定
REGIONS="${INPUT_DIR}/regions.tsv"

# ---------- 可选：群体分组（目前未提供，注释掉） ----------
# GROUP="${INPUT_DIR}/sample_group.tsv"

# ---------- 输出目录 ----------
OUTDIR="${OUT_DIR}/EasyHap_results"

# ---------- 分析参数 ----------
MODE="genotype"            # genotype=相位无关多位点基因型(自交默认)；copy=需 phasing 的染色体拷贝单倍型
MAX_MISSING=0.2            # 位点最大缺失率（0~1）
MIN_MAF=0.05               # 最小次等位基因频率
THREADS=4                  # 并行线程数
VCF_BACKEND="bcftools"     # VCF 读取后端：auto/bcftools/cyvcf2/pysam/plain
PLOT_FORMAT="pdf,png"      # 输出图格式
CELL_TEXT="auto"           # 热图碱基标注：auto/always/never
