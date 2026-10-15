#!/bin/bash
# ===========================================================================
# 脚本名称: 01_run_easyhap.sh
# 功能描述: EasyHap 单倍型分析完整流程
#   1) VCF 完整性校验 + BGZF 检查 + tabix 索引
#   2) traits 表去重合并（同一样本多次测量取平均）
#   3) 运行 easyhap analyze
#
# 使用方法:
#   登录节点直接跑（单/少量区间）:  bash 01_run_easyhap.sh
#   批量任务提交到集群:            bash run_all.sh
#
# 依赖: easyhap / bcftools / tabix（均已在 ~/easyhap_env）
# ===========================================================================

set -e
trap 'echo "[错误] 第 $LINENO 行失败"; exit 1' ERR

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SDIR}/00_config.sh"

mkdir -p "${INPUT_DIR}" "${OUTDIR}"

# ---------- 1. VCF 预处理 ----------
log "步骤 1/3: VCF 校验与索引"
if [ ! -f "${VCF}" ]; then
    log "[错误] 找不到 VCF: ${VCF}"; exit 1
fi

# 1.1 完整性校验（大文件必须，损坏会导致后续全失败）
gzip -t "${VCF}" && log "VCF gzip 完整性 OK" || { log "[错误] VCF 损坏，请重新拷贝"; exit 1; }

# 1.2 BGZF 格式检查（EasyHap 随机访问要求 BGZF，普通 gzip 无法建 tabix 索引）
MAGIC=$(head -c 4 "${VCF}" | od -An -tx1 | tr -d ' \n')
if [ "${MAGIC}" = "1f8b0804" ]; then
    log "VCF 为 BGZF 格式，符合要求"
else
    log "[错误] VCF 不是 BGZF 格式（普通 gzip），需用 bgzip 重新压缩"; exit 1
fi

# 1.3 建 tabix 索引（已存在则跳过）
if [ ! -f "${VCF}.tbi" ]; then
    log "建 tabix 索引..."
    bcftools index -f -t "${VCF}"
else
    log "索引已存在，跳过"
fi

# ---------- 2. traits 表去重合并 ----------
log "步骤 2/3: traits 表去重（同一样本多次测量取平均）"
if [ -f "${TRAITS_RAW}" ]; then
    # 输出表头
    awk -F'\t' 'NR==1{print; exit}' "${TRAITS_RAW}" > "${TRAITS}"
    # 按样本 ID 合并，对数值列求平均
    awk -F'\t' '
        NR==1 {ncol=NF; next}
        {
            id=$1
            cnt[id]++
            for(i=2;i<=NF;i++) sum[id,i]+=$i
        }
        END {
            for(id in cnt) {
                printf "%s", id
                for(i=2;i<=ncol;i++) printf "\t%.9g", sum[id,i]/cnt[id]
                printf "\n"
            }
        }' "${TRAITS_RAW}" | sort >> "${TRAITS}"
    log "traits 样本行数: $(tail -n +2 "${TRAITS}" | wc -l)（原始 $(tail -n +2 "${TRAITS_RAW}" | wc -l)）"
    log "重复 ID 数: $(awk -F'\t' 'NR>1{print $1}' "${TRAITS}" | sort | uniq -d | wc -l)（应为 0）"
else
    log "[警告] 未找到 traits 文件，本次不分析性状关联"
    TRAITS=""
fi

# ---------- 3. 运行 easyhap analyze ----------
log "步骤 3/3: 运行 EasyHap"

# 拼接可选参数
EXTRA=""
[ -n "${GTF}" ] && EXTRA="${EXTRA} --gff ${GTF}"
[ -n "${TRAITS}" ] && EXTRA="${EXTRA} --traits ${TRAITS}"
[ -n "${GROUP:-}" ] && EXTRA="${EXTRA} --group ${GROUP}"

# 区间参数：优先用 regions 文件，否则单区间
if [ -n "${REGIONS:-}" ] && [ -f "${REGIONS}" ]; then
    REGION_ARG="--region-file ${REGIONS}"
elif [ -n "${REGION:-}" ]; then
    REGION_ARG="--region ${REGION}"
else
    log "[错误] 未指定分析区间（REGIONS 或 REGION）"; exit 1
fi

log "运行命令:"
log "  easyhap analyze --vcf ${VCF} ${REGION_ARG} --mode ${MODE} ${EXTRA} --outdir ${OUTDIR}"

easyhap analyze \
    --vcf "${VCF}" \
    ${REGION_ARG} \
    --mode "${MODE}" \
    --max-missing "${MAX_MISSING}" \
    --min-maf "${MIN_MAF}" \
    --threads "${THREADS}" \
    --vcf-backend "${VCF_BACKEND}" \
    --plot --plot-format "${PLOT_FORMAT}" --cell-text "${CELL_TEXT}" \
    ${EXTRA} \
    --outdir "${OUTDIR}"

log "===== 完成 ====="
log "结果目录: ${OUTDIR}"
log "日志: ${OUTDIR}/EasyHap.log"
log "区间状态: ${OUTDIR}/EasyHap.RegionStatus.tsv"
exit 0
