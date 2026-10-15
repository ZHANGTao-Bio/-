#!/bin/bash
# ===========================================================================
# 脚本名称: 00_check_vcf.sh
# 功能描述: VCF 完整性校验 + BGZF 检查 + tabix 索引（不含分析）
#
# 使用方法:
#   登录节点/计算节点:  bash 00_check_vcf.sh
#   提交到集群:         bash run_all.sh check
# ===========================================================================

set -e
trap 'echo "[错误] 第 $LINENO 行失败"; exit 1' ERR

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SDIR}/00_config.sh"

log "VCF: ${VCF}"

# 1. 大小
log "文件大小: $(stat -c '%s' "${VCF}") 字节"

# 2. 完整性校验（9GB 约 1~3 分钟）
log "gzip 完整性校验..."
gzip -t "${VCF}" && log "VCF gzip 完整性 OK" || { log "[错误] VCF 损坏"; exit 1; }

# 3. BGZF 格式检查
MAGIC=$(head -c 4 "${VCF}" | od -An -tx1 | tr -d ' \n')
if [ "${MAGIC}" = "1f8b0804" ]; then
    log "VCF 为 BGZF 格式，符合要求"
else
    log "[错误] VCF 不是 BGZF 格式（普通 gzip），需用 bgzip 重新压缩"; exit 1
fi

# 4. 建 tabix 索引
if [ ! -f "${VCF}.tbi" ]; then
    log "建 tabix 索引..."
    bcftools index -f -t "${VCF}"
else
    log "索引已存在，跳过"
fi
ls -lh "${VCF}.tbi"

# 5. 样本数与染色体统计
log "样本数: $(bcftools query -l "${VCF}" | wc -l)"
log "染色体/contig 数: $(bcftools index -s "${VCF}" | wc -l)"

log "===== 完成 ====="
exit 0
