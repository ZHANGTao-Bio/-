#!/bin/bash
# ===========================================================================
# 脚本名称: run_all.sh
# 功能描述: EasyHap 主控脚本 —— 通过 csub 提交到集群计算节点
#
# 使用方法:
#   bash run_all.sh            # 显示帮助
#   bash run_all.sh easyhap    # 提交 EasyHap 分析到集群
#   bash run_all.sh check      # 提交 VCF 校验+索引到集群
#
# 说明:
#   - 输出/日志必须写在 GPFS 共享盘（不要写 /tmp，那是节点本地盘）
#   - 队列 c02 排队少，适合轻量任务；大批量用 c01 或 c01_yzwl
# ===========================================================================

set -e

SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SDIR}/logs"
mkdir -p "${LOG_DIR}"

QUEUE="c02"       # 提交队列
THREADS=4         # 线程数
JOBNAME="easyhap"

show_help() {
    echo ""
    echo "EasyHap 主控脚本"
    echo "================="
    echo "用法: bash run_all.sh <步骤>"
    echo ""
    echo "  easyhap   提交完整 EasyHap 分析（校验+索引+去重+分析）"
    echo "  check     只提交 VCF 校验+索引"
    echo ""
    echo "提交示例（等价手动命令）:"
    echo "  csub -q ${QUEUE} -n ${THREADS} -J ${JOBNAME} \\"
    echo "    -o ${LOG_DIR}/${JOBNAME}.out -e ${LOG_DIR}/${JOBNAME}.err << 'EOF'"
    echo "  export PATH=\$HOME/easyhap_env/bin:\$PATH"
    echo "  bash ${SDIR}/01_run_easyhap.sh"
    echo "  EOF"
    echo ""
}

# ---------- 通过 csub 提交到集群 ----------
csub_job() {
    local jobname=$1
    local script=$2
    shift 2

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 提交 ${jobname} 到集群 (队列:${QUEUE}, CPU:${THREADS})"
    echo "  脚本: ${script}"
    echo "  日志: ${LOG_DIR}/${jobname}.out / ${jobname}.err"

    csub -q "${QUEUE}" -n "${THREADS}" -J "${jobname}" \
        -o "${LOG_DIR}/${jobname}.out" \
        -e "${LOG_DIR}/${jobname}.err" \
        << EOF
export PATH=\$HOME/easyhap_env/bin:\$PATH
bash "${script}" "$@"
EOF
}

# ---------- 主入口 ----------
STEP="${1:-help}"
shift 2>/dev/null || true

case "${STEP}" in
    help|--help|-h)
        show_help
        ;;
    easyhap)
        csub_job "easyhap" "${SDIR}/01_run_easyhap.sh" "$@"
        ;;
    check)
        csub_job "vcf_check" "${SDIR}/00_check_vcf.sh" "$@"
        ;;
    *)
        echo "[错误] 未知步骤: ${STEP}"
        show_help
        exit 1
        ;;
esac

echo ""
echo "提交后可查看:"
echo "  状态: cjobs <返回的 jobid>"
echo "  输出: cat ${LOG_DIR}/${STEP}.out"
exit 0
