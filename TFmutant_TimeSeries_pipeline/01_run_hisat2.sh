#!/bin/bash
# script_name: 01_run_hisat2.sh

# ================= 防御性设置 =================
# set -e: 遇到任何错误立即退出
# set -o pipefail: 管道中只要有一环失败（如 samtools 挂了），整个命令就算失败
set -e
set -o pipefail

# 加载全局配置
# ========================================================
# 自动获取脚本所在目录并加载配置
# ========================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "❌ Error: project_config.sh not found in ${SCRIPT_DIR}"
    exit 1
fi
# 输出目录
OUT_DIR="${BASE_DIR}/01_assembly"
mkdir -p "${OUT_DIR}"

echo ">>> 启动 Pipeline 1: HISAT2 比对"
echo ">>> 参考基因组: ${INDEX_PREFIX}"
echo ">>> 并行任务数: ${PARALLEL_JOBS}"

# 定义比对函数 (供 xargs 调用)
run_alignment() {
    # 在子 Shell/函数中同样开启严格模式，防止管道错误被忽略
    set -e
    set -o pipefail

    local sample_id=$1
    local fq1=$2
    local fq2=$3
    local out_dir=$4
    local index=$5
    local hisat2=$6
    local samtools=$7
    local threads=$8
    
    echo "[$(date +'%T')] 开始处理: ${sample_id}"
    
    local bam_file="${out_dir}/${sample_id}.sort.bam"
    local log_file="${out_dir}/${sample_id}.hisat2.log"
    
    # ---------------------------------------------------------
    # 修复点 1: 修正变量名 $bam -> $bam_file
    # 修复点 2: 增加 || exit 255 确保错误能被 xargs 捕获
    # ---------------------------------------------------------
    $hisat2 -p $threads --dta -x "$index" -1 "$fq1" -2 "$fq2" 2> "$log_file" | \
    $samtools sort -@ 2 -o "$bam_file" - 
    
    # 建立索引
    $samtools index "$bam_file"
    
    echo "[$(date +'%T')] 完成: ${sample_id}"
}

export -f run_alignment

# 从 metadata 读取样本信息并构建任务列表
# 跳过标题行，提取 SampleID, FQ1, FQ2
# ---------------------------------------------------------
# 修复点 3: 如果 xargs 中的任何一个任务失败，脚本必须捕获并在主线程报错
# ---------------------------------------------------------
tail -n +2 "${METADATA_FILE}" | cut -f 1,5,6 | \
xargs -P ${PARALLEL_JOBS} -n 3 -I {} bash -c '
    # 解析参数
    args=({})
    sample_id=${args[0]}
    fq1=${args[1]}
    fq2=${args[2]}
    
    # 调用函数
    run_alignment "$sample_id" "$fq1" "$fq2" \
    "'"${OUT_DIR}"'" "'"${INDEX_PREFIX}"'" "'"${HISAT2_BIN}"'" \
    "'"${SAMTOOLS_BIN}"'" "'"${HISAT2_THREADS}"'"
'

# ---------------------------------------------------------
# 修复点 4: 显式检查 xargs 的退出状态
# ---------------------------------------------------------
if [ $? -ne 0 ]; then
    echo "❌ [FATAL ERROR] Pipeline 1 在并行处理阶段检测到失败！"
    echo "   请检查 ${OUT_DIR} 下的日志文件寻找具体原因。"
    # 强制以非 0 状态退出，通知 00_run_all_pipelines.sh 停止后续步骤
    exit 1
fi

# 生成简单的统计报告
echo "Sample ID,Total Reads,Alignment Rate" > "${OUT_DIR}/alignment_stats.csv"
for log in ${OUT_DIR}/*.hisat2.log; do
    # 确保日志文件存在且不为空
    if [ -s "$log" ]; then
        sname=$(basename "$log" .hisat2.log)
        total=$(grep "reads; of these:" "$log" | awk '{print $1}')
        rate=$(grep "overall alignment rate" "$log" | awk '{print $1}')
        echo "${sname},${total},${rate}" >> "${OUT_DIR}/alignment_stats.csv"
    fi
done

echo ">>> Pipeline 1 完成。输出位于: ${OUT_DIR}"