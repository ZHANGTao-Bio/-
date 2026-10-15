# EasyHap 单倍型分析流程（以 JMJ13a 区间为例）

> 工具：EasyHap v1.4.0（区域单倍型分析、群体比较、性状关联、LD 分析与可视化）
> 参考：EasyHap Linux 版使用教程（微信公众号）；官方 GitHub：https://github.com/Guangqi-He/EasyHap

## 1. 分析原理

EasyHap 从 phased VCF/VCF.GZ/BCF 出发，围绕指定基因/基因组区间，一步完成：

```
变异提取 → 单倍型构建 → 单倍型聚类 → 群体比较 → 性状关联 → LD 分析 → 可视化
```

支持两种模式（v1.4.0）：
- `--mode genotype`（默认，= 旧版 `inbred`）：不依赖 phasing 的多位点基因型类别，适合大豆等自交/纯合材料；
- `--mode copy`（= 旧版 `hybrid`）：需要 phase 信息，重建单条染色体拷贝的单倍型。

## 2. 目录结构

```
~/01.project/haplo_analysis/
├── 01.rawdata/          # 原始输入（VCF / GTF / traits）
│   ├── cohort.filtered.vcf.gz      # 9GB，360 样本，GATK 联合基因分型（BGZF）
│   ├── cohort.filtered.vcf.gz.tbi  # tabix 索引（bcftools index 生成）
│   ├── gene.a4.gtf                 # Phytozome Wm82.a4.v1 注释（Gm01 命名）
│   ├── traits.txt                  # 原始性状表（240 行，含 12 个重复样本 ID）
│   └── traits.txt.clean            # 去重合并后的性状表（脚本自动生成，228 行）
├── 02.inputs/           # 分析输入表
│   ├── genes.txt                    # 候选基因 ID 列表（每行一个）
│   ├── regions.tsv                  # 生成的区间（10 个：5 基因 × body/promoter）
│   └── regions.gene_map.tsv         # 区间 <-> 基因/类型 对照表
├── 03.scripts/          # 脚本
│   ├── 00_config.sh                 # 主配置（路径、参数）
│   ├── 00_check_vcf.sh              # VCF 校验+索引
│   ├── 01_run_easyhap.sh            # 完整分析流程
│   ├── generate_regions.py          # 从 GTF + 基因列表生成 regions.tsv
│   └── run_all.sh                   # csub 提交主控
└── 04.results/          # 输出结果
    ├── JMJ13a_test/                 # 第一次测试结果（已跑通）
    └── EasyHap_results/             # 批量分析结果（候选基因）
```

## 3. 环境准备（已完成）

```bash
export PATH=$HOME/easyhap_env/bin:$PATH
easyhap --version    # EasyHap 1.4.0
bcftools --version   # bcftools 1.24
```

## 4. 输入文件格式

### 4.1 VCF（必需）

- **BGZF 压缩**（不是普通 gzip）+ tabix 索引（随机访问必需）。检查 BGZF 魔数：`head -c 4 x.vcf.gz | od -An -tx1` 应显示 `1f 8b 08 04`（普通 gzip 是 `1f 8b 08 00`）。
- 大文件务必先 `gzip -t` 做完整性校验（本次曾发现原文件 42% 处 CRC 损坏）。
- 染色体命名必须与 GTF 一致（本例均为 `Gm01`…`Gm20` + scaffold）。
- 样本名必须与 traits/group 文件一致。

### 4.2 regions.tsv（批量分析区间，TAB 分隔、无表头）

```
Gm01	1428163	1430162
```

- 坐标 1-based、包含两端；
- 单区间也可直接用 `--region Gm01:1428163-1430162`；
- 变异数少于 `--min-variants`（默认 2）的区间会被自动跳过。

### 4.3 traits.txt（可选，性状关联）

```
Accession	Plant_height	Shoot_fresh_weight	Root_dry_weight
23NQFD106	0.865501829	1.324754599	0.805605787
...
```

- 首列必须是样本/accession ID，与 VCF 一致；
- **样本 ID 必须唯一**：本次发现 12 个样本各出现两次（两次测量），EasyHap 会报 `empty or duplicate sample IDs`，需先合并（取平均）。脚本 `01_run_easyhap.sh` 已内置该步骤，自动生成 `traits.txt.clean`。

### 4.4 sample_group.tsv（可选，群体比较；本例第一次测试暂不加）

```
Sample01	Cultivar
Sample02	Landrace
```

不提供时所有样本归入 `All` 组。

### 4.5 gene.a4.gtf（可选，基因结构图 + 变异注释筛选）

- GTF/GFF3 均可；染色体命名与 VCF 一致；
- 配合 `--gene-feature all|exon|intron|cds|utr` 可只取某类 feature 的变异。

## 5. VCF 预处理（重要）

```bash
cd ~/01.project/haplo_analysis/01.rawdata
export PATH=$HOME/easyhap_env/bin:$PATH

# 1) 完整性校验（大文件务必先做，损坏会导致后续全失败）
gzip -t cohort.filtered.vcf.gz && echo "VCF OK" || echo "VCF 损坏"

# 2) BGZF 格式检查（前 4 字节应为 1f 8b 08 04）
head -c 4 cohort.filtered.vcf.gz | od -An -tx1

# 3) 建 tabix 索引
bcftools index -f -t cohort.filtered.vcf.gz

# 4) 用索引随机读取目标区间验证
bcftools query -r Gm01:1428163-1430162 -f '%POS\t%REF\t%ALT\n' cohort.filtered.vcf.gz
```

> 大文件校验/索引建议用 `csub` 提交到计算节点（见《集群任务提交与使用说明.md》），或直接用 `bash 03.scripts/00_check_vcf.sh`。

## 6. 运行命令（JMJ13a 区间示例）

**方式一：直接用脚本（推荐）**

```bash
cd ~/01.project/haplo_analysis/03.scripts
bash 01_run_easyhap.sh              # 登录节点直接跑（单/少量区间）
# 或提交到集群：
bash run_all.sh easyhap             # csub 提交
```

**方式二：手动命令（等价）**

```bash
HA=/share/org/YZWL/yzwl_zhangtao/01.project/haplo_analysis
export PATH=$HOME/easyhap_env/bin:$PATH
cd $HA

easyhap analyze \
  --vcf 01.rawdata/cohort.filtered.vcf.gz \
  --region-file 02.inputs/regions.tsv \
  --traits 01.rawdata/traits.txt.clean \
  --gff 01.rawdata/gene.a4.gtf \
  --mode genotype \
  --max-missing 0.2 \
  --min-maf 0.05 \
  --threads 4 \
  --vcf-backend bcftools \
  --plot --plot-format pdf,png --cell-text auto \
  --outdir 04.results/EasyHap_results
```

- 单区间测试直接在登录节点跑即可（实测 48 秒）；
- 大批量候选基因（regions.tsv 多行）用 `run_all.sh easyhap` 提交到计算节点，并适当调大 `THREADS`。

## 7. 输出文件说明

| 文件 | 内容 |
|---|---|
| `*.HapSummary.tsv` | 每个 haplotype/cluster/变异状态及包含的样本（核心汇总） |
| `*.HapGroup.tsv` | 每个样本对应的 haplotype、cluster、group、trait |
| `*.HaplotypeFrequency.tsv` | 各群体中的单倍型频率 |
| `*.HaplotypeDiversity.tsv` | 各群体单倍型多样性 Hd |
| `*.PrivateHaplotypes.tsv` | 群体特有单倍型 |
| `*.TraitHaplotypeSummary.tsv` | 各单倍型的性状均值/效应 |
| `*.TraitAssociationTests.tsv` | Kruskal–Wallis、Mann–Whitney U、P 值与 BH 校正 |
| `*.LD_r2_matrix.tsv` | 区域内变异间 LD r² 矩阵 |
| `*.LD_allele_r2_matrix.tsv` | 等位基因水平的 LD r² 矩阵 |
| `*.VariantQCFilter.tsv` | 被 QC 过滤的变异记录 |
| `*.ProcessedVariants.tsv` | 处理后保留的变异及等位状态 |
| `*.SuperiorHaplotypeCandidates.tsv` | 优势单倍型候选 |
| `*.Haplotype.fa / .phy / .nex` | 单倍型序列导出（FASTA/PHYLIP/NEXUS） |
| `EasyHap.log` | 每个区间的完成/跳过/失败状态 |
| `EasyHap.RegionStatus.tsv` | 各区间状态汇总 |
| `*.pdf / *.png / *.svg` | 单倍型热图、群体饼图/堆叠柱图、性状箱线图、LD 热图、单倍型网络图等 |

## 8. 关键参数速查

| 参数 | 说明 | 本例取值 |
|---|---|---|
| `--mode` | genotype（自交）/ copy（需 phasing） | genotype |
| `--vcf-backend` | auto/bcftools/cyvcf2/pysam/plain | bcftools |
| `--max-missing` | 位点最大缺失率（0~1） | 0.2 |
| `--min-maf` | 最小次等位基因频率 | 0.05 |
| `--min-variants` | 少于该变异数则跳过区间 | 默认 2 |
| `--gene-feature` | all/exon/intron/cds/utr | all |
| `--cluster-threshold` | 单倍型聚类距离阈值 | 默认 0.15 |
| `--plot-format` | 图格式 pdf/svg/png | pdf,png |
| `--cell-text` | 热图碱基标注 auto/always/never | auto |
| `--no-ld` | 关闭 LD 分析 | 不加 |

## 9. 常见问题

1. **`Trait file ... has empty or duplicate sample IDs`**：traits 表里有重复样本 ID，需合并（脚本已内置取平均合并）。
2. **group/traits 里的样本没识别**：先用 `bcftools query -l cohort.filtered.vcf.gz` 核对样本名完全一致。
3. **某区间没结果**：查 `EasyHap.log`，多半是过滤后变异数少于 `--min-variants`，可放宽 `--min-maf`/`--max-missing`。
4. **没有 Gene+Haplotype 图**：确认提供了 `--gff`、VCF 与 GTF 染色体命名一致、region 覆盖目标 gene、GTF 含可识别的 gene/transcript/CDS/UTR。
5. **`bcftools index` 报 `Inflate operation failed`**：VCF 不是 BGZF 或文件损坏，先 `gzip -t` 校验、`head -c 4 | od` 查 BGZF 魔数。
6. **几十 GB 的大 VCF**：务必 `bcftools index` 建索引，EasyHap 才能随机访问加速。

## 11. 候选基因批量分析

**当前候选基因（5 个，每个做 body + promoter 两个区间）：**

| 基因 | 链 | 基因本体（body） | 启动子（ATG 上游 2000bp） |
|---|---|---|---|
| Glyma_13G107900 | - | Gm13:21251787-21256409 | Gm13:21253545-21255544 |
| Glyma_12G221500 | - | Gm12:39521788-39523926 | Gm12:39523792-39525791 |
| Glyma_04G208300 | - | Gm04:46881791-46883682 | Gm04:46883520-46885519 |
| Glyma_01G014700 | + | Gm01:1430164-1436702 | Gm01:1428587-1430586 |
| Glyma_09G207400 | - | Gm09:44060006-44066534 | Gm09:44065981-44067980 |

> 注：Glyma_01G014700 的启动子（1428587-1430586）与之前测试的 JMJ13a 区间（1428163-1430162）部分重叠。
> 对于 4 个负链基因，「ATG 上游 2000bp」会包含部分 5'UTR 并越过 TSS 继续向上游延伸；若想改成「TSS 上游 2000bp」，用 `generate_regions.py` 调整即可。

**生成 regions.tsv（本地已生成，或集群上重跑）：**

```bash
cd ~/01.project/haplo_analysis
export PATH=$HOME/easyhap_env/bin:$PATH

# 候选基因列表（每行一个 gene_name）
cat > 02.inputs/genes.txt << 'EOF'
Glyma_13G107900
Glyma_12G221500
Glyma_04G208300
Glyma_01G014700
Glyma_09G207400
EOF

# 生成 regions.tsv（基因本体 + ATG 上游 2000bp 启动子）
python 03.scripts/generate_regions.py \
  --gtf 01.rawdata/gene.a4.gtf \
  --genes 02.inputs/genes.txt \
  --upstream 2000 \
  --out 02.inputs/regions.tsv \
  --map 02.inputs/regions.gene_map.tsv
```

**运行前检查各区间的变异数（避免某些区间因变异太少被跳过）：**

```bash
cd ~/01.project/haplo_analysis
export PATH=$HOME/easyhap_env/bin:$PATH
while IFS=$'\t' read -r c s e; do
  [ -z "$c" ] && continue
  n=$(bcftools query -r "${c}:${s}-${e}" -f '%POS\n' 01.rawdata/cohort.filtered.vcf.gz | wc -l)
  printf "%-28s %5d 个变异\n" "${c}:${s}-${e}" "$n"
done < 02.inputs/regions.tsv
```

**批量运行（提交到集群）：**

```bash
cd ~/01.project/haplo_analysis/03.scripts
bash run_all.sh easyhap      # csub 提交到 c02 队列
```

## 10. 已完成 & 下一步

**已完成（JMJ13a 区间验证）**
- [x] 建 conda 环境 `~/easyhap_env`（Python 3.10.22）+ 装 bcftools/tabix/EasyHap 1.4.0
- [x] 从原服务器重新拷贝完整 VCF（原文件 42% 处 CRC 损坏）
- [x] `gzip -t` 校验 + BGZF 检查 + `bcftools index` 建索引
- [x] traits 表去重合并（240 → 228 样本）
- [x] 跑通 `Gm01:1428163-1430162`（10 变异位点，48 种单倍型，Hd=0.69，性状关联不显著）

**下一步（进行中/可选）**
- [x] 候选基因列表（5 个基因，见第 11 节）
- [x] 生成 regions.tsv（5 基因 × body/promoter = 10 区间）
- [ ] 上传 regions.tsv 到集群，检查各区间变异数
- [ ] 批量运行 EasyHap（`bash run_all.sh easyhap`）
- [ ] 解读各基因单倍型/性状关联/优势单倍型
- [ ] （可选）准备群体分组 `sample_group.tsv`（按种质「来源地」）做群体比较
