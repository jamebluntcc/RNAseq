% 结果文件说明


# fastqc 目录
>存放 fastqc 的分析结果以及对结果整理的图表，包含 fastqc_results, reads_quality_plot, gc_plot 三个子目录以及 fastqc 分析一览表 (fastqc_general_stats.txt, fastqc_general_stats.xls)。

## fastqc_results 子目录
>存放样品 fastqc 分析结果。

### 文件列表
- sample.clean_fastqc : fastqc 结果文件
- sample.clean_fastqc.html : fastqc 网页报告

## reads_quality_plot 子目录
>存放样品 reads 质量值分布图。

### 文件列表
- sample.reads_quality.png : 质量值分布图 pdf 格式
- sample.reads_quality.pdf : 质量值分布图 png 格式

### 图例
>x轴为质量分数，y轴为该质量值 reads 比例。


## gc_plot 子目录
>存放样品 reads GC 分布图。

### 文件列表
- sample.gc.png : GC 分布图 pdf 格式
- sample.gc.pdf : GC 分布图 png 格式

### 图例
>x轴代表4种碱基以及N在 read1, read2 的位置 (因为将read1, read2 放在同一坐标轴展示，因此碱基在 read2 的位置为 (x-150)bp)，y轴代表4种碱基以及N的比例。

## fastqc 分析一览表 (fastqc_general_stats )
>整合了 fastqc 的主要结果

### 表头说明

|表头     | 说明                |
|:--------|:--------------------|
|Sample_ID|样品名称|
|Reads_number(M)|reads 数目，单位为兆|
|Reads_length(bp)|reads 长度，单位为bp|
|Data_size(G)|数据量，单位为G|
|Q30(%)|质量值大于30的 reads 的比例，百分数|
|GC(%)|reads GC 含量，百分数|

# quantification 目录
>存放定量分析以及差异分析的结果。包括 kallisto, differential_analysis, expression_summary 三个子目录。

## kallisto 子目录
>存放 kallisto 定量的原始结果，包括 tsv, h5 两种格式。

### 文件列表
- sample/abundance.tsv : kallisto 定量结果 tsv 格式
- sample/abundance.h5 : kallisto 定量结果 hdf5 格式

### 表头说明
|表头       |说明                |
|:----------|:-------------------|
|target_id|转录本 ID|
|length|转录本序列长度|
|eff_length|Effective length, read pair 在转录本中可能的起始位点的总和|
|est_counts|kallisto 定量 read count 值|
|tpm|Transcripts per million (TPM), kallisto 定量 tpm 值|

## differential_analysis 子目录
>存放 edgeR 分析各比较组差异分析结果，包括差异分析

### 文件列表

- group1_vs_group2/group1_vs_group2.edgeR.DE_results.txt
- group1_vs_group2.edgeR.DE_results.xlsx
- group1_vs_group2.ALL.edgeR.DE_results.diffgenes.txt
- group1_vs_group2.group1-UP.edgeR.DE_results.txt
- group1_vs_group2.group1-UP.edgeR.DE_results.xlsx
- group1_vs_group2.group1-UP.edgeR.DE_results.diffgenes.txt
- group1_vs_group2.group2-UP.edgeR.DE_results.txt
- group1_vs_group2.group2-UP.edgeR.DE_results.xlsx
- group1_vs_group2.group2-UP.edgeR.DE_results.diffgenes.txt


|表头       |说明                |
|:----------|:-------------------|
|||
