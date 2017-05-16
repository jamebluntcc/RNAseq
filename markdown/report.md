
## mapping

## 比对分析

数据质控合格的样品，我们使用 STAR 软件，将测序的 reads 比对到基因组，比对得到的包含 reads 在基因组位置信息的 bam 文件将作为后续 SNP，可变剪切等分析的重要输入文件。

各样品测序 reads mapping 到参考基因组的情况如下图所示。

Dobin A, Davis C A, Schlesinger F, et al. Bioinformatics, 2013, 29(1): 15-21.

## SNP 分析

根据比对分析的结果，我们使用 GATK 软件对各样品转录区域进行 SNP 分析。分析结果如下表所示（只展示部分结果，完整结果见结果文件）。

我们使用 bcftools 对分析结果进行统计，统计结果如下表所示：
