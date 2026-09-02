# 研究输入文件

本包已包含以下九个处理后定量数据文件。请保持文件名和目录结构不变。`INPUT_MD5SUMS.tsv` 用于核对文件是否与本包使用的输入一致，运行 `validate_setup.R` 时会自动检查。

蛋白组文件放在 `input/protein_pairwise/`：

```text
FL_vs_FM.result.xlsx
FM_vs_FH.result.xlsx
FL_vs_FH.result.xlsx
ML_vs_MM.result.xlsx
MM_vs_MH.result.xlsx
ML_vs_MH.result.xlsx
```

每个文件读取 `ALL` 工作表，必须包含 `Accession`、`PG.Genes`、`PG.ProteinDescriptions` 及12个样本定量列。定量值为已处理的 log2 数值；跨文件的同一样本和蛋白值必须一致。

代谢组文件放在 `input/metabolite/`：

```text
heart_data_unfiltered.xlsx
serum_data_unfiltered.xlsx
urine_data_unfiltered.xlsx
```

每个文件读取 `FL-FM`、`FM-FH`、`FL-FH`、`ML-MM`、`MM-MH`、`ML-MH` 六个工作表，必须包含 `peak_name`、`mz`、`rt`、`name`、`id_kegg`、`Ion mode` 及每张表的12个样本定量列。

样本名必须保留原有组别和动物编号，以便重建每层36个样本，并按动物编号、性别和年龄组进行蛋白与代谢物配对。输入工作簿内的已有 P 值、VIP、FC 和差异标签不作为新模型的统计结果。

这些工作簿为处理后的蛋白和代谢特征定量数据，不是原始质谱文件。蛋白定量行与代谢特征定量行分别用于对应数据层的分析，不能直接等同于肽段或修饰位点。
