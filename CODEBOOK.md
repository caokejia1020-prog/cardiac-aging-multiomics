# 输出字段说明

## 分组映射

源文件沿用 L/M/H 年龄编码；正式展示统一为 E/M/L：源 L = Early，源 M = Middle，源 H = Late。源组 FL、FM、FH 分别对应 Female E、M、L；ML、MM、MH 分别对应 Male E、M、L。

## 差异表核心字段

- `Feature_ID`：分析层内唯一特征 ID；蛋白为 Accession，代谢特征为 `peak_name_IonMode`。
- `First_group`、`Second_group`：展示方向中的前组与后组。
- `N_first`、`N_second`：实际进入该特征模型的有限观测数。
- `Mean_first_log2`、`Mean_second_log2`：两组 log2 平均值。
- `logFC`：limma contrast `First - Second`；等于两组 log2 平均值之差。
- `FC_first_over_second`：`2^logFC`。
- `AveExpr`、`t`、`B`：limma 标准输出。
- `P.Value`：limma moderated t test 双侧 P 值。
- `adj.P.Val`：同一层级、同一比较全部可检验特征内的 BH-FDR。
- `Testable`：是否满足每组最少有限观测数。
- `Exclusion_reason`：不可检验原因。
- `Formal_call`：正式 BH-FDR 和 fold-change 双阈值判定，取 Up、Down 或 Not_significant。

## UpSet 输入字段

- 每行一个正式差异 Feature_ID。
- 三个逻辑列对应同一性别的 E–M、M–L 和 E–L。
- `Pattern` 依次编码三列 membership；例如 101 表示属于 E–M 和 E–L，不属于 M–L。

## 富集字段

- `Query_hits`、`Query_size`：正式差异 query 中命中该 term 的数量及 query 总数。
- `Background_term_size`、`Background_size`：研究背景内该 term 的成员数及背景总数。
- `P.Value`：右尾超几何检验。
- `adj.P.Val`：完整 eligible term universe 内的 BH-FDR，包括零命中 term。
- `Matching_IDs`：命中的基因、ENTREZID 或 KEGG compound ID。
- `Formal_enrichment`：配置的富集 FDR 阈值判定。

## 网络字段

- `Complete_pairs`：按动物 ID、性别和年龄组匹配后的完整样本对数。
- `Correlation`：Pearson r。
- `P.Value`：相关检验 P 值。
- `adj.P.Val`：相同代谢层级、性别和目标蛋白内全部候选相关的 BH-FDR。
- `Formal_edge`：同时满足 BH-FDR 与绝对相关阈值的正式边。
