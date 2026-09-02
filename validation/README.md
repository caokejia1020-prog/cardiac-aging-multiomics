# 运行环境记录

本目录记录2026年8月28日运行本代码时使用的软件环境：R 4.4.3、Bioconductor 3.20、limma 3.62.1，操作系统为 Ubuntu 26.04。

- `sessionInfo.txt`：R、操作系统及已加载软件包的信息。
- `package_versions_actual.csv`：实际安装的软件包版本。
- `runtime_version_verification.csv`：目标版本与实际版本的逐项核对结果。

目标版本清单见 `config/package_versions_target.tsv`。重新运行时，程序会在 `outputs/reproducibility/` 记录当前使用的环境，便于比较不同运行之间的软件版本。
