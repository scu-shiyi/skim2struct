# Skim2Struct

**Skim2Struct** 是一个面向非模式植物的全流程分析工具，支持从**浅层测序数据（genome skimming）**出发，快速实现**系统发育分析**、**进化保守性计算**与**结构预测和分子对接分析**。本工具特别适用于缺乏参考基因组的系统发育与功能研究场景。

---

## 📦 安装方式

### 推荐方式：使用 Conda 安装依赖

```bash
# 创建环境（可选）
conda create -n skim2struct python=3.11
conda activate skim2struct

# 安装依赖（建议在 Skim2StructProject 目录下执行）
pip install -e .
```

⚠️ 注意：某些模块（如 `DockingModule`）依赖外部工具（AutoDock、MGLTools、Open Babel），需根据文档手动安装。

---

## 📁 项目结构

```
Skim2StructProject/
├── skim2struct/               # 主程序目录
│   ├── TreeConservationModule/
│   ├── EvoDnDsModule/
│   ├── DockingModule/
│   └── __main__.py            # CLI 主入口
├── setup.py
└── README.md
```

---

## 🚀 快速使用指南

### 一、构建系统发育与保守性热图（TreeConservationModule）

```bash
skim2struct tree \
    --aligned example/tree/aligned.fasta \
    --tree example/tree/tree.nwk \
    --output tree_heatmap.pdf
```

📌 支持自动识别保守性分数并与系统发育树对齐展示。

---

### 二、计算 dN/dS 进化速率（EvoDnDsModule）

```bash
skim2struct evodnds \
    --fasta example/evo/alignment.fasta \
    --tree example/evo/species_tree.nwk \
    --output results/evo_summary.tsv
```

📌 内部调用 HyPhy 工具，支持多种模型选择与样本批处理。

---

### 三、分子结构预测与对接（DockingModule）

```bash
skim2struct docking \
    --receptor data/UGT1.pdb \
    --ligand data/UDP-Glc.mol2 \
    --output docking_results/
```

📌 支持 AlphaFold 蛋白结构输入、UDP-Glc 等配体分子自动处理，对接区域自动定位。

---

## 🔧 外部依赖

以下工具不包含在 pip 安装中，需用户自行安装：

| 工具 | 功能 | 安装建议 |
|------|------|----------|
| [MGLTools](http://mgltools.scripps.edu/) | AutoDock 对接支持工具 | 请参考官网或使用预编译版本 |
| [AutoDock Vina](http://vina.scripps.edu/) | 分子对接核心工具 | 推荐通过 conda 安装 |
| [Open Babel](http://openbabel.org/) | 分子格式转换 | `conda install -c conda-forge openbabel` |
| [HyPhy](https://github.com/veg/hyphy) | dN/dS 分析 | 推荐源码安装或预构建版本 |

---

## 📄 引用与协议

本工具使用 [MIT License](LICENSE)。如在您的研究中使用 Skim2Struct，请引用我们的项目链接：

📌 GitHub: [https://github.com/scu-shiyi/skim2struct](https://github.com/scu-shiyi/skim2struct)
