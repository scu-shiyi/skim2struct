import pandas as pd
import numpy as np
from skim2struct.utils.TreeFunction import load_tree, compute_leaf_positions
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from Bio import Phylo
import pandas as pd
import numpy as np
from typing import Optional

def compute_universal_activity_matrix(
    csv_path: str, 
    out_csv_path: Optional[str] = None,
    activity_method: str = "ratio",
    standardize: bool = "False"
) -> pd.DataFrame:
    """
    计算通用催化活性矩阵
    
    参数:
    csv_path -- 输入CSV文件路径
    out_csv_path -- 输出CSV文件路径(可选)
    activity_method -- 活性计算方法: 
        "ratio" (默认): 产物/底物结合能比值
        "delta": 产物-底物结合能差值
        "exp_delta": exp(βΔG) 热力学指标
    standardize -- 是否按基因进行Z-score标准化
    
    返回:
    催化活性矩阵DataFrame
    """
    # 物理常数
    R_kcal = 1.987e-3      # kcal·mol⁻¹·K⁻¹
    T = 298                # K
    RT = R_kcal * T        # ≈ 0.592 kcal/mol
    
    # ---------- 数据读取与预处理 ----------
    df_raw = pd.read_csv(csv_path, header=None)
    gene_row = df_raw.iloc[0, 1:]
    role_row = df_raw.iloc[1, 1:].str.lower()
    ligand_row = df_raw.iloc[2, 1:]
    
    # 创建列名 (基因_角色_配体)
    col_names = [f"{g}_{r}_{l}" for g, r, l in zip(gene_row, role_row, ligand_row)]
    
    # 处理数据行
    df_energy = df_raw.iloc[3:].copy()
    df_energy.columns = ["Sample"] + col_names
    df_energy.set_index("Sample", inplace=True)
    
    # ---------- 整理为长表格式 ----------
    tidy = (
        df_energy.stack()
                 .reset_index()
                 .rename(columns={0: "ΔG", "level_1": "GeneRole"})
    )
    tidy[["Gene", "Role", "Ligand"]] = tidy["GeneRole"].str.split("_", n=2, expand=True)
    tidy["ΔG"] = pd.to_numeric(tidy["ΔG"], errors="coerce")  # 转换为数值
    
    # 过滤无效值(缺失值或正值表示不结合)
    tidy = tidy[(tidy["ΔG"].notna()) & (tidy["ΔG"] < 0)]
    
    # ---------- 计算催化活性指标 ----------
    results = []
    for (sample, gene), sub_df in tidy.groupby(["Sample", "Gene"]):
        prod_energy = sub_df[sub_df["Role"] == "product"]["ΔG"].values
        subs_energy = sub_df[sub_df["Role"] == "substrate"]["ΔG"].values

        # 加入严格过滤条件：任意为正值或缺失 → 跳过该样本-基因组合
        if (
            len(prod_energy) == 0 or len(subs_energy) == 0 or
            np.any(prod_energy >= 0) or np.any(subs_energy >= 0)
        ):
            continue

        mean_prod = np.mean(prod_energy)
        mean_subs = np.mean(subs_energy)

        if activity_method == "ratio":
            activity = abs(mean_prod / mean_subs)
        elif activity_method == "delta":
            activity = mean_subs - mean_prod
        elif activity_method == "exp_delta":
            delta_g = mean_subs - mean_prod
            activity = np.exp(delta_g / RT)
        else:
            raise ValueError(f"未知的activity_method: {activity_method}")

        results.append({
            "Sample": sample,
            "Gene": gene,
            "Activity": activity,
            "Mean_Substrate": mean_subs,
            "Mean_Product": mean_prod
        })
    
    # 创建结果DataFrame
    df_activity = pd.DataFrame(results)
    
    # 标准化处理
    if standardize and not df_activity.empty:
        # 按基因分组标准化
        df_activity["Activity_Z"] = df_activity.groupby("Gene")["Activity"].transform(
            lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
        )
    
    # ---------- 转换为矩阵形式 ----------
    value_col = "Activity_Z" if standardize else "Activity"
    df_matrix = df_activity.pivot(index="Sample", columns="Gene", values=value_col)
    
    # 重置索引使Sample成为列
    df_matrix = df_matrix.reset_index()
    
    # ---------- 输出结果 ----------
    if out_csv_path:
        df_matrix.to_csv(out_csv_path, index=False, float_format="%.6g")
    
    return df_matrix

def plot_clade(clade, ax, depths, max_depth, leaf_positions, tick=0.25):
    x0 = depths[clade]
    if clade.is_terminal():
        y = leaf_positions[clade.name]
        ax.hlines(y, x0, max_depth, colors='gray', linestyles='--', linewidth=0.8)
        ax.vlines(max_depth, y - tick, y + tick, colors='black', linewidth=0.8)
        return y

    ys = []
    for child in clade.clades:
        y_child = plot_clade(child, ax, depths, max_depth, leaf_positions, tick)
        ys.append(y_child)
        ax.hlines(y_child, x0, depths[child], colors='black', linestyles='-', linewidth=1)

    y_min, y_max = min(ys), max(ys)
    ax.vlines(x0, y_min, y_max, colors='black', linestyles='-', linewidth=1)
    return 0.5 * (y_min + y_max)

def draw_tree_and_heatmap(
    tree, depths, max_depth,
    leaf_positions, heatmap_df,
    picture_path: str,
    name_limit: int = 20,
    cmap: str = "RdBu_r",
    norm_label: str = "kd-based (norm)",
):
    fig = plt.figure(figsize=(25, 0.4 * len(leaf_positions) + 2), dpi=150)
    gs = GridSpec(1, 2, width_ratios=[1, 6], wspace=0.15)

    ax_tree = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1])

    plot_clade(tree.root, ax_tree, depths, max_depth, leaf_positions)
    ax_tree.set_ylim(-0.5, len(leaf_positions) - 0.5)
    ax_tree.set_xlim(0, max_depth * 1.02)
    ax_tree.invert_yaxis()
    ax_tree.axis("off")

    for name, y in leaf_positions.items():
        ax_tree.text(max_depth * 1.005, y, name[:name_limit], va='center', ha='left', fontsize=8)

    sns.heatmap(
        heatmap_df.values,
        ax=ax_heat,
        cmap=cmap,
        cbar_kws={"label": norm_label},
        linewidths=0.5,
        linecolor="white",
    )
    ax_heat.set_yticklabels([])
    ax_heat.tick_params(axis="y", which="both", labelleft=False, left=True)
    ax_heat.yaxis.set_ticks_position("left")
    ax_heat.set_yticks(np.arange(len(heatmap_df)) + 0.5)

    ax_heat.set_xticks(np.arange(len(heatmap_df.columns)) + 0.5)
    ax_heat.set_xticklabels(heatmap_df.columns, rotation=45, ha='right', fontsize=8)
    ax_heat.set_xlabel("Gene", fontsize=10)

    plt.tight_layout()
    plt.savefig(picture_path)
    plt.close(fig)



def plot(csv_path, out_csv_path, tree_path, pic_path):
    tree, depths, max_depths = load_tree(tree_path)
    leaf_position = compute_leaf_positions(tree)
    leaf_position_lower = {k.lower(): v for k, v in leaf_position.items()}
    kd_martrix = compute_universal_activity_matrix(csv_path, out_csv_path)
    df_heatmap = kd_martrix.set_index("Sample")
    df_heatmap = df_heatmap.loc[leaf_position_lower.keys()]

    draw_tree_and_heatmap(tree, depths, max_depths, leaf_position, df_heatmap, pic_path)



# def draw_tree_and_heatmap(
#     tree, leaf_order,
#     heatmap_df,
#     picture_path: str,
#     name_limit: int = 20,
#     cmap: str = "RdBu_r",
#     norm_label: str = "kd-based (norm)",
# ):
#     fig = plt.figure(figsize=(26, 0.4 * len(leaf_order) + 2), dpi=150)
#     gs = GridSpec(1, 2, width_ratios=[1.5, 6], wspace=0.1)

#     ax_tree = fig.add_subplot(gs[0])
#     ax_heat = fig.add_subplot(gs[1])

#     Phylo.draw(tree, do_show=False, axes=ax_tree, show_confidence=False)
#     ax_tree.invert_yaxis()
#     ax_tree.axis("off")

#     sns.heatmap(
#         heatmap_df.values,
#         ax=ax_heat,
#         cmap=cmap,
#         cbar_kws={"label": norm_label},
#         linewidths=0.5,
#         linecolor="white",
#         xticklabels=heatmap_df.columns,
#         yticklabels=heatmap_df.index
#     )
#     ax_heat.tick_params(axis="y", labelsize=8)
#     ax_heat.tick_params(axis="x", labelrotation=90)
#     ax_heat.set_ylabel("Sample")
#     ax_heat.set_xlabel("Gene")

#     plt.tight_layout()
#     plt.savefig(picture_path)
#     plt.close(fig)

# def plot(csv_path, out_csv_path, tree_path, pic_path):
#     tree, depths, max_depths = load_tree(tree_path)
#     leaf_position = compute_leaf_positions(tree)
#     leaf_position_lower = {k.lower(): v for k, v in leaf_position.items()}
#     kd_martrix = compute_universal_activity_matrix(csv_path, out_csv_path)
#     df_heatmap = kd_martrix.set_index("Sample")
#     df_heatmap = df_heatmap.loc[leaf_position_lower.keys()]

#     draw_tree_and_heatmap(tree, list(df_heatmap.index), df_heatmap, pic_path)

if __name__ == '__main__':
    plot('/home/shiyi/Skim2StructProject/PART3_OUT/energy.csv',
         '/home/shiyi/Skim2StructProject/PART3_OUT/energy_process.csv',
         '/home/shiyi/Skim2StructProject/example_data/test3/consensus.fasta.treefile', 
         '/home/shiyi/Skim2StructProject/PART3_OUT/picture.png')
