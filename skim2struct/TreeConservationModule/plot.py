# skim2struct/TreeConservationModule/plot.py

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.gridspec import GridSpec
from pathlib import Path
import io




def plot_clade(clade, ax, depths, max_depth, leaf_positions):
    """递归绘制进化树枝条结构"""
    x0 = depths[clade]
    if clade.is_terminal():
        y = leaf_positions[clade.name]
        ax.hlines(y, x0, max_depth, colors='gray', linestyles='--', linewidth=0.8)
        return y

    ys = []
    for child in clade.clades:
        y_child = plot_clade(child, ax, depths, max_depth, leaf_positions)
        ys.append(y_child)
        ax.hlines(y_child, x0, depths[child], colors='black', linestyles='-', linewidth=1)

    y_min, y_max = min(ys), max(ys)
    ax.vlines(x0, y_min, y_max, colors='black', linestyles='-', linewidth=1)
    return 0.5 * (y_min + y_max)


def draw_tree_and_heatmap(
    tree, depths, max_depth,
    leaf_positions, heatmap_df,
    output_dir: str, gene_name: str,
    name_limit: int = 20
) -> str:
    """
    绘制左树右热图组合图，直接保存为 PNG 文件，返回保存路径
    """

    fig = plt.figure(figsize=(25, 10), dpi=150)
    gs = GridSpec(1, 2, width_ratios=[1, 6], wspace=0.15)

    ax_tree = fig.add_subplot(gs[0])
    ax_heatmap = fig.add_subplot(gs[1])

    # 左侧绘制进化树
    plot_clade(tree.root, ax_tree, depths, max_depth, leaf_positions)
    ax_tree.set_ylim(-0.5, len(leaf_positions) - 0.5)
    ax_tree.set_xlim(0, max_depth * 1.02)
    ax_tree.invert_yaxis()
    ax_tree.axis('off')

    for name, y in leaf_positions.items():
        ax_tree.text(max_depth * 1.005, y, name[:name_limit], va='center', fontsize=8)

    # 右侧绘制热图
    heatmap_data = heatmap_df.values
    sns.heatmap(
        heatmap_data,
        ax=ax_heatmap,
        cmap="YlGnBu",
        cbar_kws={'label': 'Entropy'}
    )
    ax_heatmap.set_yticklabels([])
    ax_heatmap.tick_params(axis='y', which='both', labelleft=False, left=True)
    ax_heatmap.yaxis.set_ticks_position('left')
    ax_heatmap.set_yticks(np.arange(len(heatmap_data)) + 0.5)

    # 保存图像到本地
    plt.tight_layout()
    output_path = Path(output_dir) / f"{gene_name}.png"
    plt.savefig(output_path)
    plt.close(fig)

    return str(output_path)