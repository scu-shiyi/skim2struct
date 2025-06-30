# skim2struct/TreeConservationModule/core_logic.py

from pathlib import Path
import os
from skim2struct.utils import TreeFunction
from skim2struct.TreeConservationModule import conservation_calcul
from skim2struct.TreeConservationModule.plot import draw_tree_and_heatmap

def RunTreeConservation(
    fasta_path: str,
    output_dir: str,
    tree_path: str,
    name_limit: int = 20,
    heatmap_path: str = None
) -> str:
    # 正式的时候： heatmap_path应该删除，现在做调试
    """
    主函数：根据输入FASTA，构建进化树和保守性热图，并保存结果图像。

    参数：
        fasta_path: 输入的FASTA对齐文件
        output_dir: 输出目录
        name_limit: 每个叶节点名称显示的最大字符数
        tree_path: 可选的进化树文件（Newick 格式）
        heatmap_path: 可选的已计算的熵值CSV文件

    返回：
        保存图像的路径
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    

    gene_name = Path(fasta_path).stem

    # 1. 构建或读取进化树
    if tree_path is None:
        temp_tree = output_dir / "tree_temp"
        temp_tree.mkdir(exist_ok=True)
        tree_path = TreeFunction.build_tree(fasta_path, str(temp_tree))

    tree, depths, max_depth = TreeFunction.load_tree(tree_path)
    leaf_positions = TreeFunction.compute_leaf_positions(tree)

    # 2. 构建或读取热图数据
    if heatmap_path is None:
        heatmap_path = conservation_calcul.compute_entropy_matrix(fasta_path, str(output_dir))  
    
    heatmap_df = conservation_calcul.load_heatmap_data(heatmap_path, leaf_positions.keys())

    # 3. 绘图保存
    output_path = draw_tree_and_heatmap(
        tree, depths, max_depth,
        leaf_positions, heatmap_df,
        output_dir=str(output_dir),
        gene_name=gene_name,
        name_limit=name_limit
    )

    return output_path