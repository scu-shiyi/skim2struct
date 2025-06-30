# skim2struct/TreeConservationModule/heatmap_utils.py

import pandas as pd
import ast
from skim2struct.utils.EvoScoring import EvoScoring

def compute_entropy_matrix(fasta_path: str, output_dir: str) -> str:
    """
    调用 EvoScoring 计算位置熵矩阵，返回 CSV 路径。
    """
    model = EvoScoring(fasta_path, output_dir)
    return model.position_entropy()

def load_heatmap_data(heatmap_path: str, leaf_labels: list[str]) -> pd.DataFrame:
    """
    加载 CSV，排序并转换为数值矩阵。
    """
    heatmap_df = pd.read_csv(heatmap_path, index_col=0, header=0)
    heatmap_df = heatmap_df.map(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    ordered_heatmap_df = heatmap_df.reindex(leaf_labels)

    return ordered_heatmap_df