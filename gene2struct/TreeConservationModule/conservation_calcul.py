# gene2struct/TreeConservationModule/heatmap_utils.py

import pandas as pd
import ast
from gene2struct.utils.EvoScoring import position_entropy

def compute_entropy_matrix(fasta_path: str, output_dir: str) -> str:
    return position_entropy(fasta_path, output_dir)

def load_heatmap_data(heatmap_path: str, leaf_labels: list[str], outgroups=None) -> pd.DataFrame:
    heatmap_df = pd.read_csv(heatmap_path, index_col=0, header=0)
    heatmap_df = heatmap_df.map(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    if outgroups:
        for og in outgroups:
            if og in heatmap_df.index:
                heatmap_df = heatmap_df.drop(index=og)
                
    ordered_heatmap_df = heatmap_df.reindex(leaf_labels)

    return ordered_heatmap_df
