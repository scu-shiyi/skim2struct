"""
skim2struct/EvoDnDsModule/batch.py
批量运行 dN/dS + EvoScore 并汇总结果
"""

import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from skim2struct.utils.TreeFunction import build_tree
from skim2struct.utils.EvoScoring import EvoScoring
from skim2struct.EvoDnDsModule.calcul_dnds import CalculDnDs


# ----------------------------------------------------------------------
def _process_one_dnds(fasta_path: str, output_dir: str, tree_path: str) -> dict | None:
    """单基因 dN/dS 计算包装（子线程）"""
    gene = Path(fasta_path).stem
    if tree_path is None:
        temp_tree_dir  = os.path.join(output_dir, 'dnds_dir', 'temp_tree')
        os.makedirs(temp_tree_dir, exist_ok=True)
        tree_path = build_tree(fasta_path, temp_tree_dir)

    try:
        dnds_val = CalculDnDs(fasta_path, output_dir, tree_file=tree_path).run_single_return()  # outputdir 传入总体输出目录
        return {"Gene": gene, "DnDs_Score": dnds_val}
    except Exception as exc:
        print(f"⚠ dN/dS 计算失败: {gene} -> {exc}")
        return None


# ----------------------------------------------------------------------
def batch_run(
    processed_paths: list[str],
    output_dir: str,
    tree_path: str | None = None,
) -> pd.DataFrame:
    """
    参数
    processed_paths : List[str] 已经做完过滤 / 清洁命名的 FASTA 文件路径列表
    output_dir : str 总输出目录
    tree_path : str | None 可选的 Newick 树文件；若为 None，CalculDnDs 内部可自行构树
    返回
    ----
    pandas.DataFrame
        合并后的结果表，含 ['Gene', 'EvoMean', 'EvoStd', 'EvoValues', 'DnDs_Score']
    """
    
    # 1. evo评分计算
    evo_model = EvoScoring(processed_paths, output_dir)
    evo_score_file = evo_model.scoring()
    evo_df = pd.read_csv(evo_score_file)

    # 2. dN/dS 并行计算
    dnds_scores: list[dict] = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(_process_one_dnds, fasta_path, output_dir, tree_path): fasta_path
            for fasta_path in processed_paths
        }
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                dnds_scores.append(res)

    dnds_df = (
        pd.DataFrame(dnds_scores)
        if dnds_scores
        else pd.DataFrame(columns=["Gene", "DnDs_Score"])
    )
    dnds_df.to_csv(os.path.join(output_dir, "dnds_dir", "dnds_scores.csv"), index=False)
    # 3. 汇总 Evo Score：均值、标准差、原始向量
    print("准备合并两个得分文件...")
    genes = [col for col in evo_df.columns if col != "name"]
    evo_summary = []
    for gene in genes:
        vals = evo_df[gene].dropna().values.astype(float)
        evo_summary.append(
            {"Gene": gene, "EvoMean": vals.mean(), "EvoStd": vals.std(), "EvoValues": vals}
        )
    evo_df_summary = pd.DataFrame(evo_summary)


    # 4. 合并 dN/dS 与 EvoScore
    combined_df = pd.merge(evo_df_summary, dnds_df, on="Gene", how="left")
    return combined_df



