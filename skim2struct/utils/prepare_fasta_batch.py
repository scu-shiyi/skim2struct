# skim2struct/utils/prepare_fasta_batch.py
# 对序列预处理，修订文件名称/去除外类群
import os
import shutil
import re
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO

# ----------------- 工具函数 ----------------- #
def clean_sample_id(file_name: str) -> str:
    """将文件名标准化为简单 ID"""
    basename = Path(file_name).stem
    cleaned = re.sub(r'[\s\-]+', '_', basename)
    cleaned = re.sub(r'[^A-Za-z0-9_]+', '_', cleaned)
    cleaned = re.sub(r'_+', '_', cleaned).strip('_')
    return cleaned


def _filter_fasta(in_path: str, out_path: str, outgroups: List[str]) -> bool:
    """读取 FASTA，剔除外类群序列后保存；若剔除后为空则返回 False"""
    records = list(SeqIO.parse(in_path, "fasta"))
    if outgroups:
        records = [r for r in records if all(og not in r.id for og in outgroups)]
    if not records:          # 全部被剔除
        return False
    SeqIO.write(records, out_path, "fasta")
    return True


# ------------- 主要接口函数 ----------------- #
def preprocess_fasta_dir(
        fasta_files: List[str],
        output_dir: str,
        outgroups: List[str] | None = None,
) -> Tuple[List[str], str]:
    """
    扫描目录 -> (可选)去外类群 -> 清洁命名 -> 复制/过滤到输出目录
    返回 (处理后 FASTA 列表, 输出目录)
    """
    outgroups = outgroups or []

    processed_dir = os.path.join(output_dir, "processed_fasta")
    os.makedirs(processed_dir, exist_ok=True)


    processed_files = []
    for fasta_path in fasta_files:
        clean_gene_name = clean_sample_id(fasta_path)
        dst_path = os.path.join(processed_dir, f"{clean_gene_name}.fa")

        ok = _filter_fasta(str(fasta_path), str(dst_path), outgroups)
        if ok:
            processed_files.append(str(dst_path))

    processed_files = list(dict.fromkeys(processed_files)) # 去重复
    return processed_files, str(processed_dir)

