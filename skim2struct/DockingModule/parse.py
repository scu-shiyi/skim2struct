import re, os, sys, logging
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def _split_ligands(raw: str):
    """将 'Tyrosol;UDPG' 拆分为 ['Tyrosol', 'UDPG']，兼容 ; , | / 等"""
    if not raw or pd.isna(raw):
        return []
    return [x.strip().replace(" ", "-") for x in re.split(r"[;,|/]+", str(raw)) if x.strip()]


def load_mapping(csv_path: str):
    """解析 ligand 映射表（多配体拆开）"""
    df = pd.read_csv(csv_path)
    df.columns = [c.strip().capitalize() for c in df.columns]
    if not {"Gene", "Substrate", "Product"}.issubset(df.columns):
        raise ValueError("CSV 缺必要列 Gene / Substrate / Product")

    mapping = {}
    for _, row in df.iterrows():
        gene = str(row["Gene"]).strip()
        subs = _split_ligands(row["Substrate"])
        prods = _split_ligands(row["Product"])
        mapping[gene] = {
            "substrate": subs,
            "product": prods
        }
    return mapping


_aff_pat = re.compile(r"^\s*1\s+(-?\d+\.\d+)")


def best_affinity(log_file: Path) -> float | None:
    """从 vina log 中提取第一条打分（rank=1）"""
    try:
        with open(log_file) as f:
            for line in f:
                m = _aff_pat.match(line)
                if m:
                    return float(m.group(1))
    except Exception as e:
        logger.warning(f"读取失败 {log_file}: {e}")
    return None


def build_table(docking_dir: str, mapping_csv: str, out_csv: str | None = None):
    mapping = load_mapping(mapping_csv)

    # 生成所有 Gene / Role / Ligand 三元组列名
    columns = [(gene, role, lig)
               for gene, rmap in mapping.items()
               for role, ligs in rmap.items()
               for lig in ligs]
    col_map = {col: idx for idx, col in enumerate(columns)}  # 提高索引效率
    rows_data, species_seen = {}, set()

    for gene_dir in Path(docking_dir).iterdir():
        if not gene_dir.is_dir():
            continue
        gene = gene_dir.name
        if gene not in mapping:
            logger.warning(f"Gene {gene} 不在 CSV 映射中，忽略")
            continue

        for run_dir in gene_dir.iterdir():
            if not run_dir.is_dir():
                continue
            parts = run_dir.name.split("__", 2)
            if len(parts) != 3:
                logger.warning(f"目录 {run_dir} 命名不符合规则 gene__species__ligand，跳过")
                continue
            g, species, ligand = parts
            if g != gene:
                continue
            species_seen.add(species)

            # 判断 ligand 属于哪个 role
            role = None
            for r, lst in mapping[gene].items():
                if ligand.lower() in (l.lower() for l in lst):
                    role = r
                    break
            if role is None:
                logger.debug(f"{gene} 配体 {ligand} 不在 substrate/product 映射中，跳过")
                continue

            log_file = next(run_dir.glob("*.log"), None) or next(run_dir.glob("*.txt"), None)
            if not log_file:
                logger.warning(f"{run_dir} 无 vina log 文件")
                continue
            aff = best_affinity(log_file)
            if aff is None:
                logger.warning(f"{log_file} 未成功解析打分")
                continue

            col_key = (gene, role, ligand)
            if col_key not in col_map:
                logger.warning(f"{col_key} 未出现在预定义列中，跳过")
                continue
            col_idx = col_map[col_key]
            rows_data.setdefault(species, [None] * len(columns))
            rows_data[species][col_idx] = aff

    # 构建 DataFrame
    df = pd.DataFrame.from_dict(
        rows_data, orient="index",
        columns=pd.MultiIndex.from_tuples(columns, names=["Gene", "Role", "Ligand"])
    ).sort_index(axis=1).sort_index(axis=0)

    if out_csv is not None and out_csv.strip() != "":
        df.to_csv(out_csv, float_format="%.6g")
        logger.info(f"保存对接表至 {out_csv}")
    return out_csv

# ---------- 4) CLI（可与 click 或 argparse 集成，这里简单示例） ----------
if __name__ == "__main__":
    # print(load_mapping("/home/shiyi/Skim2StructProject/example_data/test3/test.csv"))
    table = build_table('/home/shiyi/Skim2StructProject/PART3_OUT/docking', 
                        '/home/shiyi/Skim2StructProject/example_data/test3/test.csv',
                        '/home/shiyi/Skim2StructProject/PART3_OUT/energy.csv')
    # print("✅ 汇总完成，存到", '/home/shiyi/Skim2StructProject/PART3_OUT/energy.csv')
    # print(table.head())