import re, os, sys, logging
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def _split_ligands(raw: str):
    if not raw or pd.isna(raw):
        return []
    return [x.strip().replace(" ", "-") for x in re.split(r"[;,|/]+", str(raw)) if x.strip()]


def load_mapping(csv_path: str):
    df = pd.read_csv(csv_path)
    df.columns = [c.strip().capitalize() for c in df.columns]
    if not {"Gene", "Substrate", "Product"}.issubset(df.columns):
        raise ValueError("CSV file must contain columns: Gene / Substrate / Product")

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
    try:
        with open(log_file) as f:
            for line in f:
                m = _aff_pat.match(line)
                if m:
                    return float(m.group(1))
    except Exception as e:
        logger.warning(f"Failed to read {log_file}: {e}")
    return None


def build_table(docking_dir: str, mapping_csv: str, out_csv: str | None = None):
    mapping = load_mapping(mapping_csv)

    columns = [(gene, role, lig)
               for gene, rmap in mapping.items()
               for role, ligs in rmap.items()
               for lig in ligs]
    col_map = {col: idx for idx, col in enumerate(columns)}  
    rows_data, species_seen = {}, set()

    for gene_dir in Path(docking_dir).iterdir():
        if not gene_dir.is_dir():
            continue
        gene = gene_dir.name
        if gene not in mapping:
            logger.warning(f"Gene {gene} not in mapping file, skipped.")
            continue

        for run_dir in gene_dir.iterdir():
            if not run_dir.is_dir():
                continue
            parts = run_dir.name.split("__", 2)
            if len(parts) != 3:
                logger.warning(f"Directory {run_dir} does not follow naming rule: gene__species__ligand, skipped.")
                continue
            g, species, ligand = parts
            if g != gene:
                continue
            species_seen.add(species)

            role = None
            for r, lst in mapping[gene].items():
                if ligand.lower() in (l.lower() for l in lst):
                    role = r
                    break
            if role is None:
                logger.debug(f"Ligand {ligand} of {gene} not found in substrate/product mapping, skipped.")
                continue

            log_file = next(run_dir.glob("*.log"), None) or next(run_dir.glob("*.txt"), None)
            if not log_file:
                logger.warning(f"No vina log file found in {run_dir}")
                continue
            aff = best_affinity(log_file)
            if aff is None:
                logger.warning(f"Failed to parse affinity from {log_file}")
                continue

            col_key = (gene, role, ligand)
            if col_key not in col_map:
                logger.warning(f"{col_key} not in predefined columns, skipped.")
                continue
            col_idx = col_map[col_key]
            rows_data.setdefault(species, [None] * len(columns))
            rows_data[species][col_idx] = aff

    df = pd.DataFrame.from_dict(
        rows_data, orient="index",
        columns=pd.MultiIndex.from_tuples(columns, names=["Gene", "Role", "Ligand"])
    ).sort_index(axis=1).sort_index(axis=0)

    if out_csv is not None and out_csv.strip() != "":
        df.to_csv(out_csv, float_format="%.6g")
        logger.info(f"Docking summary table saved to {out_csv}")
    return out_csv

