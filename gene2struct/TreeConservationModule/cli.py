# skim2struct/TreeConservationModule/tree_conservation_cli.py
# 只做调试使用

import argparse
import os
from pathlib import Path
from skim2struct.TreeConservationModule.core import run

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="skim2struct siteview",
        description="Skim2Struct | SiteView module: Visualize site-wise heterogeneity (positive selection and entropy) combined with a phylogenetic tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    parser.add_argument("-f", "--fasta_path", required=True,
                        metavar="FASTA_FILE",
                        help="Path to a single aligned FASTA file (e.g., -f gene1.aln.fasta)")
    
    parser.add_argument("-o", "--output_dir", required=True,
                    metavar="OUT_DIR",
                    help="Directory to save results")
    
    parser.add_argument("-t", "--tree_path", 
                        default=None,
                        metavar="TREE_FILE",
                        help="Optional: Newick-format phylogenetic tree file (.treefile). If not provided, a tree will be inferred automatically.")
    
    parser.add_argument("--site",
                        action="store_true",
                        help="Enable site-level conservation calculation (default: disabled).")
    
    parser.add_argument("-n", "--name-limit", type=int, 
                        metavar="CHAR_LIMIT",
                        help="Maximum number of characters shown per leaf label on the tree.")
    
    parser.add_argument("-g", "--outgroups",
                        nargs="+",
                        default=None,
                        metavar="SPECIES",
                        help="Optional: One or more outgroup species names for tree rooting (e.g., -g Sp1 Sp2 Sp3).")
    
    parser.add_argument("--thr",
                        type=float,
                        default=1.4,
                        metavar="THRESHOLD",
                        help="Threshold for distinguishing conserved vs. poorly conserved sites (default: 1.4).")

    parser.add_argument("--heatmap_path",
                        default=None,
                        metavar="ENTROPY_CSV",
                        help="Optional: Precomputed entropy matrix CSV (for debugging only).")
    
    parser.add_argument("--mlc_path",
                        default=None,
                        metavar="MLC_FILE",
                        help="Optional: Path to PAML .mlc file for BEB site information.")
    

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    run(args)


if __name__ == "__main__":
    main()
