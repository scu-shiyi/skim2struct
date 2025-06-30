# skim2struct/TreeConservationModule/tree_conservation_cli.py
# 只做调试使用

import argparse
import os
from pathlib import Path
from skim2struct.TreeConservationModule.core import run

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="skim2struct treeconservation",
        description="Skim2Struct | TreeConservation module: Visualize phylogenetic tree and conservation heatmap from a single aligned FASTA file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    parser.add_argument("-f", "--fasta_path", required=True,
                        metavar="FASTA_FILE",
                        help="Path to a single aligned FASTA file (e.g., -f gene1.aln.fasta)")
    
    parser.add_argument("-t", "--tree_path", 
                        default=None,
                        metavar="TREE_FILE",
                        help="Optional: Newick-format phylogenetic tree file (.treefile). If not provided, a tree will be inferred automatically.")
    
    parser.add_argument("-n", "--name-limit", type=int, 
                        metavar="CHAR_LIMIT",
                        help="Maximum number of characters shown per leaf label on the tree.")
    
    parser.add_argument("-o", "--output_dir",
                        default=os.path.join(os.getcwd(), "PhyloEntropyAnalysis"),
                        metavar="OUT_DIR",
                        help="Directory to save results (default: ./PhyloEntropyAnalysis)")

    parser.add_argument("--heatmap_path",
                        default=None,
                        metavar="ENTROPY_CSV",
                        help="Optional: Precomputed entropy matrix CSV (for debugging only).")
    
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    run(args)


if __name__ == "__main__":
    main()