# skim2struct/cli.py
import argparse
import os
from skim2struct.TreeConservationModule.core import run as run_tree
from skim2struct.EvoDnDsModule.core import run as run_dnds
from skim2struct.DockingModule.core import run as run_docking

def main():
    parser = argparse.ArgumentParser(
        prog="skim2struct",
        description="Skim2Struct: A modular toolkit for gene evolution and structure-guided functional analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # TreeConservation 子命令
    parser_tree = subparsers.add_parser("treeconservation", help="Run phylogenetic conservation analysis.")
    parser_tree.add_argument("-f", "--fasta_path", required=True, metavar="FASTA_FILE",
                             help="Path to a single aligned FASTA file (e.g., -f gene1.aln.fasta)")
    parser_tree.add_argument("-t", "--tree_path", default=None, metavar="TREE_FILE",
                             help="Optional: Newick-format phylogenetic tree file (.treefile). If not provided, a tree will be inferred automatically.")
    parser_tree.add_argument("-n", "--name-limit", type=int, default=20, metavar="CHAR_LIMIT",
                             help="Maximum number of characters shown per leaf label on the tree.")
    parser_tree.add_argument("-o", "--output_dir", default=os.path.join(os.getcwd(), "PhyloEntropyAnalysis"),
                             metavar="OUT_DIR",
                             help="Directory to save results (default: ./PhyloEntropyAnalysis)")
    parser_tree.add_argument("--heatmap_path", default=None, metavar="ENTROPY_CSV",
                             help="Optional: Precomputed entropy matrix CSV (for debugging only).")

    # EvoDnDs 子命令
    parser_dnds = subparsers.add_parser("evodnds", help="Run dN/dS estimation analysis.")
    parser_dnds.add_argument("-d", "--fasta_dir", required=True, metavar="FASTA_DIR",
                             help="Directory containing multiple aligned FASTA files (one per gene).")
    parser_dnds.add_argument("-t", "--tree_path", required=True, metavar="TREE_FILE",
                             help="Phylogenetic tree in Newick format.")
    parser_dnds.add_argument("-g", "--outgroup", default=None, metavar="SPECIES_NAME",
                             help="Optional: Outgroup species name for tree rooting.")
    parser_dnds.add_argument("-o", "--output_dir", default=os.path.join(os.getcwd(), "EvoSelectionAnalysis"),
                             metavar="OUT_DIR",
                             help="Directory to save results (default: ./EvoSelectionAnalysis)")

    # Docking 子命令
    parser_dock = subparsers.add_parser("docking", help="Run molecular docking and enzyme prediction.")
    parser_dock.add_argument("-p", "--protein_dir", required=True, metavar="PROTEIN_DIR",
                             help="Directory containing receptor structure subfolders.")
    parser_dock.add_argument("-m", "--mapping_csv", required=True, metavar="MAPPING_CSV",
                             help="CSV file mapping gene to substrate/product CIDs.")
    parser_dock.add_argument("-t", "--tree_path", required=True, metavar="TREE_FILE",
                             help="Phylogenetic tree in Newick format.")
    parser_dock.add_argument("-o", "--output_dir", default=os.path.join(os.getcwd(), "DockingActivityAnalysis"),
                             metavar="OUT_DIR",
                             help="Directory to save results (default: ./DockingActivityAnalysis)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.command == "treeconservation":
        run_tree(args)
    elif args.command == "evodnds":
        run_dnds(args)
    elif args.command == "docking":
        run_docking(args)


if __name__ == "__main__":
    main()