# skim2struct/cli.py
import argparse
import os
from skim2struct.TreeConservationModule.core import run as run_tree
from skim2struct.EvoDnDsModule.core import run as run_dnds
from skim2struct.DockingModule.core import run as run_docking
import urllib.request

# Constants for MGLTools runtime installation
MGL_URL      = "https://ccsb.scripps.edu/download/532/"
MGL_ARCHIVE  = "mgltools_x86_64Linux2.1.5.7.tar.gz"
MGL_DIR_NAME = "mgltools_x86_64Linux2.1.5.7"

def ensure_mgltools():
    """Download, extract and patch MGLTools into skim2struct/utils if not already present."""
    base_dir   = os.path.abspath(os.path.dirname(__file__))
    utils_dir  = os.path.join(base_dir, 'utils')
    target_dir = os.path.join(utils_dir, MGL_DIR_NAME)

    if not os.path.isdir(target_dir):
        os.makedirs(utils_dir, exist_ok=True)
        archive_path = os.path.join(utils_dir, MGL_ARCHIVE)
        print(f"Downloading MGLTools from {MGL_URL} …")
        urllib.request.urlretrieve(MGL_URL, archive_path)
        print("Extracting MGLTools …")
        with tarfile.open(archive_path) as tar:
            tar.extractall(path=utils_dir)

        # Patch prepare_ligand4.py
        patch_path = os.path.join(
            target_dir,
            'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_ligand4.py'
        )
        print(f"Patching {patch_path} …")
        with open(patch_path, 'r') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            if 'ligand_filename = os.path.basename(a)' in line and not line.strip().startswith('#'):
                new_lines.append('# ' + line)
            elif 'ligand_filename = a' in line and line.strip().startswith('#'):
                new_lines.append(line.lstrip('# '))
            else:
                new_lines.append(line)
        with open(patch_path, 'w') as f:
            f.writelines(new_lines)

        print("MGLTools installed and patched under utils.")
    else:
        print("MGLTools already present, skipping download and patch.")

    return target_dir


def main():
    ensure_mgltools()
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