# skim2struct/EvoDnDsModule/cli.py

import argparse
import os
from skim2struct.EvoDnDsModule.core import run

def main():
    parser = argparse.ArgumentParser(
        prog="skim2struct evoselect",
        description="Skim2Struct | EvoDnDs: Batch dN/dS and evolutionary scoring module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-d", "--fasta_dir", required=True,
                        metavar="FASTA_DIR",
                        help="Directory containing multiple aligned FASTA files (one per gene).")

    parser.add_argument("-t", "--tree_path",
                        metavar="TREE_FILE",
                        default=None,
                        help="Optional: Phylogenetic tree in Newick format. If omitted, gene trees will be built individually.")

    parser.add_argument("-g", "--outgroup",
                        metavar="SPECIES_NAME",
                        default=None,
                        help="Optional: Outgroup species name for tree rooting (used in tree construction).")
    
    parser.add_argument("-o", "--output_dir",
                        metavar="OUT_DIR",
                        default=os.path.join(os.getcwd(), "EvoSelectionAnalysis"),
                        help="Directory to save results (default: ./EvoSelectionAnalysis)")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    run(args)

if __name__ == "__main__":
    main()
