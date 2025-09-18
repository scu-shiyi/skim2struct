
import argparse
from gene2struct.DockingModule.core import run
import os
def main():
    parser = argparse.ArgumentParser(
        prog="gene2struct docking",
        description="Skim2Struct | DockingModule: Predict enzyme activity via molecular docking and phylogenetic context.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-p", "--protein_dir", required=True,
                        metavar="PROTEIN_DIR",
                        help="Directory containing receptor structure subfolders.")
    
    parser.add_argument("-m", "--mapping_csv", required=True,
                        metavar="MAPPING_CSV",
                        help="CSV file mapping gene to substrate/product CIDs.")
    
    parser.add_argument("-t", "--tree_path", required=True,
                        metavar="TREE_FILE",
                        help="Phylogenetic tree in Newick format (used for result visualization).")
    
    parser.add_argument("-o", "--output_dir", 
                        default=os.path.join(os.getcwd(), 'DockingActivityAnalysis'),
                        metavar="OUT_DIR",
                        help="Directory to save results (default: ./DockingActivityAnalysis)")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    run(args)

if __name__ == "__main__":
    main()
