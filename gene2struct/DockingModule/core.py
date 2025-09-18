from pathlib import Path
from typing import Optional
from gene2struct.DockingModule.AutoDocking import AutoDocking



def run(args):
    """Thin wrapper around pipeline.AutoDocking so CLI / API share the same entry."""
    input_protein_dir = str(Path(args.protein_dir).resolve())
    mapping_csv       = str(Path(args.mapping_csv).resolve())
    tree_path         = str(Path(args.tree_path).resolve())
    output_dir        = str(Path(args.output_dir).resolve())


    AutoDocking(
        input_protein_dir=input_protein_dir,
        mapping_csv=mapping_csv,
        tree_path=tree_path,
        output_dir=output_dir,)

    return None
