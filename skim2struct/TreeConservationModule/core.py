# TreeConservationModule/core.py
from skim2struct.TreeConservationModule.RunTreeConsevation import RunTreeConservation
from pathlib import Path


def run(args):
    fasta_path = str(Path(args.fasta_path).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    tree_path = str(Path(args.tree_path).resolve()) if args.tree else None
    heatmap_path = str(Path(args.heatmap_path).resolve()) if hasattr(args, "heatmap") and args.heatmap else None
    name_limit=args.name_limit if args.name_limit is not None else 20

    output_path = RunTreeConservation(
        fasta_path=fasta_path,
        output_dir=output_dir,
        name_limit=name_limit,
        tree_path=tree_path,
        heatmap_path=heatmap_path
    )
    
    return output_path