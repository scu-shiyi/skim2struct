# TreeConservationModule/core.py
from skim2struct.TreeConservationModule.RunTreeConservation import RunTreeConservation
from pathlib import Path


def run(args):
    fasta_path = str(Path(args.fasta_path).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    tree_path = str(Path(args.tree_path).resolve()) if args.tree_path else None
    is_calculate_site = args.site
    name_limit=args.name_limit if args.name_limit is not None else 20
    outgroups    = args.og if args.og else []
    thr = args.thr if args.thr is not None else 1.4
    mlc_path     = str(Path(args.mlc_path).resolve()) if args.mlc_path else None
    heatmap_path = str(Path(args.heatmap_path).resolve()) if args.heatmap_path else None


    output_path = RunTreeConservation(
        fasta_path=fasta_path,
        output_dir=output_dir,
        tree_path=tree_path,
        is_calculate_site=is_calculate_site,
        name_limit=name_limit,
        mlc_path=mlc_path,
        outgroups=outgroups,
        heatmap_path=heatmap_path,
        thr=thr
    )
    
    return output_path
