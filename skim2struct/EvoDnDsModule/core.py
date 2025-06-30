# skim2struct/EvoDnDsModule/core.py
from pathlib import Path
from skim2struct.EvoDnDsModule.RunEvoDnDs import RunEvoDnDs

def run(args):
    fasta_dir = str(Path(args.fasta_dir).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    tree_path = str(Path(args.tree_path).resolve()) if args.tree_path else None
    outgroup = args.outgroup

    output_files = RunEvoDnDs(
        fasta_dir=fasta_dir,
        output_dir=output_dir,
        tree_path=tree_path,
        outgroup=outgroup)
    
    return output_files