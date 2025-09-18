# gene2struct/EvoDnDsModule/core.py
from pathlib import Path
from gene2struct.EvoDnDsModule.RunEvoDnDs import RunEvoDnDs

def run(args):
    fasta_input = str(Path(args.fasta_input).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    fasta_tree_map = str(Path(args.tree_map).resolve()) if args.tree_map else None
    outgroup = args.og

    out_png = RunEvoDnDs(
        fasta_input=fasta_input,
        output_dir=output_dir,
        fasta_tree_map=fasta_tree_map,
        outgroup=outgroup)
    
    return out_png
