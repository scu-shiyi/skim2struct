# skim2struct/EvoDnDsModule/RunEvoDnDs.py
import os
from pathlib import Path
import pandas as pd
from skim2struct.utils.prepare_fasta_batch import preprocess_fasta_dir
from skim2struct.EvoDnDsModule.calcul_dnds import CalculDnDs
from skim2struct.utils import EvoScoring
from skim2struct.EvoDnDsModule.plot import generate_summary_plot
from skim2struct.EvoDnDsModule.prepare_data import batch_run                # 复用你已有的批量函数





def RunEvoDnDs(fasta_dir, output_dir, tree_path, outgroup=None):
    os.makedirs(output_dir, exist_ok=True)
    fasta_paths = [str(p) for p in Path(fasta_dir).rglob("*") if p.suffix.lower() in [".fa", ".fasta", ".fas", ".txt"]]
    processed_paths, processed_dir = preprocess_fasta_dir(fasta_paths, output_dir, outgroup)

    combined_df = batch_run(processed_paths, str(output_dir), tree_path)

    summary_csv = os.path.join(output_dir, "EvoDnDs_summary.csv")
    combined_df.to_csv(summary_csv, index=False)
    print(f"Summary CSV saved to {summary_csv}")

    plot_path = os.path.join(output_dir, "EvoDnDs_Result.png")
    generate_summary_plot(combined_df, plot_path)

    return summary_csv, plot_path