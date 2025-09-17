#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

from Bio import SeqIO
from Bio import Phylo

from skim2struct.TreeConservationModule import conservation_calcul
from skim2struct.utils import TreeFunction
from skim2struct.TreeConservationModule.plot import draw_tree_and_heatmap
from skim2struct.utils.Phylip_Prepare import prepare_paml_input2
from skim2struct.utils.site import run_pair_model

def remove_outgroups(fasta_path, out_path, outgroups=None):

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if outgroups:
        records = [record for record in records if record.id not in outgroups]
    if len(records) < 3:
        raise ValueError(f"After removing outgroups, fewer than 3 sequences remain. Remaining: {len(records)}.")
    SeqIO.write(records, out_path, "fasta")
    return out_path


def RunTreeConservation(
    fasta_path: str,
    output_dir: str,
    tree_path: str = None,
    is_calculate_site: bool = False,
    name_limit: int = 20,
    thr: float = 1.4,
    outgroups: list[str] = None,
    mlc_path=None,
    heatmap_path: str = None
) -> str:

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_name = Path(fasta_path).stem

    if not is_calculate_site:

        file_input = output_dir / gene_name / "file_input"
        file_input.mkdir(parents=True, exist_ok=True)
        fasta_path = remove_outgroups(fasta_path, str(file_input / f"{gene_name}.fasta"), outgroups)


        if tree_path is None:
            tree_path = TreeFunction.build_tree1(fasta_path, str(file_input))


        for suffix in ["model.gz", "log", "iqtree", "ckp.gz", "mldist", "bionj"]:
            p = file_input / f"{gene_name}.{suffix}"
            if p.exists():
                p.unlink()

        tree, depths, max_depth = TreeFunction.load_tree(tree_path, outgroups)
        leaf_positions = TreeFunction.compute_leaf_positions(tree)


        if heatmap_path is None:
            evo_output = output_dir / gene_name / "evo_output"
            evo_output.mkdir(parents=True, exist_ok=True)
            heatmap_path = conservation_calcul.compute_entropy_matrix(fasta_path, str(evo_output))

        heatmap_df = conservation_calcul.load_heatmap_data(heatmap_path, leaf_positions.keys(), outgroups)


        output_path = draw_tree_and_heatmap(
            tree, depths, max_depth,
            leaf_positions, heatmap_df,
            output_dir=str(output_dir / gene_name),
            gene_name=gene_name,
            mlc_path=mlc_path,
            name_limit=name_limit,
            thr=thr,
            site_model=False

        )

        return output_path
    else:
        paml_output = output_dir / gene_name / "paml_output"
        paml_output.mkdir(parents=True, exist_ok=True)
        work_dir, phylip_file, tree_file, species = prepare_paml_input2(fasta_path, output_dir, tree_path=tree_path, outgroups=outgroups)
        tree, depths, max_depth = TreeFunction.load_tree(tree_file, outgroups)
        leaf_positions = TreeFunction.compute_leaf_positions(tree)

        m0m3_result = run_pair_model(phylip_file, tree_file, str(paml_output), "M0M3", only_branch=False,gene_name=gene_name)

        if m0m3_result["p"] < 0.05:
            if m0m3_result["p"] < 0.01:
                print("LRT (M0 vs M3) is significant at the 0.01 level: strong evidence of site-specific ω heterogeneity.)
            else:
                print("LRT (M0 vs M3) is significant at the 0.05 level: evidence of site-specific ω heterogeneity.")
            m7m8_result = run_pair_model(phylip_file, tree_file, str(paml_output), "M7M8", only_branch=False,gene_name=gene_name)
            if m7m8_result["p"] < 0.05:
                if m7m8_result["p"] < 0.01:
                    print("LRT (M7 vs M8) is significant at the 0.01 level: positive selection sites detected.")
                else:
                    print("LRT (M7 vs M8) is significant at the 0.05 level: positive selection sites detected.")
            else:
                print("LRT (M7 vs M8) is not significant: no positive selection sites detected.")
        mlc_path = paml_output / gene_name / paml_output / "M7M8" / "result" / "M8.mlc"

        if heatmap_path is None:
            evo_output = output_dir / gene_name / "evo_output"
            evo_output.mkdir(parents=True, exist_ok=True)
            heatmap_path = conservation_calcul.compute_entropy_matrix(fasta_path, str(evo_output))

        heatmap_df = conservation_calcul.load_heatmap_data(heatmap_path, leaf_positions.keys(), outgroups)


        output_path = draw_tree_and_heatmap(
            tree, depths, max_depth,
            leaf_positions, heatmap_df,
            output_dir=str(output_dir / gene_name),
            gene_name=gene_name,
            mlc_path=mlc_path,
            name_limit=name_limit,
            thr=thr,
            site_model=True
        )

        return output_path
        



