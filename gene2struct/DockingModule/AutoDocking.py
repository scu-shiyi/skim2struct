
from pathlib import Path
from typing import Dict, Tuple, Set
import os
import shutil
import pandas as pd
from gene2struct.utils.PreLigand import process_ligand
from gene2struct.DockingModule.DockingExecutor import DockingExecutor
from gene2struct.DockingModule.plot import plot
from gene2struct.utils.PreReceptor import process_receptors, convert_cif_to_pdb,clean_filename
import re
from gene2struct.DockingModule.parse import build_table

def cif2pdb(orign_dir, receptor_pdb_dir):
    for gene in os.listdir(orign_dir):
        gene_path = os.path.join(orign_dir, gene)
        if not os.path.isdir(gene_path):
            continue

        out_dir = os.path.join(receptor_pdb_dir, gene)
        os.makedirs(out_dir, exist_ok=True)  

        for file in os.listdir(gene_path):
            src_path = os.path.join(gene_path, file)
            if file.lower().endswith('.pdb'):

                base, ext = os.path.splitext(file)
                basename = clean_filename(base, gene)
                dst_path = os.path.join(out_dir, f"{basename}{ext}")
                shutil.copy2(src_path, dst_path)

            elif file.lower().endswith('.cif'):
                convert_cif_to_pdb(src_path, out_dir, gene)


def parse_mapping(mapping_file: str) -> Tuple[Dict[str, list[str]], list[str], list[str]]:

    def split_multiple(val):
        val = str(val).strip()
        for sep in [",", ";", "|", "/"]:
            if sep in val:
                return [v.strip().replace(" ", "-") for v in val.split(sep) if v.strip()]
        return [val.replace(" ", "-")] if val else []


    try:
        df = pd.read_csv(mapping_file, sep=None, engine="python",encoding='utf-8-sig')

    except Exception as e:
        raise ValueError(f"Failed to parse ligand mapping file: {e}")


    df.columns = [col.strip().capitalize() for col in df.columns]
    print(df.columns)
    if not {'Gene', 'Substrate', 'Product'}.issubset(df.columns):
        raise ValueError("Missing required columns: Gene, Substrate, Product")

    gene_ligand_map = {}
    ligand_set = set()

    for _, row in df.iterrows():
        gene = str(row.get("Gene", "")).strip()
        if not gene:
            continue
        substrates = split_multiple(row.get("Substrate", ""))
        products = split_multiple(row.get("Product", ""))
        ligands = list(set(substrates + products))
        gene_ligand_map[gene] = ligands
        ligand_set.update(ligands)

    genes = list(gene_ligand_map.keys())
    return gene_ligand_map, list(ligand_set), genes



def AutoDocking(input_protein_dir, mapping_csv, output_dir, tree_path):
    pdb_dir = os.path.join(output_dir, 'pdb')
    receptor_pdb = os.path.join(pdb_dir, 'receptor_pdb')
    ligand_pdb = os.path.join(pdb_dir, 'ligand_pdb')
    os.makedirs(receptor_pdb, exist_ok=True)
    os.makedirs(ligand_pdb, exist_ok=True)

    temp_dir = os.path.join(output_dir, "temp")
    temp_ligand = os.path.join(temp_dir, 'ligand')
    temp_center = os.path.join(temp_dir, 'center')
    os.makedirs(temp_ligand, exist_ok=True)
    os.makedirs(temp_center, exist_ok=True)
    
    pdbqt_dir = os.path.join(output_dir, 'pdbqt')
    receptor_pdbqt = os.path.join(pdbqt_dir, "receptor_pdbqt")
    ligand_pdbqt = os.path.join(pdbqt_dir, "ligand_pdbqt")
    os.makedirs(receptor_pdbqt, exist_ok=True)
    os.makedirs(ligand_pdbqt, exist_ok=True)

    docking_dir = os.path.join(output_dir, 'docking')
    cif2pdb(input_protein_dir, receptor_pdb)
    gene_list = [] 
    for gene in sorted(os.listdir(receptor_pdb)):
        gene_path = os.path.join(receptor_pdb, gene)
        if os.path.isdir(gene_path):
            valid = any(f.endswith(('.pdb', '.cif')) for f in os.listdir(gene_path))
            if valid:
                gene_list.append(gene)

    gene_ligand_map, ligand_set, genes = parse_mapping(mapping_csv)

    if set(gene_list) == set(genes):
        process_ligand(ligand_set, temp_ligand, ligand_pdb, ligand_pdbqt)
    else:
        print(set(gene_list))
        print(set(genes))
        raise ValueError("Gene names in protein folder and mapping CSV must match exactly.")


    for root,_,files in os.walk(receptor_pdb):
        for file in files:
            if file.endswith(".pdb"):
                src_file = os.path.join(root, file)
                relative_path = os.path.relpath(root, receptor_pdb)
                out_dir = os.path.join(receptor_pdbqt, relative_path)
                os.makedirs(out_dir, exist_ok=True)
                process_receptors(src_file, out_dir)


    executor = DockingExecutor(
        receptor_pdb_dir = receptor_pdb,
        receptor_pdbqt_dir = receptor_pdbqt,
        ligand_pdbqt_dir = ligand_pdbqt,
        docking_dir= docking_dir,
        temp_center=temp_center,
        gene_ligand_map = gene_ligand_map,
    )
    executor.run_all()

    row_binding_energy = os.path.join(output_dir, 'raw_binding_energy.csv')
    catalytic_activity_matrix = os.path.join(output_dir,'catalytic_activity_matrix.csv')
    phylo_activity_heatmap = os.path.join(output_dir, 'phylo_activity_heatmap.png')
    build_table(docking_dir, mapping_csv, row_binding_energy)

    plot(row_binding_energy, catalytic_activity_matrix, tree_path, phylo_activity_heatmap)

