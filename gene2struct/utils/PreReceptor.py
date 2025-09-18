import os
from gene2struct.utils.PreparePDBQT import receptor2pdbqt  # 你已有的函数
from typing import Union
from openbabel import pybel
from pathlib import Path
import shutil
import subprocess
import numpy as np
import tempfile
import re


def clean_filename(basename, gene_name):
    if gene_name:
        pattern = re.compile(rf"[_-]?{re.escape(gene_name)}[_-]?", flags=re.I)
        basename = pattern.sub("", basename)
    for tag in ['fold', 'model', 'structure','protein', 'gene']:
        pattern = re.compile(rf"[_-]?{tag}(?:[_-]?\d+)?", flags=re.I)
        basename = pattern.sub("", basename)

    basename = re.sub(r"[_-]{2,}", "_", basename).strip("_-")
    if not basename:
        raise ValueError("Invalid receptor naming inside the protein directory.")
    return basename


def convert_cif_to_pdb(cif_file, receptor_pdb,gene_name) -> Union[str, None]:
    """Convert CIF file to PDB file using OpenBabel."""
    try:
        mol = next(pybel.readfile("cif", cif_file))
    except StopIteration:
        print(f"Invalid CIF file: {cif_file}")
        return None
    base = os.path.splitext(os.path.basename(cif_file))[0]
    basename = clean_filename(base, gene_name)

    pdb_file = os.path.join(receptor_pdb, f"{basename}.pdb")
    mol.write("pdb", pdb_file, overwrite=True)
    print(f"Receptor PDB file saved to {receptor_pdb}")
    return pdb_file

def process_receptors(receptor_input, receptor_pdbqt) -> str:

    pdbqt_file = receptor2pdbqt(receptor_input, receptor_pdbqt)
    return pdbqt_file

def run_fpocket(receptor_pdb, temp_receptor_center, gene_name):
    FPOCKET_BIN = shutil.which("fpocket")
    receptor_pdb = Path(receptor_pdb)
    receptor_name = receptor_pdb.stem
    dest = Path(temp_receptor_center) / gene_name / f"{receptor_name}_out"

    if dest.exists():
        return dest

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tmp_receptor = tmpdir / receptor_pdb.name
        shutil.copy2(receptor_pdb, tmp_receptor)

        command = [FPOCKET_BIN, "-f", str(tmp_receptor)]
        try:
            subprocess.run(command, check=True, cwd=tmpdir)
        except subprocess.CalledProcessError as e:
            print(f"fpocket failed for {receptor_pdb.name}: {e}")
            return None
        
        out_dir = tmpdir / f"{receptor_name}_out"
        if out_dir.exists():
            try:
                shutil.move(str(out_dir), str(dest))
            except shutil.Error as e:
                print(f"shutil.move failed across devices, using copytree: {e}")
                shutil.copytree(str(out_dir), str(dest))
                shutil.rmtree(out_dir, ignore_errors=True)
            return dest
        else:
            print(f"fpocket completed but no output directory generated: {out_dir}")
            return None




def parse_fpocket_pqr(output_file):

    pqr_file = os.path.join(output_file, 'pockets/pocket1_vert.pqr')
    if not os.path.exists(pqr_file):
        print("No pockets.pqr file found. Please check fpocket results.")
        return None, None
    coords = []
    with open(pqr_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
                coords.append([x, y, z])
    if len(coords) == 0:
        print("PQR file contains no atom data. Cannot calculate pocket center/size.")
        return None, None

    coords = np.array(coords)
    center = np.mean(coords, axis=0)  
    size = np.ptp(coords, axis=0) + 10  

    return center, size
