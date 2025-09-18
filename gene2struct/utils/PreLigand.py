import requests
import os
from openbabel import pybel
from gene2struct.utils.PreparePDBQT import ligand2pdbqt, receptor2pdbqt
from typing import List, Union
import time



def is_cid(item: str) -> bool:
    return item.isdigit()

def get_cid_from_name(name: str,retries=5) -> Union[str, None]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    for attempt in range(retries):
        try:
            res = requests.get(url, timeout=10)
            res.raise_for_status()
            data = res.json()
            return name.replace(' ','-'), str(data["IdentifierList"]["CID"][0])
        except Exception as e:
            print(f"Attempt {attempt+1} failed: {e}")
            time.sleep(3)
    print(f"Network error. Please check if the URL is accessible:\n{url}")
    return None

def download_sdf_by_cid(cid: str, temp_ligand, label=None) -> Union[str, None]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    if label is None:
        ligand_file = os.path.join(temp_ligand, f"{cid}.sdf")
    elif label:
        ligand_file = os.path.join(temp_ligand, f"{label}.sdf")
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        with open(ligand_file, "wb") as f:
            f.write(response.content)
        if os.path.getsize(ligand_file) < 100:
            print(f"CID {cid}: downloaded SDF file seems invalid, deleted.")
            os.remove(ligand_file)
            return None
        print(f"CID {cid} downloaded successfully.")
        return ligand_file
    except Exception as e:
        print(f"Failed to download CID {cid}: {e}")
        return None

def convert_sdf_to_pdb(sdf_file, ligand_pdb) -> Union[str, None]:
    name = os.path.splitext(os.path.basename(sdf_file))[0]
    try:
        mol = next(pybel.readfile("sdf", sdf_file))
    except StopIteration:
        print(f"Invalid SDF file: {sdf_file}")
        return None
    pdb_file = os.path.join(ligand_pdb, f"{name}.pdb")
    mol.write("pdb", pdb_file, overwrite=True)

    # return ligand_pdb

def process_ligand(inputs,temp_ligand,ligand_pdb,ligand_pdbqt):
    for item in inputs:
        print(f"\nProcessing ligand: {item}")
        if item:
            if is_cid(item):
                cid = item
                sdf_file = os.path.join(temp_ligand, f'{cid}.sdf')
                if os.path.exists(sdf_file):
                    pdb_file = os.path.join(ligand_pdb, f'{cid}.pdb')
                    if not os.path.exists(pdb_file):
                        convert_sdf_to_pdb(sdf_file, ligand_pdb)
                else:
                    sdf_file = download_sdf_by_cid(cid, temp_ligand)
                    convert_sdf_to_pdb(sdf_file, ligand_pdb)
            else: 
                name, cid = get_cid_from_name(item)
                sdf_file = os.path.join(temp_ligand, f'{name}.sdf')
                if os.path.exists(sdf_file):
                    pdb_file = os.path.join(ligand_pdb, f'{name}.pdb')
                    if not os.path.exists(pdb_file):
                        convert_sdf_to_pdb(sdf_file, ligand_pdb)
                else:
                    sdf_file = download_sdf_by_cid(cid, temp_ligand, label=name)
                    convert_sdf_to_pdb(sdf_file, ligand_pdb)

    ligand2pdbqt(ligand_pdb, ligand_pdbqt)
