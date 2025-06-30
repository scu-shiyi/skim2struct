import requests
import os
from openbabel import pybel
from skim2struct.utils.PreparePDBQT import ligand2pdbqt, receptor2pdbqt
from typing import List, Union
import time
#
# class LigandDownloaderConverter:
#     def __init__(self, output_dir):
#         self.dock_dir = os.path.join(output_dir, 'autodock_dir')
#         os.makedirs(self.dock_dir, exist_ok=True)
#         self.temp_ligand = os.path.join(self.dock_dir, 'temp_ligand')
#         os.makedirs(self.temp_ligand, exist_ok=True) # 储存ligand过程文件----初始下载、pdb
#         # self.temp_receptor = os.path.join(self.dock_dir,'temp_receptor')
#         # os.makedirs(self.temp_receptor,  exist_ok=True)
#
#         self.ligand_pdbqt = os.path.join(self.dock_dir, "ligand_pdbqt")
#         os.makedirs(self.ligand_pdbqt, exist_ok=True)


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
            print(f"⚠️ 第 {attempt+1} 次尝试失败: {e}")
            time.sleep(3)
    print(f"网络异常，检查网址是否可访问\n{url}")
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
            print(f"CID {cid} 下载的 SDF 文件可能无效，已删除")
            os.remove(ligand_file)
            return None
        print(f"CID {cid} 下载成功")
        return ligand_file
    except Exception as e:
        print(f"下载 CID {cid} 失败：{e}")
        return None

def convert_sdf_to_pdb(sdf_file, ligand_pdb) -> Union[str, None]:
    name = os.path.splitext(os.path.basename(sdf_file))[0]
    try:
        mol = next(pybel.readfile("sdf", sdf_file))
    except StopIteration:
        print(f"SDF 文件无效: {sdf_file}")
        return None
    pdb_file = os.path.join(ligand_pdb, f"{name}.pdb")
    mol.write("pdb", pdb_file, overwrite=True)

    # return ligand_pdb

def process_ligand(inputs,temp_ligand,ligand_pdb,ligand_pdbqt):
    for item in inputs:
        print(f"\n🧬 处理配体: {item}")
        if item:
            if is_cid(item): # 如果输入cid
                cid = item
                sdf_file = os.path.join(temp_ligand, f'{cid}.sdf')
                if os.path.exists(sdf_file):
                    pdb_file = os.path.join(ligand_pdb, f'{cid}.pdb')
                    if not os.path.exists(pdb_file):
                        convert_sdf_to_pdb(sdf_file, ligand_pdb)
                else:
                    sdf_file = download_sdf_by_cid(cid, temp_ligand)
                    convert_sdf_to_pdb(sdf_file, ligand_pdb)
            else: # 输入name
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