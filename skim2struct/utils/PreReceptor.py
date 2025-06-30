import os
from skim2struct.utils.PreparePDBQT import receptor2pdbqt  # 你已有的函数
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
        raise ValueError('蛋白质文件夹内命名不规范')
    return basename


def convert_cif_to_pdb(cif_file, receptor_pdb,gene_name) -> Union[str, None]:
    """将cif文件转化为pdb文件"""
    try:
        mol = next(pybel.readfile("cif", cif_file))
    except StopIteration:
        print(f"cif 文件无效: {cif_file}")
        return None
    base = os.path.splitext(os.path.basename(cif_file))[0]
    basename = clean_filename(base, gene_name)

    pdb_file = os.path.join(receptor_pdb, f"{basename}.pdb")
    mol.write("pdb", pdb_file, overwrite=True)
    print(f"受体pdb文件保存在{receptor_pdb}")
    return pdb_file

def process_receptors(receptor_input, receptor_pdbqt) -> str:
    """
    支持输入单个文件或文件夹，文件格式包括cif或pdb
    """
    #    # 创建中间/输出目录
    # receptor_pdb = os.path.join(output_dir, 'autodock_dir/temp_receptor/receptor_pdb')
    # receptor_pdbqt = os.path.join(output_dir, 'autodock_dir/receptor_pdbqt')
    # # if os.path.exists(receptor_pdb):
    # #     shutil.rmtree(receptor_pdb)
    # # if os.path.exists(receptor_pdbqt):
    # #     shutil.rmtree(receptor_pdbqt)
    # os.makedirs(receptor_pdb, exist_ok=True)
    # os.makedirs(receptor_pdbqt, exist_ok=True)



    # if Path(receptor_input).suffix.lower() == '.cif':
    #     pdb_file = convert_cif_to_pdb(receptor_input, receptor_pdb)
    #     pdbqt_file = receptor2pdbqt(pdb_file, receptor_pdbqt)
    # elif Path(receptor_input).suffix.lower() == '.pdb':
    #     dst_file = Path(receptor_pdb) / Path(receptor_input).name
    #     if not dst_file.exists():
    #         shutil.copy(receptor_input, receptor_pdb)
    pdbqt_file = receptor2pdbqt(receptor_input, receptor_pdbqt)
    return pdbqt_file

def run_fpocket(receptor_pdb, temp_receptor_center, gene_name):
    FPOCKET_BIN = shutil.which("fpocket")
    receptor_pdb = Path(receptor_pdb)
    receptor_name = receptor_pdb.stem
    dest = Path(temp_receptor_center) / gene_name / f"{receptor_name}_out"

    if dest.exists():
        return dest

    # 创建一个完全独立的临时目录（避免fpocket内部冲突）
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tmp_receptor = tmpdir / receptor_pdb.name
        shutil.copy2(receptor_pdb, tmp_receptor)

        command = [FPOCKET_BIN, "-f", str(tmp_receptor)]
        try:
            subprocess.run(command, check=True, cwd=tmpdir)
        except subprocess.CalledProcessError as e:
            print(f"❌ fpocket 执行失败: {receptor_pdb.name} → {e}")
            return None
        
        out_dir = tmpdir / f"{receptor_name}_out"
        if out_dir.exists():
            try:
                shutil.move(str(out_dir), str(dest))
            except shutil.Error as e:
                # 跨设备移动失败 fallback
                print(f"⚠️ shutil.move 跨设备失败，fallback copytree: {e}")
                shutil.copytree(str(out_dir), str(dest))
                shutil.rmtree(out_dir, ignore_errors=True)
            return dest
        else:
            print(f"❌ fpocket 成功运行但未生成输出目录: {out_dir}")
            return None




def parse_fpocket_pqr(output_file):
    """
    在 fpocket 结果目录中自动查找 pockets.pqr 文件，计算结合口袋的中心和大小
    """
    pqr_file = os.path.join(output_file, 'pockets/pocket1_vert.pqr')
    if not os.path.exists(pqr_file):
        print("❌ 未找到 pockets.pqr 文件，请检查 fpocket 运行结果！")
        return None, None
    coords = []
    with open(pqr_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
                coords.append([x, y, z])
    if len(coords) == 0:
        print("❌ PQR 文件无原子数据，无法计算中心和大小！")
        return None, None

    coords = np.array(coords)
    center = np.mean(coords, axis=0)  # 计算几何中心
    size = np.ptp(coords, axis=0) + 10  # 计算尺寸，增加10Å缓冲
    # print(f"✅ 口袋中心: {center}")
    # print(f"✅ 口袋大小: {size}")
    return center, size
# if __name__ == "__main__":
