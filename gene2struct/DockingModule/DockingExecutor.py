import os
import subprocess
from multiprocessing import Pool, cpu_count
from pathlib import Path
from gene2struct.utils.PreReceptor import run_fpocket, parse_fpocket_pqr
from gene2struct.utils.PreReceptor import convert_cif_to_pdb  # 兼容非标准文件
from gene2struct.utils.PreReceptor import process_receptors
import shutil

class DockingExecutor:
    def __init__(self, receptor_pdb_dir, receptor_pdbqt_dir, ligand_pdbqt_dir, docking_dir,temp_center, gene_ligand_map):
        self.receptor_pdbqt_dir = receptor_pdbqt_dir
        self.ligand_pdbqt_dir = ligand_pdbqt_dir
        self.docking_dir = docking_dir
        self.gene_ligand_map = gene_ligand_map
        self.receptor_pdb_dir =receptor_pdb_dir
        self.temp_center = temp_center

        self.VINA_BIN = shutil.which('vina')
        if self.VINA_BIN is None:
            raise FileNotFoundError("vina not found in PATH. Please check your Conda environment.")
        
    def collect_tasks(self):
        tasks = []
        for root, _, files in os.walk(self.receptor_pdb_dir):
            for file in files:
                if file.endswith('.pdb'):
                    receptor_pdb_file = os.path.join(root, file)
                    gene_name = Path(root).name  

                    if gene_name not in self.gene_ligand_map:
                        continue


                    relative_subfolder = os.path.relpath(root, self.receptor_pdb_dir)
                    receptor_pdbqt_file = os.path.join(self.receptor_pdbqt_dir, relative_subfolder, file.replace(".pdb", ".pdbqt"))
                    if not Path(receptor_pdbqt_file).exists():
                        print(f"[ERROR] Missing receptor PDBQT: {receptor_pdbqt_file}")
                        continue

                    ligands = self.gene_ligand_map[gene_name]
                    for ligand_name in ligands:
                        ligand_pdbqt_file = os.path.join(self.ligand_pdbqt_dir, f"{ligand_name}.pdbqt")
                        if Path(ligand_pdbqt_file).exists():
                            tasks.append((receptor_pdb_file, receptor_pdbqt_file, ligand_pdbqt_file, gene_name, ligand_name))
                        else:
                            print(f"[WARN] Missing ligand PDBQT: {ligand_pdbqt_file}")
        return tasks

    def run_task(self, args):
        receptor_pdb_file, receptor_pdbqt_file, ligand_pdbqt_file, gene_name, ligand_name = args
        receptor_name = Path(receptor_pdb_file).stem

        fpocket_dir = os.path.join(self.temp_center, gene_name, f"{receptor_name}_out")
        if Path(fpocket_dir).exists():
            fpocket_out = fpocket_dir
            center, size = parse_fpocket_pqr(fpocket_out)
            if center is None:
                fpocket_out = run_fpocket(str(receptor_pdb_file), str(self.temp_center))
                if not fpocket_out:
                    print(f"[ERROR] fpocket failed: {receptor_name}")
                    return
                center, size = parse_fpocket_pqr(fpocket_out)
                if center is None:
                    print(f"[ERROR] fpocket parsing failed: {receptor_name}")
                    return
        else:
            fpocket_out = run_fpocket(str(receptor_pdb_file), str(self.temp_center), gene_name)
            if not fpocket_out:
                print(f"[ERROR] fpocket failed: {receptor_name}")
                return
            center, size = parse_fpocket_pqr(fpocket_out)
            if center is None:
                print(f"[ERROR] fpocket parsing failed: {receptor_name}")
                return

        out_dir = os.path.join(self.docking_dir, gene_name, f"{gene_name}__{receptor_name}__{ligand_name}")
        os.makedirs(out_dir, exist_ok=True)

        vina_config = os.path.join(out_dir, "vina_config.txt")
        with open(vina_config, 'w') as f:
            f.write(f"receptor = {receptor_pdbqt_file}\n")
            f.write(f"ligand = {ligand_pdbqt_file}\n")
            f.write(f"center_x = {center[0]:.3f}\ncenter_y = {center[1]:.3f}\ncenter_z = {center[2]:.3f}\n")
            f.write(f"size_x = {size[0]:.3f}\nsize_y = {size[1]:.3f}\nsize_z = {size[2]:.3f}\n")
            f.write(f"out = {out_dir}/result.pdbqt\nlog = {out_dir}/result.log\n")
            f.write(f"exhaustiveness = 32\nnum_modes = 20\n")

        try:
            subprocess.run([self.VINA_BIN, "--config", vina_config], check=True)
            print(f"[OK] Docking completed: {gene_name} x {ligand_name}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] vina failed: {e}")

    def run_all(self):
        for root, _, files in os.walk(self.receptor_pdb_dir):
            for file in files:
                if file.endswith('.pdb'):
                    receptor_pdb_file = os.path.join(root, file)
                    receptor_name = Path(file).stem
                    gene_name = os.path.basename(root)
                    fpocket_dir = os.path.join(self.temp_center, f"{receptor_name}_out")
                    if not Path(fpocket_dir).exists():
                        print(f"[INFO] Running fpocket: {receptor_name}")
                        run_fpocket(receptor_pdb_file, self.temp_center, gene_name)

        tasks = self.collect_tasks()
        with Pool(cpu_count()) as pool:
            pool.map(self.run_task, tasks)
