import subprocess
import os


# 动态获取项目路径
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MGLTOOLS_DIR = os.path.join(BASE_DIR, 'utils/mgltools_x86_64Linux2_1.5.7')
PYTHONSH = os.path.join(MGLTOOLS_DIR, 'bin/pythonsh')


def run_subprocess(command):
    try:
        subprocess.run(command, check=True)
        print(f"转换成功: {command[-1]}")
    except subprocess.CalledProcessError as e:
        print(f"转换失败: {e}")


def ligand2pdbqt(ligand_input, ligand_dir):
    prepare_ligand_script = os.path.join(MGLTOOLS_DIR, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py")
    ligand_files = [ligand_input] if os.path.isfile(ligand_input) else [
        os.path.join(ligand_input, f) for f in os.listdir(ligand_input)
        if f.endswith(('.pdb'))
    ]
    for ligand_file in ligand_files:
        file_name = os.path.splitext(os.path.basename(ligand_file))[0] + ".pdbqt"
        output_file = os.path.join(ligand_dir, file_name)
        command = [PYTHONSH, prepare_ligand_script,
                   "-l", ligand_file,
                   "-o", output_file,
                   "-U", "waters",
                   "-A", "hydrogens"
                   ]
        run_subprocess(command)
    # print(f'所有pbdqt格式的配体文件储存在{ligand_dir}')
    print(f'{file_name}转化成功，储存在{output_file}')
    return ligand_dir

def receptor2pdbqt(receptor_file, receptor_dir):
    """只支持输入单个文件,由于计算活性中心需要对pdb文件处理"""
    """
    :param receptor_file: 输入pdb文件完整路径
    :param output_dir: 输出目录（自动在目录内创建同名pdbqt文件）
    :return: pdbqt文件完整路径
    """
    prepare_receptor_script = os.path.join(MGLTOOLS_DIR, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py")
    file_name = os.path.splitext(os.path.basename(receptor_file))[0] + ".pdbqt"
    output_file = os.path.join(receptor_dir, file_name)
    if os.path.exists(output_file):
        print(f"pass: {file_name}")
        return output_file
    command = [PYTHONSH, prepare_receptor_script,
               "-r", receptor_file,
               "-o", output_file,
               "-U", "waters",
               "-A", "hydrogens"]
    run_subprocess(command)
    print(f'{output_file}转化成功!')
    return output_file

# if __name__ == "__main__":
#     ligand2pdbqt('/home/shiyi/Autodock/test_data/UDP-Glc.pdb', '/home/shiyi/Autodock/Outputs/ligand_pdbqt')
#     receptor2pdbqt('/home/shiyi/Autodock/test_data/1.pdb', '/home/shiyi/Autodock/Outputs/receptor_pdbqt')
# # class PreparePDBQT:
#     """将配体和受体由pdb转化为pdbqt文件"""
#     def __init__(self,receptor_input, ligand_input, output_dir ):
#         self.mgltools_path = "/home/shiyi/mgltools_x86_64Linux2_1.5.7"
#         self.receptor_input = receptor_input
#         self.ligand_input = ligand_input
#         # 配体文件夹和受体文件夹
#         self.receptor_dir = os.path.join(output_dir, 'receptor_pdbqt')
#         self.ligand_dir = os.path.join(output_dir, 'ligand_pdbqt')
#         os.makedirs(self.receptor_dir, exist_ok=True)
#         os.makedirs(self.ligand_dir, exist_ok=True)
#
#         self.python_exec = os.path.join(self.mgltools_path, "bin/pythonsh")
#
#     def _run_subprocess(self, command):
#         try:
#             subprocess.run(command, check=True)
#             print(f"转换成功: {command[-1]}")
#         except subprocess.CalledProcessError as e:
#             print(f"转换失败: {e}")
#
#     def ligand_to_pdbqt(self):
#         """使用 MGLTools 的 prepare_ligand4.py 将配体文件转换为 PDBQT 格式。"""
#         prepare_ligand_script = os.path.join(self.mgltools_path, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py")
#         ligand_files = [self.ligand_input] if os.path.isfile(self.ligand_input) else [
#             os.path.join(self.ligand_input, f) for f in os.listdir(self.ligand_input)
#             if f.endswith(('.pdb', '.mol2'))
#         ]
#         for ligand_file in ligand_files:
#             file_name = os.path.splitext(os.path.basename(ligand_file))[0] + ".pdbqt"
#             output_file = os.path.join(self.ligand_dir, file_name)
#             command = [self.python_exec, prepare_ligand_script,
#                 "-l", ligand_file,
#                 "-o", output_file,
#                 "-U", "waters",
#                 "-A", "hydrogens"
#             ]
#             self._run_subprocess(command)
#
#     def receptor_to_pdbqt(self):
#         prepare_receptor_script = os.path.join(self.mgltools_path, "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py")
#
#         receptor_files = [self.receptor_input] if os.path.isfile(self.receptor_input) else [
#             os.path.join(self.receptor_input, f) for f in os.listdir(self.receptor_input)
#             if f.endswith(('.pdb', '.mol2'))
#         ]
#
#         for receptor_file in receptor_files:
#             file_name = os.path.splitext(os.path.basename(receptor_file))[0] + ".pdbqt"
#             output_file = os.path.join(self.receptor_dir, file_name)
#             command = [self.python_exec, prepare_receptor_script,
#                        "-r", receptor_file,
#                        "-o", output_file,
#                        "-U", "waters",
#                        "-A", "hydrogens"]
#             self._run_subprocess(command)


# if __name__ == "__main__":
    # pdbqt_converter = PreparePDBQT('/home/shiyi/Autodock/test_data/1.pdb', "/home/shiyi/Autodock/test_data/UDP-Glc.pdb", '/home/shiyi/Autodock/Outputs')
    # pdbqt_converter.ligand_to_pdbqt()
    # pdbqt_converter.receptor_to_pdbqt()




