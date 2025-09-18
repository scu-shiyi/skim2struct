import subprocess
import os



BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MGLTOOLS_DIR = os.path.join(BASE_DIR, 'utils/mgltools_x86_64Linux2_1.5.7')
PYTHONSH = os.path.join(MGLTOOLS_DIR, 'bin/pythonsh')


def run_subprocess(command):
    """Run a subprocess command with error handling."""
    try:
        subprocess.run(command, check=True)
        print(f"Conversion successful: {command[-1]}")
    except subprocess.CalledProcessError as e:
        print(f"Conversion failed: {e}")



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
        print(f"{file_name} was successfully converted and saved to {output_file}")
    return ligand_dir

def receptor2pdbqt(receptor_file, receptor_dir):

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
    print(f"Receptor was successfully converted and saved to {output_file}")
    return output_file





