from setuptools import setup, find_packages
import os
import urllib.request
import tarfile
from setuptools.command.install import install

MGL_URL = "https://ccsb.scripps.edu/download/532/"
MGL_ARCHIVE = "mgltools_x86_64Linux2_1.5.7.tar.gz"
MGL_DIR_NAME = "mgltools_x86_64Linux2_1.5.7"

class PostInstallMGL(install):
    """After standard install, download, extract, and patch MGLTools into utils."""
    def run(self):
        # Run standard install
        super().run()

        # Prepare paths
        base_dir = os.path.abspath(os.path.dirname(__file__))
        utils_dir = os.path.join(base_dir, 'skim2struct', 'utils')
        target_dir = os.path.join(utils_dir, MGL_DIR_NAME)

        # If not already installed, download and extract
        if not os.path.isdir(target_dir):
            os.makedirs(utils_dir, exist_ok=True)
            archive_path = os.path.join(utils_dir, MGL_ARCHIVE)
            print(f"Downloading MGLTools from {MGL_URL}...")
            urllib.request.urlretrieve(MGL_URL, archive_path)
            print("Extracting MGLTools...")
            with tarfile.open(archive_path) as tar:
                tar.extractall(path=utils_dir)

            # Patch prepare_ligand4.py
            patch_path = os.path.join(
                target_dir,
                'MGLToolsPckgs',
                'AutoDockTools',
                'Utilities24',
                'prepare_ligand4.py'
            )
            print(f"Patching {patch_path}...")
            with open(patch_path, 'r') as f:
                lines = f.readlines()
            new_lines = []
            for line in lines:
                if 'ligand_filename = os.path.basename(a)' in line and not line.strip().startswith('#'):
                    new_lines.append('# ' + line)
                elif 'ligand_filename = a' in line and line.strip().startswith('#'):
                    new_lines.append(line.lstrip('# '))
                else:
                    new_lines.append(line)
            with open(patch_path, 'w') as f:
                f.writelines(new_lines)

            print("MGLTools installed and patched under utils.")
        else:
            print("MGLTools already present, skipping download and patch.")


setup(
    name="skim2struct",
    version="1.0.0",
    description="Skim2Struct: A modular toolkit for conservation scoring, dN/dS inference, and molecular docking.",
    author="Shishi",
    author_email="2024222040131@stu.scu.edu.cn",
    url="https://github.com/Shishi/skim2struct",  # 可选
    packages=find_packages(),  # 自动查找 skim2struct 下所有子包
    python_requires='>=3.8',
    install_requires=[
        "biopython",
        "matplotlib",
        "numpy",
        "pandas",
        "scipy",
        "seaborn",
        "psutil",
        # 视你的模块依赖还可加入其他项
    ],
    entry_points={
        'console_scripts': [
            'skim2struct = skim2struct.cli:main',  # 👈 cli.py 中的 main 函数作为命令行入口
        ]
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
