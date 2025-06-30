from setuptools import setup, find_packages

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