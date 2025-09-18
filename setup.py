from setuptools import setup, find_packages


setup(
    name="gene2struct",
    version="1.1.0",
    description="Gene2Struct: A modular toolkit for conservation scoring, dN/dS inference, and molecular docking.",
    author="Shishi",
    author_email="2024222040131@stu.scu.edu.cn",
    url="https://github.com/Shishi/gene2struct",  
    packages=find_packages(), 
    python_requires='>=3.8',
    install_requires=[
        "biopython",
        "matplotlib",
        "numpy",
        "pandas",
        "scipy",
        "seaborn",
        "psutil",
        "ete3",
        "requests"],
    entry_points={
        'console_scripts': [
            'gene2struct = gene2struct.cli:main',  
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
