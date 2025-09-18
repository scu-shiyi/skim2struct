# Gene2Struct

**Gene2Struct** is a comprehensive toolkit designed for functional gene analysis in non-model organisms.  
It enables researchers to extract coding sequences (CDS) from shallow sequencing data, rapidly perform phylogenetic analysis, evaluate conservation and positive selection at both site and branch levels, predict protein structures, and conduct large-scale molecular docking.  
This pipeline bridges phylogenomics and structural bioinformatics, providing an end-to-end solution from **gene** to **structure**.

---

## Features

- **CDS extraction** from shallow sequencing data
- **Phylogenetic analysis** with site- and branch-level selection tests
- **Conservation analysis** across genes and clades
- **Protein structure prediction** (compatible with external structure predictors)
- **Automated molecular docking** (batch docking of ligands and receptors)
- Designed for **non-model organisms**

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/scu-shiyi/Gene2Struct.git
cd Gene2Struct

2. Install MGLTools (required for docking preparation)
	•	Download mgltools_x86_64Linux2_1.5.7 from the official website.
	•	Place the extracted folder under the utils/ directory:
Gene2Struct/utils/mgltools_x86_64Linux2_1.5.7
You may also include a screenshot of this step for clarity:
![MGLTools installation screenshot](docs/images/mgltools_install.png)
