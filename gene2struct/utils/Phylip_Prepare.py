import os
from pdb import run
import sys
import shutil
import warnings
from pathlib import Path
from Bio import SeqIO
from Bio import BiopythonDeprecationWarning
from ete3 import Tree
from skim2struct.utils.TreeLoad import TreeLoad  # 你已有的模块
import subprocess
from Bio.SeqRecord import SeqRecord

from Bio import SeqIO

def check_cds(fasta_file, outgroups=None):
    """
    Check CDS file:
      - Length must be a multiple of 3
      - No premature stop codons allowed
      - Remove terminal stop codon if present
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError(f"{fasta_file} contains no valid sequences")

    cleaned_records = []
    stops = {"TAA", "TAG", "TGA"}

    for rec in records:
        seq = str(rec.seq).upper().replace("-", "")
        if len(seq) % 3 != 0:
            raise ValueError(f"{rec.id}: sequence length is not a multiple of 3")

        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]

        # Check premature stops
        if any(codon in stops for codon in codons[:-1]):
            raise ValueError(f"{rec.id}: contains a premature stop codon")

        # Remove terminal stop codon if present
        if codons[-1] in stops:
            print(f"[Warning] {rec.id}: terminal stop codon removed")
            seq = seq[:-3]

        rec.seq = rec.seq.__class__(seq)  # keep Seq type
                # ---- 简单规范化 ID ----
        header = rec.description if rec.description else rec.id
        parts = header.split('|')
        species = parts[0].strip().replace(' ', '_')  # 取 | 前物种，空格→_
        pid = None
        for p in parts[1:]:
            if p.startswith('protein_id:'):
                pid = p.split(':', 1)[1].strip()
                break
        new_id = f"{species}_{pid}" if pid else species

        rec.id = rec.name = new_id
        rec.description = ""  # 关键：清空描述，避免空格再次混入
        cleaned_records.append(rec)
    
    # Remove outgroups if specified
    if outgroups:
        cleaned_records = [rec for rec in cleaned_records if rec.id not in outgroups]
    return cleaned_records

def run_pal2nal_revised(cds_fasta, output_prefix):
    """
    Generate codon alignment using MUSCLE (AA alignment) + PAL2NAL (back-translation).

    Args:
        cds_fasta: Input CDS FASTA file (DNA, coding sequences, aligned or unaligned)
        output_prefix: Prefix for output files

    Returns:
        (codon_alignment_path, aa_aln_path): paths to FASTA files
    """
    cds_path = Path(cds_fasta)
    out_prefix = Path(output_prefix)
    aa_fasta_path = out_prefix.with_suffix(".aa.fasta")
    aa_aln_path = out_prefix.with_suffix(".aa.aln.fasta")
    codon_alignment_path = out_prefix.with_suffix(".codon.aln.fasta")

    aa_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        # 1. Translate CDS → Protein
        new_records = []
        for rec in SeqIO.parse(cds_path, "fasta"):
            aa_seq = rec.seq.translate(to_stop=True)
            new_rec = SeqRecord(aa_seq, id=rec.id, description="")
            new_records.append(new_rec)
        SeqIO.write(new_records, aa_fasta_path, "fasta")


        # 2. Muscle alignment (protein)
        muscle_path = shutil.which("muscle")
        if not muscle_path:
            raise FileNotFoundError(
            "MUSCLE not found. Please install MUSCLE and make sure it is on your PATH, "
            "or set the MUSCLE environment variable to the absolute path of the executable.\n"
            "Examples:\n"
            "  conda install -c bioconda muscle\n"
            "  export MUSCLE=/path/to/muscle"
    )
        cmd_mafft = [muscle_path, "-align", str(aa_fasta_path), '-output', aa_aln_path]
        try:
            subprocess.run(cmd_mafft, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"MUSCLE failed with exit code {e.returncode}.\n"
                f"Command: {' '.join(cmd_mafft)}\n"
                f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
            ) from e
        
        if not aa_aln_path.exists() or aa_aln_path.stat().st_size == 0:
            raise RuntimeError("MUSCLE produced no output alignment file.")
        # 3. PAL2NAL: map back to codon alignment
        pal2nal_path = shutil.which("pal2nal.pl")
        if not pal2nal_path:
            raise FileNotFoundError(
                "PAL2NAL not found. Please install PAL2NAL and make sure it is on your PATH, "
                "or set the PAL2NAL environment variable to the absolute path of the executable.\n"
                "Examples:\n"
                "  conda install -c bioconda pal2nal\n"
            )
        with open(codon_alignment_path, "w") as fout:
            cmd_p2n = [pal2nal_path, str(aa_aln_path), str(cds_path), "-output", "fasta"]
            try:
                subprocess.run(cmd_p2n,stdout=fout,check=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"PAL2NAL failed. Exit code {e.returncode}.\n"
                    f"Command: {' '.join(cmd_p2n)}\n"
                    f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
                ) from e
        if not codon_alignment_path.exists() or codon_alignment_path.stat().st_size == 0:
            raise RuntimeError("PAL2NAL produced no output codon alignment file.")
        return str(codon_alignment_path), str(aa_aln_path)
    finally:
        # Clean up temporary files
        if aa_fasta_path.exists():
            aa_fasta_path.unlink()


def run_trimal(input_fasta, output_phy, mode='automated1'):
    trimal_exe = shutil.which('trimal')
    if not trimal_exe:
        raise FileNotFoundError(
            "trimalnot found. Please install trimAl and ensure it is on your PATH, "
            "or set the TRIMAL environment variable to the absolute path of the executable.\n"
            "Examples:\n"
            "  conda install -c bioconda trimal\n"
            "  export TRIMAL=/path/to/trimal"
        )
    if mode == "automated1":
        cmd = [trimal_exe, "-in", str(input_fasta), "-out", str(output_phy), "-phylip_paml", "-automated1"]
    else:
        cmd = [trimal_exe, "-in", str(input_fasta), "-out", str(output_phy), "-phylip_paml", "-gt", str(mode)]
    print(f"[trimAl] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return Path(output_phy)


def check_species_match_before_alignment(fasta_file, species_from_tree, outgroups=None):
    """
    Check if sequence IDs in fasta match species names in tree (before alignment).
    Allows removing outgroups first.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError(f"{fasta_file} contains no sequences")

    seq_ids = {rec.id for rec in records}

    # remove outgroups if specified
    if outgroups:
        seq_ids = seq_ids - set(outgroups)

    tree_species = set(species_from_tree)

    if seq_ids != tree_species:
        missing_in_seq = tree_species - seq_ids
        missing_in_tree = seq_ids - tree_species
        msg = []
        if missing_in_seq:
            msg.append(f"Missing in fasta: {missing_in_seq}")
        if missing_in_tree:
            msg.append(f"Missing in tree: {missing_in_tree}")
        raise ValueError("Tree and fasta species do not match!\n" + "\n".join(msg))
    else:
        print("[Check] Sequence IDs match tree species.")


def prepare_paml_input1(fasta_path, output_dir, tree_path=None, outgroups=None):

    fasta_name = Path(fasta_path).stem if os.path.isfile(fasta_path) else Path(fasta_path).name
    work_dir = Path(output_dir) / "file_input" /fasta_name
    work_dir.mkdir(parents=True, exist_ok=True)
    # if work_dir.exists():
    #     shutil.rmtree(work_dir)
    # os.makedirs(work_dir, exist_ok=True)

    tmp_clean_fasta = work_dir / "cleaned.fasta"
    cleaned_records = check_cds(fasta_path, outgroups)
    SeqIO.write(cleaned_records, tmp_clean_fasta, "fasta")

    codon_alignment_path, aa_aln_path = run_pal2nal_revised(tmp_clean_fasta, work_dir / fasta_name)
    fasta_phy_path = run_trimal(codon_alignment_path, work_dir / f"{fasta_name}.phy", mode="automated1")
    if tmp_clean_fasta.exists():
        tmp_clean_fasta.unlink()

    # --- tree ---
    if tree_path:
        treeid = work_dir / Path(tree_path).with_suffix('.paml.tree').name
        treeloader = TreeLoad(tree_path)
        species, newick_tree = treeloader.parse_tree()

        if outgroups:
            outgroup_leaves = []
            without_outgroup_tree = Tree(newick_tree, format=1)
            for outgroup in outgroups:
                if outgroup in species:
                    species.remove(outgroup)
                leaves = without_outgroup_tree.get_leaves_by_name(outgroup)
                if leaves:
                    outgroup_leaves.extend(leaves)
            if outgroup_leaves:
                remaining_leaves = [leaf for leaf in without_outgroup_tree.get_leaves() if leaf not in outgroup_leaves]
                without_outgroup_tree.prune(remaining_leaves)
            without_outgroup_tree.unroot()
            n_species = len(species)
            header = f"   {n_species}  1\n"
            with open(treeid, 'w') as f:
                f.write(header)
                f.write(without_outgroup_tree.write(format=1))
        else:
            n_species = len(species)
            header = f"   {n_species}  1\n"
            with open(treeid, 'w') as f:
                f.write(header)
                f.write(newick_tree)
        check_species_match_before_alignment(fasta_phy_path, species, outgroups)
    else:
        treeid = work_dir / f"{fasta_name}.paml.tree"
        reused = False
        for cand in (treeid, work_dir / f"{fasta_name}.treefile"):
            if cand.exists() and cand.stat().st_size > 0:
                treeloader = TreeLoad(cand)
                species, newick_tree = treeloader.parse_tree()  


                n_species = len(species)
                header = f"   {n_species}  1\n"
                with open(treeid, 'w') as f:
                    f.write(header)
                    f.write(newick_tree)
                reused = True
                break

        if not reused:
            cmd = [
                "iqtree2",
                "-s", str(aa_aln_path),
                "-m", "MFP",
                "-T", '5',
                "-fast",
                '-pre', str(work_dir / fasta_name)
            ]
            subprocess.run(cmd, cwd=work_dir, check=True)


            raw_tree = work_dir / f"{fasta_name}.treefile"
            treeloader = TreeLoad(raw_tree)
            species, newick_tree = treeloader.parse_tree()

            n_species = len(species)
            header = f"   {n_species}  1\n"
            with open(treeid, 'w') as f:
                f.write(header)
                f.write(newick_tree)

            for suffix in ["model.gz", "log", "iqtree", "ckp.gz", "mldist", "bionj"]:
                f = work_dir / f"{fasta_name}.{suffix}"
                if f.exists():
                    f.unlink()


    return work_dir, fasta_phy_path, treeid, species


def prepare_paml_input2(fasta_path, output_dir, tree_path=None, outgroups=None):
    """
    Prepare input for PAML:
      - check CDS
      - codon alignment by MACSE
      - trim with trimAl and convert to PHYLIP
      - prepare tree file
    """

    fasta_name = Path(fasta_path).stem if os.path.isfile(fasta_path) else Path(fasta_path).name
    work_dir = Path(output_dir) / fasta_name / "file_input"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)


    tmp_clean_fasta = work_dir / "cleaned.fasta"
    cleaned_records = check_cds(fasta_path, outgroups)
    SeqIO.write(cleaned_records, tmp_clean_fasta, "fasta")

    codon_alignment_path, aa_aln_path = run_pal2nal_revised(tmp_clean_fasta, work_dir / fasta_name)
    fasta_phy_path = run_trimal(codon_alignment_path, work_dir / f"{fasta_name}.phy", mode="automated1")
    if tmp_clean_fasta.exists():
        tmp_clean_fasta.unlink()


    if tree_path:
        treeid = work_dir / Path(tree_path).with_suffix('.paml.tree').name
        treeloader = TreeLoad(tree_path)
        species, newick_tree = treeloader.parse_tree()

        if outgroups:
            outgroup_leaves = []
            without_outgroup_tree = Tree(newick_tree, format=1)
            for outgroup in outgroups:
                if outgroup in species:
                    species.remove(outgroup)
                leaves = without_outgroup_tree.get_leaves_by_name(outgroup)
                if leaves:
                    outgroup_leaves.extend(leaves)
            if outgroup_leaves:
                remaining_leaves = [leaf for leaf in without_outgroup_tree.get_leaves() if leaf not in outgroup_leaves]
                without_outgroup_tree.prune(remaining_leaves)
            without_outgroup_tree.unroot()
            n_species = len(species)
            header = f"   {n_species}  1\n"
            with open(treeid, 'w') as f:
                f.write(header)
                f.write(without_outgroup_tree.write(format=1))
        else:
            n_species = len(species)
            header = f"   {n_species}  1\n"
            with open(treeid, 'w') as f:
                f.write(header)
                f.write(newick_tree)
        check_species_match_before_alignment(fasta_phy_path, species, outgroups)
    else:
        treeid = work_dir / f"{fasta_name}.paml.tree"


        cmd = [
            "iqtree2",
            "-s", str(aa_aln_path),
            "-m", "MFP",
            "-T", '5',
            "-fast",
            '-pre', str(work_dir / fasta_name)
        ]
        subprocess.run(cmd, cwd=work_dir, check=True)


        raw_tree = work_dir / f"{fasta_name}.treefile"
        treeloader = TreeLoad(raw_tree)
        species, newick_tree = treeloader.parse_tree()

        n_species = len(species)
        header = f"   {n_species}  1\n"
        with open(treeid, 'w') as f:
            f.write(header)
            f.write(newick_tree)

        for suffix in ["model.gz", "log", "iqtree", "ckp.gz", "mldist", "bionj"]:
            f = work_dir / f"{fasta_name}.{suffix}"
            if f.exists():
                f.unlink()


    return work_dir, fasta_phy_path, treeid, species



