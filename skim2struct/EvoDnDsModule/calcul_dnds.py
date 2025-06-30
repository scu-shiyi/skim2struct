# 通过hyphy的slac模式计算dn/ds评分
import os
import subprocess
import re
import pandas as pd
import psutil
from pathlib import Path
from Bio import SeqIO, Phylo
import shutil

class CalculDnDs:
    def __init__(self, fasta_file, output_dir, tree_file, outgroups=None):
        self.fasta_file = fasta_file
        self.output_dir = output_dir
        self.tree_file = tree_file
        self.outgroups = outgroups or []

        self.temp_hyphy_dir = os.path.join(self.output_dir, "dnds_dir", "temp_hyphy")
        os.makedirs(self.temp_hyphy_dir, exist_ok=True)
        self.dnds_dir = os.path.join(self.output_dir, "dnds_dir")
        os.makedirs(self.dnds_dir, exist_ok=True)

        self.hyphy_exec = shutil.which("hyphy")

    def run_single_return(self):
        label = Path(self.fasta_file).stem
        dnds_value = self._calculate_dnds(self.fasta_file, self.tree_file)
        return dnds_value

    def _calculate_dnds(self, fasta_file, tree_file):
        label = Path(fasta_file).stem
        log_file = os.path.join(self.temp_hyphy_dir, label + '.log')
        json_file = os.path.join(self.temp_hyphy_dir, label + '.json')
        physical_cores = psutil.cpu_count(logical=False)
        threads = int(physical_cores * 1.2)
        if os.path.isfile(log_file):
            try:
                existing_value = self._parse_dnds(log_file)
                if existing_value is not None:
                    print(f"已检测到旧的 log 文件，直接读取 {label}基因的dN/dS")
                    return existing_value
            except Exception:
                pass
        try:
            subprocess.run([
                self.hyphy_exec, "slac",
                "--alignment", fasta_file,
                "--tree", tree_file,
                "--threads", str(threads),
                "--output", json_file
            ], stdout=open(log_file, "w"), stderr=subprocess.STDOUT, check=True)

            return self._parse_dnds(log_file)
        except subprocess.CalledProcessError:
            alignment_count, tree_count = self._check_alignment_tree_consistency(fasta_file, tree_file)
            if alignment_count != tree_count:
                print(f"⚠ 警告：序列文件中物种数量 = {alignment_count}，进化树中物种数量 = {tree_count}，二者不一致！")
            else:
                print(f"⚠ {label}HyPhy失败，但序列与树的物种数一致，可能为其他参数错误")
            return None

    def _check_alignment_tree_consistency(self, fasta_file, tree_file):
        alignment_count = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        tree = Phylo.read(tree_file, "newick")
        tree_count = len(tree.get_terminals())
        return alignment_count, tree_count

    def _parse_dnds(self, logfile):
        pattern = re.compile(r'non-synonymous/synonymous rate ratio.*?=\s*([-+]?\d*\.?\d+)', re.IGNORECASE)
        with open(logfile) as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    return float(match.group(1))
        return None