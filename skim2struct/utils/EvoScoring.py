import requests
import pandas as pd
import numpy as np
from Bio import SeqIO
import os
from pathlib import Path

class EvoScoring:

    def __init__(self, fasta_input, output_dir=None, outgroups=None):
        self.fasta_files = self._collect_fasta_files(fasta_input)
        self.output_dir = output_dir
        self.outgroups = outgroups or []

    def _collect_fasta_files(self, fasta_input):
        if isinstance(fasta_input, str):
            if os.path.isdir(fasta_input):
                return [str(Path(fasta_input) / f) for f in os.listdir(fasta_input)
                        if f.endswith(('.fa', '.fasta', '.fas', '.txt'))]
            else:
                return fasta_input.split(';')
        elif isinstance(fasta_input, list):
            return fasta_input
        else:
            raise ValueError("Invalid fasta input")

    def remove_outgroups(self, fasta_file, preserve_gap=False):
        records = list(SeqIO.parse(fasta_file, 'fasta'))
        filtered_records = [r for r in records if r.id not in self.outgroups]
        names = [r.id for r in filtered_records]

        if preserve_gap:
            seqs = [str(r.seq).upper() for r in filtered_records]
        else:
            seqs = [str(r.seq).replace('-', '').replace('?', '').upper() for r in filtered_records]

        return names, seqs, filtered_records

    def scoring(self):
        score_output = pd.DataFrame()
        out_file = None
        if self.output_dir:
            out_path = os.path.join(self.output_dir, "evo_dir")
            os.makedirs(out_path, exist_ok=True)
            out_file = os.path.join(out_path, "NNL_scores.csv")
            if os.path.exists(out_file):
                print("检测到已经计算好的evo文件")
                return out_file

        api_url = "https://3820-171-213-150-230.ngrok-free.app/score_only"

        for fasta_path in self.fasta_files:
            names, seqs, _ = self.remove_outgroups(fasta_path, preserve_gap=False)
            gene_name = Path(fasta_path).stem

            request_data = {'seq_contents': seqs}
            response = requests.post(api_url, json=request_data)

            if response.status_code == 200:
                scores = response.json()['results']
                new_data = pd.DataFrame({'name': names, gene_name: [-round(score, 2) for score in scores]})
                if score_output.empty:
                    score_output = new_data
                else:
                    score_output = pd.merge(score_output, new_data, on='name', how='outer')
            else:
                raise ValueError(f"API Error {response.status_code}: {response.text}")

        # 返回保存路径或DataFrame
        if out_file:
            score_output.to_csv(out_file, index=False)
            return out_file
        else:
            return score_output


    def align_scores_to_msa(self, msa_seq, score_list):
        result = []
        idx = 0
        for aa in msa_seq:
            if aa == "-":
                result.append(np.nan)
            else:
                result.append(score_list[idx])
                idx += 1
        return result

    def position_entropy(self):
        api_url = "https://3820-171-213-150-230.ngrok-free.app/entropy_only"
        if len(self.fasta_files) != 1:
            raise ValueError("position_entropy() 仅支持单文件处理。")

        fasta_path = self.fasta_files[0]
        names, seqs, records = self.remove_outgroups(fasta_path, preserve_gap=False)

        request_data = {'seq_contents': seqs}
        response = requests.post(api_url, json=request_data)

        if response.status_code == 200:
            entropies = response.json()['entropies']
            conservation = [[round(2 - entropy, 4) for entropy in entropy_sequence] for entropy_sequence in entropies]

            if not self.output_dir:
                return conservation
            else:
                seqs = [str(record.seq).upper() for record in records]
                aligned = {name: self.align_scores_to_msa(seq, score)
                           for name, seq, score in zip(names, seqs, conservation)}

                df_heatmap = pd.DataFrame.from_dict(aligned, orient="index")
                out_path = os.path.join(self.output_dir, "evo_dir")
                os.makedirs(out_path, exist_ok=True)
                out_file = os.path.join(out_path, "Conservation_score.csv")
                df_heatmap.to_csv(out_file, index=True)
                return out_file
        else:
            raise ValueError(f"API Error {response.status_code}: {response.text}")