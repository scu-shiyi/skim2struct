import pandas as pd
import numpy as np
from gene2struct.utils.TreeFunction import load_tree, compute_leaf_positions
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from Bio import Phylo
import pandas as pd
import numpy as np
from typing import Optional

def compute_universal_activity_matrix(
    csv_path: str, 
    out_csv_path: Optional[str] = None,
    activity_method: str = "ratio",
    standardize: bool = "False"
) -> pd.DataFrame:


    R_kcal = 1.987e-3      # kcal·mol⁻¹·K⁻¹
    T = 298                # K
    RT = R_kcal * T        # ≈ 0.592 kcal/mol
    
    df_raw = pd.read_csv(csv_path, header=None)
    gene_row = df_raw.iloc[0, 1:]
    role_row = df_raw.iloc[1, 1:].str.lower()
    ligand_row = df_raw.iloc[2, 1:]
    
    col_names = [f"{g}_{r}_{l}" for g, r, l in zip(gene_row, role_row, ligand_row)]
    
    df_energy = df_raw.iloc[3:].copy()
    df_energy.columns = ["Sample"] + col_names
    df_energy.set_index("Sample", inplace=True)
    
    tidy = (
        df_energy.stack()
                 .reset_index()
                 .rename(columns={0: "ΔG", "level_1": "GeneRole"})
    )
    tidy[["Gene", "Role", "Ligand"]] = tidy["GeneRole"].str.split("_", n=2, expand=True)
    tidy["ΔG"] = pd.to_numeric(tidy["ΔG"], errors="coerce") 
    tidy = tidy[(tidy["ΔG"].notna()) & (tidy["ΔG"] < 0)]
    
    results = []
    for (sample, gene), sub_df in tidy.groupby(["Sample", "Gene"]):
        prod_energy = sub_df[sub_df["Role"] == "product"]["ΔG"].values
        subs_energy = sub_df[sub_df["Role"] == "substrate"]["ΔG"].values

        if (
            len(prod_energy) == 0 or len(subs_energy) == 0 or
            np.any(prod_energy >= 0) or np.any(subs_energy >= 0)
        ):
            continue

        mean_prod = np.mean(prod_energy)
        mean_subs = np.mean(subs_energy)

        if activity_method == "ratio":
            activity = abs(mean_prod / mean_subs)
        elif activity_method == "delta":
            activity = mean_subs - mean_prod
        elif activity_method == "exp_delta":
            delta_g = mean_subs - mean_prod
            activity = np.exp(delta_g / RT)
        else:
            raise ValueError(f"Unknown activity_method: {activity_method}")

        results.append({
            "Sample": sample,
            "Gene": gene,
            "Activity": activity,
            "Mean_Substrate": mean_subs,
            "Mean_Product": mean_prod
        })
    
    df_activity = pd.DataFrame(results)
    
    if standardize and not df_activity.empty:
        df_activity["Activity_Z"] = df_activity.groupby("Gene")["Activity"].transform(
            lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
        )
    
    value_col = "Activity_Z" if standardize else "Activity"
    df_matrix = df_activity.pivot(index="Sample", columns="Gene", values=value_col)
    

    df_matrix = df_matrix.reset_index()
    

    if out_csv_path:
        df_matrix.to_csv(out_csv_path, index=False, float_format="%.6g")
    
    return df_matrix

def plot_clade(clade, ax, depths, max_depth, leaf_positions, tick=0.25):
    x0 = depths[clade]
    if clade.is_terminal():
        y = leaf_positions[clade.name]
        ax.hlines(y, x0, max_depth, colors='gray', linestyles='--', linewidth=0.8)
        ax.vlines(max_depth, y - tick, y + tick, colors='black', linewidth=0.8)
        return y

    ys = []
    for child in clade.clades:
        y_child = plot_clade(child, ax, depths, max_depth, leaf_positions, tick)
        ys.append(y_child)
        ax.hlines(y_child, x0, depths[child], colors='black', linestyles='-', linewidth=1)

    y_min, y_max = min(ys), max(ys)
    ax.vlines(x0, y_min, y_max, colors='black', linestyles='-', linewidth=1)
    return 0.5 * (y_min + y_max)

def draw_tree_and_heatmap(
    tree, depths, max_depth,
    leaf_positions, heatmap_df,
    picture_path: str,
    name_limit: int = 20,
    cmap: str = "RdBu_r",
    norm_label: str = "kd-based (norm)",
):
    fig = plt.figure(figsize=(25, 0.4 * len(leaf_positions) + 2), dpi=150)
    gs = GridSpec(1, 2, width_ratios=[1, 6], wspace=0.15)

    ax_tree = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1])

    plot_clade(tree.root, ax_tree, depths, max_depth, leaf_positions)
    ax_tree.set_ylim(-0.5, len(leaf_positions) - 0.5)
    ax_tree.set_xlim(0, max_depth * 1.02)
    ax_tree.invert_yaxis()
    ax_tree.axis("off")

    for name, y in leaf_positions.items():
        ax_tree.text(max_depth * 1.005, y, name[:name_limit], va='center', ha='left', fontsize=8)

    sns.heatmap(
        heatmap_df.values,
        ax=ax_heat,
        cmap=cmap,
        cbar_kws={"label": norm_label},
        linewidths=0.5,
        linecolor="white",
    )
    ax_heat.set_yticklabels([])
    ax_heat.tick_params(axis="y", which="both", labelleft=False, left=True)
    ax_heat.yaxis.set_ticks_position("left")
    ax_heat.set_yticks(np.arange(len(heatmap_df)) + 0.5)

    ax_heat.set_xticks(np.arange(len(heatmap_df.columns)) + 0.5)
    ax_heat.set_xticklabels(heatmap_df.columns, rotation=45, ha='right', fontsize=8)
    ax_heat.set_xlabel("Gene", fontsize=10)

    plt.tight_layout()
    plt.savefig(picture_path)
    plt.close(fig)



def plot(csv_path, out_csv_path, tree_path, pic_path):
    tree, depths, max_depths = load_tree(tree_path)
    leaf_position = compute_leaf_positions(tree)
    leaf_position_lower = {k.lower(): v for k, v in leaf_position.items()}
    kd_martrix = compute_universal_activity_matrix(csv_path, out_csv_path)
    df_heatmap = kd_martrix.set_index("Sample")
    df_heatmap = df_heatmap.loc[leaf_position_lower.keys()]

    draw_tree_and_heatmap(tree, depths, max_depths, leaf_position, df_heatmap, pic_path)
