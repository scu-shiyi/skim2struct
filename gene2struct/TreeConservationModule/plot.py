# gene2struct/TreeConservationModule/plot.py

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.gridspec import GridSpec
from pathlib import Path
import io
import re
import pandas as pd
from matplotlib.patches import Patch
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



def plot_clade(clade, ax, depths, max_depth, leaf_positions):
    x0 = depths[clade]
    if clade.is_terminal():
        y = leaf_positions[clade.name]
        ax.hlines(y, x0, max_depth, colors='gray', linestyles='--', linewidth=0.8)
        return y

    ys = []
    for child in clade.clades:
        y_child = plot_clade(child, ax, depths, max_depth, leaf_positions)
        ys.append(y_child)
        ax.hlines(y_child, x0, depths[child], colors='black', linestyles='-', linewidth=1)

    y_min, y_max = min(ys), max(ys)
    ax.vlines(x0, y_min, y_max, colors='black', linestyles='-', linewidth=1)
    return 0.5 * (y_min + y_max)



def _collapse_mean_by_k_cols(arr, k=3, strict=True):

    nrows, ncols = arr.shape
    if strict and (ncols % k != 0):
        raise ValueError(f"Number of columns {ncols} is not divisible by {k}, cannot collapse by codons.")
    pad = (-ncols) % k
    if pad:
        arr = np.concatenate([arr, np.full((nrows, pad), np.nan)], axis=1)
    g = arr.reshape(nrows, -1, k)          # (rows, groups, k)
    return np.nanmean(g, axis=2)           # (rows, groups)


def parse_mlc_beb_all(mlc_path: str) -> pd.DataFrame:

    text = Path(mlc_path).read_text(encoding="utf-8", errors="ignore")
    m = re.search(
        r'Bayes Empirical Bayes.*?Positively selected sites.*?(.*?)\n\s*The grid',
        text, flags=re.S | re.I
    )
    beb_block = m.group(1) if m else None
    line_pat = re.compile(
        r'^\s*(\d+)\s+([A-Za-z\*\-])\s+([01](?:\.\d+)?)'  
        r'(\*{1,2})?'                                      
        r'(?:\s+([0-9.]+)\s*\+\-\s*([0-9.]+))?'            
        r'\s*(\*{1,2})?\s*$',                             
        flags=re.M
    )

    rows = []

    def _parse_lines(block: str):
        for s in block.splitlines():
            s = s.rstrip("\r\n")
            mm = line_pat.match(s)
            if mm:
                pos   = int(mm.group(1))
                aa    = mm.group(2)
                prob  = float(mm.group(3))
                stars = (mm.group(4) or mm.group(7) or "")
                wmean = float(mm.group(5)) if mm.group(5) else None
                wse   = float(mm.group(6)) if mm.group(6) else None
                rows.append((pos, aa, prob, wmean, wse, stars))

    if beb_block is not None:
        _parse_lines(beb_block)
    else:
        in_beb = False
        for s in text.splitlines():
            if not in_beb and "Bayes Empirical Bayes" in s:
                in_beb = True
                continue
            if in_beb and "The grid" in s:
                break
            if not in_beb:
                continue
            if (("Positively selected sites" in s) or
                ("amino acids refer" in s) or
                ("Pr(w>1)" in s)):
                continue
            _parse_lines(s)

    df = pd.DataFrame(rows, columns=["pos","aa","prob","w_mean","w_se","stars"])
    if not df.empty:
        df.sort_values("pos", inplace=True, ignore_index=True)
        df["is_p95"] = df["prob"] >= 0.95
        df["is_p99"] = df["prob"] >= 0.99
    return df

def add_beb_track_all(ax_heatmap,
                      beb_df,
                      ncols,
                      bar_color='#93c5fd',
                      hi_color='#ef4444',
                      dash_color='#94a3b8',
                      gap=0.012,            
                      star_fs=8,
                      star_offset_pt=2):          

    fig = ax_heatmap.figure
    box = ax_heatmap.get_position()


    h = box.height * 0.08
    ax_b = fig.add_axes([box.x0,
                         box.y0 - h - gap,  
                         box.width,
                         h],
                        sharex=ax_heatmap)


    pr = np.zeros(int(ncols), dtype=float)
    if beb_df is not None and not beb_df.empty:
        for pos, prob in zip(beb_df['pos'].astype(int),
                             beb_df['prob'].astype(float)):
            i = pos - 1  # 1-based → 0-based
            if 0 <= i < ncols:
                pr[i] = prob

    x = np.arange(ncols)


    ax_b.bar(x, pr, width=0.9, color=bar_color, edgecolor='none',
             align='center', zorder=1)
    hi95_idx = np.where((pr >= 0.95) & (pr < 0.99))[0]
    hi99_idx = np.where(pr >= 0.99)[0]
    if hi95_idx.size:
        ax_b.bar(hi95_idx, pr[hi95_idx], width=0.9, color=hi_color,
                 edgecolor='none', align='center', zorder=2)
    if hi99_idx.size:
        ax_b.bar(hi99_idx, pr[hi99_idx], width=0.9, color=hi_color,
                 edgecolor='none', align='center', zorder=2)


    def _star_y(y):
        return min(1.001, y + 0.005) 
    ax_b.set_ylim(0, 1.10)
    for i in hi95_idx:
        ax_b.annotate(
        "*", xy=(i, pr[i]), xytext=(0, 0.2),
        textcoords="offset points", ha="center", va="bottom",
        fontsize=star_fs, color="k", clip_on=False, zorder=4
        )
    for i in hi99_idx:
        ax_b.annotate(
            "**", xy=(i, pr[i]), xytext=(0, 0.2), 
            textcoords="offset points", ha="center", va="bottom",
            fontsize=star_fs, color="k", clip_on=False, zorder=4
        )


    ax_b.axhline(0.95, ls='--', lw=0.8, color=dash_color)
    ax_b.axhline(0.99, ls='--', lw=0.8, color=dash_color)


    ax_b.set_ylim(0, 1.05)
    ax_b.set_ylabel("Pr(ω>1)", fontsize=9)
    ax_b.grid(axis='y', ls=':', lw=0.5, color='#e2e8f0')


    ax_b.tick_params(axis='x', bottom=True, labelbottom=True, labelrotation=60, pad=2)
    plt.setp(ax_b.get_xticklabels(), ha='right', rotation_mode='anchor')


    for s in ("right", "top"):
        ax_b.spines[s].set_visible(False)
def draw_tree_and_heatmap(
    tree, depths, max_depth,
    leaf_positions, heatmap_df,
    output_dir: str, gene_name: str,
    site_model,
    mlc_path = None,
    name_limit: int = 20,
    thr: float = 1.4
) -> str:

    
    fig = plt.figure(figsize=(25, 10), dpi=150)
    gs = GridSpec(1, 2, width_ratios=[1, 6], wspace=0.15)

    ax_tree = fig.add_subplot(gs[0])
    ax_heatmap = fig.add_subplot(gs[1])


    plot_clade(tree.root, ax_tree, depths, max_depth, leaf_positions)
    ax_tree.set_ylim(-0.5, len(leaf_positions) - 0.5)
    ax_tree.set_xlim(0, max_depth * 1.02)
    # ax_tree.invert_yaxis()
    ax_tree.axis('off')

    for name, y in leaf_positions.items():
        ax_tree.text(max_depth * 1.005, y, name[:name_limit], va='center', fontsize=8)


    from matplotlib.colors import ListedColormap, BoundaryNorm
    raw = heatmap_df.values.astype(float)
    agg = _collapse_mean_by_k_cols(raw, k=3, strict=True)  


    cat = np.full_like(agg, np.nan, dtype=float)
    mask = ~np.isnan(agg)
    cat[mask] = (agg[mask] >= thr).astype(int)
    # cmap = ListedColormap(["#a9be7b", '#b04552'])  
    cmap = ListedColormap(["#a9be7b", '#d24735'])  
    cmap.set_bad('#e6e6e6')  
    norm = BoundaryNorm([0, 0.5, 1.5], cmap.N)


    ax_heatmap.set_facecolor('white')
    im = ax_heatmap.imshow(cat, aspect='auto', cmap=cmap, norm=norm,
                        interpolation='nearest', origin='lower')

    
    def add_grid_gaps(ax, nrows, ncols,
                    col_every=1, row_every=1,
                    col_lw=0.3, row_lw=0.3,
                    color='#f2f2f2'):
        if col_every and col_every > 0:
            xs = np.arange(col_every, ncols, col_every) - 0.5
            for x in xs:
                ax.axvline(x, color=color, lw=col_lw, zorder=3)


        if row_every and row_every > 0:
            ys = np.arange(row_every, nrows, row_every) - 0.5
            for y in ys:
                ax.axhline(y, color=color, lw=row_lw, zorder=3)


    add_grid_gaps(ax_heatmap, nrows=cat.shape[0], ncols=cat.shape[1],
              col_every=1, row_every=1,     
              col_lw=0.3, row_lw=0.3,
              color='#f2f2f2')


    legend_elements = [
        Patch(facecolor="#a9be7b", edgecolor='k', label=f'Highly conserved (<{thr})'),
        Patch(facecolor="#d24735", edgecolor='k', label=f'Poorly conserved (≥{thr})')
    ]
    ax_leg = inset_axes(
        ax_heatmap,
        width="18%", height="22%",             
        bbox_to_anchor=(1.01, 0.78, 0.1, 0.22),
        bbox_transform=ax_heatmap.transAxes,    
        loc="upper left", borderpad=0.0
    )
    ax_leg.axis("off")
    leg = ax_leg.legend(
        handles=legend_elements,
        title="  Amino acid site  \nconservation (bits)",
        loc="upper left", frameon=True, fancybox=True,
        framealpha=0.9, edgecolor="0.3",
        fontsize=10, title_fontsize=11,
        handlelength=1.4, handleheight=1.0, borderaxespad=0.0
    )
    if site_model:
        step = max(1, cat.shape[1] // 60)
        ticks = np.arange(0, cat.shape[1], step)
        ax_heatmap.set_xticks(ticks)
        ax_heatmap.set_xticklabels(ticks, rotation=60, fontsize=6)
        ax_heatmap.tick_params(
        axis='x',
        top=False, labeltop=False,      
        bottom=False, labelbottom=False, 
        pad=2                          
        )
    else:
        ncols = cat.shape[1]
        ticks = np.arange(0, ncols, 30)
        ax_heatmap.set_xticks(ticks)
        ax_heatmap.tick_params(
        axis='x',
        top=False, labeltop=False,
        bottom=True, labelbottom=True,
    )

    from matplotlib.ticker import FixedLocator, NullFormatter

    nrows = cat.shape[0]                      
    ax_heatmap.set_ylim(-0.5, nrows - 0.5)     
    ax_heatmap.yaxis.set_major_locator(FixedLocator(np.arange(nrows))) 
    ax_heatmap.yaxis.set_major_formatter(NullFormatter())               
    ax_heatmap.tick_params(axis='y', which='major',
                        left=True, right=False, length=2.5, width=0.8,
                        labelleft=False)
    
    if mlc_path:
        beb_df = parse_mlc_beb_all(mlc_path)      
        ncols = cat.shape[1] 
        add_beb_track_all(ax_heatmap, beb_df, ncols=ncols, gap=0.005, star_fs=8)
    plt.subplots_adjust(bottom=0.12) 
        
    # plt.tight_layout()
    output_path = Path(output_dir) / f"{gene_name}.png"
    plt.savefig(output_path)
    plt.close(fig)



    return str(output_path)
