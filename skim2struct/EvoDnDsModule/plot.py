# skim2struct/EvoDnDsModule/plot.py
import matplotlib.pyplot as plt
import numpy as np

def generate_summary_plot(combined_df, out_png_path):
    """
    输入：包含 'Gene', 'EvoMean', 'EvoStd', 'EvoValues', 'DnDs_Score' 的 DataFrame
    输出：将图像保存为 PNG 文件到 out_png_path
    """
    if combined_df is None or combined_df.empty:
        print("⚠ 无数据可绘图")
        return

    fig, ax = plt.subplots(figsize=(16, 6))
    bar_width = 0.35
    dnds_bar_width = bar_width / 1.5
    x = np.arange(len(combined_df))

    ax.bar(
        x + bar_width / 2,
        combined_df['EvoMean'],
        bar_width,
        label='Evo',
        edgecolor='red',
        fill=False,
        yerr=combined_df['EvoStd'],
        capsize=5,
        error_kw={'alpha': 1, 'ecolor': 'red', 'elinewidth': 1}
    )

    if 'DnDs_Score' in combined_df.columns and combined_df['DnDs_Score'].notnull().any():
        ax.bar(
            x + bar_width / 2,
            combined_df['DnDs_Score'],
            dnds_bar_width,
            label='dN/dS',
            color='#6699CC'
        )

    for i, values in enumerate(combined_df['EvoValues']):
        x_scatter = np.random.normal(x[i] + bar_width / 2, 0.05, len(values))
        ax.scatter(x_scatter, values, color='black', s=3)

    ax.set_ylabel('Value')
    ax.set_xlabel('Gene Name')
    ax.set_title('dN/dS Scores and Evo Scores per Gene')
    ax.set_xticks(x + bar_width / 2)
    ax.set_xticklabels(combined_df['Gene'], rotation=45, ha="right")

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())
    fig.tight_layout()

    plt.savefig(out_png_path)
    plt.close(fig)
    print(f"✅ 图像已保存至：{out_png_path}")