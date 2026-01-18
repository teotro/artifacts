import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# === SETTINGS ===
layouts = ["sparse", "compact", "intermediate"]
benchmarks = [f"random_{i}" for i in range(1, 21)]
base_path = "targets"
num_perts = 30

# === PATTERN to extract DFS time ===
dfs_time_pattern = re.compile(r"Time taken for DFS:\s+([0-9.]+)\s+seconds")

# === STORE records ===
records = []

# === MAIN LOOP ===
for layout in layouts:
    for bench in benchmarks:
        for i in range(1, num_perts + 1):
            dir_path = os.path.join(base_path, layout, bench)
            txt_file = os.path.join(dir_path, f"pert_{i}.txt")
            err_file = os.path.join(dir_path, f"pert_{i}.txt.err")

            if os.path.isfile(err_file) or not os.path.isfile(txt_file):
                continue

            try:
                with open(txt_file, "r") as f:
                    content = f.read()
                match = dfs_time_pattern.search(content)
                if match:
                    time_taken = float(match.group(1))
                    records.append({
                        "Layout": layout,
                        "Benchmark": bench,
                        "DFS Time (s)": time_taken
                    })
            except Exception:
                continue

# === CONVERT TO DATAFRAME ===
df = pd.DataFrame(records)
df["Benchmark"] = df["Benchmark"].str.replace("random_", "").astype(int)

# === PLOT WITH BROKEN Y AXIS ===
sns.set(style="whitegrid", font_scale=1.35)
fig, (ax_low, ax_high) = plt.subplots(2, 1, sharex=True, figsize=(20, 8), gridspec_kw={'height_ratios': [1, 3]})

# === Plot upper (outlier) range ===
sns.boxplot(
    data=df,
    x="Benchmark",
    y="DFS Time (s)",
    hue="Layout",
    palette="Set2",
    width=0.5,
    ax=ax_low
)
sns.stripplot(
    data=df,
    x="Benchmark",
    y="DFS Time (s)",
    hue="Layout",
    dodge=True,
    palette="Set2",
    size=3,
    alpha=0.5,
    jitter=0.2,
    ax=ax_low
)
ax_low.set_ylim(5, 250)
ax_low.spines['bottom'].set_visible(False)
ax_low.tick_params(labelbottom=False, labelsize=22)
ax_low.set_ylabel("")  # ← removes y-axis label for top plot
ax_low.tick_params(axis='y', labelsize=22)

# === Plot lower (main) range ===
sns.boxplot(
    data=df,
    x="Benchmark",
    y="DFS Time (s)",
    hue="Layout",
    palette="Set2",
    width=0.5,
    ax=ax_high
)
sns.stripplot(
    data=df,
    x="Benchmark",
    y="DFS Time (s)",
    hue="Layout",
    dodge=True,
    palette="Set2",
    size=3,
    alpha=0.5,
    jitter=0.2,
    ax=ax_high
)
ax_high.set_ylim(0, 5)
ax_high.spines['top'].set_visible(False)
ax_high.set_title("DFS Execution Time per Benchmark and Layout (Broken Y-Axis)", fontsize=22, pad=14)
ax_high.set_xlabel("Synthetic Benchmark", fontsize=22, labelpad=10)
ax_high.set_ylabel("DFS Time (s)", fontsize=22, labelpad=10)
ax_high.tick_params(axis='x', labelsize=22)
ax_high.tick_params(axis='y', labelsize=22)

# === Break marks ===
d = .008
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax_low.plot([0, 1], [0, 0], transform=ax_low.transAxes, **kwargs)
ax_high.plot([0, 1], [1, 1], transform=ax_high.transAxes, **kwargs)

# === Fix legends ===
handles, labels = ax_high.get_legend_handles_labels()
ax_high.legend(handles[:3], labels[:3], title="Layout", fontsize=22, title_fontsize=24)
ax_low.get_legend().remove()

# === Final layout ===
plt.tight_layout()
plt.savefig("dfs_time_boxplot_broken_axis_fig13.pdf", dpi=300)
plt.show()
