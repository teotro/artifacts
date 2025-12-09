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

# === PATTERN to extract 'A' or 'U' from ambig column ===
ambig_pattern = re.compile(r"\|\s.*\|\s+\d+\s+\|\s+\d+\s+\|\s+([AU])\s+\|")

# === STORE full perturbation-level ambiguity values ===
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
                    lines = f.readlines()
                ambig_flags = [m.group(1) for line in lines if (m := ambig_pattern.search(line))]
                if ambig_flags:
                    ambiguous_count = ambig_flags.count("A")
                    total_count = len(ambig_flags)
                    percentage = 100 * ambiguous_count / total_count
                    records.append({
                        "Layout": layout,
                        "Benchmark": bench,
                        "Ambiguity (%)": percentage
                    })
            except Exception:
                continue

# === CONVERT TO DATAFRAME ===
df = pd.DataFrame(records)
df["Benchmark"] = df["Benchmark"].str.replace("random_", "").astype(int)

# === PLOT BOXPLOT ===
plt.figure(figsize=(20, 7))
sns.set(style="whitegrid", font_scale=1.35)

# Create boxplot
ax = sns.boxplot(
    data=df,
    x="Benchmark",
    y="Ambiguity (%)",
    hue="Layout",
    palette="Set2",
    width=0.5
)

# Optional: add jittered data points on top of boxes
sns.stripplot(
    data=df,
    x="Benchmark",
    y="Ambiguity (%)",
    hue="Layout",
    dodge=True,
    palette="Set2",
    size=3,
    alpha=0.5,
    jitter=0.2,
    ax=ax
)

# Fix legend duplication from stripplot
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:3], labels[:3], title="Layout", fontsize=26, title_fontsize=28)

# === STYLING ===
plt.title("Ambiguity Distribution per Benchmark and Layout (Box Plot)", fontsize=35, pad=12)
plt.xlabel("Synthetic Benchmark", fontsize=35, labelpad=12)
plt.ylabel("Ambiguity Percentage (%)", fontsize=35, labelpad=12)
plt.xticks(rotation=0, fontsize=35)
plt.yticks(fontsize=35)
plt.tight_layout()

# === SAVE AND SHOW ===
plt.savefig("ambiguity_boxplot_fig12.pdf", dpi=300)
plt.show()