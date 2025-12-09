import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Settings
layouts = ["sparse", "compact", "intermediate"]
benchmarks = [f"random_{i}" for i in range(1, 21)]
display_benchmarks = [f"random_mix_{i}" for i in range(1, 21)]
base_dir = "targets"
total_perts = 30

# Collect success rate data
data = {}

for layout in layouts:
    rates = []
    for bench in benchmarks:
        path = os.path.join(base_dir, layout, bench)
        if not os.path.isdir(path):
            rates.append(np.nan)
            continue
        err_count = len([
            f for f in os.listdir(path)
            if f.endswith(".err") and os.path.isfile(os.path.join(path, f))
        ])
        success = 100 - (100 * err_count / total_perts)
        rates.append(round(success, 2))
    data[layout] = rates

# Create DataFrame with updated y-axis labels
df = pd.DataFrame(data, index=display_benchmarks)

# Plot
plt.figure(figsize=(11, 13))
sns.set(font_scale=1.25)

ax = sns.heatmap(
    df,
    annot=True,
    fmt=".1f",
    cmap="cividis",  # grayscale-friendly
    linewidths=0.75,
    linecolor="gray",
    vmin=0,
    vmax=100,
    annot_kws={"size": 23},  # ← Increase annotation (cell text) font size here
    cbar_kws={"label": None}
)

# Consistent font styling
label_fontsize = 23

# Set title and axis labels
plt.title("Endpoint Pairing Success", fontsize=label_fontsize, pad=20)
ax.set_xlabel("Layout", fontsize=label_fontsize, labelpad=15)
ax.set_ylabel("Benchmark", fontsize=label_fontsize, labelpad=20)

# Colorbar label formatting
cbar = ax.collections[0].colorbar
cbar.set_label("Success Rate (%)", fontsize=20, labelpad=1)
cbar.ax.tick_params(labelsize=23)

# Tick styling
plt.xticks(rotation=0, fontsize=23)
plt.yticks(rotation=0, fontsize=23)
plt.tight_layout()

# Save
plt.savefig("success_heatmap_styled_fig14.pdf", dpi=300)
plt.show()

