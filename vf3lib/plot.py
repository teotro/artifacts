import os
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

# === CONFIGURATION ===
timeouts = [1, 10, 60, 600, 3600]
layouts = ["compact", "sparse", "intermediate"]
base_path = "."  # Base path where results_*s directories are

# === PARSE LOGS AND COUNT MATCHES ===
match_counts = {layout: {t: 0 for t in timeouts} for layout in layouts}

for layout in layouts:
    seen = set()  # avoid duplicate (bench, pert, query) triplets
    for timeout in timeouts:
        log_dir = os.path.join(base_path, f"results_{timeout}s")
        if not os.path.isdir(log_dir):
            continue
        for fname in os.listdir(log_dir):
            if not fname.endswith(f"{layout}.log"):
                continue
            fpath = os.path.join(log_dir, fname)
            with open(fpath) as f:
                bench = fname.split(f"_{layout}")[0]
                pert, query = None, None
                for line in f:
                    line = line.strip()
                    if line.startswith("Matching:"):
                        parts = line.split()
                        query = parts[1]
                        pert = parts[3]
                    elif "✅ Match found" in line and query and pert:
                        key = (bench, pert, query)
                        if key not in seen:
                            match_counts[layout][timeout] += 1
                            seen.add(key)
                        query, pert = None, None

# === CONVERT TO DATAFRAME ===
df = pd.DataFrame(match_counts).T
df = df[timeouts]  # ensure timeout order
df_cumulative = df.cumsum(axis=1)

# === PLOTTING ===
x_vals = list(range(len(timeouts)))
x_labels = [f"{t}s" for t in timeouts]

fig, ax = plt.subplots(figsize=(12, 8))
colors = {"compact": "#1f77b4", "sparse": "#2ca02c", "intermediate": "#d62728"}

for layout in layouts:
    y_vals = df_cumulative.loc[layout].values
    ax.plot(x_vals, y_vals, marker='o', linewidth=2,
            color=colors[layout], label=layout.capitalize())

ax.set_title("Cumulative Subroutines Identified", fontsize=28)
ax.set_xlabel("Timeout", fontsize=26)
ax.set_ylabel("Total Identified Subroutines", fontsize=26)
ax.set_xticks(x_vals)
ax.set_xticklabels(x_labels, fontsize=26)
ax.tick_params(axis='y', labelsize=26)
ax.legend(title="Layout", fontsize=20, title_fontsize=20)
ax.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# === SAVE FIGURES ===
# plt.savefig("cumulative_subroutines_identified.png", dpi=300)
plt.savefig("cumulative_subroutines_identified_fig17.pdf")
plt.show()
