import os
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
from pathlib import Path
import numpy as np

# === Load data from layout_results/*.txt ===
TXT_DIR = "layout_results"
records = []

for fpath in Path(TXT_DIR).glob("*.txt"):
    layout = fpath.stem
    with open(fpath) as f:
        for line in f:
            bench, sub, succ, fail, timeout, total = line.strip().split()
            succ, fail, timeout, total = map(int, [succ, fail, timeout, total])
            records.append({
                "Layout": layout,
                "Subroutine": sub,
                "Success": succ,
                "Fail": fail,
                "Timeout": timeout,
                "Total": total
            })

# Create DataFrame
df = pd.DataFrame(records)
df["Success %"] = 100 * df["Success"] / df["Total"]
df["Fail %"] = 100 * df["Fail"] / df["Total"]
df["Timeout %"] = 100 * df["Timeout"] / df["Total"]

# === Plot settings ===
layouts = ["compact", "sparse", "intermediate"]
metrics = ["Success %", "Fail %", "Timeout %"]
metric_colors = {
    "Success %": "#4CAF50",     # green
    "Fail %": "#F44336",        # red
    "Timeout %": "#FF9800"      # orange
}
metric_hatches = {
    "Success %": "",
    "Fail %": "///",
    "Timeout %": "..."
}

subroutines = sorted(df["Subroutine"].unique())
x = np.arange(len(subroutines))
bar_width = 0.6

# === Generate one figure per layout ===
for layout in layouts:
    layout_df = (
        df[df["Layout"] == layout]
        .groupby("Subroutine")[["Success %", "Fail %", "Timeout %"]]
        .mean()
        .reindex(subroutines)
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=(21, 11))
    bottoms = np.zeros(len(subroutines))

    for metric in metrics:
        values = layout_df[metric].values
        ax.bar(
            x,
            values,
            bar_width,
            bottom=bottoms,
            label=metric,
            color=metric_colors[metric],
            hatch=metric_hatches[metric],
            edgecolor='black',
            linewidth=0.6
        )
        bottoms += values

    # Formatting
    ax.set_xticks(x)
    ax.set_xticklabels(subroutines, rotation=45, ha="right", fontsize=36)
    ax.set_ylabel("Recovery Success Rate", fontsize=42)
    ax.set_xlabel("Subroutine Instances", fontsize=42)
    ax.set_title(f"Stacked Match Result Breakdown - {layout.capitalize()} Layout", fontsize=42)
    ax.tick_params(axis='y', labelsize=36)
    ax.set_ylim(0, 100)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0, 0),
        fontsize=36,
        frameon=True
    )
    plt.tight_layout()

    # Save per layout
    outname = f"layout_{layout}_fig16.pdf"
    plt.savefig(outname, dpi=300)
    plt.close()

# import os
# import pandas as pd
# from pathlib import Path

# # === Load data from layout_results/*.txt ===
# TXT_DIR = "layout_results"
# records = []

# for fpath in Path(TXT_DIR).glob("*.txt"):
#     layout = fpath.stem
#     with open(fpath) as f:
#         for line in f:
#             bench, sub, succ, fail, timeout, total = line.strip().split()
#             succ, fail, timeout, total = map(int, [succ, fail, timeout, total])
#             records.append({
#                 "Layout": layout,
#                 "Subroutine": sub,
#                 "Success": succ,
#                 "Fail": fail,
#                 "Timeout": timeout,
#                 "Total": total  # kept for completeness, not used
#             })

# # Create DataFrame
# df = pd.DataFrame(records)

# # Normalize always by 30
# df["Success %"] = 100 * df["Success"] / 30
# df["Fail %"] = 100 * df["Fail"] / 30
# df["Timeout %"] = 100 * df["Timeout"] / 30

# # Print tables per layout + average success rate
# layouts = ["compact", "sparse", "intermediate"]
# metrics = ["Success %", "Fail %", "Timeout %"]

# for layout in layouts:
#     layout_df = (
#         df[df["Layout"] == layout]
#         .groupby("Subroutine")[metrics]
#         .mean()
#         .fillna(0)
#         .sort_index()
#     )

#     avg_success = layout_df["Success %"].mean()

#     print(f"\n=== {layout.upper()} Layout ===")
#     print(f"{'Subroutine':<20} {'Success %':>10} {'Fail %':>10} {'Timeout %':>12}")
#     print("-" * 56)
#     for sub, row in layout_df.iterrows():
#         print(f"{sub:<20} {row['Success %']:>10.1f} {row['Fail %']:>10.1f} {row['Timeout %']:>12.1f}")
    
#     print(f"\nAverage Success Rate for {layout}: {avg_success:.2f}%\n")


# import os
# import pandas as pd
# from pathlib import Path

# # === Load data from layout_results/*.txt ===
# TXT_DIR = "layout_results"
# records = []

# for fpath in Path(TXT_DIR).glob("*.txt"):
#     layout = fpath.stem
#     with open(fpath) as f:
#         for line in f:
#             bench, sub, succ, fail, timeout, total = line.strip().split()
#             succ = int(succ)
#             records.append({
#                 "Layout": layout,
#                 "Subroutine": sub,
#                 "Success": succ
#             })

# # Create DataFrame
# df = pd.DataFrame(records)

# # Normalize success by 30; fail is 100 - success; timeout is always 0
# df["Success %"] = 100 * df["Success"] / 30
# df["Fail %"] = 100 - df["Success %"]
# df["Timeout %"] = 0.0

# # Print tables per layout + average success rate
# layouts = ["compact", "sparse", "intermediate"]
# metrics = ["Success %", "Fail %", "Timeout %"]

# for layout in layouts:
#     layout_df = (
#         df[df["Layout"] == layout]
#         .groupby("Subroutine")[metrics]
#         .mean()
#         .fillna(0)
#         .sort_index()
#     )

#     avg_success = layout_df["Success %"].mean()

#     print(f"\n=== {layout.upper()} Layout ===")
#     print(f"{'Subroutine':<20} {'Success %':>10} {'Fail %':>10} {'Timeout %':>12}")
#     print("-" * 56)
#     for sub, row in layout_df.iterrows():
#         print(f"{sub:<20} {row['Success %']:>10.1f} {row['Fail %']:>10.1f} {row['Timeout %']:>12.1f}")
    
#     print(f"\nAverage Success Rate for {layout}: {avg_success:.2f}%\n")
