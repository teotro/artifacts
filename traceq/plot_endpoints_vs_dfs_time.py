import os
import re
import json
import matplotlib.pyplot as plt

# Base directory where layout folders are located
base_dir = "targets"
layouts = ["compact", "sparse", "intermediate"]
colors = {"compact": "blue", "sparse": "green", "intermediate": "red"}

# Dictionary to store data per layout
layout_points = {layout: {"x": [], "y": []} for layout in layouts}

for layout in layouts:
    layout_path = os.path.join(base_dir, layout)
    if not os.path.isdir(layout_path):
        print(f"Layout folder not found: {layout_path}")
        continue

    for benchmark in os.listdir(layout_path):
        bench_path = os.path.join(layout_path, benchmark)
        if not os.path.isdir(bench_path):
            continue

        for i in range(1, 31):  # pert_1 to pert_30
            txt_path = os.path.join(bench_path, f"pert_{i}.txt")
            json_path = os.path.join(bench_path, f"pert_{i}.json")

            if not os.path.exists(txt_path) or not os.path.exists(json_path):
                continue

            # Extract total logical-unknown cells
            total_unknown_cells = 0
            with open(txt_path, "r") as f_txt:
                for line in f_txt:
                    match = re.search(r"\((\d+)\s+logical-unknown cells\)", line)
                    if match:
                        count = int(match.group(1))
                        total_unknown_cells += count

            # Extract dfs_time_sec
            try:
                with open(json_path, "r") as f_json:
                    stats = json.load(f_json)
                    dfs_time = stats.get("dfs_time_sec", 0.0)
            except (json.JSONDecodeError, FileNotFoundError):
                continue

            layout_points[layout]["x"].append(total_unknown_cells)
            layout_points[layout]["y"].append(dfs_time)

# Generate separate plots per layout
for layout in layouts:
    x_vals = layout_points[layout]["x"]
    y_vals = layout_points[layout]["y"]

    plt.figure(figsize=(12, 9))
    plt.scatter(x_vals, y_vals, color=colors[layout], alpha=0.6)
    plt.xlabel("Total Qubit Endpoints with Ambiguous Paths", fontsize=26)
    plt.ylabel("DFS Time (s)", fontsize=25)
    plt.title(f"DFS Time vs Total Qubit Endpoints with Ambiguous Paths", fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"dfs_vs_unknown_cells_{layout}_fig15.pdf", dpi=300)
    plt.show()
