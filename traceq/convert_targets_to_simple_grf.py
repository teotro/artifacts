import os
import pickle
import networkx as nx

# === SIMPLE .grf FORMAT FUNCTION ===
def save_graph(g, fn):
    g = nx.convert_node_labels_to_integers(g)
    with open(fn, "w") as f:
        f.write(f"{len(g)}\n")
        for u in g.nodes:
            f.write(f"{u} 0\n")
        for u in g.nodes:
            succs = list(g.successors(u))
            f.write(f"{len(succs)}\n")
            for v in succs:
                f.write(f"{u} {v}\n")

# === SETTINGS ===
layouts = ["sparse", "compact", "intermediate"]
benchmarks = [f"random_{i}" for i in range(1, 21)]
base_path = "targets"
num_perts = 30

# === MAIN CONVERSION LOOP ===
for layout in layouts:
    for bench in benchmarks:
        for i in range(1, num_perts + 1):
            bench_dir = os.path.join(base_path, layout, bench)
            pkl_file = os.path.join(bench_dir, f"pert_{i}.pkl")
            grf_file = os.path.join(bench_dir, f"pert_{i}.grf")
            err_file = os.path.join(bench_dir, f"pert_{i}.txt.err")

            if os.path.isfile(err_file) or not os.path.isfile(pkl_file):
                continue

            # Remove existing .grf if present
            if os.path.isfile(grf_file):
                try:
                    os.remove(grf_file)
                    print(f"🗑 Removed existing: {grf_file}")
                except Exception as e:
                    print(f"⚠ Could not remove {grf_file}: {e}")
                    continue

            # Load and convert
            try:
                with open(pkl_file, "rb") as f:
                    g_nx = pickle.load(f)

                save_graph(g_nx, grf_file)
                print(f"✔ Saved: {grf_file}")

            except Exception as e:
                print(f"❌ Failed to convert {pkl_file}: {e}")
