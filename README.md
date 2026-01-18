# 🚀 TraceQ: FTQC Access Trace Analysis Framework

**TraceQ** is a Python-based software framework designed for the analysis of space-time activity traces, defined as access traces, from lattice surgery in **Fault-Tolerant Quantum Computing (FTQC)**.

By processing these "access traces" and reconstructing the program's **Directed Acyclic Graph (DAG)**, we can effectively identify common quantum subroutines using advanced subgraph matching techniques. The framework incorporates heuristics to manage the inherent ambiguity found in these traces and has been validated on a variety of synthetic FTQC benchmarks across three distinct surface-code layouts: **Square Sparse**, **Compact**, and **Intermediate**.

## 📖 General Workflow

TraceQ implements a two-stage approach:

1.  **Trace Reconstruction:** The input traces are processed to resolve ambiguities and reconstruct the quantum program's complete Dataflow DAG.
2.  **Subroutine Identification:** The reconstructed DAGs are searched for known quantum subroutines using the **VF3 subgraph matching algorithm**.

For complete details on the methodology, please refer to the paper:

> **TraceQ: Trace-Based Reconstruction of Quantum Circuit Dataflow in Surface-Code Fault-Toelrant Quantum Computing**
> *[https://arxiv.org/pdf/2508.14533](https://arxiv.org/pdf/2508.14533)*

-----

## 🛠️ Reproducibility and Usage

### Hardware Recommendation

To minimize execution time by running processes in parallel, we **strongly recommend** using a machine with a high core count. Our experiments were run on an **Intel(R) Xeon(R) Gold 6438N Processor** with **256GB of RAM**. You will need at least 600GB of disk space.

### Step 1: Setup Environment

Use the Conda environment file to ensure you have all necessary dependencies installed:

```bash
# 1. Create the 'traceq' environment
conda env create -f environment.yml

# 2. Activate the environment
conda activate traceq
```

### Step 2: Generate Access Traces and Reconstruct DAGs

This phase runs the core trace reconstruction algorithm (`tracemaker`) on the synthetic benchmarks, testing robustness against 30 different perturbations across three layouts.

**⚠️ IMPORTANT:** Before running, navigate to the `traceq` directory and **update all necessary paths** in `run_tracemaker.sh` and `rerun_failed.sh`.

```bash
# Navigate to the main directory
cd traceq
```

Run the reconstruction for all three layouts in parallel:

```bash
# Run reconstruction for the 'compact' layout
./run_tracemaker.sh compact --parallel True --num_cores 30

# Run reconstruction for the 'intermediate' layout
./run_tracemaker.sh intermediate --parallel True --num_cores 30

# Run reconstruction for the 'sparse' layout
./run_tracemaker.sh sparse --parallel True --num_cores 30
```

**Note:** Each run has an **1-hour timeout limit**. Successful reconstructions are saved as `.pkl` files in the `targets` folder. You can check the execution status in the corresponding `.txt` (success) or `.err` (failure) files.

-----

### Step 3: Rerunning Failed Cases

We allow up to **10 extra attempts** to resolve failed perturbations (due to timeouts or odd-ambiguity errors):

```bash
./rerun_failed.sh compact --num_cores 30
./rerun_failed.sh intermediate --num_cores 30
./rerun_failed.sh sparse --num_cores 30
```

-----

## 📈 Generating Figures (Paper Results)

Once the DAGs are reconstructed, run the following scripts to generate the figures presented in the paper.

### A. Ambiguity and Performance (Figures 12, 13, 14, 15)

| Figure | Description | Command |
| :--- | :--- | :--- |
| **Figure 12** | Box plot of **ambiguity percentages** across synthetic benchmarks and layouts. | `python plot_ambiguity_by_layout.py` |
| **Figure 13** | Box plot of **DFS execution times** across benchmarks and layouts (Y-axis broken at 5s). | `python dfs_time_boxplot.py` |
| **Figure 14** | **Heatmap of success rates** (out of 30 perturbations) after post-DFS processing. | `python generate_success_heatmap.py` |
| **Figure 15** | Correlation between **DFS time and ambiguous qubit endpoints** by layout. | `python plot_endpoints_vs_dfs_time.py` |

### B. Subroutine Matching (Figures 16, 17)

This step uses the VF3 exact subgraph matching algorithm to identify subroutines.

**i. Preparation and Initial Matching (1 second timeout):**

1.  **Convert DAGs:** Convert the successful `.pkl` DAG targets into the `.grf` format required by `vf3lib`. (Queries are pre-provided in `.grf` format).
    ```bash
    python convert_targets_to_simple_grf.py
    ```
    `cd ../vf3lib && python convert_targets_to_simple_grf.py`

2.  **Run Initial Matching:** Execute the matching process with a short timeout.
    ```bash
    ./matching_all.sh compact 1
    ./matching_all.sh intermediate 1
    ./matching_all.sh sparse 1
    ```

**ii. Rerunning with Longer Timeouts:**

We progressively increase the timeout for the remaining unmatched benchmarks to gather full data:

| New Timeout (seconds) | Previous Timeout (seconds) | Command (Example for Compact) |
| :--- | :--- | :--- |
| **10** | 1 | `./rerun_timeouts.sh compact 10 1` |
| **60** | 10 | `./rerun_timeouts.sh compact 60 10` |
| **600** | 60 | `./rerun_timeouts.sh compact 600 60` |
| **3600** | 600 | `./rerun_timeouts.sh compact 3600 600` |

*Repeat the `rerun_timeouts.sh` commands above for the `intermediate` and `sparse` layouts.*

**iii. Final Figure Generation:**

1.  **Aggregate Results:**
    ```bash
    ./analyze_layouts.sh compact
    ./analyze_layouts.sh intermediate
    ./analyze_layouts.sh sparse
    ```
2.  **Figure 16:** Generates the Stacked Success Rate plot.
    ```bash
    python analyze_success_rate.py
    ```
3.  **Figure 17:** Generates the Cumulative Subroutine Identification plot.
    ```bash
    python plot.py
    ```

-----

## 📚 References

If you use TraceQ in your research, please cite our paper:

**TraceQ: Trace-Based Reconstruction of Quantum Circuit Dataflow in Surface-Code Fault-Toelrant Quantum Computing**
*Theodoros Trochatos, Christopher Kang, Andrew Wang, Frederic T Chong, Jakub Szefer*
**IEEE International Symposium on High-Performance Computer Architecture (HPCA) - 2026**