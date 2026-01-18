if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore', category=UserWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    import os
    import argparse
    import networkx as nx
    import matplotlib.pyplot as plt
    import numpy as np
    from qiskit import QuantumCircuit
    from collections import OrderedDict, deque, defaultdict
    from qiskit.converters import circuit_to_dag
    import tracemaker
    import l3_trace_dag
    import heuristics
    import utils
    import helper, helper_2
    import pickle

    def repr_6bit_matrix(matrix):
        # print matrix where all elements are 6bit
        return "\n".join(" ".join(f"{x:07b}" for x in row) for row in matrix)
    from collections import OrderedDict
    import numpy as np

    def normalize_known_endpoints(recovered_l3_traces):
        """
        Normalize all traces by promoting cells to 0b1110000 if:
        - Their (x, y) coordinate is globally identified as a known qubit endpoint
        - And their current value is exactly 0b1010000 (a wire tile)

        Logs each promotion and prints:
        - Global endpoint set before and after normalization
        - A summary list of all promoted locations
        """
        UNK_FLAG    = 0b1110000
        WIRE_FLAG   = 0b1010000

        # Step 1: Gather global set of known qubit endpoint coordinates (before)
        global_endpoints_before = set()
        for entry in recovered_l3_traces.values():
            mat = entry["trace"]
            for x in range(mat.shape[0]):
                for y in range(mat.shape[1]):
                    val = mat[x, y]
                    msb = (val >> 4) & 0b111
                    if val == UNK_FLAG or msb in {0b010, 0b011}:
                        global_endpoints_before.add((x, y))

        print(f"\nGlobal endpoints BEFORE normalization (total = {len(global_endpoints_before)}):")
        print(global_endpoints_before)

        # Step 2: Apply promotion per trace
        refined = OrderedDict()
        all_promoted_locations = set()

        for trace_id, entry in recovered_l3_traces.items():
            mat = entry["trace"].copy()
            promoted_cells = []

            for (x, y) in global_endpoints_before:
                if mat[x, y] == WIRE_FLAG:
                    mat[x, y] = UNK_FLAG
                    promoted_cells.append((x, y))
                    all_promoted_locations.add((x, y))

            if promoted_cells:
                print(f"[Trace {trace_id}] Promoted {len(promoted_cells)} cell(s) to 0b1110000 at: {promoted_cells}")

            refined[trace_id] = {
                "trace": mat,
                "paths": entry["paths"]
            }

        # Step 3: Recompute global endpoints after promotion
        global_endpoints_after = set()
        for entry in refined.values():
            mat = entry["trace"]
            for x in range(mat.shape[0]):
                for y in range(mat.shape[1]):
                    val = mat[x, y]
                    msb = (val >> 4) & 0b111
                    if val == UNK_FLAG or msb in {0b010, 0b011}:
                        global_endpoints_after.add((x, y))

        print(f"\nGlobal endpoints AFTER normalization (total = {len(global_endpoints_after)}):")
        print(global_endpoints_after)

        # Final: Print summary of all promoted locations
        if all_promoted_locations:
            print(f"\nTotal promoted locations across all traces: {len(all_promoted_locations)}")
            print(sorted(all_promoted_locations))
        else:
            print("\nThere were no locations to promote.")

        return refined


    # -----------------------------------------------
    # Parse command-line arguments
    # -----------------------------------------------
    parser = argparse.ArgumentParser(description="TraceQ trace-based reconstruction")

    parser.add_argument("--path", type=str, required=True, help="Path to QASM benchmark")
    parser.add_argument("--num_qubits", type=int, required=True, help="Number of logical qubits")
    parser.add_argument("--benchmark_pert", type=str, required=True, help="Benchmark Perturbation (e.g., pert_1)")
    parser.add_argument("--benchmark_name", type=str, required=True, help="Benchmark Name(e.g., random_1)")

    layout_group = parser.add_mutually_exclusive_group(required=True)
    layout_group.add_argument("--compact", action='store_true', help="Use compact layout")
    layout_group.add_argument("--sparse", action='store_true', help="Use sparse layout")
    layout_group.add_argument("--intermediate", action='store_true', help="Use intermediate layout")

    args = parser.parse_args()

    # -----------------------------------------------
    # Set up variables from args
    # -----------------------------------------------
    path = args.path
    num_qubits = args.num_qubits
    benchmark_name=args.benchmark_name
    benchmark_pert= args.benchmark_pert

    if args.compact:
        layout_type = "compact"
        initial_layout_program = utils.generate_compact_grid(num_qubits)
    elif args.sparse:
        layout_type = "sparse"
        initial_layout_program = utils.generate_fixed_square_sparse_matrix_v8(num_qubits)
    else:
        layout_type = "intermediate"
        initial_layout_program = utils.generate_intermediate_grid(num_qubits)

    qc = QuantumCircuit.from_qasm_file(path)

    # -----------------------------------------------
    # Generate traces and reconstruct DAG
    # -----------------------------------------------
    l3_traces = l3_trace_dag.save_traces_reconstruct_dag(num_qubits, path, initial_layout_program)

    for key, mat in list(l3_traces.items()):
        l3_traces[key] = {"trace": mat}

    l1_traces = tracemaker.reduce_l3_to_l1(l3_traces)
    recovered_l3_traces = heuristics.l1_traces_to_l3(l1_traces)

    for TRACE_IDX in l3_traces.keys():
        print("Trace ID:", TRACE_IDX, "\n")
        print(repr_6bit_matrix(l3_traces[TRACE_IDX]["trace"]))

    for TRACE_IDX in l3_traces.keys():
        print("Trace ID:", TRACE_IDX, "\n")
        print(repr_6bit_matrix(np.array(recovered_l3_traces[TRACE_IDX]["trace"], dtype=int)))
    
    refined_l3_traces = normalize_known_endpoints(recovered_l3_traces)

    for TRACE_IDX in l3_traces.keys():
        print("Trace ID:", TRACE_IDX, "\n")
        print(repr_6bit_matrix(np.array(refined_l3_traces[TRACE_IDX]["trace"], dtype=int)))

    # program_dag = helper.reconstruct_dag_from_l3_traces(recovered_l3_traces, layout_type, benchmark_name, benchmark_pert)
    # program_dag = helper.reconstruct_dag_from_l3_traces(recovered_l3_traces, layout_type, benchmark_name)
    program_dag = helper_2.reconstruct_dag_from_l3_traces(refined_l3_traces, layout_type, benchmark_name, benchmark_pert)
    # program_dag = helper_2.reconstruct_dag_from_l3_traces(refined_l3_traces, layout_type, benchmark_name)


    # -----------------------------------------------
    # Save result to pickle
    # -----------------------------------------------
    output_dir = os.path.join("targets", layout_type, benchmark_name)
    # output_dir = os.path.join("queries_3", layout_type)
    os.makedirs(output_dir, exist_ok=True)
    output_pkl = os.path.join(output_dir, f"{benchmark_pert}.pkl")
    # output_pkl = os.path.join(output_dir, f"{benchmark_name}.pkl")

    with open(output_pkl, "wb") as f:
        pickle.dump(program_dag, f)

    print(f"✅ DAG for {benchmark_name} saved to {output_pkl}")
