import sys
import re
import random
import numpy as np
sys.path.append("../")
# import qasm_reader
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso
from pathlib import Path

# Let's put all functions here

def generate_fixed_square_sparse_matrix_v8(n):
    """
    Generate a fixed-layout square sparse matrix based on k^2 thresholds:
    - For n <= 2^2, uses 5x5 grid placing 4 ones.
    - For 2^2 < n <= 3^2, uses 7x7 grid placing 9 ones.
    - For 3^2 < n <= 4^2, uses 9x9 grid placing 16 ones.
    and so on (size = 2*k + 1, ones = k^2).
    """
    # Find smallest k: k^2 >= n, starting at k=2
    k = 2
    while n > k**2:
        k += 1

    max_ones = k**2
    size = 2 * k + 1  # ensures odd dimension with single-zero border

    matrix = np.zeros((size, size), dtype=int)
    placed = 0

    # Fill in fixed pattern at odd indices
    for i in range(1, size, 2):
        for j in range(1, size, 2):
            if placed >= max_ones:
                break
            matrix[i, j] = 1
            placed += 1
        if placed >= max_ones:
            break

    return matrix

def generate_compact_grid(num_qubits):
    """
    Generate a compact 3-row grid where qubits (1s) are placed in top and bottom rows 
    at alternating positions, while the middle row is filled with 0s.
    Example output for num_qubits = 9:
    [1 0 1 0 1 0 1 0 1]
    [0 0 0 0 0 0 0 0 0]
    [1 0 1 0 1 0 1 0 0]
    """
    # Ensure enough columns to alternate 1s and 0s and fit all qubits
    cols = (num_qubits + 1) // 2 * 2 - 1
    grid = np.zeros((3, cols), dtype=int)

    placed = 0
    for i in range(0, cols, 2):
        if placed < num_qubits:
            grid[0, i] = 1
            placed += 1
        if placed < num_qubits:
            grid[2, i] = 1
            placed += 1

    return grid

def generate_intermediate_grid(num_qubits):
    """
    Generate a 2-row intermediate grid where all qubits (1s) are placed in the first row,
    and all ancillas (0s) in the second row.
    Example for num_qubits = 6:
    [1 1 1 1 1 1]
    [0 0 0 0 0 0]
    """
    grid = np.zeros((2, num_qubits), dtype=int)
    grid[0, :] = 1
    return grid

def compile_random_qasm_program(qasm_file_paths, output_path="../../bench_from_chris/random_assembled.qasm"):
    """
    Combine multiple QASM programs with randomly relabeled qubits to avoid collisions.

    Parameters:
    - qasm_file_paths: list of paths to QASM files (strings or Path objects)
    - output_path: where to write the assembled QASM output

    Returns:
    - A string containing the new assembled QASM program
    """
    assembled_qasm = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
    ]

    global_qubit_counter = 0
    global_creg_counter = 0
    used_qubits = []
    used_cregs = []
    body = []

    for idx, file_path in enumerate(qasm_file_paths):
        with open(file_path, "r") as f:
            lines = f.readlines()

        q_lines = []
        c_lines = []
        instr_lines = []
        for line in lines:
            if line.startswith("qreg"):
                qreg_match = re.match(r"qreg\s+q\[(\d+)\];", line.strip())
                if qreg_match:
                    local_qcount = int(qreg_match.group(1))
                    new_qubits = [f"q{global_qubit_counter + i}" for i in range(local_qcount)]
                    used_qubits += new_qubits
                    q_lines.append(f"qreg q[{global_qubit_counter + local_qcount}];")
                    global_qubit_counter += local_qcount
            elif line.startswith("creg"):
                creg_match = re.match(r"creg\s+c\[(\d+)\];", line.strip())
                if creg_match:
                    local_ccount = int(creg_match.group(1))
                    new_cregs = [f"c{global_creg_counter + i}" for i in range(local_ccount)]
                    used_cregs += new_cregs
                    c_lines.append(f"creg c[{global_creg_counter + local_ccount}];")
                    global_creg_counter += local_ccount
            elif line.strip() and not line.startswith("//") and not line.startswith("OPENQASM") and not line.startswith("include"):
                instr_lines.append(line.strip())

        # Random remapping of qubit indices in this subprogram
        mapping = {}
        all_indices = sorted(set(int(idx) for match in re.findall(r"q\[(\d+)\]", "".join(instr_lines)) for idx in [match]))
        shuffled = random.sample(range(global_qubit_counter - len(all_indices), global_qubit_counter), len(all_indices))
        for old, new in zip(all_indices, shuffled):
            mapping[old] = new

        # Apply mapping to instructions
        remapped_instrs = []
        for instr in instr_lines:
            remapped = re.sub(r"q\[(\d+)\]", lambda m: f"q[{mapping[int(m.group(1))]}]", instr)
            remapped_instrs.append(remapped)

        body.extend(remapped_instrs)

    # Write final assembled QASM
    assembled_qasm += [f"qreg q[{global_qubit_counter}];"]
    if global_creg_counter > 0:
        assembled_qasm += [f"creg c[{global_creg_counter}];"]
    assembled_qasm += body

    assembled_text = "\n".join(assembled_qasm)

    # Optionally save to file
    with open(output_path, "w") as f:
        f.write(assembled_text)

    return assembled_text


def extract_2q_gates_from_qasm(path):
    qasm_lines = qasm_reader.read_QASM_file(path)[8:]
    gates_2q = []
    qubits = set()
    for line in qasm_lines:
        out = qasm_reader.read_QASM_line(line)
        if out:
            gate, targs = out
            # if len(targs) > 2:
            #     raise ValueError('Encountered a gate with more than 2 qubits')
            if len(targs) == 2:
                gates_2q.append((gate, targs))
                qubits.update(targs)

    # For debugging purposes
    print("Gates for qasm-derived dag", gates_2q)

    return gates_2q, qubits

def qasm_to_dependency_dag(gates_2q, qubits):
    dag = nx.DiGraph()
    qubit_status = {q: None for q in qubits}
    num_gates = len(gates_2q)
    for idx, (gate, targs) in enumerate(gates_2q[::-1]):
        true_idx = num_gates - idx - 1
        dag.add_node(true_idx, type='gate')
        dag.nodes[true_idx]["gate"] = gate
        dag.nodes[true_idx]["targs"] = targs

        future_dependencies = [qubit_status[q] for q in targs]

        for dep in future_dependencies:
            if dep is not None:
                dag.add_edge(true_idx, dep)
        for q in targs:
            qubit_status[q] = true_idx
    return dag


def subgraph_match(full_dag: nx.DiGraph, subroutine_dag: nx.DiGraph):
    """
    Attempt to match subroutine_dag as a subgraph inside full_dag,
    allowing endpoints (control/target) to match in any order,
    and handling qubit labels like (3,1), np.int64, 'Q3_1', etc.
    Returns (True, mapping) or (False, None).
    """

    def normalize_qubit(q):
        # Accepts: 'Q(3, 1)', 'Q3_1', (3,1), etc.
        if isinstance(q, str):
            if q.startswith('Q('):
                q = q[2:-1]  # drop 'Q(' and ')'
                nums = [int(x) for x in q.split(',')]
                return f"Q{nums[0]}_{nums[1]}"
            if q.startswith('Q') and '_' in q:
                return q
            if q.startswith('('):  # just tuple, not Q
                nums = [int(x) for x in q[1:-1].split(',')]
                return f"Q{nums[0]}_{nums[1]}"
        if isinstance(q, tuple) and len(q) == 2:
            return f"Q{int(q[0])}_{int(q[1])}"
        return str(q)


    def node_match(n1_attrs, n2_attrs):
        # Gate type must match
        if n1_attrs.get('gate', '').lower() != n2_attrs.get('gate', '').lower():
            return False
        q1 = {normalize_qubit(n1_attrs.get('control')), normalize_qubit(n1_attrs.get('target'))}
        q2 = {normalize_qubit(n2_attrs.get('control')), normalize_qubit(n2_attrs.get('target'))}
        return q1 == q2

    def edge_match(e1_attrs, e2_attrs):
        return True

    GM = nx.isomorphism.DiGraphMatcher(
        full_dag,
        subroutine_dag,
        node_match=node_match,
        edge_match=edge_match
    )

    found = GM.subgraph_is_isomorphic()
    if found:
        mapping = next(GM.subgraph_isomorphisms_iter())
        print("✅ Subroutine detected!")
        return True, mapping
    else:
        print("❌ Subroutine NOT detected.")
        return False, None




    

def visualize_dag(G, title="DAG Visualization", node_size= 200, font_size=8):
    pos = nx.kamada_kawai_layout(G)
    # plt.figure(figsize=(10, 7))
    nx.draw(G, pos, with_labels=True, node_size=node_size, font_size=font_size)
    # labels = {n: G.nodes[n].get('gate', str(n)) for n in G.nodes()}
    # nx.draw_networkx_labels(G, pos, labels=labels, font_size=font_size)
    plt.title(title)
    plt.axis('off')
    plt.show()

def detect_subroutine_from_qasms(program_qasm_file: str,
                                subroutine_qasm_file: str):
    """
    Load two QASM files, build their dependency DAGs, and check
    whether subroutine_dag is an exact subgraph of program_dag.

    Args:
        program_qasm_file (str): Path to the main program QASM.
        subroutine_qasm_file (str): Path to the subroutine QASM.

    Returns:
        (bool, dict or None):
          - True and the mapping if the subroutine is found.
          - False and None otherwise.
    """
    # Build the program DAG
    gates_prog, qubits_prog = extract_2q_gates_from_qasm(program_qasm_file)
    program_dag = qasm_to_dependency_dag(gates_prog, qubits_prog)

    # Build the subroutine DAG
    gates_sub, qubits_sub = extract_2q_gates_from_qasm(subroutine_qasm_file)
    subroutine_dag = qasm_to_dependency_dag(gates_sub, qubits_sub)
    visualize_dag(subroutine_dag, title="Subroutine DAG")

    # Exact subgraph matching
    found, mapping = subgraph_match(program_dag, subroutine_dag)

    if found:
        print("✅ Subroutine detected in program!")
    else:
        print("❌ Subroutine NOT detected.")

    return program_dag, subroutine_dag


def detect_subroutine_from_traces(program_dag, subroutine_dag):
    """
    Detects a subroutine in a program DAG via exact subgraph matching,
    by delegating to the shared `subgraph_match` function.

    Args:
        program_dag (nx.DiGraph or tuple): The full program DAG (or (dag, …)).
        subroutine_dag (nx.DiGraph or tuple): The subroutine DAG (or (dag, …)).

    Returns:
        found (bool): Whether the subroutine was detected.
        mapping (dict or None): The first mapping if found, else None.
    """
    # Unpack if a tuple was passed in
    if isinstance(program_dag, tuple):
        program_dag = program_dag[0]
    if isinstance(subroutine_dag, tuple):
        subroutine_dag = subroutine_dag[0]

    # Type checks
    if not isinstance(program_dag, nx.DiGraph):
        raise TypeError(f"Expected nx.DiGraph for program_dag, got {type(program_dag)}")
    if not isinstance(subroutine_dag, nx.DiGraph):
        raise TypeError(f"Expected nx.DiGraph for subroutine_dag, got {type(subroutine_dag)}")
    
    found, mapping = subgraph_match(program_dag, subroutine_dag)

    if found:
        print("✅ Subroutine detected in program!")
    else:
        print("❌ Subroutine NOT detected.")

    return found, mapping
