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

# def prune_dag(dag, prune_fraction):
#     """
#     Prunes nodes from the DAG based on low degree only.
#     Nodes with the lowest degree are prioritized for pruning.

#     Args:
#         dag (nx.DiGraph): The DAG to prune.
#         prune_fraction (float): Fraction of nodes to prune.

#     Returns:
#         pruned_dag (nx.DiGraph): The pruned DAG.
#     """
#     # Debugging Output
#     # print(f"Type of 'dag': {type(dag)}")
#     # print(f"Content of 'dag': {dag}")

#     # Ensure the input is a DiGraph and not a tuple
#     if isinstance(dag, tuple):
#         dag = dag[0]

#     if not isinstance(dag, nx.DiGraph):
#         raise TypeError(f"Expected nx.DiGraph, got {type(dag)} instead.")

#     # Calculate node degrees
#     degrees = dict(dag.degree())

#     # Determine the number of nodes to prune
#     num_nodes_to_prune = int(len(dag.nodes) * prune_fraction)
#     if num_nodes_to_prune == 0:
#         # print("No nodes to prune.")
#         return dag.copy()

#     # Rank nodes based on degree in ascending order
#     sorted_nodes = sorted(dag.nodes(), key=lambda node: degrees[node])

#     # Select the nodes to prune
#     nodes_to_prune = sorted_nodes[:num_nodes_to_prune]

#     # Create a pruned copy of the DAG
#     pruned_dag = dag.copy()
#     pruned_dag.remove_nodes_from(nodes_to_prune)

#     print(f"Pruned {len(nodes_to_prune)} nodes out of {len(dag.nodes)} ({len(nodes_to_prune) / len(dag.nodes):.2%})")

#     return pruned_dag


# def approximate_subgraph_match(full_dag, subroutine_dag, sensitivity, prune_fraction):
#     """
#     Attempt to approximately match subroutine_dag as a subgraph inside full_dag
#     after controlled pruning based on node degree and centrality.

#     Args:
#         full_dag (nx.DiGraph): The full program DAG.
#         subroutine_dag (nx.DiGraph): The subroutine DAG.
#         sensitivity (float): Fraction of nodes that must match.
#         prune_fraction (float): Fraction of nodes to prune before matching.

#     Returns:
#         (bool, dict): Whether a match was found and the best match mapping.
#     """
#     # Prune the subroutine DAG before matching
#     pruned_subroutine = prune_dag(subroutine_dag, prune_fraction=prune_fraction)

#     # Node matcher: match nodes by gate type
#     def node_match(n1_attrs, n2_attrs):
#         return n1_attrs.get('gate', '').lower() == n2_attrs.get('gate', '').lower()

#     # Edge matcher: trivial (edges have no attributes)
#     def edge_match(e1_attrs, e2_attrs):
#         return True

#     GM = nx.isomorphism.DiGraphMatcher(
#         full_dag,
#         pruned_subroutine,
#         node_match=node_match,
#         edge_match=edge_match
#     )

#     matches = list(GM.subgraph_isomorphisms_iter())

#     if not matches:
#         return False, None

#     # Pick the match that covers the most nodes
#     best_match = max(matches, key=lambda m: len(m))
#     match_fraction = len(best_match) / pruned_subroutine.number_of_nodes()

#     print(f"Best match covers {match_fraction:.2f} of the pruned subroutine DAG.")

#     return match_fraction >= sensitivity, best_match

# def subgraph_match(full_dag: nx.DiGraph, subroutine_dag: nx.DiGraph):
#     """
#     Attempt to match subroutine_dag exactly as a subgraph inside full_dag.

#     Args:
#         full_dag (nx.DiGraph): The full program DAG.
#         subroutine_dag (nx.DiGraph): The subroutine DAG to find.

#     Returns:
#         (bool, dict): 
#           - True and a mapping dict if subroutine_dag is found inside full_dag.
#           - False and None otherwise.
#     """
#     # Node matcher: match nodes by gate type (case-insensitive)
#     def node_match(n1_attrs, n2_attrs):
#         return n1_attrs.get('gate', '').lower() == n2_attrs.get('gate', '').lower()

#     # Edge matcher: no attributes, so always True
#     def edge_match(e1_attrs, e2_attrs):
#         return True

#     GM = nx.isomorphism.DiGraphMatcher(
#         full_dag,
#         subroutine_dag,
#         node_match=node_match,
#         edge_match=edge_match
#     )

#     # Check if subroutine_dag is isomorphic to some subgraph of full_dag
#     if GM.subgraph_is_isomorphic():
#         # Grab the first mapping found
#         mapping = next(GM.subgraph_isomorphisms_iter())
#         return True, mapping
#     else:
#         return False, None


# def subgraph_match(full_dag: nx.DiGraph, subroutine_dag: nx.DiGraph):
#     """
#     Attempt to match subroutine_dag as a subgraph inside full_dag,
#     allowing endpoints (control/target) to match in any order.
#     Returns (True, mapping) or (False, None).
#     """

#     # def qubit_label(q):
#     #     # Standardize qubit label: tuple or string
#     #     if isinstance(q, tuple):
#     #         return f"Q{q[0]}_{q[1]}"
#     #     return str(q)
#     # def qubit_label(q):
#     #     # Handles both (row,col) tuples and already-formatted strings
#     #     if isinstance(q, tuple) and len(q) == 2:
#     #         return f"Q{q[0]}_{q[1]}"
#     #     if isinstance(q, str):
#     #         return q
#     #     return str(q)

#     # def qubit_label(q):
#     #     # Handles (3, 1), 'Q3_1', etc.
#     #     if isinstance(q, tuple):
#     #         return f"Q{q[0]}_{q[1]}"
#     #     if isinstance(q, str):
#     #         if q.startswith('Q') and '_' in q:
#     #             return q
#     #         # Handles cases like "Q(3, 1)" or "(3, 1)"
#     #         if q.startswith('Q(') and ',' in q:
#     #             parts = q[2:-1].split(',')
#     #             return f"Q{int(parts[0].strip())}_{int(parts[1].strip())}"
#     #         if q.startswith('(') and ',' in q:
#     #             parts = q[1:-1].split(',')
#     #             return f"Q{int(parts[0].strip())}_{int(parts[1].strip())}"
#     #     # Last resort: try to eval and format
#     #     try:
#     #         t = eval(q)
#     #         if isinstance(t, tuple) and len(t) == 2:
#     #             return f"Q{t[0]}_{t[1]}"
#     #     except:
#     #         pass
#     #     return str(q)

#     # def qubit_label(q):
#     #     """
#     #     Standardize qubit labels: supports (row, col) tuples with int, np.int64, etc.,
#     #     and string labels in 'Qrow_col' or with various parens formats.
#     #     """
#     #     import numpy as np
#     #     # Tuple: handle ints, np.int64, etc.
#     #     if isinstance(q, tuple) and len(q) == 2:
#     #         row = int(q[0])  # convert np.int64, etc. to plain int
#     #         col = int(q[1])
#     #         return f"Q{row}_{col}"
#     #     if isinstance(q, str):
#     #         if q.startswith('Q') and '_' in q:
#     #             return q
#     #         # e.g. "Q(3, 1)" or "(3, 1)"
#     #         if ('(' in q and ',' in q and ')' in q):
#     #             nums = q.replace('Q','').replace('(','').replace(')','').split(',')
#     #             if len(nums) == 2:
#     #                 try:
#     #                     row, col = int(nums[0]), int(nums[1])
#     #                     return f"Q{row}_{col}"
#     #                 except Exception:
#     #                     pass
#     #     # If still not matched, maybe it's a numpy tuple
#     #     try:
#     #         # np.int64 or weird tuple repr
#     #         if hasattr(q, '__iter__') and len(q) == 2:
#     #             row, col = int(q[0]), int(q[1])
#     #             return f"Q{row}_{col}"
#     #     except Exception:
#     #         pass
#     #     return str(q)

#     def qubit_label(q):
#         """
#         Standardize qubit label: supports (row, col) tuples with int, np.int64, etc.,
#         or strings like 'Qrow_col'. Always returns 'Qrow_col'.
#         """
#         import numpy as np

#         # Handle tuple of any numeric types (including numpy ints)
#         if isinstance(q, tuple) and len(q) == 2:
#             row = int(q[0])
#             col = int(q[1])
#             return f"Q{row}_{col}"
#         # Handle string of form 'Qrow_col'
#         if isinstance(q, str):
#             # 'Q(np.int64(3), np.int64(5))'
#             if q.startswith('Q('):
#                 # Extract inside the parentheses
#                 inside = q[2:-1]  # Remove 'Q(' and ')'
#                 items = inside.split(',')
#                 if len(items) == 2:
#                     # Remove type wrappers and spaces
#                     row = int(items[0].replace('np.int64', '').replace('(', '').replace(')', '').strip())
#                     col = int(items[1].replace('np.int64', '').replace('(', '').replace(')', '').strip())
#                     return f"Q{row}_{col}"
#             # Already in 'Qrow_col' format
#             if q.startswith('Q') and '_' in q:
#                 return q
#             # Try to parse from '(row, col)' string
#             if q.startswith('(') and q.endswith(')') and ',' in q:
#                 items = q[1:-1].split(',')
#                 if len(items) == 2:
#                     row = int(items[0].replace('np.int64', '').replace('(', '').replace(')', '').strip())
#                     col = int(items[1].replace('np.int64', '').replace('(', '').replace(')', '').strip())
#                     return f"Q{row}_{col}"
#         # Fallback: try to unpack as iterable
#         try:
#             if hasattr(q, '__iter__') and len(q) == 2:
#                 row = int(q[0])
#                 col = int(q[1])
#                 return f"Q{row}_{col}"
#         except Exception:
#             pass
#         # Otherwise just return string
#         return str(q)





#     # def node_match(n1_attrs, n2_attrs):
#     #     # Gate type must match (e.g., 'cx', 'cz', etc.)
#     #     if n1_attrs.get('gate', '').lower() != n2_attrs.get('gate', '').lower():
#     #         return False
#     #     # Unordered set of endpoints
#     #     qubits1 = {qubit_label(n1_attrs.get('control')), qubit_label(n1_attrs.get('target'))}
#     #     qubits2 = {qubit_label(n2_attrs.get('control')), qubit_label(n2_attrs.get('target'))}
#     #     # print(f"Comparing {qubits1} to {qubits2}")
#     #     return qubits1 == qubits2
#     # def node_match(n1_attrs, n2_attrs):
#     # # Gate type must match (e.g., 'cx', 'cz', etc.)
#     #     if n1_attrs.get('gate', '').lower() != n2_attrs.get('gate', '').lower():
#     #         return False
#     #     # Unordered set of endpoints: allow (q1, q2) == (q2, q1)
#     #     qubits1 = {qubit_label(n1_attrs.get('control')), qubit_label(n1_attrs.get('target'))}
#     #     qubits2 = {qubit_label(n2_attrs.get('control')), qubit_label(n2_attrs.get('target'))}
#     #     return qubits1 == qubits2
#     def node_match(n1_attrs, n2_attrs):
#         if n1_attrs.get('gate', '').lower() != n2_attrs.get('gate', '').lower():
#             return False
#         q1a = qubit_label(n1_attrs.get('control'))
#         q1b = qubit_label(n1_attrs.get('target'))
#         q2a = qubit_label(n2_attrs.get('control'))
#         q2b = qubit_label(n2_attrs.get('target'))
#         print(f"Comparing: {{'{q1a}', '{q1b}'}} <-> {{'{q2a}', '{q2b}'}}")
#         return {q1a, q1b} == {q2a, q2b}



#     def edge_match(e1_attrs, e2_attrs):
#         return True  # No edge attributes to match

#     GM = nx.isomorphism.DiGraphMatcher(
#         full_dag,
#         subroutine_dag,
#         node_match=node_match,
#         edge_match=edge_match
#     )

#     if GM.subgraph_is_isomorphic():
#         mapping = next(GM.subgraph_isomorphisms_iter())
#         return True, mapping
#     else:
#         return False, None

# def subgraph_match(full_dag: nx.DiGraph, subroutine_dag: nx.DiGraph):
#     """
#     Attempt to match subroutine_dag as a subgraph inside full_dag,
#     allowing endpoints (control/target) to match in any order,
#     and handling qubit labels like (3,1), np.int64, 'Q3_1', etc.
#     Returns (True, mapping) or (False, None).
#     """
#     def normalize_qubit(q):
#         import re
#         # Already standardized label
#         if isinstance(q, str):
#             # Q3_1 style
#             if q.startswith('Q') and '_' in q[1:]:
#                 return q
#             # Handles Q(np.int64(3), np.int64(1)), Q(3,1), (3,1), (np.int64(3), np.int64(1))
#             match = re.match(r"Q?\(?\s*([a-zA-Z0-9\._]+)\s*,\s*([a-zA-Z0-9\._]+)\s*\)?", q)
#             if match:
#                 def extract_num(val):
#                     val = val.replace('np.int64', '').replace('(', '').replace(')', '').replace(' ', '')
#                     return int(val)
#                 try:
#                     n1 = extract_num(match.group(1))
#                     n2 = extract_num(match.group(2))
#                     return f"Q{n1}_{n2}"
#                 except Exception:
#                     return q
#             # Try to convert to int (for cases like '3')
#             try:
#                 n = int(q)
#                 return f"Q{n}"
#             except Exception:
#                 return q
#         elif hasattr(q, '__iter__') and not isinstance(q, str) and len(q) == 2:
#             # tuple, list (with or without numpy types)
#             try:
#                 return f"Q{int(q[0])}_{int(q[1])}"
#             except Exception:
#                 return str(q)
#         elif q is not None:
#             try:
#                 return f"Q{int(q)}"
#             except Exception:
#                 return str(q)
#         return None


#     def node_match(n1_attrs, n2_attrs):
#         # Collect and normalize, skipping None
#         qubits1 = set(filter(lambda x: x is not None, [
#             n1_attrs.get('control'), n1_attrs.get('target')
#         ]))
#         qubits2 = set(filter(lambda x: x is not None, [
#             n2_attrs.get('control'), n2_attrs.get('target')
#         ]))
#         norm1 = set(normalize_qubit(q) for q in qubits1)
#         norm2 = set(normalize_qubit(q) for q in qubits2)
#         print(f"Comparing: {norm1} <-> {norm2}")  # For debugging
#         return norm1 == norm2

#     def edge_match(e1_attrs, e2_attrs):
#         return True  # No attributes

#     GM = nx.isomorphism.DiGraphMatcher(
#         full_dag,
#         subroutine_dag,
#         node_match=node_match,
#         edge_match=edge_match
#     )

#     if GM.subgraph_is_isomorphic():
#         mapping = next(GM.subgraph_isomorphisms_iter())
#         return True, mapping
#     else:
#         return False, None

def subgraph_match(full_dag: nx.DiGraph, subroutine_dag: nx.DiGraph):
    """
    Attempt to match subroutine_dag as a subgraph inside full_dag,
    allowing endpoints (control/target) to match in any order,
    and handling qubit labels like (3,1), np.int64, 'Q3_1', etc.
    Returns (True, mapping) or (False, None).
    """
    # def normalize_qubit(q):
    #     if isinstance(q, str):
    #         if q.startswith('Q') and '_' in q[1:]:
    #             return q
    #         if q.startswith('(') and ',' in q:
    #             try:
    #                 nums = [int(x) for x in q.strip('()').replace('np.int64', '').replace(' ', '').split(',')]
    #                 return f"Q{nums[0]}_{nums[1]}"
    #             except Exception:
    #                 return q
    #         return q
    #     elif hasattr(q, '__iter__') and len(q) == 2:
    #         try:
    #             return f"Q{int(q[0])}_{int(q[1])}"
    #         except Exception:
    #             return str(q)
    #     elif q is not None:
    #         try:
    #             return f"Q{int(q)}"
    #         except Exception:
    #             return str(q)
    #     return None
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

    # def node_match(n1_attrs, n2_attrs):
    #     # Only compare unordered endpoints, ignore gate type if you want
    #     qubits1 = {normalize_qubit(n1_attrs.get('control')), normalize_qubit(n1_attrs.get('target'))}
    #     qubits2 = {normalize_qubit(n2_attrs.get('control')), normalize_qubit(n2_attrs.get('target'))}
    #     print(f"Comparing: {qubits1} <-> {qubits2}")
    #     return qubits1 == qubits2
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

# def detect_subroutine_from_qasms(program_qasm_file, subroutine_qasm_file, sensitivity, prune_fraction=0.0):
#     # Program DAG
#     gates_2q_prog, qubits_prog = extract_2q_gates_from_qasm(program_qasm_file)
#     program_dag = qasm_to_dependency_dag(gates_2q_prog, qubits_prog)
#     # visualize_dag(program_dag, title="Program DAG")

#     # Subroutine DAG ### For testing purposes
#     gates_2q_sub, qubits_sub = extract_2q_gates_from_qasm(subroutine_qasm_file)
#     subroutine_dag = qasm_to_dependency_dag(gates_2q_sub, qubits_sub)
#     # visualize_dag(subroutine_dag, title="Subroutine DAG")
#     # print_dag_structure(subroutine_dag)

#     # Matching
#     found, match = approximate_subgraph_match(program_dag, subroutine_dag, sensitivity, prune_fraction=0.0)

#     if found:
#         print("✅ Subroutine detected in program!")
#     # else:
#     #     print("❌ Subroutine NOT detected.")

#     return subroutine_dag
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


# def detect_subroutine_from_dags(program_dag, subroutine_dag, sensitivity, prune_fraction):
#     """
#     Detects a subroutine in a program DAG using approximate subgraph matching.

#     Args:
#         program_dag (nx.DiGraph): The full program DAG.
#         subroutine_dag (nx.DiGraph): The subroutine DAG.
#         sensitivity (float): Fraction of nodes that must match.
#         prune_fraction (float): Fraction of nodes to prune before matching.

#     Returns:
#         found (bool): Whether the subroutine was detected.
#         match (dict): The best match mapping.
#     """

#     # Ensure inputs are nx.DiGraph
#     if isinstance(program_dag, tuple):
#         program_dag = program_dag[0]
#     if isinstance(subroutine_dag, tuple):
#         subroutine_dag = subroutine_dag[0]

#     if not isinstance(program_dag, nx.DiGraph):
#         raise TypeError(f"Expected nx.DiGraph for program_dag, got {type(program_dag)}")
#     if not isinstance(subroutine_dag, nx.DiGraph):
#         raise TypeError(f"Expected nx.DiGraph for subroutine_dag, got {type(subroutine_dag)}")

#     print(f"Detecting subroutine with prune_fraction={prune_fraction}")

#     # Prune the subroutine DAG before matching
#     pruned_subroutine = prune_dag(subroutine_dag, prune_fraction=prune_fraction)

#     # Node matcher: match nodes by gate type
#     def node_match(n1_attrs, n2_attrs):
#         return n1_attrs.get('gate', '').lower() == n2_attrs.get('gate', '').lower()

#     # Edge matcher: trivial (edges have no attributes)
#     def edge_match(e1_attrs, e2_attrs):
#         return True

#     # Apply subgraph matching
#     GM = nx.isomorphism.DiGraphMatcher(
#         program_dag,
#         pruned_subroutine,
#         node_match=node_match,
#         edge_match=edge_match
#     )

#     matches = list(GM.subgraph_isomorphisms_iter())

#     if not matches:
#         print("❌ Subroutine NOT detected.")
#         return False, None

#     # Pick the match that covers the most nodes
#     best_match = max(matches, key=lambda m: len(m))
#     match_fraction = len(best_match) / pruned_subroutine.number_of_nodes()

#     print(f"Best match covers {match_fraction:.2f} of the pruned subroutine DAG.")

#     return match_fraction >= sensitivity, best_match

# def detect_subroutine_from_traces(program_dag, subroutine_dag):
#     """
#     Detects a subroutine in a program DAG via exact subgraph matching,
#     by delegating to the shared `subgraph_match` function.

#     Args:
#         program_dag (nx.DiGraph or tuple): The full program DAG (or (dag, …)).
#         subroutine_dag (nx.DiGraph or tuple): The subroutine DAG (or (dag, …)).

#     Returns:
#         found (bool): Whether the subroutine was detected.
#         mapping (dict or None): The first mapping if found, else None.
#     """
#     # Unpack if a tuple was passed in
#     if isinstance(program_dag, tuple):
#         program_dag = program_dag[0]
#     if isinstance(subroutine_dag, tuple):
#         subroutine_dag = subroutine_dag[0]

#     # Type checks
#     if not isinstance(program_dag, nx.DiGraph):
#         raise TypeError(f"Expected nx.DiGraph for program_dag, got {type(program_dag)}")
#     if not isinstance(subroutine_dag, nx.DiGraph):
#         raise TypeError(f"Expected nx.DiGraph for subroutine_dag, got {type(subroutine_dag)}")

#     # Delegate to the exact subgraph matcher
#     found, mapping = subgraph_match(program_dag, subroutine_dag)

#     if found:
#         print("✅ Subroutine detected in program!")
#     else:
#         print("❌ Subroutine NOT detected.")

#     return found, mapping

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
    
    # print("Program DAG nodes:")
    # for n, d in program_dag.nodes(data=True):
    #     print(n, d)
    # print("\nSubroutine DAG nodes:")
    # for n, d in subroutine_dag.nodes(data=True):
    #     print(n, d)

    # Delegate to the exact subgraph matcher
    found, mapping = subgraph_match(program_dag, subroutine_dag)

    if found:
        print("✅ Subroutine detected in program!")
    else:
        print("❌ Subroutine NOT detected.")

    return found, mapping
