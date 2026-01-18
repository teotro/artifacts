import warnings
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning) 
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit
from heapq import heappop, heappush
from qiskit.converters import circuit_to_dag
from collections import OrderedDict
from grakel import graph_from_networkx
from grakel.kernels import RandomWalk
import os
import random

# from utils import generate_fixed_square_sparse_matrix_v8

import json

def extract_qubit_mapping(initial_matrix):
    """
    Extracts the qubit mapping from the initial matrix.
    Logical qubits are marked by '1' in the matrix.

    Parameters:
        initial_matrix (list of lists): 2D matrix with logical qubit locations.

    Returns:
        dict: Mapping of logical qubits to (row, col) coordinates.
    """
    qubit_mapping = {}
    logical_qubit = 0

    for row in range(len(initial_matrix)):
        for col in range(len(initial_matrix[0])):
            if initial_matrix[row][col] == 1:
                qubit_mapping[logical_qubit] = (row, col)
                logical_qubit += 1

    return qubit_mapping

# ✅ Function to save matrices as JSON
def save_traces_as_json(traces, filename='l3_traces.json'):
    # Convert matrices to 6-bit binary strings and then to lists for JSON compatibility
    json_traces = {
        key: [[f"{cell:06b}" for cell in row] for row in matrix]
        for key, matrix in traces.items()
    }
    with open(filename, 'w') as f:
        json.dump(json_traces, f, indent=4)
    print(f"✅ Traces saved to {filename} in JSON format")

# ✅ Function to load matrices from JSON
def load_traces_from_json(filename='l3_traces.json'):
    with open(filename, 'r') as f:
        traces = json.load(f)
    print(f"✅ Traces loaded from {filename}")
    return traces

def number_of_qubits(qc):
    return qc.num_qubits

def generate_random_matrix(n, qc):
    matrix = np.zeros((n, n), dtype=int)
    ones_positions = np.random.choice(n * n, number_of_qubits(qc), replace=False)
    for pos in ones_positions:
        row, col = divmod(pos, n)
        matrix[row, col] = 1
    return matrix

def get_two_qubit_gates(qc):
    two_qubit_gates = [
        (qc.find_bit(inst[1][0])[0], qc.find_bit(inst[1][1])[0])
        for inst in qc.data if len(inst[1]) == 2
    ]
    return two_qubit_gates

def create_graph(A):
    G = nx.Graph()
    rows, cols = len(A), len(A[0])
    for r in range(rows):
        for c in range(cols):
            if A[r][c] in {0, 1}:
                if r > 0: G.add_edge((r, c), (r - 1, c))
                if r < rows - 1: G.add_edge((r, c), (r + 1, c))
                if c > 0: G.add_edge((r, c), (r, c - 1))
                if c < cols - 1: G.add_edge((r, c), (r, c + 1))
    return G

def find_parallel_groups_from_dag(qc):
    """
    Find parallel groups of CX gates from the levels of the DAG.
    """
    # Convert circuit to DAG
    dag = circuit_to_dag(qc)

    # Extract logical layers using DAG layers
    parallel_groups = []
    for i, layer in enumerate(dag.layers()):
        cx_gates_in_layer = [
            node for node in layer['graph'].nodes()
            if hasattr(node, 'op') and (node.op.name == 'cx' or node.op.name == 'cz' or node.op.name == 'cu1')
        ]
        if cx_gates_in_layer:
            parallel_groups.append(cx_gates_in_layer)

    return parallel_groups

def find_logical_qubits(A):
    return [(r, c) for r in range(len(A)) for c in range(len(A[0])) if A[r][c] == 1]

def heuristic(node, end):
    """Manhattan distance heuristic for A* search."""
    return abs(node[0] - end[0]) + abs(node[1] - end[1])

def find_path_astar(G, start, end, blocked_nodes):
    """Find the shortest path using A* algorithm with randomized tie-breaking.
    
    The function randomly shuffles the order of neighbors and adds a small random
    perturbation to the priority to encourage exploration of alternative paths.
    """
    if start in blocked_nodes or end in blocked_nodes:
        return None
    
    frontier = []
    heappush(frontier, (0, start))
    came_from = {start: None}
    cost_so_far = {start: 0}
    
    while frontier:
        _, current = heappop(frontier)
        
        if current == end:
            # Reconstruct the path.
            path = []
            while current is not None:
                path.append(current)
                current = came_from[current]
            return path[::-1]
        
        # Get neighbors as a list and shuffle them for random exploration.
        neighbors = list(G.neighbors(current))
        random.shuffle(neighbors)
        
        for neighbor in neighbors:
            if neighbor in blocked_nodes:
                continue
            new_cost = cost_so_far[current] + 1
            if neighbor not in cost_so_far or new_cost < cost_so_far[neighbor]:
                cost_so_far[neighbor] = new_cost
                # Add a small random perturbation to the priority.
                priority = new_cost + heuristic(neighbor, end) + random.uniform(0, 0.015)
                # priority = new_cost + heuristic(neighbor, end)
                heappush(frontier, (priority, neighbor))
                came_from[neighbor] = current
                
    return None

def is_path_valid(new_path, existing_paths):
    """
    Checks if a path is valid (does not intersect existing paths).
    """
    for path in existing_paths:
        if set(new_path).intersection(set(path)):
            return False
    return True

from itertools import permutations


def iterative_adjust_paths_alternating(current_set, new_gate, G, blocked_nodes, qubit_to_matrix_mapping, logical_qubits, max_iter=100, max_attempts=200):
    """
    Adjusts paths for all gates in the current set plus the new gate.
    For any conflict between two paths, first tries to re-route the later (second) gate
    up to 'max_attempts' times before attempting to re-route the earlier (first) gate.
    Every re-routing takes into account the global blocked nodes (inactive logical qubit positions)
    as well as nodes used by the other paths.
    If local iterative adjustments fail (i.e. conflicts persist for max_iter iterations),
    then a global optimization is attempted by trying all permutations of the gates.
    Returns a list of tuples (control, target, path) if successful, or None if not.
    """
    gates = current_set + [new_gate]
    gate_paths = {}

    for idx, (control, target, _) in enumerate(gates):
        start = qubit_to_matrix_mapping[control]
        end = qubit_to_matrix_mapping[target]
        bn = set(logical_qubits) - {start, end}

        path = find_path_astar(G, start, end, bn)
        if not path:
            return None
        gate_paths[idx] = path

    for iteration in range(max_iter):
        conflict_pairs = []

        for i in range(len(gates)):
            for j in range(i + 1, len(gates)):
                if set(gate_paths[i][1:-1]).intersection(set(gate_paths[j][1:-1])):
                    conflict_pairs.append((i, j))

        if not conflict_pairs:
            return [(gates[i][0], gates[i][1], gate_paths[i]) for i in range(len(gates))]

        for i, j in conflict_pairs:
            control_j, target_j, _ = gates[j]
            start_j = qubit_to_matrix_mapping[control_j]
            end_j = qubit_to_matrix_mapping[target_j]
            bn_j = set(logical_qubits) - {start_j, end_j}
            blocked_for_j = bn_j.union(*(gate_paths[k] for k in range(len(gates)) if k != j))

            new_path_j = None
            for _ in range(max_attempts):
                new_path_j = find_path_astar(G, start_j, end_j, blocked_for_j)
                if new_path_j:
                    break

            if new_path_j:
                gate_paths[j] = new_path_j
                continue

            control_i, target_i, _ = gates[i]
            start_i = qubit_to_matrix_mapping[control_i]
            end_i = qubit_to_matrix_mapping[target_i]
            bn_i = set(logical_qubits) - {start_i, end_i}
            blocked_for_i = bn_i.union(*(gate_paths[k] for k in range(len(gates)) if k != i))

            new_path_i = find_path_astar(G, start_i, end_i, blocked_for_i)
            if new_path_i:
                gate_paths[i] = new_path_i
            else:
                return None

    return None




def update_matrix_for_path_l3(matrix, path, control, target):
    """
    Updates the matrix in 6-bit format:
    - 2 MSBs:
        - 00: Idle qubit.
        - 01: Connection qubit (part of the path).
        - 10: Control qubit (applies the two-qubit gate).
        - 11: Target qubit (receives the two-qubit gate).
    - 4 LSBs:
        - North, South, West, East edges.

    Args:
        matrix: The 2D array representing the qubit layout.
        path: The path between control and target qubits.
        control: The (row, col) location of the control qubit.
        target: The (row, col) location of the target qubit.

    Returns:
        The updated matrix in 6-bit format.
    """
    rows, cols = matrix.shape[:2]

    # Update control qubit
    matrix[control[0], control[1]] = (matrix[control[0], control[1]] & 0b001111) | 0b100000

    # Update target qubit
    matrix[target[0], target[1]] = (matrix[target[0], target[1]] & 0b001111) | 0b110000

    # Update all connection qubits in the path (except control and target)
    for (r1, c1), (r2, c2) in zip(path, path[1:]):
        # Skip the control and target qubits
        if (r1, c1) != control and (r1, c1) != target:
            matrix[r1, c1] = (matrix[r1, c1] & 0b001111) | 0b010000  # Set as connection qubit
        if (r2, c2) != control and (r2, c2) != target:
            matrix[r2, c2] = (matrix[r2, c2] & 0b001111) | 0b010000  # Set as connection qubit

        # Determine the edge direction and update both qubits
        if r1 == r2:  # Same row, horizontal movement
            if c1 < c2:  # Moving east
                matrix[r1, c1] |= 0b00001  # Set east edge active
                matrix[r2, c2] |= 0b00010  # Set west edge active
            else:  # Moving west
                matrix[r1, c1] |= 0b00010  # Set west edge active
                matrix[r2, c2] |= 0b00001  # Set east edge active
        elif c1 == c2:  # Same column, vertical movement
            if r1 < r2:  # Moving south
                matrix[r1, c1] |= 0b00100  # Set south edge active
                matrix[r2, c2] |= 0b01000  # Set north edge active
            else:  # Moving north
                matrix[r1, c1] |= 0b01000  # Set north edge active
                matrix[r2, c2] |= 0b00100  # Set south edge active

    # Reset connection qubits without active edges (010000) to 000000
    for r in range(rows):
        for c in range(cols):
            if matrix[r, c] == 0b010000:
                matrix[r, c] = 0b000000

    return matrix



def draw_graph(G, A, qubit_to_matrix_mapping, gate_paths, step_counter):
    pos = {(r, c): (c, -r) for r in range(A.shape[0]) for c in range(A.shape[1])}
    labels = {
        (row, col): f"Q{qubit_idx}" for qubit_idx, (row, col) in qubit_to_matrix_mapping.items()
    }
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, with_labels=True, labels=labels, node_size=500, node_color="lightblue", font_size=10)
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    for idx, path in enumerate(gate_paths):
        edges = list(zip(path, path[1:]))
        nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=colors[idx % len(colors)], width=2)
    # # plt.title("Graph with Logical Qubit Connections for Parallel CNOT Gates")
    # plt.savefig(f'/home/george/FTQC_Security/targets_5/{step_counter}.pdf')
    # plt.show()
    # plt.close()
    out_dir = '/home/george/FTQC_Security/targets_6/'
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(f"{out_dir}/{step_counter}.pdf")
    plt.show()



def build_connectivity_graph(two_qubit_gates):
    """
    Builds a connectivity graph based on the frequency of interaction between logical qubits.
    Vertices: Logical qubits (labeled by indices, e.g., Q1, Q2)
    Edges: Frequency of interaction between two logical qubits
    """
    connectivity_graph = nx.Graph()

    for control, target in two_qubit_gates:
        control_qubit = f"Q{control}"  # Use logical qubit index as label
        target_qubit = f"Q{target}"  # Use logical qubit index as label

        # Add edge or update weight
        if connectivity_graph.has_edge(control_qubit, target_qubit):
            connectivity_graph[control_qubit][target_qubit]['weight'] += 1
        else:
            connectivity_graph.add_edge(control_qubit, target_qubit, weight=1)

    return connectivity_graph

def draw_connectivity_graph(connectivity_graph):
    """
    Draws the connectivity graph with weights displayed as edge labels.
    Vertices are labeled with logical qubit indices.
    """
    pos = nx.kamada_kawai_layout(connectivity_graph)  # Position nodes for a visually appealing layout. Other choices circular_layout, kamada_kawai_layout
    plt.figure(figsize=(8, 8))
    nx.draw(
        connectivity_graph, pos, with_labels=True, node_size=500, node_color="lightblue", font_size=10
    )
    edge_labels = nx.get_edge_attributes(connectivity_graph, 'weight')
    nx.draw_networkx_edge_labels(connectivity_graph, pos, edge_labels=edge_labels)
    plt.title("Connectivity Graph")
    plt.show()

from collections import deque

def reconstruct_dag_from_traces(filename='l3_traces.json', qubit_mapping=None):
    """
    Reconstruct the DAG from traces by correctly pairing controls
    and targets even when multiple gates occur in parallel.
    """
    # Load the traces
    with open(filename, 'r') as f:
        traces = json.load(f)

    dag = nx.DiGraph()
    all_gates = []

    # 1) Flatten each trace step into (step, control, target, gate_type)
    for step_str, matrix in sorted(traces.items(), key=lambda x: int(x[0])):
        step = int(step_str)
        rows, cols = len(matrix), len(matrix[0])
        control_qubits = []
        target_qubits  = []

        # extract raw control & target coordinates
        for r in range(rows):
            for c in range(cols):
                cell = int(matrix[r][c], 2)
                role = (cell >> 4) & 0b11
                if role == 0b10:
                    control_qubits.append((r, c))
                elif role == 0b11:
                    target_qubits.append((r, c))

        # build the connection graph from the 4 LSBs
        G_conn = nx.Graph()
        for r in range(rows):
            for c in range(cols):
                cell = int(matrix[r][c], 2)
                north = (cell >> 3) & 1
                south = (cell >> 2) & 1
                west  = (cell >> 1) & 1
                east  = (cell >> 0) & 1
                if north and r > 0:
                    G_conn.add_edge((r, c), (r - 1, c))
                if south and r < rows - 1:
                    G_conn.add_edge((r, c), (r + 1, c))
                if west and c > 0:
                    G_conn.add_edge((r, c), (r, c - 1))
                if east and c < cols - 1:
                    G_conn.add_edge((r, c), (r, c + 1))

        # BFS‐pair each control to its true target
        pairs = []
        for ctrl in control_qubits:
            queue = deque([ctrl])
            seen  = {ctrl}
            match = None
            while queue and match is None:
                u = queue.popleft()
                if u in target_qubits:
                    match = u
                    break
                for nbr in G_conn.neighbors(u):
                    if nbr not in seen:
                        seen.add(nbr)
                        queue.append(nbr)
            if match is None:
                raise RuntimeError(f"No matching target for control {ctrl} at step {step}")
            pairs.append((ctrl, match))
            target_qubits.remove(match)

        # record each paired gate
        for control, target in pairs:
            # map to labels
            control_label = f"Q{control[0]}_{control[1]}"
            target_label  = f"Q{target[0]}_{target[1]}"
            if qubit_mapping:
                control_label = qubit_mapping.get(control, control_label)
                target_label  = qubit_mapping.get(target,  target_label)
            gate_type = 'cx'  # or 'cz' as appropriate

            all_gates.append((step, control_label, target_label, gate_type))

    # 2) Reverse-iterate to build the dependency DAG
    dag = nx.DiGraph()
    qubit_status = {}
    num_gates = len(all_gates)

    for idx, (step, control, target, gate_type) in enumerate(all_gates[::-1]):
        true_idx = num_gates - idx - 1
        gate_id  = f"G{true_idx}"
        dag.add_node(gate_id,
                     control=control,
                     target=target,
                     step=step,
                     gate=gate_type)

        # link to next use of each qubit
        for dep in (qubit_status.get(control), qubit_status.get(target)):
            if dep:
                dag.add_edge(gate_id, dep)

        # update last‐seen pointers
        qubit_status[control] = gate_id
        qubit_status[target]  = gate_id

    return dag


def visualize_dag(G, title="Reconstructed DAG Visualization", node_size= 200, font_size=8):
    pos = nx.kamada_kawai_layout(G)
    # plt.figure(figsize=(10, 7))
    nx.draw(G, pos, with_labels=True, node_size=node_size, font_size=font_size)
    # labels = {n: G.nodes[n].get('gate', str(n)) for n in G.nodes()}
    # nx.draw_networkx_labels(G, pos, labels=labels, font_size=font_size)
    # plt.title(title)
    out_dir = '/home/george/FTQC_Security/targets_6/'
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(f"{out_dir}/figure.pdf")
    plt.axis('off')
    plt.show()

def save_traces_reconstruct_dag(n, path, initial_matrix):
    """
    Function for extracting traces and reconstructing the DAG.

    Args:
        n (int): Number of logical qubits.
        path (str): Path to the QASM file.

    Returns:
        nx.DiGraph: The reconstructed DAG.
    """
    # Extract benchmark name from the path
    benchmark_name = path.split('/')[-1].replace('.qasm', '')

    # Load the Quantum Circuit
    qc = QuantumCircuit.from_qasm_file(path)

    # Generate random matrix with logical qubit locations
    import random

    def generate_random_matrix(n, size=(30, 30)):
        matrix = np.zeros(size, dtype=int)
        if n > size[0] * size[1]:
            raise ValueError("Number of 1s exceeds matrix capacity")

        # Generate unique random positions
        positions = random.sample(range(size[0] * size[1]), n)
        
        for pos in positions:
            row, col = divmod(pos, size[1])
            matrix[row, col] = 1

        return matrix

    # initial_matrix = generate_random_matrix(n)
    # initial_matrix = generate_fixed_square_sparse_matrix_v8(n)
    qubit_mapping = extract_qubit_mapping(initial_matrix)

    # Transform the matrix to 5-bit binary format
    transformed_matrix = np.array([
        [0b100000 if cell == 1 else 0b000000 for cell in row]
        for row in initial_matrix
    ], dtype=int)

    # Create graph and identify logical qubits
    G = create_graph(initial_matrix)
    logical_qubits = find_logical_qubits(initial_matrix)
    two_qubit_gates = get_two_qubit_gates(qc)
    qubit_to_matrix_mapping = {i: logical_qubits[i] for i in range(number_of_qubits(qc))}

    # Step 1: Extract parallel groups using the DAG layers
    parallel_groups = find_parallel_groups_from_dag(qc)

    parallel_gate_sets = []  # Initialize the list to store all parallel gate sets
    traces = OrderedDict()
    step_counter = 0

    # Step 2: Process each parallel group
    for group in parallel_groups:
        current_set = []
        for node in group:
            control, target = [qc.find_bit(q).index for q in node.qargs]
            start = qubit_to_matrix_mapping[control]
            end = qubit_to_matrix_mapping[target]
            blocked_nodes = set(logical_qubits) - {start, end}

            # Adjust paths iteratively
            new_gate = (control, target, None)
            updated_set = iterative_adjust_paths_alternating(
                current_set, 
                new_gate, 
                G, 
                blocked_nodes, 
                qubit_to_matrix_mapping, 
                logical_qubits, 
                max_iter=20
            )

            if updated_set:
                current_set = updated_set
            else:
                parallel_gate_sets.append(current_set)
                current_set = [(control, target, find_path_astar(G, start, end, blocked_nodes))]

        if current_set:
            parallel_gate_sets.append(current_set)


    # Step 3: Process and update each parallel gate set
    for gate_set in parallel_gate_sets:
        gate_paths = []
        # print("\nParallel Gate Set:")
        updated_matrix = transformed_matrix.copy()

        for control, target, path in gate_set:
            # print(f"Control: Q{control}, Target: Q{target}, Path: {path}")
            gate_paths.append(path)
            # Update the matrix in 6-bit L3 format
            updated_matrix = update_matrix_for_path_l3(updated_matrix, path, qubit_to_matrix_mapping[control], qubit_to_matrix_mapping[target])

        # Construct a unique key for each step
        # matrix_key = f"Level {step_counter}: {str(updated_matrix)}"

        # --- zero out any pure-logical-qubit flags (0b100000) ---
        updated_matrix[updated_matrix == 0b100000] = 0b000000
        
        traces[step_counter] = updated_matrix
        step_counter += 1

        # draw_graph(G, initial_matrix, qubit_to_matrix_mapping, gate_paths, step_counter)
        # visualize_dag(reconstructed_dag, title=f"Reconstructed DAG - {benchmark_name}")

    # Compute and print average parallelism degree
    # num_parallel_steps = len(parallel_groups)
    # if num_parallel_steps > 0:
    #     total_gates = sum(len(group) for group in parallel_groups)
    #     avg_parallelism = total_gates / num_parallel_steps
    #     print(f"Average parallelism degree: {avg_parallelism:.2f}")
    # else:
    #     print("No parallel two-qubit gates found.")


        # Display the updated transformed matrix in binary format
        # print("Updated Transformed Matrix (6-bit Binary):")
        # print(np.vectorize(lambda x: f"{x:06b}")(updated_matrix))  # Display in binary format

        # draw_graph(G, initial_matrix, qubit_to_matrix_mapping, gate_paths, step_counter)

        # Save traces to a JSON file with dynamic naming
        # trace_filename = f"l3_traces_{benchmark_name}.json"
        # save_traces_as_json(traces, trace_filename)

        # Reconstruct the DAG from the saved traces
        # reconstructed_dag = reconstruct_dag_from_traces(filename=trace_filename, qubit_mapping=qubit_mapping)

        # Visualize the DAG
        # visualize_dag(reconstructed_dag, title=f"Reconstructed DAG - {benchmark_name}")

    # # Return the DAG object
    # return reconstructed_dag
    return traces