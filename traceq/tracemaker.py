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

import json


# ✅ Function to save matrices as JSON
def save_traces_as_json(traces, filename='l3_traces.json'):
    """
    Save the traces as JSON format. 6 bit binary format is used in each cell.
    
    Args:
        traces (dict): A dictionary of traces with keys
            representing the level and values as matrices.
        filename (str): The name of the JSON file to save.
    """

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
    """
    Load traces from a JSON file. The matrices are in 6-bit binary format.
    
    Args:
        filename (str): The name of the JSON file to load.
        
    Returns:
        dict: A dictionary of traces with keys representing the level and values as matrices.
    """

    with open(filename, 'r') as f:
        traces = json.load(f)
    print(f"✅ Traces loaded from {filename}")
    return traces

def generate_random_matrix(n, qc):
    # generate docstring
    """
    Generate a random matrix with 1s at random positions.
    Use qc.num_qubits number of 1s.
    
    Args:
        n (int): The size of the matrix.
        qc (QuantumCircuit): The quantum circuit.
        
    Returns:
        np.ndarray: A random matrix with 1s at random positions.
    """
    matrix = np.zeros((n, n), dtype=int)
    ones_positions = np.random.choice(n * n, qc.num_qubits, replace=False)
    for pos in ones_positions:
        row, col = divmod(pos, n)
        matrix[row, col] = 1
    return matrix

def get_two_qubit_gates(qc):
    """Get 2 qubit gates from a quantum circuit.

    Args:
        qc (QuantumCircuit): quantum circuit

    Returns:
        List of tuples: List of 2qubit gates
    """
    two_qubit_gates = [
        (qc.find_bit(inst[1][0])[0], qc.find_bit(inst[1][1])[0])
        for inst in qc.data if len(inst[1]) == 2
    ]
    return two_qubit_gates

def create_graph(A):
    """Builds connectivity graph of 2D array

    Args:
        A (2D array): Surface code plane

    Returns:
        nx.Graph: Networkx Connectivity Graph
    """
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

    Args:
        qc (QuantumCircuit): Circuit

    Returns:
        List of lists: List of phases where all gates may be scheduled in parallel.
    """
    # Convert circuit to DAG
    dag = circuit_to_dag(qc)

    # Extract logical layers using DAG layers
    parallel_groups = []
    for i, layer in enumerate(dag.layers()):
        cx_gates_in_layer = [
            node for node in layer['graph'].nodes()
            if hasattr(node, 'op') and node.op.name == 'cx'
        ]
        if cx_gates_in_layer:
            parallel_groups.append(cx_gates_in_layer)

    return parallel_groups

def find_logical_qubits(A):
    """
    Find logical qubits from the 2D array.
    
    Args:
        A (np.ndarray): 2D array representing the surface code plane.
        If A[r][c] == 1, it represents a logical qubit.
        
    Returns:
        List of tuples: List of logical qubits as (row, col) tuples.
    """
    return [(r, c) for r in range(len(A)) for c in range(len(A[0])) if A[r][c] == 1]

def heuristic(node, end):
    """Manhattan distance heuristic for A* search."""
    return abs(node[0] - end[0]) + abs(node[1] - end[1])

def find_shortest_path(G, start, end, blocked_nodes):
    """Finds shortest path on graph G without using any blocked nodes.

    Args:
        G (nx.Graph): Connectivity graph of underlying device
        start (node): Start node on G
        end (node): End node on G
        blocked_nodes (List): List of blocked nodes
    """
    
    if start in blocked_nodes or end in blocked_nodes:
        return None
    
    G_restricted = G.copy()
    G_restricted.remove_nodes_from(blocked_nodes)
    
    try:
        return nx.shortest_path(G_restricted, start, end)
    except nx.NetworkXNoPath:
        return None


def find_path_astar(G, start, end, blocked_nodes):
    """Find the shortest path using A* algorithm."""
    if start in blocked_nodes or end in blocked_nodes:
        return None
    
    frontier = []
    heappush(frontier, (0, start))
    came_from = {start: None}
    cost_so_far = {start: 0}
    
    while frontier:
        _, current = heappop(frontier)
        
        if current == end:
            # Reconstruct path
            path = []
            while current:
                path.append(current)
                current = came_from[current]
            return path[::-1]
        
        for neighbor in G.neighbors(current):
            if neighbor in blocked_nodes:
                continue
            new_cost = cost_so_far[current] + 1
            if neighbor not in cost_so_far or new_cost < cost_so_far[neighbor]:
                cost_so_far[neighbor] = new_cost
                priority = new_cost + heuristic(neighbor, end)
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

def iterative_adjust_paths_with_reevaluation(current_set, new_gate, G, blocked_nodes):
    """
    Iteratively adjusts paths for all gates in the current set and the new gate to achieve better parallelization.
    Reevaluates the last path when conflicts occur.
    """
    gates = current_set + [new_gate]
    gate_paths = {}

    # Initialize paths for all gates
    for idx, (control, target, _) in enumerate(gates):
        start = qubit_to_matrix_mapping[control]
        end = qubit_to_matrix_mapping[target]
        blocked_nodes = set(logical_qubits) - {start, end}  # Recalculate blocked nodes
        path = find_path_astar(G, start, end, blocked_nodes)

        if not path:
            return None  # Parallelization fails if any gate cannot find a valid path

        gate_paths[idx] = path

    # Iterative conflict resolution
    for iteration in range(len(gates) * 2):  # Limit iterations to prevent infinite loops
        conflicts = []

        # Check for conflicts
        for i, path in gate_paths.items():
            for j, other_path in gate_paths.items():
                if i != j and not is_path_valid(path, [other_path]):
                    conflicts.append(i)
                    break

        if not conflicts:
            # No conflicts; all paths are valid
            return [(gates[i][0], gates[i][1], gate_paths[i]) for i in range(len(gates))]

        # Adjust conflicting paths, starting from the last path
        for idx in reversed(conflicts):
            control, target, _ = gates[idx]
            start = qubit_to_matrix_mapping[control]
            end = qubit_to_matrix_mapping[target]

            # Block all other paths except the current one
            blocked_for_adjustment = set()
            for i, p in gate_paths.items():
                if i != idx:
                    blocked_for_adjustment.update(p)

            # Recalculate path with blocking
            new_path = find_path_astar(G, start, end, blocked_for_adjustment)

            if not new_path:
                return None  # Fail parallelization if no alternative path is found

            gate_paths[idx] = new_path

    # Return None if valid paths are not achieved within the limit
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
                matrix[r1, c1] |= 0b000001  # Set east edge active
                matrix[r2, c2] |= 0b000010  # Set west edge active
            else:  # Moving west
                matrix[r1, c1] |= 0b000010  # Set west edge active
                matrix[r2, c2] |= 0b000001  # Set east edge active
        elif c1 == c2:  # Same column, vertical movement
            if r1 < r2:  # Moving south
                matrix[r1, c1] |= 0b000100  # Set south edge active
                matrix[r2, c2] |= 0b001000  # Set north edge active
            else:  # Moving north
                matrix[r1, c1] |= 0b001000  # Set north edge active
                matrix[r2, c2] |= 0b000100  # Set south edge active

    # Reset connection qubits without active edges (010000) to 000000
    for r in range(rows):
        for c in range(cols):
            if matrix[r, c] == 0b010000:
                matrix[r, c] = 0b000000

    return matrix



def draw_graph(G, A, qubit_to_matrix_mapping, gate_paths):
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
    plt.title("Graph with Logical Qubit Connections for Parallel CNOT Gates")
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




def make_l3_trace(qc, initial_layout):
    """Given a quantum circuit and initial layout, generate L3 traces.

    Args:
        qc (QuantumCircuit): Circuit
        initial_layout (Matrix): Initial layout of surface code plane

    Returns:
        OrderedDict: Traces constrained by physical interactions
    """

    # initialize the assignment of qubits to the surface code plane
    G = create_graph(initial_layout)
    logical_qubits = find_logical_qubits(initial_layout)
    qubit_to_matrix_mapping = {i: logical_qubits[i] for i in range(qc.num_qubits)}
    
    # Extract parallel groups using the DAG layers
    # There may still be conflicts coming from scheduling!
    parallel_groups = find_parallel_groups_from_dag(qc)

    parallel_gate_sets = []  # Initialize the list to store all parallel gate sets

    traces = OrderedDict()  # Use OrderedDict to maintain the order of updates
    step_counter = 0  # Counter for unique trace keys


    # Step 2: Build trace for each parallel group
    for parallel_pass in parallel_groups:
        # TODO: this should instead be implemented as a while loop with heap for optimal
        # this is because of the following: during the scheduling, one can imagine that there are conflicts that prevent one edge from being scheduled but then can be scheduled 
        # e.g. cnot 1, 2 | cnot 3, 4 | cnot 5, 6 | cnot 2, 3
        # say that cnot 5,6 cannot be scheduled due to physical constraints
        # then, we should implement two passes:
        #  cnot 1, 2 | cnot 3, 4 // cnot 5, 6 | cnot 2, 3
        # but our current program would do:
        # cnot 1, 2 | cnot 3, 4 // cnot 5, 6 // cnot 2, 3
        current_set = []  # Reset the current set for each layer
        blocked_nodes = set(logical_qubits)
        
        for gate_node in parallel_pass:
            control, target = [qc.find_bit(q).index for q in gate_node.qargs]
            start = qubit_to_matrix_mapping[control]
            end = qubit_to_matrix_mapping[target]
            
            # attempt to find path
            
            proposed_path = find_shortest_path(G, start, end, blocked_nodes - {start, end})
            if proposed_path:
                current_set.append((control, target, proposed_path))
            else:
                # Finalize current set and start a new one
                parallel_gate_sets.append(current_set)
                
                blocked_nodes = set(logical_qubits)
                proposed_path = find_shortest_path(G, start, end, blocked_nodes - {start, end})
                
                if not proposed_path:
                    raise ValueError("No path found for gate")
                
                current_set = [(control, target, proposed_path)]
                
            blocked_nodes.update(proposed_path)

        if current_set:
            parallel_gate_sets.append(current_set)

    # let's do an assertion to check validity
    for parallel_set in parallel_gate_sets:
        involved_qubits = set()
        for control, target, path in parallel_set:
            assert path is not None, "Invalid path found in parallel set"
            assert control not in involved_qubits, "Control qubit already involved in parallel set"
            assert target not in involved_qubits, "Target qubit already involved in parallel set"
            # check that the path isn't in the involved qubits
            for node in path:
                assert node not in involved_qubits, "Node already involved in parallel set"
            involved_qubits.update(path)

    # Step 3: Convert each parallel set to be a trace
    transformed_matrix = np.array(initial_layout, dtype=int)
    
    for idx, parallel_set in enumerate(parallel_gate_sets):
        print("\nParallel Gate Set:")
        updated_matrix = transformed_matrix.copy()

        for control, target, path in parallel_set:
            print(f"Control: Q{control}, Target: Q{target}, Path: {path}")
            
            updated_matrix = update_matrix_for_path_l3(updated_matrix, path, qubit_to_matrix_mapping[control], qubit_to_matrix_mapping[target])
            
        traces[idx] = {"trace": updated_matrix, "parallel_set": parallel_set}

    # Draw the connectivity graph
    return traces

def _reduce_l3_to_l1(trace_mat):
    # given matrix with 6 bitstring set, reduce to free/busy matrix
    # Recall format:
    #     Updates the matrix in 6-bit format:
    # - 2 MSBs:
    #     - 00: Idle qubit.
    #     - 01: Connection qubit (part of the path).
    #     - 10: Control qubit (applies the two-qubit gate).
    #     - 11: Target qubit (receives the two-qubit gate).
    # - 4 LSBs:
    #     - North, South, West, East edges.
    
    # I.e.: if MSBs are not 00, then the qubit is busy
    # i.e., if # > 0b010000 = 16, then the qubit is busy
    
    return np.array(trace_mat > 0b010000, dtype=int)

def reduce_l3_to_l1(traces):
    # generate docstring
    """
    Given a set of traces, reduce to L1 traces.
    
    Args:
        traces (dict): A dictionary of traces with keys
            representing the level and values as matrices.
            
    Returns:
        dict: A dictionary of traces with keys representing the level
            and values as L1 matrices.
    """
    out = OrderedDict()
    for key, trace in traces.items():
        out[key] = {"trace": _reduce_l3_to_l1(trace["trace"])}
    
    return out


# def visualize_clustering_coefficients(G):
#     """Visualize clustering coefficients."""
#     local_cc = nx.clustering(G)
#     avg_cc = nx.average_clustering(G)
#     global_cc = nx.transitivity(G)

#     plt.figure(figsize=(12, 8))
#     # 1. Draw the graph
#     plt.subplot(1, 2, 1)
#     pos = nx.spring_layout(G)
#     nx.draw(
#         G, pos, with_labels=True, node_color='lightblue',
#         edge_color='gray', node_size=500, font_size=10
#     )
#     plt.title("Connectivity Graph")

#     # 2. Plot clustering coefficients
#     plt.subplot(1, 2, 2)
#     local_cc_values = list(local_cc.values())
#     plt.hist(local_cc_values, bins=10, color='skyblue', edgecolor='black', alpha=0.7)
#     plt.axvline(avg_cc, color='red', linestyle='--', label=f"Average CC: {avg_cc:.2f}")
#     plt.axvline(global_cc, color='green', linestyle='--', label=f"Global CC: {global_cc:.2f}")
#     plt.xlabel("Clustering Coefficient")
#     plt.ylabel("Frequency")
#     plt.title("Distribution of Clustering Coefficients")
#     plt.legend()
#     plt.tight_layout()
#     plt.show()

# def visualize_betweenness_centrality(G):
#     """Visualize betweenness centrality."""
#     centrality = nx.betweenness_centrality(G, weight='weight')
#     pos = nx.spring_layout(G)
#     node_sizes = [500 + 3000 * centrality[node] for node in G.nodes]

#     plt.figure(figsize=(10, 8))
#     nx.draw(
#         G, pos, with_labels=True, node_color='lightblue',
#         edge_color='gray', node_size=node_sizes, font_size=10
#     )
#     plt.title("Betweenness Centrality (Node Size Represents Centrality)")
#     plt.show()

# def visualize_communities(G):
#     """Visualize communities using the Girvan-Newman algorithm."""
#     from networkx.algorithms.community import girvan_newman
#     communities = next(girvan_newman(G))
#     pos = nx.spring_layout(G)
#     color_map = ['red', 'blue', 'green', 'orange', 'purple']
#     community_colors = {}
#     for i, community in enumerate(communities):
#         for node in community:
#             community_colors[node] = color_map[i % len(color_map)]

#     node_colors = [community_colors[node] for node in G.nodes]

#     plt.figure(figsize=(10, 8))
#     nx.draw(
#         G, pos, with_labels=True, node_color=node_colors,
#         edge_color='gray', node_size=500, font_size=10
#     )
#     plt.title("Communities in Connectivity Graph")
#     plt.show()

# def visualize_edge_frequency_distribution(G):
#     """Visualize the frequency distribution of edges in the connectivity graph."""
#     weights = nx.get_edge_attributes(G, 'weight').values()
#     plt.figure(figsize=(10, 6))
#     plt.hist(weights, bins=10, color='lightgreen', edgecolor='black', alpha=0.7)
#     plt.xlabel("Edge Weight (Frequency of Interaction)")
#     plt.ylabel("Count")
#     plt.title("Edge Weight Distribution")
#     plt.tight_layout()
#     plt.show()






    
