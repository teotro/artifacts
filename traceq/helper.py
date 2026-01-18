import warnings
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit
from collections import OrderedDict, deque, defaultdict
from qiskit.converters import circuit_to_dag

import matplotlib.pyplot as plt
from matplotlib import cm, colors

import time

from collections import OrderedDict
import numpy as np





def repr_6bit_matrix(matrix):
    # print matrix where all elements are 6bit
    return "\n".join(" ".join(f"{x:07b}" for x in row) for row in matrix)


def perfect_matchings(nodes):
    """
    Recursively generate all perfect matchings of an even‐sized list 'nodes'.
    Each matching is a list of 2‐tuples.
    """
    if not nodes:
        yield []
        return
    u = nodes[0]
    for v in nodes[1:]:
        rest = [x for x in nodes if x not in (u, v)]
        for pm in perfect_matchings(rest):
            yield [(u, v)] + pm

def filter_non_overlapping_paths(solutions):
    """
    Filters out solutions where any of the paths in a solution overlap in tiles.
    Assumes solutions are in the form: [ [ ((q1, q2), path), ((q3, q4), path2), ... ], ... ]
    """
    valid = []
    for sol in solutions:
        seen = set()
        overlap = False
        for _, path in sol:
            for tile in path:
                if tile in seen:
                    overlap = True
                    break
                seen.add(tile)
            if overlap:
                break
        if not overlap:
            valid.append(sol)
    return valid

def enumerate_all_connections_backtrack_dfs(trace_matrix):
    """
    Backtracking search that finds all non-crossing ways to pair up
    the logical-unknown qubits (0b1110000) by paths over the wire-tiles
    (0b1010000), using an on-the-fly DFS to generate only those paths
    that obey “no crossing” and “include at least one wire-tile.”
    """
    n_rows, n_cols = trace_matrix.shape

    # 1) Identify qubit endpoints and wire-tiles
    qubits = [tuple(rc) for rc in np.argwhere(trace_matrix == 0b1110000)]
    wires  = {tuple(rc)
              for rc in np.argwhere((trace_matrix & 0b1110000) == 0b1010000)}
    target_cover = set(wires)

    # 2) Build subgraph on qubits + wires
    G = nx.grid_2d_graph(n_rows, n_cols)
    keep = set(qubits) | wires
    G.remove_nodes_from([n for n in G if n not in keep])

    # ####Comment this out
    # G = nx.grid_2d_graph(n_rows, n_cols)
    # keep = set(qubits) | wires
    # G.remove_nodes_from([n for n in G if n not in keep])

    # # ✅ Add direct edges between adjacent qubit endpoints (if not removed above)
    # for u in qubits:
    #     for v in qubits:
    #         if u != v and v in G[u]:
    #             G.add_edge(u, v)


    solutions = []
    visited_fails = set()  # memo of (remaining_qubits, used) states that fail

    def constrained_paths(u, v, used):
        """
        DFS generator of simple paths from u to v in G that:
        - never revisit nodes
        - never step into any in `used`
        - include at least one wire-tile
        """
        def dfs(path, seen_wires):
            node = path[-1]
            if node == v:
                if seen_wires:
                    yield list(path)
                return
            # if node == v:
            #     if seen_wires or len(path) == 2:
            #         yield list(path)
            #     return
            for nbr in G.neighbors(node):
                if nbr in path or nbr in used:
                    continue
                yield from dfs(path + [nbr], seen_wires or (nbr in wires))

        # start with u; mark if u itself is a wire
        yield from dfs([u], u in wires)

    def backtrack(pairs_left, used, current_solution):
        state_key = (frozenset(pairs_left), frozenset(used))
        if state_key in visited_fails:
            return

        if not pairs_left:
            if used == target_cover:
                solutions.append(list(current_solution))
            else:
                visited_fails.add(state_key)
            return

        found_any = False
        # fixed ordering: pick the first qubit available
        u = pairs_left[0]
        for v in pairs_left[1:]:
            rest = [q for q in pairs_left if q not in (u, v)]
            for path in constrained_paths(u, v, used):
                path_set = set(path)
                # must cover only within target_cover
                new_used = used | (path_set & target_cover)
                if not new_used.issubset(target_cover):
                    continue
                # recurse
                current_solution.append(((u, v), path))
                backtrack(rest, new_used, current_solution)
                current_solution.pop()
                found_any = True

        if not found_any:
            visited_fails.add(state_key)

    backtrack(qubits, set(), [])
    return solutions

def filter_same_length_duplicates(solutions):
    """
    Given a list of solutions, each
      sol = [ ((q1,q2), path1), ((q3,q4), path2), … ]
    Returns only one rep per unique set of (q1, q2) pairs as ints, ignoring route/path.
    All pairs are saved as ((int, int), (int, int)), never np.int64.
    """
    def to_int_tuple(pair):
        return (int(pair[0][0]), int(pair[0][1])), (int(pair[1][0]), int(pair[1][1]))

    seen = {}
    for sol in solutions:
        # Canonical signature with all pairs as int and sorted
        sig = frozenset(
            tuple(sorted((int(pair[0][0]), int(pair[0][1]), int(pair[1][0]), int(pair[1][1]))))
            for pair, _ in sol
        )
        if sig not in seen:
            # Store solution as ints everywhere
            sol_as_ints = [((int(pair[0][0]), int(pair[0][1])), (int(pair[1][0]), int(pair[1][1])), path) for pair, path in sol]
            # Collapse to ((q1,q2), path), where q1 and q2 are each (int, int)
            sol_as_ints = [((int(pair[0]), int(pair[1])), (int(pair2[0]), int(pair2[1])), path)
                           if isinstance(pair[0], (tuple, list)) and isinstance(pair2[0], (tuple, list)) else ((int(pair[0]), int(pair[1])), (int(pair2[0]), int(pair2[1])), path)
                           for (pair, pair2), path in [((pair[0], pair[1]), path) for pair, path in sol]]
            # Or, more generally, if your pairs are already (q1, q2), each is a tuple of ints:
            sol_as_ints = [((int(q1[0]), int(q1[1])), (int(q2[0]), int(q2[1])), path) for (q1, q2), path in sol]
            seen[sig] = sol_as_ints
    # Return solutions as [ [ ((q1,q2), path), ... ], ... ] with q1, q2 tuples of ints
    return [
        [((q1, q2), path) for (q1, q2, path) in sol]
        for sol in seen.values()
    ]

def infer_ambiguity_and_count(trace_matrix):
    """
    Returns (ambiguous: bool, logical_unknown_count: int)
    - ambiguous is True if any cell has MSB (bit-6) = 1
    - logical_unknown_count is the number of cells exactly == 0b1110000
    """
    # mask for MSB=1 (bit-6 in a 7-bit encoding)
    msb_mask = (trace_matrix & 0b1000000) != 0
    ambiguous = bool(msb_mask.any())

    # count the “pure logical‐unknown” pattern 0b1110000
    logical_unknown_count = int(np.sum(trace_matrix == 0b1110000))

    return ambiguous, logical_unknown_count

def visualize_dag(G, title="Reconstructed DAG Visualization", node_size= 200, font_size=8):
    pos = nx.kamada_kawai_layout(G)
    # plt.figure(figsize=(10, 7))
    nx.draw(G, pos, with_labels=True, node_size=node_size, font_size=font_size)
    # labels = {n: G.nodes[n].get('gate', str(n)) for n in G.nodes()}
    # nx.draw_networkx_labels(G, pos, labels=labels, font_size=font_size)
    plt.title(title)
    plt.axis('off')
    plt.show()


from multiprocessing import Pool, cpu_count

def enumerate_all_connections_backtrack_dfs_parallel(trace_matrix):
    """
    Parallel version of enumerate_all_connections_backtrack_dfs.
    Automatically adapts the number of processes and reports it.
    """
    # 1) Identify qubit endpoints and wire-tiles
    qubits = [tuple(rc) for rc in np.argwhere(trace_matrix == 0b1110000)]
    wires = {
        tuple(rc)
        for rc in np.argwhere((trace_matrix & 0b1110000) == 0b1010000)
    }

    # nothing to pair
    if len(qubits) < 2:
        return []

    u = qubits[0]
    targets = qubits[1:]
    tasks = [(trace_matrix, u, v) for v in targets]

    num_tasks = len(tasks)
    max_cores = cpu_count()
    num_procs = min(num_tasks, max_cores)

    print(f"[DFS] Spawning {num_procs} process{'es' if num_procs > 1 else ''} for {num_tasks} (u,v) task{'s' if num_tasks > 1 else ''}.")

    if num_procs <= 1:
        results = map(_worker_pair_and_backtrack, tasks)
    else:
        with Pool(processes=num_procs) as pool:
            results = pool.map(_worker_pair_and_backtrack, tasks)

    solutions = [sol for sublist in results for sol in sublist]
    return solutions


def _worker_pair_and_backtrack(args):
    """
    Worker that handles one initial match (u,v),
    then backtracks on the remaining qubits.
    """
    trace_matrix, u, v = args

    # rebuild the little graph & data
    n_rows, n_cols = trace_matrix.shape
    wires = {
        tuple(rc)
        for rc in np.argwhere((trace_matrix & 0b1110000) == 0b1010000)
    }
    target_cover = set(wires)

    qubits = [tuple(rc) for rc in np.argwhere(trace_matrix == 0b1110000)]
    G = nx.grid_2d_graph(n_rows, n_cols)
    keep = set(qubits) | wires
    G.remove_nodes_from([n for n in G if n not in keep])

    # ##### Comment this out
    # # ✅ Add direct edges between adjacent qubit endpoints
    # for q1 in qubits:
    #     for q2 in qubits:
    #         if q1 != q2 and q2 in G[q1]:
    #             G.add_edge(q1, q2)


    rest = [q for q in qubits if q not in (u, v)]
    solutions = []
    visited_fails = set()

    def constrained_paths(u0, v0, used0):
        """Same as before, but closed over G and wires."""
        def dfs(path, seen_wire):
            node = path[-1]
            if node == v0:
                if seen_wire:
                    yield list(path)
                return
            # if node == v0:
            #     if seen_wire or len(path) == 2:  # direct neighbor
            #         yield list(path)
            #     return

            for nbr in G.neighbors(node):
                if nbr in path or nbr in used0:
                    continue
                yield from dfs(path + [nbr], seen_wire or (nbr in wires))

        yield from dfs([u0], u0 in wires)

    def backtrack(pairs_left, used, curr_solution):
        key = (frozenset(pairs_left), frozenset(used))
        if key in visited_fails:
            return
        if not pairs_left:
            if used == target_cover:
                solutions.append(curr_solution.copy())
            else:
                visited_fails.add(key)
            return

        found = False
        u1 = pairs_left[0]
        for v1 in pairs_left[1:]:
            nxt = [q for q in pairs_left if q not in (u1, v1)]
            for path in constrained_paths(u1, v1, used):
                new_used = used | (set(path) & target_cover)
                if not new_used.issubset(target_cover):
                    continue
                curr_solution.append(((u1, v1), path))
                backtrack(nxt, new_used, curr_solution)
                curr_solution.pop()
                found = True
        if not found:
            visited_fails.add(key)

    # for the initial (u,v) we need to generate all prefix paths
    for prefix in constrained_paths(u, v, set()):
        used0 = set(prefix) & target_cover
        if not used0.issubset(target_cover):
            continue
        backtrack(rest, used0, [((u, v), prefix)])

    return solutions

def reconstruct_dag_from_l3_traces(recovered_l3_traces, layout_type, benchmark_name, benchmark_pert):
# def reconstruct_dag_from_l3_traces(refined_l3_traces, layout_type, benchmark_name):
    import numpy as np
    import networkx as nx
    from collections import deque, defaultdict
    from bisect import bisect_right

    total_time_dfs = 0
    dfs_times_per_step = []

    # ------------------------------------------------------------------------------
    # (1) Collect all unambiguous detections and all ambiguous‐step solutions
    # ------------------------------------------------------------------------------
    all_gates = []                    # will hold (step, u, v, 'cx', False) for unambiguous
    ambiguous_step_solutions = {}     # maps step → [ [ (u,v), … ],  … ] for each branch

    for step_str, entry in sorted(recovered_l3_traces.items(), key=lambda x: int(x[0])):
    # for step_str, entry in sorted(refined_l3_traces.items(), key=lambda x: int(x[0])):
        step = int(step_str)
        mat = entry["trace"]
        rows, cols = mat.shape
        ambiguous = bool(np.any(mat & 0b1000000))
        gates_this_step = set()

        # --- (1a) detect unambiguous qubit‐endpoints (no MSB=1) via endpoint matching ---
        endpoints = []
        for r in range(rows):
            for c in range(cols):
                cell = int(mat[r, c])
                #  (cell & 0b1000000)==0 → no “ambiguity bit” set
                #  ((cell >> 4)&0b11)==0b11 → this is a qubit‐endpoint
                if ((cell & 0b1000000) == 0) and (((cell >> 4) & 0b11) == 0b11):
                    endpoints.append((r, c))

        if endpoints:
            G_conn = nx.Graph()
            G_conn.add_nodes_from(endpoints)
            for r in range(rows):
                for c in range(cols):
                    cell = int(mat[r, c])
                    b = cell & 0b1111
                    if (b & 0b1000) and (r > 0):
                        G_conn.add_edge((r, c), (r - 1, c))
                    if (b & 0b0100) and (r < rows - 1):
                        G_conn.add_edge((r, c), (r + 1, c))
                    if (b & 0b0010) and (c > 0):
                        G_conn.add_edge((r, c), (r, c - 1))
                    if (b & 0b0001) and (c < cols - 1):
                        G_conn.add_edge((r, c), (r, c + 1))

            endpoints_set = set(endpoints)
            while endpoints_set:
                u = endpoints_set.pop()
                queue = deque([u])
                seen = {u}
                match = None
                while queue and (match is None):
                    node = queue.popleft()
                    if node in endpoints_set:
                        match = node
                        break
                    for nbr in G_conn.neighbors(node):
                        if nbr not in seen:
                            seen.add(nbr)
                            queue.append(nbr)

                if match:
                    # unambiguous gate: (step, coordinate_u, coordinate_match, 'cx', False)
                    gates_this_step.add((step, u, match, 'cx', False))
                    endpoints_set.remove(match)

        # --- (1b) if MSB=1 anywhere, do the DFS backtracking to collect ambiguous branches ---
        if ambiguous:
            ambiguous, unknown_count = infer_ambiguity_and_count(mat)

            if not ambiguous:
                print(f"Trace {step_str}: Unambiguous, skipping.")
                continue

            print(f"\nTrace {step_str}: Ambiguous! ({unknown_count} logical-unknown cells)")


            if unknown_count % 2 != 0:
                    raise ValueError(f"Trace {step_str}: found odd number of logical‐unknown cells ({unknown_count}).")
            
            qubits = [tuple(rc) for rc in np.argwhere(mat == 0b1110000)]
            wires = {tuple(rc) for rc in np.argwhere((mat & 0b1110000) == 0b1010000)}
            if qubits and wires:
                start_dfs = time.time()
                raw = enumerate_all_connections_backtrack_dfs_parallel(mat)
                end_dfs= time.time()
                total_time_dfs += (end_dfs - start_dfs)

                dfs_times_per_step.append(end_dfs - start_dfs)

                clean = filter_same_length_duplicates(raw)

                non_overlapping = filter_non_overlapping_paths(clean)

                print(f"Found {len(non_overlapping)} distinct non-overlapping connection(s):")
                for sol_idx, sol in enumerate(non_overlapping, start=1):
                    print(f"  Solution #{sol_idx}:")
                    for (q1, q2), path in sol:
                        print(f"    {q1} ↔ {q2}, path length = {len(path)}")
                if non_overlapping:
                    ambiguous_step_solutions[step] = []
                    for sol in non_overlapping:
                        branch_pairs = []
                        for (u, v), path in sol:
                            branch_pairs.append((u, v))
                        ambiguous_step_solutions[step].append(branch_pairs)

        else:
            dfs_times_per_step.append(0.0) # no DFS performed

        # add any unambiguous gates from this step to all_gates
        if gates_this_step:
            all_gates.extend(sorted(gates_this_step))

    # ------------------------------------------------------------------------------
    # (2) Build a “gate_tuples” list that includes every unambiguous gate _and_ every
    #     ambiguous branch as its own entry.  Then sort by (step asc, ambiguous last) and
    #     reverse so that higher true_idx appears first.
    # ------------------------------------------------------------------------------
    time_dag_reconstruction = 0
    start_dag_reconstruction = time.time()

    gate_tuples = []  
    # first, all unambiguous gates
    for (step, u, v, gate_type, is_ambig) in all_gates:
        gate_tuples.append((step, u, v, gate_type, None, is_ambig))

    # then, for each ambiguous step, each branch
    for step, solutions in ambiguous_step_solutions.items():
        for branch_idx, branch_pairs in enumerate(solutions):
            for (u, v) in branch_pairs:
                gate_tuples.append((step, u, v, 'cx', branch_idx, True))

    # sort by (step ascending, then unambiguous first)
    gate_tuples.sort(key=lambda x: (x[0], 0 if not x[5] else 1))
    num_gates = len(gate_tuples)

    # assign “true_idx” by reversing that sort
    gate_tuples = gate_tuples[::-1]  # now index 0 in this list has true_idx = num_gates - 1

    # ------------------------------------------------------------------------------
    # (3) Build a helper per‐qubit mapping so we can “look forward” from each gate on each qubit.
    # ------------------------------------------------------------------------------
    # We’ll store: per_qubit[q] = [ (true_idx, gate_id, trace_id, is_ambiguous), … ] sorted ascending by true_idx
    per_qubit = defaultdict(list)

    # First create a node_info list (so we can refer to “true_idx” and “gate_id” and “step”)
    node_info = []
    for idx, (step, u, v, gate_type, branch_idx, is_ambig) in enumerate(gate_tuples):
        true_idx = num_gates - idx - 1
        gate_id  = f"G{true_idx}"
        node_info.append({
            'true_idx'  : true_idx,
            'gate_id'   : gate_id,
            'trace_id'  : step,
            'u'         : u,
            'v'         : v,
            'ambiguous' : is_ambig
        })
        # record in per_qubit
        per_qubit[u].append((true_idx, gate_id, step, is_ambig))
        per_qubit[v].append((true_idx, gate_id, step, is_ambig))

    # Sort each qubit’s list by true_idx ascending (so we can bisect to find “next higher true_idx”)
    for q in per_qubit:
        per_qubit[q].sort(key=lambda tpl: tpl[0])  # tpl[0] = true_idx

    # ------------------------------------------------------------------------------
    # (4) Build the DAG edges by “looking forward” on each qubit separately.
    #     For gate G_t on qubit q, we find the minimal true_idx'>t with a different trace_id.
    #     If that next gate is unambiguous, link G_t→G_{true_idx'}.
    #     If that next gate is ambiguous (trace_id = S), collect _all_ entries for that same step S
    #     whose true_idx>S_t (still > t) and add edges G_t→each one.
    # ------------------------------------------------------------------------------
    dag = nx.DiGraph()
    for info in node_info:
        dag.add_node(info['gate_id'],
                     control = f"Q{info['u']}",
                     target  = f"Q{info['v']}",
                     step    = info['trace_id'],
                     gate    = 'cx',
                     ambiguous = info['ambiguous'],
                     branch  = None)

    # We'll record edges_added_by_gate so we can print later
    edges_added_by_gate = {info['gate_id']: [] for info in node_info}

    # Now, for each gate in node_info (any order is fine, since we look “forward” by true_idx),
    #   but let’s iterate ascending true_idx to keep it conceptually neat:
    node_info_sorted = sorted(node_info, key=lambda x: x['true_idx'])
    for info in node_info_sorted:
        t0 = info['true_idx']
        g0 = info['gate_id']
        step0 = info['trace_id']
        q_u = info['u']
        q_v = info['v']

        def find_and_link_on_qubit(q):
            """
            Look up per_qubit[q] to find the next (true_idx'>t0) with trace_id != step0.
            If that gate is unambiguous, link exactly to it.
            If it’s ambiguous (trace_id = S), collect all entries with that same S and link to each.
            """
            lst = per_qubit[q]
            # find insertion point for (t0, …) to get index of next higher‐true_idx
            idx = bisect_right(lst, (t0, "", 0, False))
            # skip forward until we find an entry with a different trace_id
            while idx < len(lst) and lst[idx][2] == step0:
                idx += 1
            if idx >= len(lst):
                return  # no further gate on that qubit

            next_true, next_id, next_step, next_amb = lst[idx]
            if not next_amb:
                # unambiguous → link only to that one
                dag.add_edge(g0, next_id)
                edges_added_by_gate[g0].append(f"{t0}->{next_true}")
            else:
                # ambiguous: collect all entries in lst with trace_id == next_step
                # starting from idx (they may not be contiguous, but we can scan)
                for (t1, gid1, s1, amb1) in lst:
                    if t1 <= t0:
                        continue
                    if s1 == next_step:
                        dag.add_edge(g0, gid1)
                        edges_added_by_gate[g0].append(f"{t0}->{t1}")
                # once we've linked all branches of that same trace_id, stop

        # look for dependency on q_u
        find_and_link_on_qubit(q_u)
        # look for dependency on q_v
        find_and_link_on_qubit(q_v)
    end_dag_reconstruction = time.time()
    # ------------------------------------------------------------------------------
    # (5) Print the full “trace_rows” table and a mini‐table
    # ------------------------------------------------------------------------------
    # We want to print in descending true_idx order; build a list of rows accordingly:
    node_info_desc = sorted(node_info, key=lambda x: x['true_idx'], reverse=True)

    trace_rows = []
    for rev_idx, info in enumerate(node_info_desc):
        t0   = info['true_idx']
        g0   = info['gate_id']
        step = info['trace_id']
        amb  = "A" if info['ambiguous'] else "U"
        u, v = info['u'], info['v']
        targs = f"({u},{v})"
        edges_added = sorted(edges_added_by_gate[g0])
        trace_rows.append({
            'rev_idx'     : rev_idx,
            'targs'       : targs,
            'true_idx'    : t0,
            'edges_added' : ", ".join(edges_added) if edges_added else "none",
            'trace_id'    : step,
            'ambig'       : amb
        })

    # --- column widths for the full table ---
    width_rev    = 8
    width_targs  = 18
    width_true   = 14
    width_edges  = 130
    width_trace  = 20
    width_ambig  = 8

    header = (
        f"|{'rev_idx':^{width_rev}}"
        f"|{'targs':^{width_targs}}"
        f"|{'true_idx':^{width_true}}"
        f"|{'edges_added':^{width_edges}}"
        f"|{'trace_id':^{width_trace}}"
        f"|{'ambig':^{width_ambig}}|"
    )
    sep = (
        f"|{'-'*width_rev}"
        f"|{'-'*width_targs}"
        f"|{'-'*width_true}"
        f"|{'-'*width_edges}"
        f"|{'-'*width_trace}"
        f"|{'-'*width_ambig}|"
    )
    print(header)
    print(sep)
    for row in trace_rows:
        print(
            f"|{str(row['rev_idx']):^{width_rev}}"
            f"|{str(row['targs']):^{width_targs}}"
            f"|{str(row['true_idx']):^{width_true}}"
            f"|{row['edges_added']:<{width_edges}}"
            f"|{str(row['trace_id']):^{width_trace}}"
            f"|{str(row['ambig']):^{width_ambig}}|"
        )

    # --- mini‐table: only (targs / true_idx / trace_id / ambig) ---
    print()
    print("Mini-table (only targs / true_idx / trace_id / ambig):")
    w_targs = 18
    w_true  = 12
    w_trace = 12
    w_ambig = 8

    mini_header = (
        f"|{'targs':^{w_targs}}"
        f"|{'true_idx':^{w_true}}"
        f"|{'trace_id':^{w_trace}}"
        f"|{'ambig':^{w_ambig}}|"
    )
    mini_sep = (
        f"|{'-'*w_targs}"
        f"|{'-'*w_true}"
        f"|{'-'*w_trace}"
        f"|{'-'*w_ambig}|"
    )
    print(mini_header)
    print(mini_sep)
    for row in trace_rows:
        print(
            f"|{row['targs']:^{w_targs}}"
            f"|{str(row['true_idx']):^{w_true}}"
            f"|{str(row['trace_id']):^{w_trace}}"
            f"|{row['ambig']:^{w_ambig}}|"
        )

    print(f"\nTotal nodes: {dag.number_of_nodes()}")

    print(f"Time taken for DFS: {total_time_dfs:.4f} seconds")
    print(f"Time taken for DAG reconstruction: {end_dag_reconstruction - start_dag_reconstruction:.4f} seconds")

    import json
    import os

    # Ensure these are set per benchmark instance
    # benchmark_name = "add_7"
    # layout_type = "sparse"  # or "compact", "intermediate"

    # Collect timing + DAG stats
    timing_data = {
        # "benchmark_name": benchmark_name,
        "layout": layout_type,
        "dfs_time_sec": total_time_dfs,
        "dag_reconstruction_time_sec": end_dag_reconstruction - start_dag_reconstruction,
        "num_nodes": dag.number_of_nodes(),
        "num_edges": dag.number_of_edges(),
        "per_step_dfs_time_sec": dfs_times_per_step
    }

    # Dynamically generate output path
    base_dir = f"/home/george/FTQC_Security/targets_2/"
    output_dir = os.path.join(base_dir, layout_type)
    os.makedirs(output_dir, exist_ok=True)

    # File name: e.g., add_7_sparse.json
    output_path = os.path.join(output_dir, f"{benchmark_name}/{benchmark_pert}.json")
    # output_path = os.path.join(output_dir, f"{benchmark_name}.json")

    # Save JSON
    with open(output_path, "w") as f:
        json.dump(timing_data, f, indent=4)

    print(f"✅ Timing data saved to: {output_path}")

    return dag

