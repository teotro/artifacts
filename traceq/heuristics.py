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

from scipy.signal import convolve2d

def repr_6bit_matrix(matrix):
    # print matrix where all elements are 6bit
    
    return "\n".join(" ".join(f"{x:06b}" for x in row) for row in matrix)
    

def find_l1_logicals(l1_trace):
    """
    Find L1 logicals in a given L1 trace
    
    Args:
    l1_trace: np.array
        A 2D array representing the L1 trace
        
    Returns:
    l1_logical: list of idxs
    """
    
    # find l1 logicals by identifying patches that have only one neighboring patch that is busy
    
    # we proceed in two steps: first, apply a convolution to find the neighborsum of all patches
    # then, filter this on patches that are themselves busy
    # this satisfies our requirements
    
    filter = np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])
    
    neighborsum = convolve2d(l1_trace, filter, mode='same', boundary="fill", fillvalue=0)
    
    l1_logical = l1_trace.copy()
    masked = l1_logical * (neighborsum == 1)
    
    return [tuple(x) for x in np.argwhere(masked == 1)]


def get_connected_components(l1_trace):
    # given an l1 trace, produce a list of connected components
    # this is actually quite straightforward using networkx
    dx, dy = l1_trace.shape
    
    G = nx.grid_2d_graph(dx, dy)
    
    # now, remove all nodes that are not 1 in the trace
    unneeded_nodes = []
    for node in G.nodes:
        if l1_trace[node] == 0:
            unneeded_nodes.append(node)
    
    G.remove_nodes_from(unneeded_nodes)
    return list(nx.connected_components(G))

def _construct_constricted_graph(connected_component, logical_qubits):
    """Creates a special constricted graph to facilitate analysis
    Given a connected component G and logical_qubits V', we want to create a special graph G' such that 
    G' = (V', E') where E' = {(u, v) | u, v connected in G and the simple path that connects u, v does not contain any nodes in V'}


    Args:
        connected_component (Set): Connected component of qubits
        logical_qubits (Set): Set of logical qubits

    Returns:
        nx.Graph: Graph with nodes that are just the logical qubits.
    """
    
    # create a graph with logical_qubits as nodes
    relevant_qubits = logical_qubits.intersection(set(connected_component.nodes))
    
    out_G = nx.Graph()
    out_G.add_nodes_from(relevant_qubits)
    
    # the test is the following: if we remove all logical_qubits
    # except v1, v2, are they still connected?
    for v1 in relevant_qubits:
        for v2 in relevant_qubits:
            if v1 == v2:
                continue
            
            modified_G = connected_component.copy()
            modified_G.remove_nodes_from(relevant_qubits - {v1, v2})
            
            if nx.has_path(modified_G, v1, v2):
                # connect them in out_G
                out_G.add_edge(v1, v2)
    
    return out_G

def _find_degree_1_node(G):
    for node in G.nodes:
        if G.degree[node] == 1:
            return node
    
    return None


def connected_component_to_paths(connected_component, logical_qubits):
    """
    Taking in a connected component, attempt to find the paths that connect the logical qubits.
    
    Create a graph of the connected components, start at a node with degree 1 (if there is not one, let's fail for now). Then, find the closest node via BFS. This is the end of the path and remove the edge. Continue until all edges are removed.

    Args:
        connected_component (Set): Set of tuples representing the connected component
        logical_qubits (List): List of tuples representing the coordinates where there are logical qubits

    Returns:
        set, list of sets of tuples: Connected component, list of paths
    """
    
    
    dim = max(max(x, y) for x, y in connected_component) + 1
    G = nx.grid_2d_graph(dim, dim)
    G.remove_nodes_from(G.nodes - set(connected_component))
    
    # now, we need to find the starting node
    labeled = []
    
    constricted_G = _construct_constricted_graph(G, logical_qubits)
    
    while True:
        # we want to remove the node with degree 1 and its neighbor
        start_node = _find_degree_1_node(constricted_G)
        
        # if there is no start node, we are done
        # either because there are no more nodes or because the remaining connected component is ambiguous
        if start_node is None:
            return set(G.nodes), labeled
        
        # otherwise, remove the degree 1 node and its neighbor
        assert len(list(constricted_G.neighbors(start_node))) == 1
        end_node = list(constricted_G.neighbors(start_node))[0]
        
        assert start_node != end_node
        
        # remove both nodes from the constricted graph
        constricted_G.remove_node(start_node)
        constricted_G.remove_node(end_node)
        
        # remove the path from the original graph
        path = set(nx.shortest_path(G, start_node, end_node))
        G.remove_nodes_from(path)
        
        labeled.append(path)
        


def _free_busy_to_bitstring(l1_trace):
    """Given a free/busy map, label an edge as active if both nodes are busy.

    Args:
        l1_trace (Matrix): 2D array

    Returns:
        Matrix: 2D array
    """
    
    filter = np.array(
        [[0, 8, 0],
         [2, 0, 1],
         [0, 4, 0]]
    )
    
    # flip filter 
    filter = np.flip(filter)
    
    return convolve2d(l1_trace, filter, mode='same', boundary="fill", fillvalue=0) * l1_trace
    


def _l1_traces_to_logical_qubits(l1_traces, expected_num_of_qubits=None):
    # phase I: get logical qubits
    logical_qubits = set()

    for trace in l1_traces.values():
        logical_qubits.update(find_l1_logicals(trace["trace"]))

    # check to see if this is enough logicals
    if expected_num_of_qubits:
        if len(logical_qubits) != expected_num_of_qubits:
            raise ValueError(f"Expected {expected_num_of_qubits} qubits, got {len(logical_qubits)}")
    
    return logical_qubits


def _get_known_unknown_paths(l1_trace, logical_qubits):

    # has list of sets of tuples

    # we will rebuild the connected components into paths
    # identify the easy paths (have only two logicals) and the hard paths (have more than two logicals)
    connected_components = get_connected_components(l1_trace)
    
    known = []
    unknown = []
    
    for component in connected_components:
        remainder, solved = connected_component_to_paths(component, logical_qubits)
        
        known.extend(solved)
        unknown.append(remainder)
    
    
    return known, unknown




def l1_traces_to_l3(l1_traces, expected_num_of_qubits=None):
    assert len(l1_traces) > 0, "No L1 traces given"
    
    reconstructed_traces = OrderedDict()
    dx, dy = l1_traces[0]["trace"].shape
    
    logical_qubits = _l1_traces_to_logical_qubits(l1_traces, expected_num_of_qubits)
        
    for idx, trace in l1_traces.items():
        # we need to reconstruct both the paths and the matrix free/busy
        # begin by identifying the known and unknown paths
        known_paths, unknown_paths = _get_known_unknown_paths(trace["trace"], logical_qubits)
        
        # label on the out_trace matrix:
        out_trace = np.zeros((dx, dy))
        
        # first, label all active qubits
        for idxs in np.argwhere(trace["trace"] == 1):
            x, y = idxs
            
            if (x, y) in logical_qubits:
                out_trace[x, y] = 0b110000 # we don't know anything beyond that it is a logical qubit
            else:
                out_trace[x, y] = 0b010000
        
        # for known paths, we can label the paths exactly
        # add a 1 for every node in the matrix that is busy
        # then perform the freebusy trick
        # print("Known Paths: ", known_paths)

        
        for path in known_paths:
            assert len(path) > 1, known_paths
            known_paths_freebusy = np.zeros((dx, dy))
            
            for path_node in path:
                x, y = path_node
                known_paths_freebusy[x, y] = 1
        
            # need to filter the freebusy trick by active nodes
            out_trace += _free_busy_to_bitstring(known_paths_freebusy) * trace["trace"]
        ##############
        
        # for the unknown paths, add a special flag
        base_unknown_trace = np.zeros((dx, dy))
        # print("Unknown Paths: ", unknown_paths)
        for path in unknown_paths:
            for path_node in path:
                x, y = path_node
                base_unknown_trace[x, y] = 0b1000000
        
        out_trace += base_unknown_trace

        ##############
        # # for the unknown paths, add a special flag—but only if
        # # that path actually contains ≥2 logical qubits
        # base_unknown_trace = np.zeros((dx, dy))
        # for path in unknown_paths:
        #     # find which logicals lie in this remainder
        #     logicals_in_path = [node for node in path if node in logical_qubits]
        #     if len(logicals_in_path) < 2:
        #         continue      # skip singletons or empty remainders
        #     for (x, y) in logicals_in_path:
        #         base_unknown_trace[x, y] = 0b1000000

        # out_trace += base_unknown_trace

        reconstructed_traces[idx] = {"trace": np.array(out_trace, dtype=int), "paths": (known_paths, unknown_paths)}
        
    return reconstructed_traces
