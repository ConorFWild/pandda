
import networkx as nx

def convert_edges_to_graph(edges, label_hash={}):
    """Convert a list of labels and distances to a graph"""
    g = nx.Graph()
    for lab1, lab2, edge_wt in edges: g.add_edge(label_hash.get(lab1,lab1), label_hash.get(lab2,lab2), weight=edge_wt)
    return g

def convert_tree_to_edges(cluster_node):
    """Takes a scipy ClusterNode object and returns network-like edges for the object"""
    # Top of the tree
    top = cluster_node
    # Nodes yet to process
    next_nodes = [cluster_node]
    # List of connections
    connections = []
    # Iterate through nodes and create edges to their children
    while True:
        # Get the next node
        node = next_nodes.pop(0)
        # Get attached nodes
        l = node.get_left()
        r = node.get_right()
        # Iterate through attached nodes
        for child in [l,r]:
            # Append to nodes if not leaves
            if not child.is_leaf(): next_nodes.append(child)
            # Add connections and distance
            connections.append((node.get_id(), child.get_id(), node.dist-child.dist))
        # Exit if not more nodes
        if not next_nodes: break
    return connections
