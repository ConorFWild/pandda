
from bamboo.stats.cluster import convert_tree_to_edges

def plot_tree_as_network(cluster_node):
    """Takes a scipy ClusterNode object and plot as a phylogenetic tree"""

    edges = convert_tree_to_edges(cluster_node)


