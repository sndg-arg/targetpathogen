
import pickle
import sys
import os.path
from networkx.algorithms.centrality.betweenness import betweenness_centrality
from networkx.algorithms.components.connected import connected_components
import matplotlib
import networkx as nx
def centrality(digraph_file_path):
    with open(digraph_file_path, 'rb') as file:
        digraph_data = pickle.load(file)
    base_name = os.path.basename(digraph_file_path)  # Get the base name of the file
    name_without_ext = os.path.splitext(base_name)[0]  # Remove the extension
    output_path = os.path.join(os.path.dirname(digraph_file_path), f"{name_without_ext}.gml")
    nx.write_gml(digraph_data, output_path)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py DiGraph.digraph")
    else:
        digraph_file_path = sys.argv[1]
        centrality(digraph_file_path)
