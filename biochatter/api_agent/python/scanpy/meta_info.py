"""
This file import all external data required by the current package.
"""

data_names = ["adata", "data"]
_dep_graph_path = "biochatter/api_agent/python/scanpy/graph.json"

import json

with open(_dep_graph_path, 'r') as f:
    dep_graph_dict = json.load(f)

def preprocess_dep_graph_dict(dep_graph_dict: dict) -> dict:
    """Preprocess the dependency graph data."""
    _nodes = dep_graph_dict["nodes"]
    nodes = []
    _edges = dep_graph_dict["edges"]
    edges = []

    for node in _nodes:
        if node.get('_deprecated', False):
            continue
        else:
            nodes.append(node)

    for edge in _edges:
        if edge.get('_deprecated', False):
            continue
        else:
            edges.append(edge)

    dep_graph_dict["nodes"] = nodes
    dep_graph_dict["edges"] = edges
    dep_graph_dict["node_index"] = {node["api"]: i for i, node in enumerate(nodes)}
    dep_graph_dict["edge_index"] = {f"{edge['source']}:{edge['target']}": i for i, edge in enumerate(edges)}

    return dep_graph_dict

# Jiahang: two dependencies representations different, one is a:b, another is (a, b), should be unified.
dep_graph_dict = preprocess_dep_graph_dict(dep_graph_dict)
api_names = [node["api"] for node in dep_graph_dict["nodes"]]
dependencies = [(edge["source"], edge["target"]) for edge in dep_graph_dict["edges"]]