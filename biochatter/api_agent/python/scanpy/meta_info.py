data_names = ["adata", "data"]
_dep_graph_path = "biochatter/api_agent/python/scanpy/graph.json"

import json

with open(_dep_graph_path, 'r') as f:
    dep_graph_data = json.load(f)

def preprocess_dep_graph_data(dep_graph_data: dict) -> dict:
    """Preprocess the dependency graph data."""
    _nodes = dep_graph_data["nodes"]
    nodes = []

    for node in _nodes:
        if node.get('_deprecated', False):
            continue
        else:
            nodes.append(node)

    dep_graph_data["nodes"] = nodes
    return dep_graph_data

dep_graph_data = preprocess_dep_graph_data(dep_graph_data)