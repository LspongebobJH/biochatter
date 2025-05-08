from .graph import DependencyGraph, BaseDependency
from .helplers import TopoSortLayers

from pydantic import BaseModel

def get_topo_sort_layers(G: DependencyGraph, reverse=False) -> list[list[BaseModel]]:
    """Topological sort each layer of the dependency graph.
    If reverse is True, the topological sort is done in reverse order.
    """
    if reverse:
        get_zero_in_or_outdeg_apis = G.get_zero_outdeg_apis
    else:
        get_zero_in_or_outdeg_apis = G.get_zero_indeg_apis
    all_layers = []
    while G.number_of_nodes() > 0:    
        all_layers.append(get_zero_in_or_outdeg_apis())
        G.remove_apis_from(all_layers[-1])

    if reverse:
        all_layers.reverse()
    return TopoSortLayers(all_layers)

