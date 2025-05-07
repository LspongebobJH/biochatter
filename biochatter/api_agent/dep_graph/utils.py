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

def get_active_in_apis(deps: list[BaseDependency], api: BaseModel) -> list[BaseDependency]:
        """Get the active dependencies of the given API."""
        active_deps = []
        for dep in deps:
            if check_active_dep_args(dep.args, api):
                active_deps.append(dep)
        return [dep.u_api for dep in active_deps]

def check_active_dep_args(dep_args: dict, api: BaseModel) -> bool:
    """Check if the edge required argument is activated by the target api."""
    e_arg_name = list(dep_args.keys())[0] # Jiahang: only one edge required arg being considered for now
    e_arg_val = dep_args[e_arg_name]
    if e_arg_name in api.model_fields.keys() and getattr(api, e_arg_name) == e_arg_val:
        return True
    return False

