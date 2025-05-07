import networkx as nx
from networkx import DiGraph
from pydantic import BaseModel, PrivateAttr, Field
import json
from typing import Any
from queue import Queue

from biochatter.api_agent.base.agent_abc import BaseAPI, BaseDependency



class DependencyGraph(DiGraph):
    """A class representing a dependency graph for API calls.

    This class extends the DiGraph class from NetworkX to represent a directed
    graph where nodes are API calls and edges represent dependencies between them.

    Not all DiGraph methods have corresponding methods in this class. Only methods used
    by other codes are implemented. The rest are inherited from DiGraph.
    """

    def __init__(self, dep_graph_path: str | None = None, api_class_dict: dict | None = None):
        super.__init__()

        self.api_dict = {}
        self.dep_dict = {}

        if dep_graph_path:
            assert api_class_dict is not None, "API class dictionary must be provided if dep_graph_path is given."
            with open(dep_graph_path, "r") as f:
                dep_graph = json.load(f)
            dep_graph = nx.node_link_graph(dep_graph, edges='edges')

            api_names = dep_graph.nodes()
            self.add_apis_from([api_class_dict[api_name]() for api_name in api_names])
            for u_api_name, v_api_name, e_data in dep_graph.edges(data=True):
                dep = BaseDependency(
                    u_api_name=u_api_name,
                    v_api_name=v_api_name,
                    args=e_data.get("args", {}),
                    arg_typs=e_data.get("arg_typs", {})
                )
                self.add_dep(dep)

    
    def add_api(self, api: BaseAPI):
        super().add_node(api._api_name)
        self.api_dict[api._api_name] = api

    def add_apis_from(self, api_list: list[BaseModel]):
        for api in api_list:
            self.add_node(api)

    def remove_api(self, api: BaseAPI):
        super().remove_node(api._api_name)
        del self.api_dict[api._api_name]

    def remove_apis_from(self, api_list: list[BaseModel]):
        for api in api_list:
            self.remove_node(api)

    def add_dep(self, dep: BaseDependency):
        u_api_name, v_api_name = dep.u_api._api_name, dep.v_api._api_name
        super().add_edge(u_api_name, v_api_name)
        self.dep_dict[(u_api_name, v_api_name)] = dep

    def add_deps_from(self, dep_list: list[BaseDependency]):
        for dep in dep_list:
            self.add_edge(dep)
        
    def remove_dep(self, u_api: BaseAPI, v_api: BaseAPI):
        super().remove_edge(u_api._api_name, v_api._api_name)
    
    def remove_deps_from(self, u_v_api_list: list[(BaseModel, BaseModel)]):
        u_v_names = [(u_api._api_name, v_api._api_name) for u_api, v_api in u_v_api_list]
        super().remove_edges_from(u_v_names)

    def in_apis(self, api: BaseAPI) -> list[BaseModel]:
        """Get the dependent APIs of the given API."""
        in_nodes = super().predecessors(api._api_name)
        return self.get_apis(list(in_nodes))
    
    def out_apis(self, api: BaseAPI) -> list[BaseModel]:
        """Get the APIs that depend on the given API."""
        out_nodes = super().successors(api._api_name)
        return self.get_apis(list(out_nodes))
    
    def in_deps(self, api: BaseAPI) -> list[BaseDependency]:
        """Get the dependencies of the given API."""
        in_edges = super().in_edges(api._api_name, data=True)
        return self.get_deps(list(in_edges))
    
    def out_deps(self, api: BaseAPI) -> list[BaseDependency]:
        """Get the dependencies that depend on the given API."""
        out_edges = super().out_edges(api._api_name, data=True)
        return self.get_deps(list(out_edges))
    
    def clear(self):
        super().clear()
        self.api_dict.clear()
        self.dep_dict.clear()

    def clear_deps(self):
        super().clear_edges()
        self.dep_dict.clear()

    def get_api(self, api_name: str) -> BaseModel:
        """Get the API object associated with the given API name."""
        return self.api_dict.get(api_name)
    
    def get_apis(self, api_names: list[str]) -> list[BaseModel]:
        """Get the API objects associated with the given API names."""
        return [self.get_api(api_name) for api_name in api_names]
    
    def get_dep(self, u_api_name: str, v_api_name: str) -> BaseDependency:
        """Get the dependency object associated with the given API names."""
        return self.dep_dict.get((u_api_name, v_api_name))
    
    def get_deps(self, u_v_api_names: list[(str, str)]) -> list[BaseDependency]:
        """Get the dependency objects associated with the given API names."""
        return [self.get_dep((u_api_name, v_api_name)) for u_api_name, v_api_name in u_v_api_names]
    
    def retrieve_sub_g(self, sub_g: DiGraph) -> DiGraph:
        """Reteieve the sub DependencyGraph given the sub latent DiGraph
        The latent DiGraph only needs the node set and egde set, where node names are indexed by API names.
        There is no need for any node or edge attributes for latent DiGraph.
        Specifically, this method retrieves sub API and BaseDependency sets.
        In principle, all DiGraph subgraph retrieval operations need to run this method afterwards.
        """
        sub_apis = self.get_apis(sub_g.nodes())
        sub_deps = self.get_deps(sub_g.edges())
        graph = DependencyGraph()
        graph.add_apis_from(sub_apis)
        graph.add_deps_from(sub_deps)
        return graph
    
    def get_zero_indeg_apis(self) -> list[BaseModel]:
        """Get apis with zero indegree in the dependency graph."""
        nodes = [node for node in self.nodes() if self.in_degree(node) == 0]
        return self.get_apis(nodes)
    
    def get_zero_outdeg_apis(self) -> list[BaseModel]:
        """Get apis with zero outdegree in the dependency graph."""
        nodes = [node for node in self.nodes() if self.out_degree(node) == 0]
        return self.get_apis(nodes)
    
    def get_actual_deps(self, api: BaseAPI) -> list[BaseDependency]:
        """Get actual dependencies of the given API."""
        deps = self.in_deps(api)
        for dep in deps:
            deps_info, deps = dep.deps_info, dep.deps
            u_api = self.get_api(dep.u_api_name)
            products = u_api._products