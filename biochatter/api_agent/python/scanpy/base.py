from biochatter.api_agent.base.agent_abc import BaseAPI, InputAPI, BaseDependency, InputDependency
from typing import Any
from biochatter.api_agent.dep_graph.utils import _str_list_to_keys_info
from .meta_info import dep_graph_data

class ScanpyAPI(BaseAPI):
    def model_post_init(self, __context: Any):
        # initialize api products keys_info
        # Jiahang: this may be duplicate with read_apis_from_graph_dict.
        # Jiahang: the best way to input keys_info into API should be
        # 1. precompile keys_info into API class
        # 2. initialize keys_info from external data file (this implementation)
        # 3. explicitly pass keys_info to API instance after initialization
        # (in graph.py implementation), which is actually incorrect because
        # each time when API class is used by langchain to parametrise API,
        # the keys_info will be re-initialized with None and thus will be lost.
        # In the future, we should first remove implementation 3.
        # Note that the precompilation is recommened since keys_info is static and
        # bound to the API class.
        #
        # Jiahang: then.. how about Dependency?
        
        node_idx = dep_graph_data['node_index'][self._api_name]
        node = dep_graph_data['nodes'][node_idx]
        input_api = InputAPI.model_validate(node)
        self._products.keys_info = _str_list_to_keys_info(input_api.products)


class ScanpyDependency(BaseDependency):
    def model_post_init(self, __context: Any):
        edge_idx = dep_graph_data['edge_index'][f"{self.u_api_name}:{self.v_api_name}"]
        edge = dep_graph_data['edges'][edge_idx]
        input_dep = InputDependency.model_validate(edge)
        self.deps.keys_info = _str_list_to_keys_info(input_dep.dependencies)
        self.args = input_dep.args
        self.arg_types = input_dep.arg_types
