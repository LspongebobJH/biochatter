from biochatter.api_agent.base.agent_abc import BaseAPI
from pydantic import PrivateAttr
from typing import Any
import json
from biochatter.api_agent.dep_graph.utils import _str_list_to_keys_info

class ScanpyAPI(BaseAPI):
    _graph_data_path: str = PrivateAttr(default='biochatter/api_agent/python/scanpy/graph_test.json')

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
        with open(self._graph_data_path, 'r') as f:
            graph_data = json.load(f)
        node_idx = graph_data['node_index'][self._api_name]
        node = graph_data['nodes'][node_idx]
        self._products.keys_info = _str_list_to_keys_info(node['products'])

