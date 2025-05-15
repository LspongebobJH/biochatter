from biochatter.api_agent.base.agent_abc import BaseAPI, BaseDependency
from .meta_info import dep_graph_dict
from pydantic import PrivateAttr, Field

class ScanpyAPI(BaseAPI):
    _dep_graph_dict: dict = PrivateAttr(default=dep_graph_dict)

class ScanpyDependency(BaseDependency):
    _dep_graph_dict: dict = PrivateAttr(default=dep_graph_dict)