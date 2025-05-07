from biochatter.api_agent.base.agent_abc import BaseAPIModel, BaseQueryBuilder, BaseFetcher, BaseInterpreter
from biochatter.llm_connect import Conversation
from biochatter.utils import get_zero_outdeg_nodes, get_zero_outdeg_nodes_if_remove_nodes
from biochatter.api_agent.dep_graph import DependencyGraph
from biochatter.api_agent.dep_graph.utils import get_topo_sort_layers, get_active_in_apis
from biochatter.api_agent.dep_graph.helplers import TopoSortLayers
from biochatter.api_agent.base.agent_abc import BaseDependency
from .meta_api import TARGET_TOOLS_DICT, TOOLS_DICT

from langchain_core.output_parsers import PydanticToolsParser
from langchain_core.language_models.chat_models import BaseChatModel
from pydantic import BaseModel
import networkx as nx
from networkx.classes.function import set_node_attributes

import json
from queue import Queue


class ScanpyQueryBuilder(BaseQueryBuilder):
    def __init__(self, 
                 conversation: Conversation,
                 dep_graph_path: str = 'biochatter/api_agent/python/scanpy/graph.json'
                 ):
        super().__init__(conversation=conversation)

        self.dep_graph = DependencyGraph(dep_graph_path=dep_graph_path, api_class_dict=TOOLS_DICT)

    def build_api_query(
        self,
        question: str,
    ) -> list[BaseModel]:

        tools = list(TARGET_TOOLS_DICT.values())

        # Jiahang: only one target API being considered for now
        # can we somehow restrict LLM to predict only one API?
        api = self._parametrise_api(question, tools)
        api: BaseModel = api[0]
        api_chain = self._trace_back(question, api)

        return api_chain
    
    def _parametrise_api(
        self,
        question: str,
        tools: list[BaseModel],
    ):
        """Parametrise the API data model using LLM."""
        llm_with_tools: BaseChatModel = self.conversation.chat.bind_tools(tools, tool_choice="required")
        parser = PydanticToolsParser(tools=tools)
        runnable = llm_with_tools | parser
        tools = runnable.invoke(question)
        return tools
    
    def _trace_back(
            self, 
            question: str,
            api: BaseModel
        ) -> list[BaseAPIModel]:
        pass
    
class ScanpyFetcher(BaseFetcher):
    def fetch_results(
        self,
        query_models: list[BaseModel],
        data: object,
        retries: int | None = 3,
    ):
        pass

class ScanpyInterpreter(BaseInterpreter):
     def summarise_results(
        self,
        question: str,
        conversation: Conversation,
        response: object,
    ) -> str:
        pass