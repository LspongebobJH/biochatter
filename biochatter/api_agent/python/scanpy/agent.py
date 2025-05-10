from biochatter.api_agent.base.agent_abc import BaseQueryBuilder, BaseFetcher, BaseInterpreter
from biochatter.llm_connect import Conversation
from biochatter.api_agent.dep_graph import DependencyGraph
from biochatter.api_agent.dep_graph.utils import is_active_dep, retrieve_products, aggregate_deps
from biochatter.api_agent.base.agent_abc import BaseAPI
from .meta_api import TARGET_TOOLS_DICT, TOOLS_DICT

from langchain_core.output_parsers import PydanticToolsParser
from langchain_core.language_models.chat_models import BaseChatModel
from pydantic import BaseModel
import networkx as nx

import json
from queue import Queue
import scanpy
from typing import Any

class ScanpyQueryBuilder(BaseQueryBuilder):
    def __init__(self, 
                 conversation: Conversation,
                 dep_graph_path: str = 'biochatter/api_agent/python/scanpy/graph.json'
                 ):
        super().__init__(conversation=conversation)

        with open(dep_graph_path, 'r') as f:
            dep_graph = json.load(f)
        self.dep_graph = DependencyGraph(dep_graph=dep_graph, api_class_dict=TOOLS_DICT)

    def build_api_query(
        self,
        question: str,
    ) -> list[BaseModel]:

        tools = list(TARGET_TOOLS_DICT.values())

        # Jiahang: only one target API being considered for now
        # can we somehow restrict LLM to predict only one API?
        api = self._parametrise_api(question, tools)
        api: BaseModel = api[0]
        execution_graph = self._trace_back(question, api)

        return [execution_graph]
    
    def _parametrise_api(
        self,
        question: str,
        tools: list[BaseAPI | type],
    ):
        """Parametrise the API data model using LLM."""
        if isinstance(tools[0], BaseAPI):
            tools = [tool.__class__ for tool in tools]
        llm_with_tools: BaseChatModel = self.conversation.chat.bind_tools(tools, tool_choice="required")
        parser = PydanticToolsParser(tools=tools)
        runnable = llm_with_tools | parser
        tools = runnable.invoke(question)
        return tools
    
    def _trace_back(
            self, 
            question: str,
            api: BaseAPI
        ) -> list[BaseAPI]:
        execution_graph = DependencyGraph()
        execution_graph.add_api(api)
        next_api_queue = Queue()
        next_api_queue.put(api)

        while not next_api_queue.empty():
            api = next_api_queue.get()
            for in_dep in self.dep_graph.in_deps(api._api_name):
                if is_active_dep(in_dep, api):
                    # Jiahang: note that the api (node) names of execution graph is a primary key.
                    # It is worth considering whether multiple predecessors of the same API could be added.
                    # That means, the same API with different arguments may occur in the final API chain.
                    # For now, we only allow a single instance of each API.
                    if in_dep.u_api_name not in execution_graph.nodes:
                        active_predecessor = self.dep_graph.get_api(in_dep.u_api_name)
                        active_predecessor = self._parametrise_api(question, [active_predecessor])[0]
                        execution_graph.add_api(active_predecessor)
                        next_api_queue.put(active_predecessor)
                    execution_graph.add_dep(in_dep)
        return execution_graph
    
class ScanpyFetcher(BaseFetcher):
    def fetch_results(
        self,
        execution_graph: list[DependencyGraph], # Jiahang: we pass a list to follow the interface. Bad practice.
        data: object,
        retries: int | None = 3,
    ) -> object:
        execution_graph: DependencyGraph = execution_graph[0]
        root = execution_graph.get_api("root")
        root._deps.data = data
        execution_graph.update_api(root)

        # topological sort
        topo_sort = list(nx.topological_sort(execution_graph))
        topo_sort = execution_graph.get_apis(topo_sort)
        for api in topo_sort:
            in_deps = execution_graph.in_deps(api._api_name)
            if len(in_deps) > 0:
                api = aggregate_deps(in_deps, api)
            api.execute(state={'sc': scanpy})
            execution_graph.update_api(api)
            out_deps = execution_graph.out_deps(api._api_name)
            for out_dep in out_deps:
                out_dep = retrieve_products(api, out_dep)
                execution_graph.update_dep(out_dep)
        return api._products.data
    
class ScanpyInterpreter(BaseInterpreter):
    def summarise_results(
        self,
        question: str,
        response: object,
    ) -> str:
        # Jiahang: no need to summarise the results for Scanpy.
        return response