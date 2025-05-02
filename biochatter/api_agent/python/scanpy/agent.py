from biochatter.api_agent.base.agent_abc import BaseAPIModel, BaseQueryBuilder, BaseFetcher, BaseInterpreter
from biochatter.api_agent.python.scanpy.pl import TOOLS_DICT
from biochatter.llm_connect import Conversation

from langchain_core.output_parsers import PydanticToolsParser
from langchain_core.language_models.chat_models import BaseChatModel
from pydantic import BaseModel

import json
import networkx as nx

class ScanpyQueryBuilder(BaseQueryBuilder):
    def __init__(self, 
                 dep_graph_path: str = 'biochatter/api_agent/python/scanpy/graph.json'):
        super().__init__()
        
        with open(dep_graph_path, "r") as f:
            dep_graph = json.load(f)
        self.dep_graph = nx.node_link_graph(dep_graph, edges='edges')
        self.data_names = dep_graph['data_names']

    def parameterise_query(
        self,
        question: str,
        conversation: Conversation,
    ) -> list[BaseModel]:

        tools = list(TOOLS_DICT.values())
        llm_with_tools: BaseChatModel = conversation.chat.bind_tools(tools)
        parser = PydanticToolsParser(tools=tools)
        runnable = llm_with_tools | parser
        api = runnable.invoke(question)
        api: BaseModel = api[0] # only target API

        api_chain = self._trace_back(api)
        return api_chain
    
    def _trace_back(self, api: BaseModel) -> list[BaseAPIModel]:
        """Trace back the dependency graph to find dependent API."""
        if api.api_name == 'root':
            return api
        
        api_name = api.api_name
        in_edges = self.dep_graph.in_edges(api_name, data=True)

        for src, dst, e_data in in_edges:
            e_arg = e_data.get('arg', None)
            if self._check_e_arg(e_arg, api):
                src_api = self.dep_graph.nodes[src]['api']
                src_api = TOOLS_DICT[src_api]
                return [src_api] + self._trace_back(src_api)

    def _check_e_arg(self, e_arg: dict, api: BaseModel) -> bool:
        """Check if the edge required argument is activated by the target api."""
        if e_arg is None:
            return False
        e_arg_name = list(e_arg.keys())[0]
        e_arg_val = e_arg[e_arg_name]
        # if the edge required arg is data itself, then this dependency is a must.
        if e_arg_name in self.data_names:
            return True
        if e_arg_name in api.model_fields_set and api.model_fields[e_arg_name] == e_arg_val:
            return True
        return False
    
                


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