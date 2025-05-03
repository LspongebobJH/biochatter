from biochatter.api_agent.base.agent_abc import BaseAPIModel, BaseQueryBuilder, BaseFetcher, BaseInterpreter
from biochatter.llm_connect import Conversation
from .meta_api import TARGET_TOOLS_DICT, TOOLS_DICT

from langchain_core.output_parsers import PydanticToolsParser
from langchain_core.language_models.chat_models import BaseChatModel
from pydantic import BaseModel

import json
import networkx as nx

class ScanpyQueryBuilder(BaseQueryBuilder):
    def __init__(self, 
                 conversation: Conversation,
                 dep_graph_path: str = 'biochatter/api_agent/python/scanpy/graph.json',
                 meta_info_path: str = 'biochatter/api_agent/python/scanpy/meta_info.json'):
        super().__init__(conversation=conversation)
        
        with open(dep_graph_path, "r") as f:
            dep_graph = json.load(f)
        with open(meta_info_path, "r") as f:
            meta_info = json.load(f)
        self.dep_graph = nx.node_link_graph(dep_graph, edges='edges')
        self.data_names = meta_info['data_names']

    def build_api_query(
        self,
        question: str,
    ) -> list[BaseModel]:

        tools = list(TARGET_TOOLS_DICT.values())
        api = self._parametrise_api(question, tools)
        # Jiahang: only one target API being considered for now
        # can we somehow restrict LLM to predict only one API?
        api: BaseModel = api[0] 

        api_subg = nx.DiGraph()
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
        """Trace back the dependency graph to find dependent API."""
        api_name = api._api_name
        if api_name == 'root':
            return api
        
        in_edges = self.dep_graph.in_edges(api_name, data=True)

        api_list = []
        dep_list = []
        for src, dst, e_data in in_edges:
            e_args = e_data["args"]
            if self._check_e_arg(e_args, api):
                src_api = self.dep_graph.nodes[src]['api']
                src_api = TOOLS_DICT[src_api]
                src_api = self._parametrise_api(question, [src_api])[0]
                api_list.append(src_api)
                dep_list.append()
            return [src_api] + self._trace_back(question, src_api)

    def _check_e_arg(self, e_args: dict, api: BaseModel) -> bool:
        """Check if the edge required argument is activated by the target api."""
        e_arg_name = list(e_args.keys())[0] # Jiahang: only one edge required arg being considered for now
        e_arg_val = e_args[e_arg_name]
        if e_arg_name in api.model_fields.keys() and getattr(api, e_arg_name) == e_arg_val:
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