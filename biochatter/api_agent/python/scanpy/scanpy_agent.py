from biochatter.api_agent.base.agent_abc import BaseAPIModel, BaseQueryBuilder, BaseFetcher, BaseInterpreter
from biochatter.api_agent.python.scanpy.scanpy_pl_full import TOOLS
from biochatter.llm_connect import Conversation

from langchain_core.output_parsers import PydanticToolsParser
from langchain_core.language_models.chat_models import BaseChatModel
from pydantic import BaseModel

class ScanpyQueryBuilder(BaseQueryBuilder):
    def parameterise_query(
        self,
        question: str,
        conversation: Conversation,
    ) -> list[BaseModel]:

        llm_with_tools: BaseChatModel = conversation.chat.bind_tools(TOOLS)
        parser = PydanticToolsParser(tools=TOOLS)
        api_query = llm_with_tools.invoke(question)
        api_query = parser.invoke(api_query)

        return api_query

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