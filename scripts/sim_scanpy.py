from biochatter.llm_connect import GptConversation
from biochatter.api_agent.base.api_agent import APIAgent
from biochatter.api_agent.python.scanpy import ScanpyQueryBuilder, ScanpyFetcher, ScanpyInterpreter

import os
import scanpy
from scanpy.datasets import pbmc3k
from dotenv import load_dotenv
load_dotenv()

scanpy.settings.datasetdir = os.environ.get("DATA")

# Create an API agent for OncoKB
conversation = GptConversation(model_name="gpt-4", prompts={})
conversation.set_api_key(os.environ.get("API_KEY"))
scanpy_agent = APIAgent(
    conversation=conversation,
    query_builder=ScanpyQueryBuilder(),
    fetcher=ScanpyFetcher(),
    interpreter=ScanpyInterpreter()
)

# Execute a query
question = "Visualize umap embedding of cells' gene counts data."
data = pbmc3k()
result = scanpy_agent.execute(question, data=data)

print(result)