from biochatter.llm_connect import GptConversation
from biochatter.api_agent.base.api_agent import APIAgent
from biochatter.api_agent.python.scanpy import ScanpyQueryBuilder, ScanpyFetcher, ScanpyInterpreter

import os
import scanpy
from scanpy.datasets import pbmc3k

scanpy.settings.datasetdir = os.environ.get("DATA")

# Create an API agent for OncoKB
query_builder_conv = GptConversation(model_name="gpt-4", prompts={})
interpreter_conv = GptConversation(model_name="gpt-4", prompts={})

scanpy_agent = APIAgent(
    query_builder=ScanpyQueryBuilder(
        conversation=query_builder_conv,
        dep_graph_path = 'biochatter/api_agent/python/scanpy/graph_test.json'),
    fetcher=ScanpyFetcher(),
    interpreter=ScanpyInterpreter(
        conversation=interpreter_conv,
    )
)

# Execute a query
question = "Visualize umap embedding of cells' gene counts data where cells are colored by leiden clustering."
data = pbmc3k()
result = scanpy_agent.execute(question, data=data)

print(result)