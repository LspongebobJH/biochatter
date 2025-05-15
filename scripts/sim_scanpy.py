from biochatter.llm_connect import GptConversation
from biochatter.api_agent.base.api_agent import APIAgent
from biochatter.api_agent.python.scanpy import ScanpyQueryBuilder, ScanpyFetcher, ScanpyInterpreter

import os
import scanpy
from scanpy.datasets import pbmc3k

scanpy.settings.datasetdir = os.environ.get("DATA")

system_prompt = """
You are a professional bioinformatician. You have access to the data object named `data`.
Please only use the provided tools. Do not use any tools that are not provided.
"""
# Create an API agent for OncoKB
query_builder_conv = GptConversation(model_name="gpt-3.5-turbo", prompts={
    "primary_model_prompts": system_prompt
})
interpreter_conv = GptConversation(model_name="gpt-3.5-turbo", prompts={
    "primary_model_prompts": system_prompt
})

scanpy_agent = APIAgent(
    query_builder=ScanpyQueryBuilder(
        conversation=query_builder_conv,
    ),
    fetcher=ScanpyFetcher(),
    interpreter=ScanpyInterpreter( # Jiahang: explain codes, args, etc. see biomania.
        conversation=interpreter_conv,
    )
)

# Execute a query
# question = "Visualize umap embedding of cells' gene counts data where cells are colored by leiden clustering with smooth coloring."
question = "Visualize umap embedding density of cells' gene counts data where cells are grouped by louvain clustering."
data = pbmc3k()
result = scanpy_agent.execute(question, data=data)

pass