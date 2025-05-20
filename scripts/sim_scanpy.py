from biochatter.llm_connect import GptConversation
from biochatter.api_agent import APIAgent, ScanpyQueryBuilder, ScanpyFetcher, ScanpyInterpreter

import os
import scanpy
from scanpy.datasets import pbmc3k, krumsiek11

scanpy.settings.datasetdir = os.environ.get("DATA")

system_prompt = """
You are a professional bioinformatician. 
1. You have access to the data object named `data`.
2. Please only use the provided tools. Do not use any tools that are not provided.
3. n_comps = 15
"""
# Create an API agent for OncoKB
query_builder_conv = GptConversation(
    model_name="gpt-3.5-turbo", 
    prompts={
        "primary_model_prompts": system_prompt
    }
)
interpreter_conv = GptConversation(
    model_name="gpt-3.5-turbo", 
    prompts={
        "primary_model_prompts": system_prompt
    }
)

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
question = "visualize stacked violin plot of gene expressions of genes Gata2, Gata1 and Fog1, where cells are clustered by leiden algorithm."
data = krumsiek11()
# data = pbmc3k()
result = scanpy_agent.execute(question, data=data)

pass