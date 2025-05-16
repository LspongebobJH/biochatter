from biochatter.llm_connect import GptConversation
from biochatter.api_agent import APIAgent, ScanpyQueryBuilder, ScanpyFetcher, ScanpyInterpreter

import json
import os
import scanpy
from scanpy.datasets import pbmc3k

scanpy.settings.datasetdir = os.environ.get("DATA", "data")

system_prompt = """
You are a professional bioinformatician. You have access to the data object named `data`.
Please only use the provided tools. Do not use any tools that are not provided.
"""
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

with open('experiment/scanpy/cases.json', 'r') as f:
    cases = json.load(f)

# Execute a query
for question in cases:
    data = pbmc3k()
    result = scanpy_agent.execute(question, data=data)

    pass

