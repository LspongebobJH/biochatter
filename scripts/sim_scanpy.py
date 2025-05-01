from biochatter.llm_connect import GptConversation
from biochatter.api_agent.base.api_agent import APIAgent
from biochatter.api_agent.python.scanpy import ScanpyQueryBuilder, ScanpyFetcher, ScanpyInterpreter

# Create an API agent for OncoKB
scanpy_agent = APIAgent(
    conversation_factory=GptConversation(model_name="gpt-4", prompts={}, correct=False),
    query_builder=ScanpyQueryBuilder(),
    fetcher=ScanpyFetcher(),
    interpreter=ScanpyInterpreter()
)

# Execute a query
question = "Visualize umap embedding of cells' gene counts data."
result = scanpy_agent.execute(question)

print(result)