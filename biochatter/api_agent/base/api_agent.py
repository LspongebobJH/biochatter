"""Base API agent module."""

from typing import Any

from biochatter.llm_connect import Conversation

from .agent_abc import BaseFetcher, BaseInterpreter, BaseQueryBuilder

### logic

# 1. User asks question req. API to find the answer
# 2. rag_agent receives the question and sends it to the api_agent // api_agent is manually selected
# 3. api_agent writes query for API specific to the question
# 3.1 question + API prompt template + BlastQuery are input into
# langchain.chains.openai_functions.create_structured_output_runnable -> returns a structured output which is our API call object
# 4. API call
# 4.1 API call is made using the structured output from 3.1
# 4.2 API returns a response which is saved
#
# 5. Read response from and uses it to answer question
# 6. answer is returned to rag_agent


## Agent class
class APIAgent:
    def __init__(
        self,
        conversation: Conversation,
        query_builder: BaseQueryBuilder,
        fetcher: BaseFetcher,
        interpreter: BaseInterpreter,
    ):
        """API agent class to interact with a tool's API for querying and fetching
        results.  The query fields have to be defined in a Pydantic model
        (`BaseModel`) and used (i.e., parameterised by the LLM) in the query
        builder. Specific API agents are defined in submodules of this directory
        (`api_agent`). The agent's logic is implemented in the `execute` method.

        Attributes
        ----------
            conversation (Conversation): BioChatter conversation

            query_builder (BaseQueryBuilder): An instance of a child of the
                BaseQueryBuilder class.

            result_fetcher (BaseFetcher): An instance of a child of the
                BaseFetcher class.

            result_interpreter (BaseInterpreter): An instance of a child of the
                BaseInterpreter class.

        """
        self.conversation = conversation
        self.query_builder = query_builder
        self.fetcher = fetcher
        self.interpreter = interpreter
        self.final_answer = None

    def execute(self, question: str, data: Any) -> str | None:
        """Wrapper that uses class methods to execute the API agent logic. Consists
        of 1) query generation, 2) query submission, 3) results fetching, and
        4) answer extraction. The final answer is stored in the final_answer
        attribute.

        Args:
        ----
            question (str): The question to be answered.

        """
        # Generate query
        try:
            query_models = self.query_builder.parameterise_query(question, self.conversation)
        except Exception as e:
            raise Exception(f"Failed to generate query: {e}")

        # Fetch results
        try:
            response = self.fetcher.fetch_results(query_models, data, 100)
        except Exception as e:
            raise Exception(f"Failed to fetch results: {e}")

        # Extract answer from results
        try:
            final_answer = self.interpreter.summarise_results(
                question=question,
                conversation=self.conversation,
                response=response,
            )
        except Exception as e:
            raise Exception(f"Failed to extract answer from results: {e}")

        self.final_answer = final_answer
        return final_answer

    def get_description(self, tool_name: str, tool_desc: str):
        return f"This API agent interacts with {tool_name}'s API for querying and fetching results. {tool_desc}"
