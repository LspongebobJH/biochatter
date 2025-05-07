"""Abstract base classes for API interaction components.

Provides base classes for query builders, fetchers, and interpreters used in
API interactions and result processing.
"""

import os
from dotenv import load_dotenv
from abc import ABC, abstractmethod

from langchain_core.prompts import ChatPromptTemplate
from pydantic import BaseModel, ConfigDict, Field, create_model, PrivateAttr
from typing import Any

from biochatter.llm_connect import Conversation

class BaseQueryBuilder(ABC):
    """An abstract base class for query builders."""
    def __init__(self, conversation: Conversation):
        """Initialise the query builder with a conversation object."""
        self.conversation = conversation

        load_dotenv()
        self.conversation.set_api_key(
            os.environ.get("API_KEY")
        )

    @property
    def structured_output_prompt(self) -> ChatPromptTemplate:
        """Define a structured output prompt template.

        This provides a default implementation for an API agent that can be
        overridden by subclasses to return a ChatPromptTemplate-compatible
        object.
        """
        return ChatPromptTemplate.from_messages(
            [
                (
                    "system",
                    "You are a world class algorithm for extracting information in structured formats.",
                ),
                (
                    "human",
                    "Use the given format to extract information from the following input: {input}",
                ),
                ("human", "Tip: Make sure to answer in the correct format"),
            ],
        )

    
    @abstractmethod
    def build_api_query(
        self,
        question: str,
    ) -> list[BaseModel]:
        """Parameterise a query object.

        Parameterises a Pydantic model with the fields of the API based on the
        given question using a BioChatter conversation instance. Must be
        implemented by subclasses.

        Args:
        ----
            question (str): The question to be answered.

        Returns:
        -------
            A list containing one or more parameterised instance(s) of the query
            object (Pydantic BaseModel).

        """


class BaseFetcher(ABC):
    """Abstract base class for fetchers.

    A fetcher is responsible for submitting queries (in systems where
    submission and fetching are separate) and fetching and saving results of
    queries. It has to implement a `fetch_results()` method, which can wrap a
    multi-step procedure to submit and retrieve. Should implement retry method to
    account for connectivity issues or processing times.
    """

    @abstractmethod
    def fetch_results(
        self,
        query_models: list[BaseModel],
        data: object,
        retries: int | None = 3,
    ):
        """Fetch results by submitting a query.

        Can implement a multi-step procedure if submitting and fetching are
        distinct processes (e.g., in the case of long processing times as in the
        case of BLAST).

        Args:
        ----
            query_models: list of Pydantic models describing the parameterised
                queries

        """


class BaseInterpreter(ABC):
    """Abstract base class for result interpreters.

    The interpreter is aware of the nature and structure of the results and can
    extract and summarise information from them.
    """
    def __init__(self, conversation: Conversation):
        """Initialise the interpreter with a conversation object."""
        self.conversation = conversation

        load_dotenv()
        self.conversation.set_api_key(
            os.environ.get("API_KEY")
        )

    @abstractmethod
    def summarise_results(
        self,
        question: str,
        response: object,
    ) -> str:
        """Summarise an answer based on the given parameters.

        Args:
        ----
            question (str): The question that was asked.

            conversation_factory (Callable): A function that creates a
                BioChatter conversation.

            response (object): The response.text returned from the request.

        Returns:
        -------
            A summary of the answer.

        Todo:
        ----
            Genericise (remove file path and n_lines parameters, and use a
            generic way to get the results). The child classes should manage the
            specifics of the results.

        """


class BaseAPIModel(BaseModel):
    """A base class for all API models.

    Includes default fields `uuid` and `api_name`.
    """

    uuid: str | None = Field(
        None,
        description="Unique identifier for the model instance",
    )
    model_config = ConfigDict(arbitrary_types_allowed=True)


class BaseTools:
    """Abstract base class for tools."""

    def make_pydantic_tools(self) -> list[BaseAPIModel]:
        """Uses pydantics create_model to create a list of pydantic tools from a dictionary of parameters"""
        tools = []
        for func_name, tool_params in self.tools_params.items():
            tools.append(create_model(func_name, **tool_params, __base__=BaseAPIModel))
        return tools

class BaseObject(BaseModel):
    """A class representing an object, such as an API, dependency, data, keys_info, etc."""
    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")
    def __hash__(self):
        members = self._hash_members()
        members = tuple(f"{k}:{v}" for k, v in members.items())
        return hash(members)
    
    def _hash_members(self) -> dict:
        """A dict of members to be hased = {member_name: member_value}"""
        return self.model_dump()

class BaseKeysInfo(BaseObject):
    """A class representing a keys info object."""
    membership: str = Field(
        default="item", 
        choices=["item", "attr", "self"],
        description="The membership method to get the data of the key"
    )
    keys: dict[str, "BaseKeysInfo"] | None = Field(
        default=None,
        description="A dictionary of keys and their keys info"
    )

class BaseData(BaseObject):
    """A class representing a data object.

    Data example: # dict form of keys_info
    keys_info: {
        "membership": "self",
        "keys": {
            "layer1": {
                "membership": "item",
                "keys": {
                    "key1": {"membership": "item", "keys": None},
                    "key2": {"membership": "attr", "keys": None}
                }
            },
            "layer2": {
                "membership": "item", 
                "keys": {
                    "key3": {"membership": "item", "keys": None},
                    "key4": {"membership": "attr", "keys": None}
                }
            }
        }
    }
    layer1, layer2, keys1, keys2, ... are keys of data object
    membership is the membership method to get the data of the key
    
    for instance, given the keys_info as above:
    Then to access object of "key1", we use data[layer1][key1] or data.__getitem__(layer1).__getitem__(key1);
    To access object of "key2", we use data[layer1].key2 or data.__getitem__(layer1).__getattribute__("key2").
    For flexibility, we use membership to specify the membership method to get the data of the key rather than using 
    [] and . to access the data.
    """
    data: Any = None
    keys_info: BaseKeysInfo = Field(default=BaseKeysInfo(), description="The keys of the data object")

    def _hash_members(self) -> dict:
        members = self.model_dump()
        members.pop('data') # data can be complex structure not hashable
        return members


class BaseAPI(BaseObject):
    """Base class for all API models.
    
    We use PrivateAttr to store api_name and products to avoid them being
    included in the argument prediction of LLM through langchain.
    """

    _api_name: str = PrivateAttr(default="")
    _products: BaseData = PrivateAttr(default=BaseData())

    def _hash_members(self):
        members = self.model_dump()
        members['_api_name'] = self._api_name
        members['_products'] = self._products._hash_members()
        return members
    
    def execute(self):
        """Execute the API call with the given arguments."""
        pass


class BaseDependency(BaseObject):
    """A class representing an edge in the dependency graph.

    This class is used to represent the dependencies between API calls in the
    dependency graph.
    """
    u_api_name: str
    v_api_name: str
    args: dict
    arg_typs: dict
    deps: BaseData

    def _hash_members(self):
        members = self.model_dump()
        members['deps'] = self.deps._hash_members()
        return members
