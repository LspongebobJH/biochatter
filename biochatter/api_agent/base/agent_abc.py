"""Abstract base classes for API interaction components.

Provides base classes for query builders, fetchers, and interpreters used in
API interactions and result processing.
"""

import os
from dotenv import load_dotenv
from abc import ABC, abstractmethod
import ast
import json
from typing import Any
from copy import deepcopy
from langchain_core.prompts import ChatPromptTemplate
from pydantic import BaseModel, ConfigDict, Field, create_model, PrivateAttr, field_validator, model_validator
from pydantic.fields import FieldInfo
from biochatter.llm_connect import Conversation
from ._python_interpreter import evaluate_python_code

def run_codes(code: str, state: dict[str, object]):
    """
    Run codes

    Parameters
    -----------
    code : str
        A single valid code snippet as a string.
    state: dict[str, object]
        A dictionary of variables to be used in the code snippet. E.g. {'sc': sc, 'adata': adata}
    """
    
    try:
        result = str(evaluate_python_code(code, state=state)[0])
    except Exception as e:
        return f"ERROR: {str(e)}", e
    return result, None


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
        # members = tuple(f"{k}:{v}" for k, v in members.items())
        members = json.dumps(members, sort_keys=True, ensure_ascii=True)
        return hash(members)
    
    def _hash_members(self) -> dict:
        """A dict of members to be hased = {member_name: member_value}"""
        return self.model_dump()
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

class BaseKeysInfo(BaseObject):
    """A class representing a keys info object."""
    membership: str = Field(
        default="self", 
        choices=["item", "attr", "self"],
        description="The membership method to get the data of the key"
    )
    keys: dict[str, "BaseKeysInfo"] = Field(
        default={},
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
        members.pop('data') # Jiahang: data can be complex structure not hashable, so we don't consider it for now.
        return members


class BaseAPI(BaseObject):
    """Base class for all API models.
    
    We use PrivateAttr to store api_name and products to avoid them being
    included in the argument prediction of LLM through langchain.

    Jiahang: in these classes, some data members should not be set in initialization.
    Instead, they should only be set dynamically during forward pass over execution graph,
    which is conducted internally. Please implement relevant validator to ensure this.

    Jiahang: we predefine the data name to be "data" for all APIs. This is because we assume
    all input data should be stored in one data object. Even if there are multiple data objects, they
    can all be stored in the same data object through ways like dict or list. This also means that 
    the variable name "data" should not be overlapped. This is one of the standards.
    """

    _api_name: str = PrivateAttr(default="")
    _products: BaseData = PrivateAttr(default=BaseData())
    # The dependencies of the API.
    # This object should NOT be set in initialization.
    # It should only be set dynamically during forward pass over execution graph,
    # which is conducted internally.
    _deps: BaseData = PrivateAttr(default=BaseData())

    # state should only be things like imported packages, environment variables, etc., which
    # will not be used in actual computation and not be modified. all computed things should be
    # stored in _deps.data and _products.data.
    # Jiahang: we need to some how check this.
    # Jiahang: these notes should be written in the docstring of the class.

    def _hash_members(self):
        members = self.model_dump()
        members['_api_name'] = self._api_name
        members['_products'] = self._products._hash_members()
        members['_deps'] = self._deps._hash_members()
        return members
    
    def _var_repr(self, var) -> str:
        if type(var) == str and var != "data":
            return f"'{var}'"
        return var
    
    def to_api_calling(self) -> str:
        """Convert a BaseAPI object to a string of api calling."""
        params = []
        for name in self.model_fields.keys():
            params.append(f"{name}={self._var_repr(self.__getattribute__(name))}")
        return f"{self._api_name}({', '.join(params)})"

    # Jiahang: be abstractmethod in the future
    def execute(self, state: dict[str, object]):
        """Execute the API call with the given arguments."""
        api_calling = self.to_api_calling()
        state["data"] = deepcopy(self._deps.data)
        results, error = run_codes(api_calling, state)
        if error:# Jiahang: error handling and multiple retry are not implemented yet.
            raise ValueError(error)
        else:
            self._products.data = state["data"]
            return results, api_calling


class BaseDependency(BaseObject):
    """A class representing an edge in the dependency graph.

    This class is used to represent the dependencies between API calls in the
    dependency graph.
    """
    u_api_name: str = Field(default="", description="The name of the source API")
    v_api_name: str = Field(default="", description="The name of the target API")
    args: dict[str, str] = Field(default={}, description="The arguments of the dependency")
    arg_types: dict[str, str] = Field(default={}, description="The argument types of the dependency")
    deps: BaseData = Field(default=BaseData(), description="The data of the dependency")

    def _hash_members(self):
        members = self.model_dump()
        members['deps'] = self.deps._hash_members()
        return members
    

class InputAPI(BaseObject):
    """A class representing an input API.
    
    InputAPI is input from dependency graph JSON structure,
    and will be converted to BaseAPI for internal use.

    This class is created to ease users' efforts to manually create dependency graph.
    But this class is not internally friendly, so will be converted to BaseAPI.
    """
    api: str = Field(default="", description="The name of the API")
    products: list[str] = Field(default=[], description="The products of the API")
    id: str = Field(default="", description="The id of the API")

    @field_validator("products", mode="after")
    @classmethod
    def _check_product(cls, products: list[str]) -> list[str]:
        for product in products:
            product = ast.parse(product)
            assert len(product.body) == 1, "Each product should be a single data object."

            # Jiahang: this assert needs to be carefully considered
            assert isinstance(product.body[0], ast.Expr) and \
                not isinstance(product.body[0].value, ast.Constant), \
                "Each product should be an variable. " \
                "Functions, classes, assignment, constants, etc. are not supported."
        return products
    
    @model_validator(mode="after")
    def _check_id(self) -> "InputAPI":
        assert self.id == self.api, "The id of the API should be the same as the api name."
        return self

class InputDependency(BaseObject):
    """A class representing an input dependency.
    
    InputDependency is input from dependency graph JSON structure,
    and will be converted to BaseDependency for internal use.

    This class is created to ease users' efforts to manually create dependency graph.
    But this class is not internally friendly, so will be converted to BaseDependency.
    """
    dependencies: list[str] = Field(default=[], description="The dependencies of the dependency")
    source: str = Field(default="", description="The source of the dependency")
    target: str = Field(default="", description="The target of the dependency")
    args: dict = Field(default={}, description="The arguments of the dependency")
    arg_types: dict = Field(default={}, description="The argument types of the dependency")

    @model_validator(mode="after")
    def _check_args(self) -> "InputDependency":
        assert len(self.args) == 1 and len(self.arg_types) == 1, "Only one activation arg is permitted for the dependency."
        return self
    
    