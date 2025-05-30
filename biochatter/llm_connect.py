"""Module for managing connections to LLM providers and handling conversations.

This module provides classes for connecting to different LLM APIs (OpenAI,
Anthropic, Ollama, etc.) and managing conversations with them, including message
history, context injection, and response correction capabilities.
"""

import json
import logging
import urllib.parse
from abc import ABC, abstractmethod
from typing import Callable  # noqa: UP035

try:
    import streamlit as st
except ImportError:
    st = None

import anthropic
import nltk
import openai
from langchain_anthropic import ChatAnthropic
from langchain_community.chat_models import ChatOllama
from langchain_community.llms.huggingface_hub import HuggingFaceHub
from langchain_core.messages import AIMessage, HumanMessage, SystemMessage
from langchain_openai import AzureChatOpenAI, ChatOpenAI

from ._image import encode_image, encode_image_from_url
from ._stats import get_stats
from .rag_agent import RagAgent
from .selector_agent import RagAgentSelector

logger = logging.getLogger(__name__)

OPENAI_MODELS = [
    "gpt-3.5-turbo",
    "gpt-3.5-turbo-16k",
    "gpt-3.5-turbo-1106",  # further updated 3.5-turbo
    "gpt-4",
    "gpt-4-32k",
    "gpt-4-1106-preview",  # gpt-4 turbo, 128k tokens
]

HUGGINGFACE_MODELS = ["bigscience/bloom"]

XINFERENCE_MODELS = ["custom-endpoint"]

TOKEN_LIMITS = {
    "gpt-3.5-turbo": 4000,
    "gpt-3.5-turbo-16k": 16000,
    "gpt-3.5-turbo-1106": 16000,
    "gpt-4": 8000,
    "gpt-4-32k": 32000,
    "gpt-4-1106-preview": 128000,
    "bigscience/bloom": 1000,
    "custom-endpoint": 1,  # Reasonable value?
}


# Jiahang: it's a bad practice to put LLM initialization in `set_api_key` method,
# since the method name is not related to LLM initialization.
#
# Jiahang: it's unnecessary to define chat, ca_chat in so many ways. We should make it
# unified and simple.
class Conversation(ABC):
    """Use this class to set up a connection to an LLM API.

    Can be used to set the user name and API key, append specific messages for
    system, user, and AI roles (if available), set up the general context as
    well as manual and tool-based data inputs, and finally to query the API
    with prompts made by the user.

    The conversation class is expected to have a `messages` attribute to store
    the conversation, and a `history` attribute, which is a list of messages in
    a specific format for logging / printing.

    """

    def __init__(
        self,
        model_name: str,
        prompts: dict,
        temperature: float = 0.0,
        correct: bool = False,
        split_correction: bool = False,
        use_ragagent_selector: bool = False,
    ) -> None:
        super().__init__()
        self.model_name = model_name
        self.prompts = prompts
        self.correct = correct
        self.split_correction = split_correction
        self.rag_agents: list[RagAgent] = []
        self.history = []
        self.messages = []
        self.ca_messages = []
        self.current_statements = []
        self._use_ragagent_selector = use_ragagent_selector
        self._chat = None
        self._ca_chat = None
        self.temperature = temperature

    @property
    def chat(self):
        """Access the chat attribute with error handling."""
        if self._chat is None:
            msg = "Chat attribute not initialized. Did you call set_api_key()?"
            logger.error(msg)
            raise AttributeError(msg)
        return self._chat

    @chat.setter
    def chat(self, value):
        """Set the chat attribute."""
        self._chat = value

    @property
    def ca_chat(self):
        """Access the correcting agent chat attribute with error handling."""
        if self._ca_chat is None:
            msg = "Correcting agent chat attribute not initialized. Did you call set_api_key()?"
            logger.error(msg)
            raise AttributeError(msg)
        return self._ca_chat

    @ca_chat.setter
    def ca_chat(self, value):
        """Set the correcting agent chat attribute."""
        self._ca_chat = value

    @property
    def use_ragagent_selector(self) -> bool:
        """Whether to use the ragagent selector."""
        return self._use_ragagent_selector

    @use_ragagent_selector.setter
    def use_ragagent_selector(self, val: bool) -> None:
        """Set the use_ragagent_selector attribute."""
        self._use_ragagent_selector = val

    def set_user_name(self, user_name: str) -> None:
        """Set the user name."""
        self.user_name = user_name

    def set_rag_agent(self, agent: RagAgent) -> None:
        """Update or insert rag_agent.

        If the rag_agent with the same mode already exists, it will be updated.
        Otherwise, the new rag_agent will be inserted.
        """
        i, _ = self.find_rag_agent(agent.mode)
        if i < 0:
            # insert
            self.rag_agents.append(agent)
        else:
            # update
            self.rag_agents[i] = agent

    def find_rag_agent(self, mode: str) -> tuple[int, RagAgent]:
        """Find the rag_agent with the given mode."""
        for i, val in enumerate(self.rag_agents):
            if val.mode == mode:
                return i, val
        return -1, None

    @abstractmethod
    def set_api_key(self, api_key: str, user: str | None = None) -> None:
        """Set the API key."""

    def get_prompts(self) -> dict:
        """Get the prompts."""
        return self.prompts

    def set_prompts(self, prompts: dict) -> None:
        """Set the prompts."""
        self.prompts = prompts

    def append_ai_message(self, message: str) -> None:
        """Add a message from the AI to the conversation.

        Args:
        ----
            message (str): The message from the AI.

        """
        self.messages.append(
            AIMessage(
                content=message,
            ),
        )

    def append_system_message(self, message: str) -> None:
        """Add a system message to the conversation.

        Args:
        ----
            message (str): The system message.

        """
        self.messages.append(
            SystemMessage(
                content=message,
            ),
        )

    def append_ca_message(self, message: str) -> None:
        """Add a message to the correcting agent conversation.

        Args:
        ----
            message (str): The message to the correcting agent.

        """
        self.ca_messages.append(
            SystemMessage(
                content=message,
            ),
        )

    def append_user_message(self, message: str) -> None:
        """Add a message from the user to the conversation.

        Args:
        ----
            message (str): The message from the user.

        """
        self.messages.append(
            HumanMessage(
                content=message,
            ),
        )

    def append_image_message(
        self,
        message: str,
        image_url: str,
        local: bool = False,
    ) -> None:
        """Add a user message with an image to the conversation.

        Also checks, in addition to the `local` flag, if the image URL is a
        local file path. If it is local, the image will be encoded as a base64
        string to be passed to the LLM.

        Args:
        ----
            message (str): The message from the user.
            image_url (str): The URL of the image.
            local (bool): Whether the image is local or not. If local, it will
                be encoded as a base64 string to be passed to the LLM.

        """
        parsed_url = urllib.parse.urlparse(image_url)
        if local or not parsed_url.netloc:
            image_url = f"data:image/jpeg;base64,{encode_image(image_url)}"
        else:
            image_url = f"data:image/jpeg;base64,{encode_image_from_url(image_url)}"

        self.messages.append(
            HumanMessage(
                content=[
                    {"type": "text", "text": message},
                    {"type": "image_url", "image_url": {"url": image_url}},
                ],
            ),
        )

    def setup(self, context: str) -> None:
        """Set up the conversation with general prompts and a context."""
        for msg in self.prompts["primary_model_prompts"]:
            if msg:
                self.append_system_message(msg)

        for msg in self.prompts["correcting_agent_prompts"]:
            if msg:
                self.append_ca_message(msg)

        self.context = context
        msg = f"The topic of the research is {context}."
        self.append_system_message(msg)

    def setup_data_input_manual(self, data_input: str) -> None:
        """Set up the data input manually."""
        self.data_input = data_input
        msg = f"The user has given information on the data input: {data_input}."
        self.append_system_message(msg)

    def setup_data_input_tool(self, df, input_file_name: str) -> None:
        """Set up the data input tool."""
        self.data_input_tool = df

        for tool_name in self.prompts["tool_prompts"]:
            if tool_name in input_file_name:
                msg = self.prompts["tool_prompts"][tool_name].format(df=df)
                self.append_system_message(msg)

    def query(
        self,
        text: str,
        image_url: str | None = None,
    ) -> tuple[str, dict | None, str | None]:
        """Query the LLM API using the user's query.

        Appends the most recent query to the conversation, optionally injects
        context from the RAG agent, and runs the primary query method of the
        child class.

        Args:
        ----
            text (str): The user query.

            image_url (str): The URL of an image to include in the conversation.
                Optional and only supported for models with vision capabilities.

        Returns:
        -------
            tuple: A tuple containing the response from the API, the token usage
                information, and the correction if necessary/desired.

        """
        if not image_url:
            self.append_user_message(text)
        else:
            self.append_image_message(text, image_url)

        self._inject_context(text)

        msg, token_usage = self._primary_query()

        if not token_usage:
            # indicates error
            return (msg, token_usage, None)

        if not self.correct:
            return (msg, token_usage, None)

        cor_msg = "Correcting (using single sentences) ..." if self.split_correction else "Correcting ..."

        if st:
            with st.spinner(cor_msg):
                corrections = self._correct_query(text)
        else:
            corrections = self._correct_query(text)

        if not corrections:
            return (msg, token_usage, None)

        correction = "\n".join(corrections)
        return (msg, token_usage, correction)

    def _correct_query(self, msg: str) -> list[str]:
        corrections = []
        if self.split_correction:
            nltk.download("punkt")
            tokenizer = nltk.data.load("tokenizers/punkt/english.pickle")
            sentences = tokenizer.tokenize(msg)
            for sentence in sentences:
                correction = self._correct_response(sentence)

                if str(correction).lower() not in ["ok", "ok."]:
                    corrections.append(correction)
        else:
            correction = self._correct_response(msg)

            if str(correction).lower() not in ["ok", "ok."]:
                corrections.append(correction)

        return corrections

    @abstractmethod
    def _primary_query(self, text: str) -> tuple[str, dict | None]:
        """Run the primary query."""

    @abstractmethod
    def _correct_response(self, msg: str) -> str:
        """Correct the response."""

    def _inject_context_by_ragagent_selector(self, text: str) -> list[str]:
        """Inject the context generated by RagAgentSelector.

        The RagAgentSelector will choose the appropriate rag agent to generate
        context according to user's question.

        Args:
        ----
            text (str): The user query to be used for choosing rag agent

        """
        rag_agents: list[RagAgent] = [agent for agent in self.rag_agents if agent.use_prompt]
        decider_agent = RagAgentSelector(
            rag_agents=rag_agents,
            conversation_factory=lambda: self,
        )
        result = decider_agent.execute(text)
        if result.tool_result is not None and len(result.tool_result) > 0:
            return result.tool_result
        # find rag agent selected
        rag_agent = next(
            [agent for agent in rag_agents if agent.mode == result.answer],
            None,
        )
        if rag_agent is None:
            return None
        return rag_agent.generate_responses(text)

    def _inject_context(self, text: str) -> None:
        """Inject the context received from the RAG agent into the prompt.

        The RAG agent will find the most similar n text fragments and add them
        to the message history object for usage in the next prompt. Uses the
        document summarisation prompt set to inject the context. The ultimate
        prompt should include the placeholder for the statements, `{statements}`
        (used for formatting the string).

        Args:
        ----
            text (str): The user query to be used for similarity search.

        """
        sim_msg = "Performing similarity search to inject fragments ..."

        if st:
            with st.spinner(sim_msg):
                statements = []
                if self.use_ragagent_selector:
                    statements = self._inject_context_by_ragagent_selector(text)
                else:
                    for agent in self.rag_agents:
                        try:
                            docs = agent.generate_responses(text)
                            statements = statements + [doc[0] for doc in docs]
                        except ValueError as e:
                            logger.warning(e)

        else:
            statements = []
            if self.use_ragagent_selector:
                statements = self._inject_context_by_ragagent_selector(text)
            else:
                for agent in self.rag_agents:
                    try:
                        docs = agent.generate_responses(text)
                        statements = statements + [doc[0] for doc in docs]
                    except ValueError as e:
                        logger.warning(e)

        if statements and len(statements) > 0:
            prompts = self.prompts["rag_agent_prompts"]
            self.current_statements = statements
            for i, prompt in enumerate(prompts):
                # if last prompt, format the statements into the prompt
                if i == len(prompts) - 1:
                    self.append_system_message(
                        prompt.format(statements=statements),
                    )
                else:
                    self.append_system_message(prompt)

    def get_last_injected_context(self) -> list[dict]:
        """Get a formatted list of the last context.

        Get the last context injected into the conversation. Contains one
        dictionary for each RAG mode.

        Returns
        -------
            List[dict]: A list of dictionaries containing the mode and context
            for each RAG agent.

        """
        return [{"mode": agent.mode, "context": agent.last_response} for agent in self.rag_agents]

    def get_msg_json(self) -> str:
        """Return a JSON representation of the conversation.

        Returns a list of dicts of the messages in the conversation in JSON
        format. The keys of the dicts are the roles, the values are the
        messages.

        Returns
        -------
            str: A JSON representation of the messages in the conversation.

        """
        d = []
        for msg in self.messages:
            if isinstance(msg, SystemMessage):
                role = "system"
            elif isinstance(msg, HumanMessage):
                role = "user"
            elif isinstance(msg, AIMessage):
                role = "ai"
            else:
                error_msg = f"Unknown message type: {type(msg)}"
                raise TypeError(error_msg)

            d.append({role: msg.content})

        return json.dumps(d)

    def reset(self) -> None:
        """Reset the conversation to the initial state."""
        self.history = []
        self.messages = []
        self.ca_messages = []
        self.current_statements = []


class WasmConversation(Conversation):
    """Conversation class for the wasm model."""

    def __init__(
        self,
        model_name: str,
        prompts: dict,
        temperature: float = 0.0,
        correct: bool = False,
        split_correction: bool = False,
    ) -> None:
        """Initialize the WasmConversation class.

        This class is used to return the complete query as a string to be used
        in the frontend running the wasm model. It does not call the API itself,
        but updates the message history similarly to the other conversation
        classes. It overrides the `query` method from the `Conversation` class
        to return a plain string that contains the entire message for the model
        as the first element of the tuple. The second and third elements are
        `None` as there is no token usage or correction for the wasm model.

        """
        super().__init__(
            model_name=model_name,
            prompts=prompts,
            correct=correct,
            split_correction=split_correction,
            temperature=temperature,
        )

    def query(self, text: str) -> tuple:
        """Return the entire message history as a single string.

        This is the message that is sent to the wasm model.

        Args:
        ----
            text (str): The user query.

        Returns:
        -------
            tuple: A tuple containing the message history as a single string,
                and `None` for the second and third elements of the tuple.

        """
        self.append_user_message(text)

        self._inject_context(text)

        return (self._primary_query(), None, None)

    def _primary_query(self):
        """Concatenate all messages in the conversation.

        Build a single string from all messages in the conversation.
        Currently discards information about roles (system, user).

        Returns
        -------
            str: A single string from all messages in the conversation.

        """
        return "\n".join([m.content for m in self.messages])

    def _correct_response(self, msg: str) -> str:
        """Do not use for the wasm model."""
        return "ok"

    def set_api_key(self, api_key: str, user: str | None = None) -> bool:
        """Do not use for the wasm model."""
        return True


class XinferenceConversation(Conversation):
    """Conversation class for the Xinference deployment."""

    def __init__(
        self,
        base_url: str,
        prompts: dict,
        temperature: float = 0.0,
        model_name: str = "auto",
        correct: bool = False,
        split_correction: bool = False,
    ) -> None:
        """Connect to an open-source LLM via the Xinference client.

        Connect to a running Xinference deployment and set up a conversation
        with the user. Also initialise a second conversational agent to
        provide corrections to the model output, if necessary.

        Args:
        ----
            base_url (str): The base URL of the Xinference instance (should not
            include the /v1 part).

            prompts (dict): A dictionary of prompts to use for the conversation.

            model_name (str): The name of the model to use. Will be mapped to
            the according uid from the list of available models. Can be set to
            "auto" to use the first available model.

            correct (bool): Whether to correct the model output.

            split_correction (bool): Whether to correct the model output by
            splitting the output into sentences and correcting each sentence
            individually.

        """
        # Shaohong: Please keep this xinference importing code here, so that,
        # we don't need to depend on xinference if we dont need it (xinference
        # is expensive to install)
        from xinference.client import Client

        super().__init__(
            model_name=model_name,
            prompts=prompts,
            correct=correct,
            split_correction=split_correction,
            temperature=temperature,
        )
        self.client = Client(base_url=base_url)

        self.models = {}
        self.load_models()

        self.ca_model_name = model_name

        self.set_api_key()

        # TODO make accessible by drop-down

    def load_models(self) -> None:
        """Load the models from the Xinference client."""
        for id, model in self.client.list_models().items():
            model["id"] = id
            self.models[model["model_name"]] = model

    def append_system_message(self, message: str) -> None:
        """Override the system message addition.

        Xinference does not accept multiple system messages. We concatenate them
        if there are multiple.

        Args:
        ----
            message (str): The message to append.

        """
        # if there is not already a system message in self.messages
        if not any(isinstance(m, SystemMessage) for m in self.messages):
            self.messages.append(
                SystemMessage(
                    content=message,
                ),
            )
        else:
            # if there is a system message, append to the last one
            for i, msg in enumerate(self.messages):
                if isinstance(msg, SystemMessage):
                    self.messages[i].content += f"\n{message}"
                    break

    def append_ca_message(self, message: str) -> None:
        """Override the system message addition for the correcting agent.

        Xinference does not accept multiple system messages. We concatenate them
        if there are multiple.

        TODO this currently assumes that the correcting agent is the same model
        as the primary one.

        Args:
        ----
            message (str): The message to append.

        """
        # if there is not already a system message in self.messages
        if not any(isinstance(m, SystemMessage) for m in self.ca_messages):
            self.ca_messages.append(
                SystemMessage(
                    content=message,
                ),
            )
        else:
            # if there is a system message, append to the last one
            for i, msg in enumerate(self.ca_messages):
                if isinstance(msg, SystemMessage):
                    self.ca_messages[i].content += f"\n{message}"
                    break

    def _primary_query(self) -> tuple:
        """Query the Xinference client API.

        Use the user's message and return the response using the message history
        (flattery system messages, prior conversation) as context. Correct the
        response if necessary.

        LLaMA2 architecture does not accept separate system messages, so we
        concatenate the system message with the user message to form the prompt.
        'LLaMA enforces a strict rule that chats should alternate
        user/assistant/user/assistant, and the system message, if present,
        should be embedded into the first user message.' (from
        https://discuss.huggingface.co/t/issue-with-llama-2-chat-template-and-out-of-date-documentation/61645/3)

        Returns
        -------
            tuple: A tuple containing the response from the Xinference API
            (formatted similarly to responses from the OpenAI API) and the token
            usage.

        """
        try:
            history = self._create_history()
            # TODO this is for LLaMA2 arch, may be different for newer models
            prompt = history.pop()
            response = self.model.chat(
                prompt=prompt["content"],
                chat_history=history,
                generate_config={"max_tokens": 2048, "temperature": self.temperature},
            )
        except (
            openai._exceptions.APIError,
            openai._exceptions.OpenAIError,
            openai._exceptions.ConflictError,
            openai._exceptions.NotFoundError,
            openai._exceptions.APIStatusError,
            openai._exceptions.RateLimitError,
            openai._exceptions.APITimeoutError,
            openai._exceptions.BadRequestError,
            openai._exceptions.APIConnectionError,
            openai._exceptions.AuthenticationError,
            openai._exceptions.InternalServerError,
            openai._exceptions.PermissionDeniedError,
            openai._exceptions.UnprocessableEntityError,
            openai._exceptions.APIResponseValidationError,
        ) as e:
            return str(e), None

        msg = response["choices"][0]["message"]["content"]
        token_usage = response["usage"]

        self._update_usage_stats(self.model_name, token_usage)

        self.append_ai_message(msg)

        return msg, token_usage

    def _create_history(self) -> list:
        """Create a history of messages from the conversation.

        Returns
        -------
            list: A list of messages from the conversation.

        """
        history = []
        # extract text components from message contents
        msg_texts = [m.content[0]["text"] if isinstance(m.content, list) else m.content for m in self.messages]

        # check if last message is an image message
        is_image_message = False
        if isinstance(self.messages[-1].content, list):
            is_image_message = self.messages[-1].content[1]["type"] == "image_url"

        # find location of last AI message (if any)
        last_ai_message = None
        for i, m in enumerate(self.messages):
            if isinstance(m, AIMessage):
                last_ai_message = i

        # concatenate all messages before the last AI message into one message
        if last_ai_message:
            history.append(
                {
                    "role": "user",
                    "content": "\n".join(
                        [m for m in msg_texts[:last_ai_message]],
                    ),
                },
            )
            # then append the last AI message
            history.append(
                {
                    "role": "assistant",
                    "content": msg_texts[last_ai_message],
                },
            )

            # then concatenate all messages after that
            # into one HumanMessage
            history.append(
                {
                    "role": "user",
                    "content": "\n".join(
                        [m for m in msg_texts[last_ai_message + 1 :]],
                    ),
                },
            )

        # if there is no AI message, concatenate all messages into one user
        # message
        else:
            history.append(
                {
                    "role": "user",
                    "content": "\n".join([m for m in msg_texts[:]]),
                },
            )

        # if the last message is an image message, add the image to the history
        if is_image_message:
            history[-1]["content"] = [
                {"type": "text", "text": history[-1]["content"]},
                {
                    "type": "image_url",
                    "image_url": {
                        "url": self.messages[-1].content[1]["image_url"]["url"],
                    },
                },
            ]
        return history

    def _correct_response(self, msg: str) -> str:
        """Correct the response from the Xinference API.

        Send the response to a secondary language model. Optionally split the
        response into single sentences and correct each sentence individually.
        Update usage stats.

        Args:
        ----
            msg (str): The response from the model.

        Returns:
        -------
            str: The corrected response (or OK if no correction necessary).

        """
        ca_messages = self.ca_messages.copy()
        ca_messages.append(
            HumanMessage(
                content=msg,
            ),
        )
        ca_messages.append(
            SystemMessage(
                content="If there is nothing to correct, please respond with just 'OK', and nothing else!",
            ),
        )
        history = []
        for m in self.messages:
            if isinstance(m, SystemMessage):
                history.append({"role": "system", "content": m.content})
            elif isinstance(m, HumanMessage):
                history.append({"role": "user", "content": m.content})
            elif isinstance(m, AIMessage):
                history.append({"role": "assistant", "content": m.content})
        prompt = history.pop()
        response = self.ca_model.chat(
            prompt=prompt["content"],
            chat_history=history,
            generate_config={"max_tokens": 2048, "temperature": self.temperature},
        )

        correction = response["choices"][0]["message"]["content"]
        token_usage = response["usage"]

        self._update_usage_stats(self.ca_model_name, token_usage)

        return correction

    def _update_usage_stats(self, model: str, token_usage: dict) -> None:
        """Update redis database with token usage statistics.

        Use the usage_stats object with the increment method.

        Args:
        ----
            model (str): The model name.

            token_usage (dict): The token usage statistics.

        """

    def set_api_key(self) -> bool:
        """Try to get the Xinference model from the client API.

        If the model is found, initialise the conversational agent. If the model
        is not found, `get_model` will raise a RuntimeError.

        Returns
        -------
            bool: True if the model is found, False otherwise.

        """
        try:
            if self.model_name is None or self.model_name == "auto":
                self.model_name = self.list_models_by_type("chat")[0]
            self.model = self.client.get_model(
                self.models[self.model_name]["id"],
            )

            if self.ca_model_name is None or self.ca_model_name == "auto":
                self.ca_model_name = self.list_models_by_type("chat")[0]
            self.ca_model = self.client.get_model(
                self.models[self.ca_model_name]["id"],
            )
            return True

        except RuntimeError:
            self._chat = None
            self._ca_chat = None
            return False

    def list_models_by_type(self, model_type: str) -> list[str]:
        """List the models by type.

        Args:
        ----
            model_type (str): The type of model to list.

        Returns:
        -------
            list[str]: A list of model names.

        """
        names = []
        if model_type in ["embed", "embedding"]:
            for name, model in self.models.items():
                if "model_ability" in model:
                    if "embed" in model["model_ability"]:
                        names.append(name)
                elif model["model_type"] == "embedding":
                    names.append(name)
            return names
        for name, model in self.models.items():
            if "model_ability" in model:
                if model_type in model["model_ability"]:
                    names.append(name)
            elif model["model_type"] == model_type:
                names.append(name)
        return names


class OllamaConversation(Conversation):
    """Conversation class for the Ollama model."""

    def set_api_key(self, api_key: str, user: str | None = None) -> bool:
        """Set the API key for the Ollama API. Not implemented.

        Args:
        ----
            api_key (str): The API key for the Ollama API.

            user (str): The user for usage statistics.

        Returns:
        -------
            bool: True if the API key is valid, False otherwise.

        """
        err = "Ollama does not require an API key."
        raise NotImplementedError(err)

    def __init__(
        self,
        base_url: str,
        prompts: dict,
        temperature: float = 0.0,
        model_name: str = "llama3",
        correct: bool = False,
        split_correction: bool = False,
    ) -> None:
        """Connect to an Ollama LLM via the Ollama/Langchain library.

        Set up a conversation with the user. Also initialise a second
        conversational agent to provide corrections to the model output, if
        necessary.

        Args:
        ----
            base_url (str): The base URL of the Ollama instance.

            prompts (dict): A dictionary of prompts to use for the conversation.

            model_name (str): The name of the model to use. Can be any model
                name available in your Ollama instance.

            correct (bool): Whether to correct the model output.

            split_correction (bool): Whether to correct the model output by
                splitting the output into sentences and correcting each sentence
                individually.

        """
        super().__init__(
            model_name=model_name,
            prompts=prompts,
            correct=correct,
            temperature=temperature,
            split_correction=split_correction,
        )
        self.model_name = model_name
        self.model = ChatOllama(
            base_url=base_url,
            model=self.model_name,
            temperature=self.temperature,
        )

        self.ca_model_name = "mixtral:latest"

        self.ca_model = ChatOllama(
            base_url=base_url,
            model_name=self.ca_model_name,
            temperature=self.temperature,
        )

    def append_system_message(self, message: str) -> None:
        """Override the system message addition.

        Ollama does not accept multiple system messages. Concatenate them if
        there are multiple.

        Args:
        ----
            message (str): The message to append.

        """
        # if there is not already a system message in self.messages
        if not any(isinstance(m, SystemMessage) for m in self.messages):
            self.messages.append(
                SystemMessage(
                    content=message,
                ),
            )
        else:
            # if there is a system message, append to the last one
            for i, msg in enumerate(self.messages):
                if isinstance(msg, SystemMessage):
                    self.messages[i].content += f"\n{message}"
                    break

    def append_ca_message(self, message: str) -> None:
        """Override the system message addition for the correcting agent.

        Ollama does not accept multiple system messages. Concatenate them if
        there are multiple.

        TODO this currently assumes that the correcting agent is the same model
        as the primary one.

        Args:
        ----
            message (str): The message to append.

        """
        # if there is not already a system message in self.messages
        if not any(isinstance(m, SystemMessage) for m in self.ca_messages):
            self.ca_messages.append(
                SystemMessage(
                    content=message,
                ),
            )
        else:
            # if there is a system message, append to the last one
            for i, msg in enumerate(self.ca_messages):
                if isinstance(msg, SystemMessage):
                    self.ca_messages[i].content += f"\n{message}"
                    break

    def _primary_query(self) -> tuple:
        """Query the Ollama client API with the user's message.

        Return the response using the message history (flattery system messages,
        prior conversation) as context. Correct the response if necessary.

        Returns
        -------
            tuple: A tuple containing the response from the Ollama API
            (formatted similarly to responses from the OpenAI API) and the token
            usage.

        """
        try:
            messages = self._create_history(self.messages)
            response = self.model.invoke(
                messages,
                # ,generate_config={"max_tokens": 2048, "temperature": 0},
            )
        except (
            openai._exceptions.APIError,
            openai._exceptions.OpenAIError,
            openai._exceptions.ConflictError,
            openai._exceptions.NotFoundError,
            openai._exceptions.APIStatusError,
            openai._exceptions.RateLimitError,
            openai._exceptions.APITimeoutError,
            openai._exceptions.BadRequestError,
            openai._exceptions.APIConnectionError,
            openai._exceptions.AuthenticationError,
            openai._exceptions.InternalServerError,
            openai._exceptions.PermissionDeniedError,
            openai._exceptions.UnprocessableEntityError,
            openai._exceptions.APIResponseValidationError,
        ) as e:
            return str(e), None
        response_dict = response.dict()
        msg = response_dict["content"]
        token_usage = response_dict["response_metadata"]["eval_count"]

        self._update_usage_stats(self.model_name, token_usage)

        self.append_ai_message(msg)

        return msg, token_usage

    def _create_history(self, messages: list) -> list:
        history = []
        for _, m in enumerate(messages):
            if isinstance(m, AIMessage):
                history.append(AIMessage(content=m.content))
            elif isinstance(m, HumanMessage):
                history.append(HumanMessage(content=m.content))
            elif isinstance(m, SystemMessage):
                history.append(SystemMessage(content=m.content))

        return history

    def _correct_response(self, msg: str) -> str:
        """Correct the response from the Ollama API.

        Send the response to a secondary language model. Optionally split the
        response into single sentences and correct each sentence individually.
        Update usage stats.

        Args:
        ----
            msg (str): The response from the model.

        Returns:
        -------
            str: The corrected response (or OK if no correction necessary).

        """
        ca_messages = self.ca_messages.copy()
        ca_messages.append(
            HumanMessage(
                content=msg,
            ),
        )
        ca_messages.append(
            SystemMessage(
                content="If there is nothing to correct, please respond with just 'OK', and nothing else!",
            ),
        )
        response = self.ca_model.invoke(
            chat_history=self._create_history(self.messages),
        ).dict()
        correction = response["content"]
        token_usage = response["eval_count"]

        self._update_usage_stats(self.ca_model_name, token_usage)

        return correction

    def _update_usage_stats(self, model: str, token_usage: dict) -> None:
        """Update redis database with token usage statistics.

        Use the usage_stats object with the increment method.

        Args:
        ----
            model (str): The model name.

            token_usage (dict): The token usage statistics.

        """


class AnthropicConversation(Conversation):
    """Conversation class for the Anthropic model."""

    def __init__(
        self,
        model_name: str,
        prompts: dict,
        temperature: float = 0.0,
        correct: bool = False,
        split_correction: bool = False,
    ) -> None:
        """Connect to Anthropic's API and set up a conversation with the user.

        Also initialise a second conversational agent to provide corrections to
        the model output, if necessary.

        Args:
        ----
            model_name (str): The name of the model to use.

            prompts (dict): A dictionary of prompts to use for the conversation.

            split_correction (bool): Whether to correct the model output by
                splitting the output into sentences and correcting each
                sentence individually.

        """
        super().__init__(
            model_name=model_name,
            prompts=prompts,
            correct=correct,
            temperature=temperature,
            split_correction=split_correction,
        )

        self.ca_model_name = "claude-3-5-sonnet-20240620"
        # TODO make accessible by drop-down

    def set_api_key(self, api_key: str, user: str | None = None) -> bool:
        """Set the API key for the Anthropic API.

        If the key is valid, initialise the conversational agent. Optionally set
        the user for usage statistics.

        Args:
        ----
            api_key (str): The API key for the Anthropic API.

            user (str, optional): The user for usage statistics. If provided and
                equals "community", will track usage stats.

        Returns:
        -------
            bool: True if the API key is valid, False otherwise.

        """
        client = anthropic.Anthropic(
            api_key=api_key,
        )
        self.user = user

        try:
            client.count_tokens("Test connection")
            self.chat = ChatAnthropic(
                model_name=self.model_name,
                temperature=self.temperature,
                api_key=api_key,
            )
            self.ca_chat = ChatAnthropic(
                model_name=self.ca_model_name,
                temperature=self.temperature,
                api_key=api_key,
            )
            if user == "community":
                self.usage_stats = get_stats(user=user)

            return True

        except anthropic._exceptions.AuthenticationError:
            self._chat = None
            self._ca_chat = None
            return False

    def _primary_query(self) -> tuple:
        """Query the Anthropic API with the user's message.

        Return the response using the message history (flattery system messages,
        prior conversation) as context. Correct the response if necessary.

        Returns
        -------
            tuple: A tuple containing the response from the Anthropic API and
                the token usage.

        """
        try:
            history = self._create_history()
            response = self.chat.generate([history])
        except (
            anthropic._exceptions.APIError,
            anthropic._exceptions.AnthropicError,
            anthropic._exceptions.ConflictError,
            anthropic._exceptions.NotFoundError,
            anthropic._exceptions.APIStatusError,
            anthropic._exceptions.RateLimitError,
            anthropic._exceptions.APITimeoutError,
            anthropic._exceptions.BadRequestError,
            anthropic._exceptions.APIConnectionError,
            anthropic._exceptions.AuthenticationError,
            anthropic._exceptions.InternalServerError,
            anthropic._exceptions.PermissionDeniedError,
            anthropic._exceptions.UnprocessableEntityError,
            anthropic._exceptions.APIResponseValidationError,
        ) as e:
            return str(e), None

        msg = response.generations[0][0].text
        token_usage = response.llm_output.get("token_usage")

        self.append_ai_message(msg)

        return msg, token_usage

    def _create_history(self) -> list:
        """Create a history of messages for the Anthropic API.

        Returns
        -------
            list: A list of messages, with the last message being the most
                recent.

        """
        history = []
        # extract text components from message contents
        msg_texts = [m.content[0]["text"] if isinstance(m.content, list) else m.content for m in self.messages]

        # check if last message is an image message
        is_image_message = False
        if isinstance(self.messages[-1].content, list):
            is_image_message = self.messages[-1].content[1]["type"] == "image_url"

        # find location of last AI message (if any)
        last_ai_message = None
        for i, m in enumerate(self.messages):
            if isinstance(m, AIMessage):
                last_ai_message = i

        # Aggregate system messages into one message at the beginning
        system_messages = [m.content for m in self.messages if isinstance(m, SystemMessage)]
        if system_messages:
            history.append(
                SystemMessage(content="\n".join(system_messages)),
            )

        # concatenate all messages before the last AI message into one message
        if last_ai_message is not None:
            history.append(
                HumanMessage(
                    content="\n".join([m for m in msg_texts[:last_ai_message]]),
                ),
            )
            # then append the last AI message
            history.append(
                AIMessage(
                    content=msg_texts[last_ai_message],
                ),
            )

            # then concatenate all messages after that
            # into one HumanMessage
            history.append(
                HumanMessage(
                    content="\n".join(
                        [m for m in msg_texts[last_ai_message + 1 :]],
                    ),
                ),
            )

        # else add human message to history (without system messages)
        else:
            last_system_message = None
            for i, m in enumerate(self.messages):
                if isinstance(m, SystemMessage):
                    last_system_message = i
            history.append(
                HumanMessage(
                    content="\n".join(
                        [m for m in msg_texts[last_system_message + 1 :]],
                    ),
                ),
            )

        # if the last message is an image message, add the image to the history
        if is_image_message:
            history[-1]["content"] = [
                {"type": "text", "text": history[-1]["content"]},
                {
                    "type": "image_url",
                    "image_url": {
                        "url": self.messages[-1].content[1]["image_url"]["url"],
                    },
                },
            ]
        return history

    def _correct_response(self, msg: str) -> str:
        """Correct the response from the Anthropic API.

        Send the response to a secondary language model. Optionally split the
        response into single sentences and correct each sentence individually.
        Update usage stats.

        Args:
        ----
            msg (str): The response from the Anthropic API.

        Returns:
        -------
            str: The corrected response (or OK if no correction necessary).

        """
        ca_messages = self.ca_messages.copy()
        ca_messages.append(
            HumanMessage(
                content=msg,
            ),
        )
        ca_messages.append(
            SystemMessage(
                content="If there is nothing to correct, please respond with just 'OK', and nothing else!",
            ),
        )

        response = self.ca_chat.generate([ca_messages])

        correction = response.generations[0][0].text
        token_usage = response.llm_output.get("token_usage")

        return correction


class GptConversation(Conversation):
    """Conversation class for the OpenAI GPT model."""

    def __init__(
        self,
        model_name: str,
        prompts: dict,
        temperature: float = 0.0,
        correct: bool = False,
        split_correction: bool = False,
        base_url: str = None,
        update_token_usage: Callable | None = None,
    ) -> None:
        """Connect to OpenAI's GPT API and set up a conversation with the user.

        Also initialise a second conversational agent to provide corrections to
        the model output, if necessary.

        Args:
        ----
            model_name (str): The name of the model to use.

            prompts (dict): A dictionary of prompts to use for the conversation.

            split_correction (bool): Whether to correct the model output by
                splitting the output into sentences and correcting each
                sentence individually.

            base_url (str): Optional OpenAI base_url value to use custom
                endpoint URL instead of default

        """
        super().__init__(
            model_name=model_name,
            prompts=prompts,
            correct=correct,
            temperature=temperature,
            split_correction=split_correction,
        )
        self.base_url = base_url
        self.ca_model_name = "gpt-3.5-turbo"
        # TODO make accessible by drop-down

        self._update_token_usage = update_token_usage

    def set_api_key(self, api_key: str, user: str | None = None) -> bool:
        """Set the API key for the OpenAI API.

        If the key is valid, initialise the conversational agent. Optionally set
        the user for usage statistics.

        Args:
        ----
            api_key (str): The API key for the OpenAI API.

            user (str, optional): The user for usage statistics. If provided and
                equals "community", will track usage stats.

        Returns:
        -------
            bool: True if the API key is valid, False otherwise.

        """
        self.user = user

        try:
            self.chat = ChatOpenAI(
                model_name=self.model_name,
                temperature=self.temperature,
                openai_api_key=api_key,
                base_url=self.base_url,
            )
            self.ca_chat = ChatOpenAI(
                model_name=self.ca_model_name,
                temperature=self.temperature,
                openai_api_key=api_key,
                base_url=self.base_url,
            )
            if user == "community":
                self.usage_stats = get_stats(user=user)

            return True

        except openai._exceptions.AuthenticationError:
            self._chat = None
            self._ca_chat = None
            return False

    def _primary_query(self) -> tuple:
        """Query the OpenAI API with the user's message.

        Return the response using the message history (flattery system messages,
        prior conversation) as context. Correct the response if necessary.

        Returns
        -------
            tuple: A tuple containing the response from the OpenAI API and the
                token usage.

        """
        try:
            response = self.chat.generate([self.messages])
        except (
            openai._exceptions.APIError,
            openai._exceptions.OpenAIError,
            openai._exceptions.ConflictError,
            openai._exceptions.NotFoundError,
            openai._exceptions.APIStatusError,
            openai._exceptions.RateLimitError,
            openai._exceptions.APITimeoutError,
            openai._exceptions.BadRequestError,
            openai._exceptions.APIConnectionError,
            openai._exceptions.AuthenticationError,
            openai._exceptions.InternalServerError,
            openai._exceptions.PermissionDeniedError,
            openai._exceptions.UnprocessableEntityError,
            openai._exceptions.APIResponseValidationError,
        ) as e:
            return str(e), None

        msg = response.generations[0][0].text
        token_usage = response.llm_output.get("token_usage")

        self._update_usage_stats(self.model_name, token_usage)

        self.append_ai_message(msg)

        return msg, token_usage

    def _correct_response(self, msg: str) -> str:
        """Correct the response from the OpenAI API.

        Send the response to a secondary language model. Optionally split the
        response into single sentences and correct each sentence individually.
        Update usage stats.

        Args:
        ----
            msg (str): The response from the OpenAI API.

        Returns:
        -------
            str: The corrected response (or OK if no correction necessary).

        """
        ca_messages = self.ca_messages.copy()
        ca_messages.append(
            HumanMessage(
                content=msg,
            ),
        )
        ca_messages.append(
            SystemMessage(
                content="If there is nothing to correct, please respond with just 'OK', and nothing else!",
            ),
        )

        response = self.ca_chat.generate([ca_messages])

        correction = response.generations[0][0].text
        token_usage = response.llm_output.get("token_usage")

        self._update_usage_stats(self.ca_model_name, token_usage)

        return correction

    def _update_usage_stats(self, model: str, token_usage: dict) -> None:
        """Update redis database with token usage statistics.

        Use the usage_stats object with the increment method.

        Args:
        ----
            model (str): The model name.

            token_usage (dict): The token usage statistics.

        """
        if self.user == "community":
            # Only process integer values
            stats_dict = {f"{k}:{model}": v for k, v in token_usage.items() if isinstance(v, int | float)}
            self.usage_stats.increment(
                "usage:[date]:[user]",
                stats_dict,
            )

        if self._update_token_usage is not None:
            self._update_token_usage(self.user, model, token_usage)


class AzureGptConversation(GptConversation):
    """Conversation class for the Azure GPT model."""

    def __init__(
        self,
        deployment_name: str,
        model_name: str,
        prompts: dict,
        temperature: float = 0.0,
        correct: bool = False,
        split_correction: bool = False,
        version: str | None = None,
        base_url: str | None = None,
        update_token_usage: Callable | None = None,
    ) -> None:
        """Connect to Azure's GPT API and set up a conversation with the user.

        Extends GptConversation.

        Args:
        ----
            deployment_name (str): The name of the Azure deployment to use.

            model_name (str): The name of the model to use. This is distinct
                from the deployment name.

            prompts (dict): A dictionary of prompts to use for the conversation.

            correct (bool): Whether to correct the model output.

            split_correction (bool): Whether to correct the model output by
                splitting the output into sentences and correcting each
                sentence individually.

            version (str): The version of the Azure API to use.

            base_url (str): The base URL of the Azure API to use.

            update_token_usage (Callable): A function to update the token usage
                statistics.

        """
        super().__init__(
            model_name=model_name,
            prompts=prompts,
            correct=correct,
            temperature=temperature,
            split_correction=split_correction,
            update_token_usage=update_token_usage,
        )

        self.version = version
        self.base_url = base_url
        self.deployment_name = deployment_name

    def set_api_key(self, api_key: str, user: str | None = None) -> bool:
        """Set the API key for the Azure API.

        If the key is valid, initialise the conversational agent. No user stats
        on Azure.

        Args:
        ----
            api_key (str): The API key for the Azure API.

            user (str, optional): The user for usage statistics.

        Returns:
        -------
            bool: True if the API key is valid, False otherwise.

        """
        try:
            self.chat = AzureChatOpenAI(
                deployment_name=self.deployment_name,
                model_name=self.model_name,
                openai_api_version=self.version,
                azure_endpoint=self.base_url,
                openai_api_key=api_key,
                temperature=self.temperature,
            )
            self.ca_chat = AzureChatOpenAI(
                deployment_name=self.deployment_name,
                model_name=self.model_name,
                openai_api_version=self.version,
                azure_endpoint=self.base_url,
                openai_api_key=api_key,
                temperature=self.temperature,
            )

            self.chat.generate([[HumanMessage(content="Hello")]])
            self.user = user if user is not None else "Azure Community"

            return True

        except openai._exceptions.AuthenticationError:
            self._chat = None
            self._ca_chat = None
            return False

    def _update_usage_stats(self, model: str, token_usage: dict) -> None:
        if self._update_token_usage is not None:
            self._update_token_usage(self.user, model, token_usage)


class BloomConversation(Conversation):
    """Conversation class for the Bloom model."""

    def __init__(
        self,
        model_name: str,
        prompts: dict,
        split_correction: bool,
        temperature: float = 0.0,
    ) -> None:
        """Initialise the BloomConversation class.

        DEPRECATED: Superceded by XinferenceConversation.
        """
        super().__init__(
            model_name=model_name,
            prompts=prompts,
            temperature=temperature,
            split_correction=split_correction,
        )

        self.messages = []

    def set_api_key(self, api_key: str, user: str | None = None) -> bool:
        """Set the API key for the HuggingFace API.

        If the key is valid, initialise the conversational agent.

        Args:
        ----
            api_key (str): The API key for the HuggingFace API.

            user (str): The user for usage statistics.

        Returns:
        -------
            bool: True if the API key is valid, False otherwise.

        """
        self.chat = HuggingFaceHub(
            repo_id=self.model_name,
            model_kwargs={"temperature": self.temperature},  # "regular sampling"
            # as per https://huggingface.co/docs/api-inference/detailed_parameters
            huggingfacehub_api_token=api_key,
        )

        try:
            self.chat.generate(["Hello, I am a biomedical researcher."])
            return True
        except ValueError:
            return False

    def _cast_messages(self, messages: list) -> str:
        """Render the different roles of the chat-based conversation."""
        cast = ""
        for m in messages:
            if isinstance(m, SystemMessage):
                cast += f"System: {m.content}\n"
            elif isinstance(m, HumanMessage):
                cast += f"Human: {m.content}\n"
            elif isinstance(m, AIMessage):
                cast += f"AI: {m.content}\n"
            else:
                error_msg = f"Unknown message type: {type(m)}"
                raise TypeError(error_msg)

        return cast

    def _primary_query(self) -> tuple:
        response = self.chat.generate([self._cast_messages(self.messages)])

        msg = response.generations[0][0].text
        token_usage = {
            "prompt_tokens": 0,
            "completion_tokens": 0,
            "total_tokens": 0,
        }

        self.append_ai_message(msg)

        return msg, token_usage

    def _correct_response(self, msg: str) -> str:
        return "ok"
