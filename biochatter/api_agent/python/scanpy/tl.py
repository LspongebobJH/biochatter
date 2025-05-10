from pydantic import PrivateAttr, Field

from biochatter.api_agent.base.agent_abc import BaseAPI
from .base import ScanpyAPI

# Jiahang: unfinished
class ScTlUmap(ScanpyAPI):
    """Embed the neighborhood graph using UMAP."""

    _api_name: str = PrivateAttr(default="sc.tl.umap")
    adata: str = Field(
        "data",
        description="Annotated data matrix",
    )
    min_dist: float = Field(
        0.5,
        description="The effective minimum distance between embedded points.",
    )
    spread: float = Field(
        1.0,
        description="The effective scale of embedded points.",
    )
    n_components: int = Field(
        2,
        description="The number of dimensions of the embedding.",
    )

     # Jiahang: there should be a super class of scanpy to import this package into state
    

# Jiahang: unfinished
class ScTlLeiden(ScanpyAPI):
    """Cluster the neighborhood graph using the Leiden algorithm."""

    _api_name: str = PrivateAttr(default="sc.tl.leiden")
    adata: str = Field(
        "data",
        description="Annotated data matrix",
    )
    resolution: float = Field(
        1.0,
        description=(
            "A parameter value controlling the coarseness of the clustering. "
            "Higher values lead to more clusters."
        ),
    )
    random_state: int = Field(
        0,
        description="Random seed for reproducibility.",
    )
    use_weights: bool = Field(
        True,
        description="Whether to use edge weights in the clustering.",
    )

    

TOOLS = [
    ScTlUmap,
    ScTlLeiden,
]

TOOLS_DICT = {tool._api_name.default: tool for tool in TOOLS}
