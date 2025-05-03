from collections.abc import Collection
from typing import Literal

from pydantic import BaseModel, PrivateAttr, Field

# Jiahang: unfinished
class ScTlUmap(BaseModel):
    """Embed the neighborhood graph using UMAP."""

    _api_name: str = PrivateAttr(default="sc.tl.umap")
    adata: str = Field(
        "adata",
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

# Jiahang: unfinished
class ScTlLeiden(BaseModel):
    """Cluster the neighborhood graph using the Leiden algorithm."""

    _api_name: str = PrivateAttr(default="sc.tl.leiden")
    adata: str = Field(
        "adata",
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
