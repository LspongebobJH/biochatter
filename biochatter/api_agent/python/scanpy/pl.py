"""Module for interacting with the `scanpy` API for plotting (`pl`)."""

import uuid
from typing import Any
from pydantic import BaseModel, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI
from .base import ScanpyAPI

class ScPlScatter(ScanpyAPI):
    """Parameters for querying the scanpy `pl.scatter` API."""

    _api_name: str = PrivateAttr(
        default="sc.pl.scatter",
    )
    adata: str = Field(default="data", description="Annotated data matrix.")
    x: str | None = Field(default=None, description="x coordinate.")
    y: str | None = Field(default=None, description="y coordinate.")
    color: str | tuple[float, ...] | list[str | tuple[float, ...]] | None = Field(
        default=None,
        description="Keys for annotations of observations/cells or variables/genes, or a hex color specification.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Whether to use raw attribute of adata. Defaults to True if .raw is present.",
    )
    layers: str | list[str] | None = Field(
        default=None,
        description="Layer(s) to use from adata's layers attribute.",
    )
    basis: str | None = Field(
        default=None,
        description="String that denotes a plotting tool that computed coordinates (e.g., 'pca', 'tsne', 'umap').",
    )
    sort_order: bool = Field(
        default=True,
        description="For continuous annotations used as color parameter, plot data points with higher values on top.",
    )
    groups: str | list[str] | None = Field(
        default=None,
        description="Restrict to specific categories in categorical observation annotation.",
    )
    projection: str = Field(
        default="2d",
        description="Projection of plot ('2d' or '3d').",
    )
    legend_loc: str | None = Field(
        default="right margin",
        description="Location of legend ('none', 'right margin', 'on data', etc.).",
    )
    size: int | float | None = Field(
        default=None,
        description="Point size. If None, automatically computed as 120000 / n_cells.",
    )
    color_map: str | None = Field(
        default=None,
        description="Color map to use for continuous variables (e.g., 'magma', 'viridis').",
    )
    show: bool | None = Field(
        default=None,
        description="Show the plot, do not return axis.",
    )
    save: str | bool | None = Field(
        default=None,
        description="If True or a str, save the figure. String is appended to default filename.",
    )


### Embeddings
class ScPlPca(ScanpyAPI):
    """Parameters for querying the scanpy `pl.pca` API."""

    _api_name: str = PrivateAttr(
        default="sc.pl.pca",
    )
    adata: str = Field(
        "data",
        description="Annotated data matrix.",
    )
    color: str | list[str] | None = Field(
        default=None,
        description="Keys for annotations of observations/cells or variables/genes.",
    )
    components: str | list[str] = Field(
        default="1,2",
        description="For example, ['1,2', '2,3']. To plot all available components use 'all'.",
    )
    projection: str = Field(
        default="2d",
        description="Projection of plot.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Location of legend.",
    )
    legend_fontsize: int | float | str | None = Field(
        default=None,
        description="Font size for legend.",
    )
    legend_fontweight: int | str | None = Field(
        default=None,
        description="Font weight for legend.",
    )
    color_map: str | None = Field(
        default=None,
        description="String denoting matplotlib color map.",
    )
    palette: str | list[str] | dict | None = Field(
        default=None,
        description="Colors to use for plotting categorical annotation groups.",
    )
    frameon: bool | None = Field(
        default=None,
        description="Draw a frame around the scatter plot.",
    )
    size: int | float | None = Field(
        default=None,
        description="Point size. If `None`, is automatically computed as 120000 / n_cells.",
    )
    show: bool | None = Field(
        default=None,
        description="Show the plot, do not return axis.",
    )
    save: str | bool | None = Field(
        default=None,
        description="If `True` or a `str`, save the figure.",
    )
    ax: str | None = Field(
        default=None,
        description="A matplotlib axes object.",
    )
    return_fig: bool = Field(
        default=False,
        description="Return the matplotlib figure object.",
    )
    marker: str | None = Field(
        default=".",
        description="Marker symbol.",
    )
    annotate_var_explained: bool = Field(
        default=False,
        description="Annotate the percentage of explained variance.",
    )


class ScPlTsne(ScanpyAPI):
    """Parameters for querying the Scanpy `pl.tsne` API."""

    _api_name: str = PrivateAttr(
        default="sc.pl.tsne",
    )
    adata: str = Field(
        "data",
        description="Annotated data matrix.",
    )
    color: str | list[str] | None = Field(
        default=None,
        description="Keys for annotations of observations/cells or variables/genes.",
    )
    gene_symbols: str | None = Field(
        default=None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
    )
    sort_order: bool = Field(
        default=True,
        description="Plot data points with higher values on top for continuous annotations.",
    )
    edges: bool = Field(
        default=False,
        description="Show edges.",
    )
    edges_width: float = Field(
        default=0.1,
        description="Width of edges.",
    )
    edges_color: str | list[float] | list[str] = Field(
        default="grey",
        description="Color of edges.",
    )
    neighbors_key: str | None = Field(
        default=None,
        description="Key for neighbors connectivities.",
    )
    arrows: bool = Field(
        default=False,
        description="Show arrows (deprecated in favor of `scvelo.pl.velocity_embedding`).",
    )
    arrows_kwds: dict[str, Any] | None = Field(
        default=None,
        description="Arguments passed to `quiver()`.",
    )
    groups: str | None = Field(
        default=None,
        description="Restrict to specific categories in categorical observation annotation.",
    )
    components: str | list[str] | None = Field(
        default=None,
        description="Components to plot, e.g., ['1,2', '2,3']. Use 'all' to plot all available components.",
    )
    projection: str = Field(
        default="2d",
        description="Projection of plot ('2d' or '3d').",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Location of legend.",
    )
    legend_fontsize: int | float | str | None = Field(
        default=None,
        description="Font size for legend.",
    )
    legend_fontweight: int | str = Field(
        default="bold",
        description="Font weight for legend.",
    )
    legend_fontoutline: int | None = Field(
        default=None,
        description="Line width of the legend font outline in pt.",
    )
    size: float | list[float] | None = Field(
        default=None,
        description="Point size. If `None`, computed as 120000 / n_cells.",
    )
    color_map: str | Any | None = Field(
        default=None,
        description="Color map for continuous variables.",
    )
    palette: str | list[str] | Any | None = Field(
        default=None,
        description="Colors for plotting categorical annotation groups.",
    )
    na_color: str | tuple[float, ...] = Field(
        default="lightgray",
        description="Color for null or masked values.",
    )
    na_in_legend: bool = Field(
        default=True,
        description="Include missing values in the legend.",
    )
    frameon: bool | None = Field(
        default=None,
        description="Draw a frame around the scatter plot.",
    )
    vmin: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="Lower limit of the color scale.",
    )
    vmax: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="Upper limit of the color scale.",
    )
    vcenter: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="Center of the color scale, useful for diverging colormaps.",
    )
    norm: Any | None = Field(
        default=None,
        description="Normalization for the colormap.",
    )
    add_outline: bool = Field(
        default=False,
        description="Add a thin border around groups of dots.",
    )
    outline_width: tuple[float, ...] = Field(
        default=(0.3, 0.05),
        description="Width of the outline as a fraction of the scatter dot size.",
    )
    outline_color: tuple[str, ...] = Field(
        default=("black", "white"),
        description="Colors for the outline: border color and gap color.",
    )
    ncols: int = Field(
        default=4,
        description="Number of panels per row.",
    )
    hspace: float = Field(
        default=0.25,
        description="Height of the space between multiple panels.",
    )
    wspace: float | None = Field(
        default=None,
        description="Width of the space between multiple panels.",
    )
    return_fig: bool | None = Field(
        default=None,
        description="Return the matplotlib figure.",
    )
    show: bool | None = Field(
        default=None,
        description="Show the plot; do not return axis.",
    )
    save: str | bool | None = Field(
        default=None,
        description="If `True` or a `str`, save the figure.",
    )
    ax: Any | None = Field(
        default=None,
        description="A matplotlib axes object.",
    )
    # kwargs: dict[str, Any] | None = Field(
    #     default=None,
    #     description="Additional arguments passed to `matplotlib.pyplot.scatter()`.",
    # )
    # Jiahang: kwargs that being sent to internal API are not supported now since it needs
    # to be carefully handled and the handling way should be a standard.


class ScPlUmap(ScanpyAPI):
    """Parameters for querying the Scanpy `pl.umap` API."""

    _api_name: str = PrivateAttr(
        default="sc.pl.umap",
    )
    adata: str = Field(
        "data",
        description="Annotated data matrix.",
    )
    color: str | list[str] | None = Field(
        default=None,
        description="Keys for annotations of observations/cells or variables/genes.",
    )
    mask_obs: str | None = Field(
        default=None,
        description="Mask for observations.",
    )
    gene_symbols: str | None = Field(
        default=None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
    )
    sort_order: bool = Field(
        default=True,
        description="Plot data points with higher values on top for continuous annotations.",
    )
    edges: bool = Field(
        default=False,
        description="Show edges.",
    )
    edges_width: float = Field(
        default=0.1,
        description="Width of edges.",
    )
    edges_color: str | list[float] | list[str] = Field(
        default="grey",
        description="Color of edges.",
    )
    neighbors_key: str | None = Field(
        default=None,
        description="Key for neighbors connectivities.",
    )
    arrows: bool = Field(
        default=False,
        description="Show arrows (deprecated in favor of `scvelo.pl.velocity_embedding`).",
    )
    arrows_kwds: dict[str, Any] | None = Field(
        default=None,
        description="Arguments passed to `quiver()`.",
    )
    groups: str | None = Field(
        default=None,
        description="Restrict to specific categories in categorical observation annotation.",
    )
    components: str | list[str] | None = Field(
        default=None,
        description="Components to plot, e.g., ['1,2', '2,3']. Use 'all' to plot all available components.",
    )
    dimensions: int | None = Field(
        default=None,
        description="Number of dimensions to plot.",
    )
    layer: str | None = Field(
        default=None,
        description="Name of the AnnData object layer to plot.",
    )
    projection: str = Field(
        default="2d",
        description="Projection of plot ('2d' or '3d').",
    )
    scale_factor: float | None = Field(
        default=None,
        description="Scale factor for the plot.",
    )
    color_map: str | Any | None = Field(
        default=None,
        description="Color map for continuous variables.",
    )
    cmap: str | Any | None = Field(
        default=None,
        description="Alias for `color_map`.",
    )
    palette: str | list[str] | Any | None = Field(
        default=None,
        description="Colors for plotting categorical annotation groups.",
    )
    na_color: str | tuple[float, ...] = Field(
        default="lightgray",
        description="Color for null or masked values.",
    )
    na_in_legend: bool = Field(
        default=True,
        description="Include missing values in the legend.",
    )
    size: float | list[float] | None = Field(
        default=None,
        description="Point size. If `None`, computed as 120000 / n_cells.",
    )
    frameon: bool | None = Field(
        default=None,
        description="Draw a frame around the scatter plot.",
    )
    legend_fontsize: int | float | str | None = Field(
        default=None,
        description="Font size for legend.",
    )
    legend_fontweight: int | str = Field(
        default="bold",
        description="Font weight for legend.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Location of legend.",
    )
    legend_fontoutline: int | None = Field(
        default=None,
        description="Line width of the legend font outline in pt.",
    )
    colorbar_loc: str = Field(
        default="right",
        description="Location of the colorbar.",
    )
    vmax: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="Upper limit of the color scale.",
    )
    vmin: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="Lower limit of the color scale.",
    )
    vcenter: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="Center of the color scale, useful for diverging colormaps.",
    )
    norm: Any | None = Field(
        default=None,
        description="Normalization for the colormap.",
    )
    add_outline: bool = Field(
        default=False,
        description="Add a thin border around groups of dots.",
    )
    outline_width: tuple[float, ...] = Field(
        default=(0.3, 0.05),
        description="Width of the outline as a fraction of the scatter dot size.",
    )
    outline_color: tuple[str, ...] = Field(
        default=("black", "white"),
        description="Colors for the outline: border color and gap color.",
    )
    ncols: int = Field(
        default=4,
        description="Number of panels per row.",
    )
    hspace: float = Field(
        default=0.25,
        description="Height of the space between multiple panels.",
    )
    wspace: float | None = Field(
        default=None,
        description="Width of the space between multiple panels.",
    )
    show: bool | None = Field(
        default=None,
        description="Show the plot; do not return axis.",
    )
    save: str | bool | None = Field(
        default=None,
        description="If `True` or a `str`, save the figure.",
    )
    ax: Any | None = Field(
        default=None,
        description="A matplotlib axes object.",
    )
    return_fig: bool | None = Field(
        default=None,
        description="Return the matplotlib figure.",
    )
    marker: str = Field(
        default=".",
        description="Marker symbol.",
    )
    # kwargs: dict[str, Any] | None = Field(
    #     default=None,
    #     description="Additional arguments passed to `matplotlib.pyplot.scatter()`.",
    # )

    


class ScPlDrawGraph(ScanpyAPI):
    """Parameters for querying the Scanpy `pl.draw_graph` API."""

    _api_name: str = PrivateAttr(
        default="sc.pl.draw_graph",
    )
    adata: str = Field(
        "data",
        description="Annotated data matrix.",
    )
    color: str | list[str] | None = Field(
        default=None,
        description="Keys for annotations of observations/cells or variables/genes.",
    )
    gene_symbols: str | None = Field(
        default=None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
    )
    sort_order: bool = Field(
        default=True,
        description=(
            "For continuous annotations used as color parameter, "
            "plot data points with higher values on top of others."
        ),
    )
    edges: bool = Field(
        default=False,
        description="Show edges.",
    )
    edges_width: float = Field(
        default=0.1,
        description="Width of edges.",
    )
    edges_color: str | list[float] | list[str] = Field(
        default="grey",
        description="Color of edges.",
    )
    neighbors_key: str | None = Field(
        default=None,
        description="Where to look for neighbors connectivities.",
    )
    arrows: bool = Field(
        default=False,
        description="Show arrows (deprecated in favor of `scvelo.pl.velocity_embedding`).",
    )
    arrows_kwds: dict[str, Any] | None = Field(
        default=None,
        description="Arguments passed to `quiver()`.",
    )
    groups: str | list[str] | None = Field(
        default=None,
        description="Restrict to a few categories in categorical observation annotation.",
    )
    components: str | list[str] | None = Field(
        default=None,
        description="For instance, ['1,2', '2,3']. To plot all available components use components='all'.",
    )
    projection: str = Field(
        default="2d",
        description="Projection of plot.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Location of legend.",
    )
    legend_fontsize: int | float | str | None = Field(
        default=None,
        description="Numeric size in pt or string describing the size.",
    )
    legend_fontweight: int | str = Field(
        default="bold",
        description="Legend font weight.",
    )
    legend_fontoutline: int | None = Field(
        default=None,
        description="Line width of the legend font outline in pt.",
    )
    colorbar_loc: str | None = Field(
        default="right",
        description="Where to place the colorbar for continuous variables.",
    )
    size: float | list[float] | None = Field(
        default=None,
        description="Point size. If None, is automatically computed as 120000 / n_cells.",
    )
    color_map: str | Any | None = Field(
        default=None,
        description="Color map to use for continuous variables.",
    )
    palette: str | list[str] | Any | None = Field(
        default=None,
        description="Colors to use for plotting categorical annotation groups.",
    )
    na_color: str | tuple[float, ...] = Field(
        default="lightgray",
        description="Color to use for null or masked values.",
    )
    na_in_legend: bool = Field(
        default=True,
        description="If there are missing values, whether they get an entry in the legend.",
    )
    frameon: bool | None = Field(
        default=None,
        description="Draw a frame around the scatter plot.",
    )
    vmin: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="The value representing the lower limit of the color scale.",
    )
    vmax: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="The value representing the upper limit of the color scale.",
    )
    vcenter: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="The value representing the center of the color scale.",
    )
    norm: Any | None = Field(
        default=None,
        description="Normalization for the colormap.",
    )
    add_outline: bool = Field(
        default=False,
        description="Add a thin border around groups of dots.",
    )
    outline_width: tuple[float, ...] = Field(
        default=(0.3, 0.05),
        description="Width of the outline as a fraction of the scatter dot size.",
    )
    outline_color: tuple[str, ...] = Field(
        default=("black", "white"),
        description="Colors for the outline: border color and gap color.",
    )
    ncols: int = Field(
        default=4,
        description="Number of panels per row.",
    )
    hspace: float = Field(
        default=0.25,
        description="Height of the space between multiple panels.",
    )
    wspace: float | None = Field(
        default=None,
        description="Width of the space between multiple panels.",
    )
    return_fig: bool | None = Field(
        default=None,
        description="Return the matplotlib figure.",
    )
    show: bool | None = Field(
        default=None,
        description="Show the plot; do not return axis.",
    )
    save: str | bool | None = Field(
        default=None,
        description="If `True` or a `str`, save the figure.",
    )
    ax: Any | None = Field(
        default=None,
        description="A matplotlib axes object.",
    )
    layout: str | None = Field(
        default=None,
        description="One of the `draw_graph()` layouts.",
    )
    # kwargs: dict[str, Any] | None = Field(
    #     default=None,
    #     description="Additional arguments passed to `matplotlib.pyplot.scatter()`.",
    # )


class ScPlSpatial(ScanpyAPI):
    """Parameters for querying the Scanpy `pl.spatial` API."""

    _api_name: str = PrivateAttr(
        default="sc.pl.spatial",
    )
    adata: str = Field(
        "data",
        description="Annotated data matrix.",
    )
    color: str | list[str] | None = Field(
        default=None,
        description="Keys for annotations of observations/cells or variables/genes.",
    )
    gene_symbols: str | None = Field(
        default=None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
    )
    layer: str | None = Field(
        default=None,
        description="Name of the AnnData object layer to plot.",
    )
    library_id: str | None = Field(
        default=None,
        description="Library ID for Visium data, e.g., key in `adata.uns['spatial']`.",
    )
    img_key: str | None = Field(
        default=None,
        description=(
            "Key for image data, used to get `img` and `scale_factor` from "
            "'images' and 'scalefactors' entries for this library."
        ),
    )
    img: Any | None = Field(
        default=None,
        description="Image data to plot, overrides `img_key`.",
    )
    scale_factor: float | None = Field(
        default=None,
        description="Scaling factor used to map from coordinate space to pixel space.",
    )
    spot_size: float | None = Field(
        default=None,
        description="Diameter of spot (in coordinate space) for each point.",
    )
    crop_coord: tuple[int, ...] | None = Field(
        default=None,
        description="Coordinates to use for cropping the image (left, right, top, bottom).",
    )
    alpha_img: float = Field(
        default=1.0,
        description="Alpha value for image.",
    )
    bw: bool = Field(
        default=False,
        description="Plot image data in grayscale.",
    )
    sort_order: bool = Field(
        default=True,
        description=(
            "For continuous annotations used as color parameter, plot data points "
            "with higher values on top of others."
        ),
    )
    groups: str | list[str] | None = Field(
        default=None,
        description="Restrict to specific categories in categorical observation annotation.",
    )
    components: str | list[str] | None = Field(
        default=None,
        description="For example, ['1,2', '2,3']. To plot all available components, use 'all'.",
    )
    projection: str = Field(
        default="2d",
        description="Projection of plot.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Location of legend.",
    )
    legend_fontsize: int | float | str | None = Field(
        default=None,
        description="Numeric size in pt or string describing the size.",
    )
    legend_fontweight: int | str = Field(
        default="bold",
        description="Legend font weight.",
    )
    legend_fontoutline: int | None = Field(
        default=None,
        description="Line width of the legend font outline in pt.",
    )
    colorbar_loc: str | None = Field(
        default="right",
        description="Where to place the colorbar for continuous variables.",
    )
    size: float = Field(
        default=1.0,
        description="Point size. If None, automatically computed as 120000 / n_cells.",
    )
    color_map: str | Any | None = Field(
        default=None,
        description="Color map to use for continuous variables.",
    )
    palette: str | list[str] | Any | None = Field(
        default=None,
        description="Colors to use for plotting categorical annotation groups.",
    )
    na_color: str | tuple[float, ...] | None = Field(
        default=None,
        description="Color to use for null or masked values.",
    )
    na_in_legend: bool = Field(
        default=True,
        description="If there are missing values, whether they get an entry in the legend.",
    )
    frameon: bool | None = Field(
        default=None,
        description="Draw a frame around the scatter plot.",
    )
    vmin: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="The value representing the lower limit of the color scale.",
    )
    vmax: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="The value representing the upper limit of the color scale.",
    )
    vcenter: str | float | Any | list[str | float | Any] | None = Field(
        default=None,
        description="The value representing the center of the color scale.",
    )
    norm: Any | None = Field(
        default=None,
        description="Normalization for the colormap.",
    )
    add_outline: bool = Field(
        default=False,
        description="Add a thin border around groups of dots.",
    )
    outline_width: tuple[float, ...] = Field(
        default=(0.3, 0.05),
        description="Width of the outline as a fraction of the scatter dot size.",
    )
    outline_color: tuple[str, ...] = Field(
        default=("black", "white"),
        description="Colors for the outline: border color and gap color.",
    )
    ncols: int = Field(
        default=4,
        description="Number of panels per row.",
    )
    hspace: float = Field(
        default=0.25,
        description="Height of the space between multiple panels.",
    )
    wspace: float | None = Field(
        default=None,
        description="Width of the space between multiple panels.",
    )
    return_fig: bool | None = Field(
        default=None,
        description="Return the matplotlib figure.",
    )
    show: bool | None = Field(
        default=None,
        description="Show the plot; do not return axis.",
    )
    save: str | bool | None = Field(
        default=None,
        description="If `True` or a `str`, save the figure.",
    )
    ax: Any | None = Field(
        default=None,
        description="A matplotlib axes object.",
    )
    # kwargs: dict[str, Any] | None = Field(
    #     default=None,
    #     description="Additional arguments passed to `matplotlib.pyplot.scatter()`.",
    # )

TOOLS = [
    ScPlScatter,
    ScPlPca,
    ScPlTsne,
    ScPlUmap,
    ScPlDrawGraph,
    ScPlSpatial,
]

TOOLS_DICT = {tool._api_name.default: tool for tool in TOOLS}