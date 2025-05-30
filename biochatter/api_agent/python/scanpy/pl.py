from __future__ import annotations

from typing import Any, Optional
from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScPlPaga(BaseAPI):
    """
    Plot the PAGA graph through thresholding low-connectivity edges. Compute a coarse-grained layout of the data. Reuse this layout for obtaining embeddings with more meaningful global topology. This method uses ForceAtlas2 or igraph's layout algorithms for most layouts.
    """

    adata: Any = Field(..., description="Annotated data matrix.", title="Adata")
    threshold: Any = Field(
        None,
        description="Defines the threshold for edge weights, below which edges are not drawn.",
        title="Threshold",
    )
    color: Any = Field(
        None,
        description="Specifies node colors based on gene names or annotations, and can be used to visualize pie charts at nodes.",
        title="Color",
    )
    layout: Any = Field(
        None,
        description="Determines the plotting layout for computing positions on the graph.",
        title="Layout",
    )
    layout_kwds: Optional[Any] = Field(
        {}, description="Keywords for the layout.", title="Layout Kwds"
    )
    init_pos: Any = Field(
        None,
        description="Initial x and y coordinates for layout initialization.",
        title="Init Pos",
    )
    root: Optional[Any] = Field(
        0, description="Defines the root node for tree layouts.", title="Root"
    )
    labels: Any = Field(
        None, description="Node labels for visualization.", title="Labels"
    )
    single_component: Optional[Any] = Field(
        False,
        description="Restricts the graph to the largest connected component.",
        title="Single Component",
    )
    solid_edges: Optional[Any] = Field(
        "connectivities",
        description="Key specifying the matrix for solid black edges.",
        title="Solid Edges",
    )
    dashed_edges: Any = Field(
        None,
        description="Key specifying the matrix for dashed grey edges.",
        title="Dashed Edges",
    )
    transitions: Any = Field(
        None, description="Key specifying the matrix for arrows.", title="Transitions"
    )
    fontsize: Any = Field(
        None, description="Font size for node labels.", title="Fontsize"
    )
    fontweight: Optional[Any] = Field(
        "bold", description="Font weight for text elements.", title="Fontweight"
    )
    fontoutline: Any = Field(
        None,
        description="Width of the white outline around fonts.",
        title="Fontoutline",
    )
    text_kwds: Optional[Any] = Field(
        {}, description="Keywords for text elements.", title="Text Kwds"
    )
    node_size_scale: Optional[Any] = Field(
        1.0, description="Scale factor for node sizes.", title="Node Size Scale"
    )
    node_size_power: Optional[Any] = Field(
        0.5,
        description="Influence of group sizes on node radius.",
        title="Node Size Power",
    )
    edge_width_scale: Optional[Any] = Field(
        1.0, description="Scale factor for edge widths.", title="Edge Width Scale"
    )
    min_edge_width: Any = Field(
        None, description="Minimum width of solid edges.", title="Min Edge Width"
    )
    max_edge_width: Any = Field(
        None,
        description="Maximum width of solid and dashed edges.",
        title="Max Edge Width",
    )
    arrowsize: Optional[Any] = Field(
        30, description="Size of arrow heads for directed graphs.", title="Arrowsize"
    )
    title: Any = Field(None, description="Title for the plot.", title="Title")
    left_margin: Optional[Any] = Field(
        0.01, description="Left margin specification.", title="Left Margin"
    )
    random_state: Optional[Any] = Field(
        0, description="Seed for random initialization.", title="Random State"
    )
    pos: Any = Field(
        None, description="Coordinates for drawing the graph.", title="Pos"
    )
    normalize_to_color: Optional[Any] = Field(
        False,
        description="Option to normalize categorical plots to color or groups.",
        title="Normalize To Color",
    )
    cmap: Any = Field(None, description="Color map for visualization.", title="Cmap")
    cax: Any = Field(
        None, description="Matplotlib axes object for colorbar.", title="Cax"
    )
    colorbar: Any = Field(
        None, description="Indicates the presence of a colorbar.", title="Colorbar"
    )
    cb_kwds: Optional[Any] = Field(
        {}, description="Keyword arguments for the colorbar.", title="Cb Kwds"
    )
    frameon: Any = Field(
        None,
        description="Option to draw a frame around the PAGA graph.",
        title="Frameon",
    )
    add_pos: Optional[Any] = Field(
        True, description="Add positions to the PAGA graph.", title="Add Pos"
    )
    export_to_gexf: Optional[Any] = Field(
        False, description="Export to the gexf format.", title="Export To Gexf"
    )
    use_raw: Optional[Any] = Field(
        True, description="Indicates whether raw data is used.", title="Use Raw"
    )
    colors: Any = Field(None, description="No description available.", title="Colors")
    groups: Any = Field(None, description="No description available.", title="Groups")
    plot: Optional[Any] = Field(
        True,
        description="Option to create the figure or only compute the layout.",
        title="Plot",
    )
    show: Any = Field(None, description="Option to display the plot.", title="Show")
    save: Any = Field(None, description="Option to save the figure.", title="Save")
    ax: Any = Field(None, description="Matplotlib axes object.", title="Ax")
    _api_name: str = PrivateAttr(default="sc.pl.paga")
    _products_original: list[str] = PrivateAttr(default=['data.uns["paga"]["pos"]'])
    _data_name: str = PrivateAttr(default="adata")


class ScPlScatter(BaseAPI):
    """
    Scatter plot along observations or variables axes. Color the plot using annotations of observations (.obs), variables (.var) or expression of genes (.var_names).
    """

    adata: Any = Field(
        ..., description="Annotated data matrix with type annotation", title="Adata"
    )
    x: Any = Field(
        None,
        description="The x coordinate with type annotation of string or None",
        title="X",
    )
    y: Any = Field(
        None,
        description="The y coordinate with type annotation of string or None",
        title="Y",
    )
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes, or a hex color specification, with possible types including str, ColorLike, Collection[str] or None",
        title="Color",
    )
    use_raw: Any = Field(
        None,
        description="Boolean indicating whether to use the raw attribute of adata, defaulting to True if .raw is present",
        title="Use Raw",
    )
    layers: Any = Field(
        None,
        description="Specification of the layer for x, y, and color using the layers attribute of adata if present, with possible types being str, Collection[str], or None",
        title="Layers",
    )
    sort_order: Optional[Any] = Field(
        True,
        description="Boolean indicating whether to plot data points with higher values on top for continuous annotations used as color parameter",
        title="Sort Order",
    )
    alpha: Any = Field(None, description="Floating point number or None", title="Alpha")
    basis: Any = Field(
        None,
        description="String denoting a plotting tool that computed coordinates or None",
        title="Basis",
    )
    groups: Any = Field(
        None,
        description="Restriction to specific categories in categorical observation annotation, with possible types including str, Iterable[str], or None",
        title="Groups",
    )
    components: Any = Field(
        None,
        description="String specifying components for plotting or 'all', with possible types being str, Collection[str], or None",
        title="Components",
    )
    projection: Optional[Any] = Field(
        "2d", description="Projection of plot, either '2d' or '3d'", title="Projection"
    )
    legend_loc: Optional[Any] = Field(
        "right margin",
        description="Location of legend, with possible values including 'on data', 'right margin', None, or a valid keyword for the loc parameter of matplotlib.legend.Legend",
        title="Legend Loc",
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size or string description for legend font size",
        title="Legend Fontsize",
    )
    legend_fontweight: Any = Field(
        None,
        description="Legend font weight as a numeric value or string, with default to 'bold' or 'normal' based on legend_loc",
        title="Legend Fontweight",
    )
    legend_fontoutline: Any = Field(
        None,
        description="Line width of the legend font outline in pt or None",
        title="Legend Fontoutline",
    )
    color_map: Any = Field(
        None,
        description="Color map to use for continuous variables, can be a name or matplotlib.colors.Colormap instance",
        title="Color Map",
    )
    palette: Any = Field(
        None,
        description="Colors to use for plotting categorical annotation groups, with various possible types including Cycler, ListedColormap, ColorLike, or Sequence[ColorLike]",
        title="Palette",
    )
    frameon: Any = Field(
        None,
        description="Boolean indicating whether to draw a frame around the scatter plot",
        title="Frameon",
    )
    right_margin: Any = Field(
        None, description="Floating point number or None", title="Right Margin"
    )
    left_margin: Any = Field(
        None, description="Floating point number or None", title="Left Margin"
    )
    size: Any = Field(
        None,
        description="Point size, automatically computed if None, or a sequence containing the size for each cell",
        title="Size",
    )
    marker: Optional[Any] = Field(
        ".", description="String or sequence of strings for marker", title="Marker"
    )
    title: Any = Field(
        None,
        description="Title for panels as a string or list of strings",
        title="Title",
    )
    show: Any = Field(
        None, description="Boolean indicating whether to show the plot", title="Show"
    )
    save: Any = Field(
        None, description="Save the figure with True, a string, or None", title="Save"
    )
    ax: Any = Field(None, description="A matplotlib axes object or None", title="Ax")
    _api_name: str = PrivateAttr(default="sc.pl.scatter")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlUmap(BaseAPI):
    """
    Scatter plot in UMAP basis.
    """

    adata: Any = Field(..., description="Annotated data matrix.", title="Adata")
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes.",
        title="Color",
    )
    mask_obs: Any = Field(
        None, description="No description available.", title="Mask Obs"
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
        title="Gene Symbols",
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
        title="Use Raw",
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points with higher values on top of others.",
        title="Sort Order",
    )
    edges: Optional[Any] = Field(False, description="Show edges.", title="Edges")
    edges_width: Optional[Any] = Field(
        0.1, description="Width of edges.", title="Edges Width"
    )
    edges_color: Optional[Any] = Field(
        "grey", description="Color of edges.", title="Edges Color"
    )
    neighbors_key: Any = Field(
        None,
        description="Where to look for neighbors connectivities.",
        title="Neighbors Key",
    )
    arrows: Optional[Any] = Field(
        False,
        description="Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).",
        title="Arrows",
    )
    arrows_kwds: Any = Field(
        None,
        description="Passed to :meth:`~matplotlib.axes.Axes.quiver`.",
        title="Arrows Kwds",
    )
    groups: Any = Field(
        None,
        description="Restrict to a few categories in categorical observation annotation.",
        title="Groups",
    )
    components: Any = Field(
        None,
        description="To plot all available components use `components='all'.",
        title="Components",
    )
    dimensions: Any = Field(
        None,
        description="0-indexed dimensions of the embedding to plot as integers.",
        title="Dimensions",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer that wants to be plotted.",
        title="Layer",
    )
    projection: Optional[Any] = Field(
        "2d", description="Projection of plot (default: `'2d'`).", title="Projection"
    )
    scale_factor: Any = Field(
        None, description="No description available.", title="Scale Factor"
    )
    color_map: Any = Field(
        None, description="Color map to use for continous variables.", title="Color Map"
    )
    cmap: Any = Field(None, description="No description available.", title="Cmap")
    palette: Any = Field(
        None,
        description="Colors to use for plotting categorical annotation groups.",
        title="Palette",
    )
    na_color: Optional[Any] = Field(
        "lightgray",
        description="Color to use for null or masked values.",
        title="Na Color",
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="If there are missing values, whether they get an entry in the legend.",
        title="Na In Legend",
    )
    size: Any = Field(None, description="Point size.", title="Size")
    frameon: Any = Field(
        None, description="Draw a frame around the scatter plot.", title="Frameon"
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string describing the size.",
        title="Legend Fontsize",
    )
    legend_fontweight: Optional[Any] = Field(
        "bold", description="Legend font weight.", title="Legend Fontweight"
    )
    legend_loc: Optional[Any] = Field(
        "right margin", description="Location of legend.", title="Legend Loc"
    )
    legend_fontoutline: Any = Field(
        None,
        description="Line width of the legend font outline in pt.",
        title="Legend Fontoutline",
    )
    colorbar_loc: Optional[Any] = Field(
        "right",
        description="Where to place the colorbar for continous variables.",
        title="Colorbar Loc",
    )
    vmax: Any = Field(
        None,
        description="The value representing the upper limit of the color scale.",
        title="Vmax",
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale.",
        title="Vmin",
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale.",
        title="Vcenter",
    )
    norm: Any = Field(None, description="No description available.", title="Norm")
    add_outline: Optional[Any] = Field(
        False,
        description="If set to True, this will add a thin border around groups of dots.",
        title="Add Outline",
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description="Tuple with two width numbers used to adjust the outline.",
        title="Outline Width",
    )
    outline_color: Optional[Any] = Field(
        ["black", "white"],
        description="Tuple with two valid color names used to adjust the add_outline.",
        title="Outline Color",
    )
    ncols: Optional[Any] = Field(
        4, description="Number of panels per row.", title="Ncols"
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.",
        title="Hspace",
    )
    wspace: Any = Field(
        None,
        description="Adjust the width of the space between multiple panels.",
        title="Wspace",
    )
    title: Any = Field(
        None,
        description="Provide title for panels either as string or list of strings.",
        title="Title",
    )
    show: Any = Field(
        None, description="Show the plot, do not return axis.", title="Show"
    )
    save: Any = Field(
        None, description="If `True` or a `str`, save the figure.", title="Save"
    )
    ax: Any = Field(None, description="A matplotlib axes object.", title="Ax")
    return_fig: Any = Field(
        None, description="Return the matplotlib figure.", title="Return Fig"
    )
    marker: Optional[Any] = Field(
        ".", description="No description available.", title="Marker"
    )
    _api_name: str = PrivateAttr(default="sc.pl.umap")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlTsne(BaseAPI):
    """
    Scatter plot in tSNE basis.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: <class 'anndata._core.anndata.AnnData'>",
        title="Adata",
    )
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes. Original type annotation: str | collections.abc.Sequence[str] | None",
        title="Color",
    )
    mask_obs: Any = Field(
        None,
        description="No description available. Original type annotation: numpy.ndarray[tuple[int, ...], numpy.dtype[numpy.bool]] | str | None",
        title="Mask Obs",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols. Original type annotation: str | None",
        title="Gene Symbols",
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression. Original type annotation: bool | None",
        title="Use Raw",
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points with higher values on top of others. Original type annotation: <class 'bool'>",
        title="Sort Order",
    )
    edges: Optional[Any] = Field(
        False,
        description="Show edges. Original type annotation: <class 'bool'>",
        title="Edges",
    )
    edges_width: Optional[Any] = Field(
        0.1,
        description="Width of edges. Original type annotation: <class 'float'>",
        title="Edges Width",
    )
    edges_color: Optional[Any] = Field(
        "grey",
        description="Color of edges. Original type annotation: str | collections.abc.Sequence[float] | collections.abc.Sequence[str]",
        title="Edges Color",
    )
    neighbors_key: Any = Field(
        None,
        description="Where to look for neighbors connectivities. Original type annotation: str | None",
        title="Neighbors Key",
    )
    arrows: Optional[Any] = Field(
        False,
        description="Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`). Original type annotation: <class 'bool'>",
        title="Arrows",
    )
    arrows_kwds: Any = Field(
        None,
        description="Passed to :meth:`~matplotlib.axes.Axes.quiver`. Original type annotation: collections.abc.Mapping[str, typing.Any] | None",
        title="Arrows Kwds",
    )
    groups: Any = Field(
        None,
        description="Restrict to a few categories in categorical observation annotation. Original type annotation: str | collections.abc.Sequence[str] | None",
        title="Groups",
    )
    components: Any = Field(
        None,
        description="For plotting specific components. Original type annotation: str | collections.abc.Sequence[str] | None",
        title="Components",
    )
    dimensions: Any = Field(
        None,
        description="0-indexed dimensions of the embedding to plot as integers. Original type annotation: tuple[int, int] | collections.abc.Sequence[tuple[int, int]] | None",
        title="Dimensions",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer that wants to be plotted. Original type annotation: str | None",
        title="Layer",
    )
    projection: Optional[Any] = Field(
        "2d",
        description="Projection of plot. Original type annotation: typing.Literal['2d', '3d']",
        title="Projection",
    )
    scale_factor: Any = Field(
        None,
        description="No description available. Original type annotation: float | None",
        title="Scale Factor",
    )
    color_map: Any = Field(
        None,
        description="Color map to use for continous variables. Original type annotation: matplotlib.colors.Colormap | str | None",
        title="Color Map",
    )
    cmap: Any = Field(
        None,
        description="No description available. Original type annotation: matplotlib.colors.Colormap | str | None",
        title="Cmap",
    )
    palette: Any = Field(
        None,
        description="Colors to use for plotting categorical annotation groups. Original type annotation: str | collections.abc.Sequence[str] | cycler.Cycler | None",
        title="Palette",
    )
    na_color: Optional[Any] = Field(
        "lightgray",
        description="Color to use for null or masked values. Original type annotation: str | tuple[float, ...]",
        title="Na Color",
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="Whether missing values get an entry in the legend. Original type annotation: <class 'bool'>",
        title="Na In Legend",
    )
    size: Any = Field(
        None,
        description="Point size. Original type annotation: float | collections.abc.Sequence[float] | None",
        title="Size",
    )
    frameon: Any = Field(
        None,
        description="Draw a frame around the scatter plot. Original type annotation: bool | None",
        title="Frameon",
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string. Original type annotation: typing.Union[float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], NoneType]",
        title="Legend Fontsize",
    )
    legend_fontweight: Optional[Any] = Field(
        "bold",
        description="Legend font weight. Original type annotation: typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]",
        title="Legend Fontweight",
    )
    legend_loc: Optional[Any] = Field(
        "right margin",
        description="Location of legend. Original type annotation: typing.Optional[typing.Literal['none', 'right margin', 'on data', 'on data export', 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center']]",
        title="Legend Loc",
    )
    legend_fontoutline: Any = Field(
        None,
        description="Line width of the legend font outline in pt. Original type annotation: int | None",
        title="Legend Fontoutline",
    )
    colorbar_loc: Optional[Any] = Field(
        "right",
        description="Location of the colorbar. Original type annotation: str | None",
        title="Colorbar Loc",
    )
    vmax: Any = Field(
        None,
        description="The value representing the upper limit of the color scale. Original type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title="Vmax",
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale. Original type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title="Vmin",
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Original type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title="Vcenter",
    )
    norm: Any = Field(
        None,
        description="No description available. Original type annotation: matplotlib.colors.Normalize | collections.abc.Sequence[matplotlib.colors.Normalize] | None",
        title="Norm",
    )
    add_outline: Optional[Any] = Field(
        False,
        description="Add a thin border around groups of dots. Original type annotation: bool | None",
        title="Add Outline",
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description="Adjust the outline width. Original type annotation: tuple[float, float]",
        title="Outline Width",
    )
    outline_color: Optional[Any] = Field(
        ["black", "white"],
        description="Adjust the outline color. Original type annotation: tuple[str, str]",
        title="Outline Color",
    )
    ncols: Optional[Any] = Field(
        4,
        description="Number of panels per row. Original type annotation: <class 'int'>",
        title="Ncols",
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels. Original type annotation: <class 'float'>",
        title="Hspace",
    )
    wspace: Any = Field(
        None,
        description="Adjust the width of the space between multiple panels. Original type annotation: float | None",
        title="Wspace",
    )
    title: Any = Field(
        None,
        description="Provide title for panels. Original type annotation: str | collections.abc.Sequence[str] | None",
        title="Title",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="Save the figure. Original type annotation: bool | str | None",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Original type annotation: matplotlib.axes._axes.Axes | None",
        title="Ax",
    )
    return_fig: Any = Field(
        None,
        description="Return the matplotlib figure. Original type annotation: bool | None",
        title="Return Fig",
    )
    marker: Optional[Any] = Field(
        ".",
        description="No description available. Original type annotation: str | collections.abc.Sequence[str]",
        title="Marker",
    )
    _api_name: str = PrivateAttr(default="sc.pl.tsne")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlHeatmap(BaseAPI):
    """
    Heatmap of the expression values of genes. If `groupby` is given, the heatmap is ordered by the respective group. For example, a list of marker genes can be plotted, ordered by clustering. If the `groupby` observation annotation is not categorical the observation annotation is turned into a categorical by binning the data into the number specified in `num_categories`.
    """

    adata: Any = Field(..., description="Annotated data matrix.", title="Adata")
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names` and can include grouping information.",
        title="Var Names",
    )
    groupby: Any = Field(
        ...,
        description="The key of the observation grouping to consider.",
        title="Groupby",
    )
    use_raw: Any = Field(
        None,
        description="Flag to indicate whether to use the `raw` attribute of `adata` if present.",
        title="Use Raw",
    )
    log: Optional[Any] = Field(
        False, description="Option to plot on a logarithmic axis.", title="Log"
    )
    num_categories: Optional[Any] = Field(
        7,
        description="Number of groups into which the groupby observation should be subdivided if not categorical.",
        title="Num Categories",
    )
    dendrogram: Optional[Any] = Field(
        False,
        description="Option to include a dendrogram based on hierarchical clustering between groupby categories.",
        title="Dendrogram",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
        title="Gene Symbols",
    )
    var_group_positions: Any = Field(
        None,
        description="Parameter to highlight groups of `var_names` using brackets or color blocks.",
        title="Var Group Positions",
    )
    var_group_labels: Any = Field(
        None,
        description="Labels for each of the `var_group_positions` that want to be highlighted.",
        title="Var Group Labels",
    )
    var_group_rotation: Any = Field(
        None, description="Rotation degrees for labels.", title="Var Group Rotation"
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer to be plotted.",
        title="Layer",
    )
    standard_scale: Any = Field(
        None,
        description="Option to standardize dimensions between 0 and 1.",
        title="Standard Scale",
    )
    swap_axes: Optional[Any] = Field(
        False, description="Option to swap x and y axes in the plot.", title="Swap Axes"
    )
    show_gene_labels: Any = Field(
        None,
        description="Toggle to show or remove gene labels based on the number of genes.",
        title="Show Gene Labels",
    )
    show: Any = Field(
        None,
        description="Option to show the plot without returning axis.",
        title="Show",
    )
    save: Any = Field(
        None,
        description="Option to save the figure with a specified filename and filetype.",
        title="Save",
    )
    figsize: Any = Field(
        None,
        description="Figure size parameter when `multi_panel=True`.",
        title="Figsize",
    )
    vmin: Any = Field(
        None, description="The lower limit of the color scale.", title="Vmin"
    )
    vmax: Any = Field(
        None, description="The upper limit of the color scale.", title="Vmax"
    )
    vcenter: Any = Field(
        None,
        description="The center of the color scale, useful for diverging colormaps.",
        title="Vcenter",
    )
    norm: Any = Field(
        None,
        description="Custom color normalization object from matplotlib.",
        title="Norm",
    )
    kwds: Any = Field(..., description="Additional keyword arguments.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.heatmap")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlDotplot(BaseAPI):
    """
    Make a dot plot of the expression values of var_names. For each var_name and each groupby category a dot is plotted. Each dot represents two values: mean expression within each category (visualized by color) and fraction of cells expressing the var_name in the category (visualized by the size of the dot). If groupby is not given, the dotplot assumes that all data belongs to a single category. A gene is considered expressed if the expression value in the adata (or adata.raw) is above the specified threshold which is zero by default. An example of dotplot usage is to visualize, for multiple marker genes, the mean value and the percentage of cells expressing the gene across multiple clusters. This function provides a convenient interface to the scanpy.pl.DotPlot class. If you need more flexibility, you should use scanpy.pl.DotPlot directly.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`. If `var_names` is a mapping, then the key is used to group values. Original type annotation: _VarNames | Mapping[str, _VarNames]",
        title="Var Names",
    )
    groupby: Any = Field(
        ...,
        description="The key of the observation grouping to consider. Original type annotation: str | Sequence[str]",
        title="Groupby",
    )
    use_raw: Any = Field(
        None,
        description="Use `raw` attribute of `adata` if present. Original type annotation: bool | None",
        title="Use Raw",
    )
    log: Optional[Any] = Field(
        False,
        description="Plot on a logarithmic axis. Original type annotation: bool",
        title="Log",
    )
    num_categories: Optional[Any] = Field(
        7,
        description="Determines the number of groups when groupby observation is not categorical. Original type annotation: int",
        title="Num Categories",
    )
    categories_order: Any = Field(
        None,
        description="Order in which to show the categories. Original type annotation: Sequence[str] | None",
        title="Categories Order",
    )
    expression_cutoff: Optional[Any] = Field(
        0.0,
        description="Threshold for gene expression to determine cell expression. Original type annotation: float",
        title="Expression Cutoff",
    )
    mean_only_expressed: Optional[Any] = Field(
        False,
        description="If True, average gene expression only over expressing cells. Original type annotation: bool",
        title="Mean Only Expressed",
    )
    standard_scale: Any = Field(
        None,
        description="Standardize dimension between 0 and 1. Original type annotation: Literal['var', 'group'] | None",
        title="Standard Scale",
    )
    title: Any = Field(
        None,
        description="Title for the figure. Original type annotation: str | None",
        title="Title",
    )
    colorbar_title: Optional[Any] = Field(
        "Mean expression\nin group",
        description="Title for the color bar. Original type annotation: str | None",
        title="Colorbar Title",
    )
    size_title: Optional[Any] = Field(
        "Fraction of cells\nin group (%)",
        description="Title for the size legend. Original type annotation: str | None",
        title="Size Title",
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`. Original type annotation: tuple[float, float] | None",
        title="Figsize",
    )
    dendrogram: Optional[Any] = Field(
        False,
        description="Add dendrogram based on hierarchical clustering between groupby categories. Original type annotation: bool | str",
        title="Dendrogram",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols. Original type annotation: str | None",
        title="Gene Symbols",
    )
    var_group_positions: Any = Field(
        None,
        description="Highlight groups of `var_names` with brackets or color blocks. Original type annotation: Sequence[tuple[int, int]] | None",
        title="Var Group Positions",
    )
    var_group_labels: Any = Field(
        None,
        description="Labels for highlighted `var_group_positions`. Original type annotation: Sequence[str] | None",
        title="Var Group Labels",
    )
    var_group_rotation: Any = Field(
        None,
        description="Label rotation degrees. Original type annotation: float | None",
        title="Var Group Rotation",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer to be plotted. Original type annotation: str | None",
        title="Layer",
    )
    swap_axes: Optional[Any] = Field(
        False,
        description="Swap x and y axes. Original type annotation: bool | None",
        title="Swap Axes",
    )
    dot_color_df: Any = Field(
        None,
        description="Dataframe for dot colors. Original type annotation: pd.DataFrame | None",
        title="Dot Color Df",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="Save the figure with option for filetype. Original type annotation: str | bool | None",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="Matplotlib axes object. Original type annotation: _AxesSubplot | None",
        title="Ax",
    )
    return_fig: Optional[Any] = Field(
        False,
        description="Return DotPlot object. Original type annotation: bool | None",
        title="Return Fig",
    )
    vmin: Any = Field(
        None,
        description="Lower limit of color scale. Original type annotation: float | None",
        title="Vmin",
    )
    vmax: Any = Field(
        None,
        description="Upper limit of color scale. Original type annotation: float | None",
        title="Vmax",
    )
    vcenter: Any = Field(
        None,
        description="Center of color scale. Original type annotation: float | None",
        title="Vcenter",
    )
    norm: Any = Field(
        None,
        description="Custom color normalization object from matplotlib. Original type annotation: Normalize | None",
        title="Norm",
    )
    cmap: Optional[Any] = Field(
        "Reds",
        description="Matplotlib color map. Original type annotation: Colormap | str | None",
        title="Cmap",
    )
    dot_max: Any = Field(
        None,
        description="Maximum dot size based on fraction value. Original type annotation: float | None",
        title="Dot Max",
    )
    dot_min: Any = Field(
        None,
        description="Minimum dot size. Original type annotation: float | None",
        title="Dot Min",
    )
    smallest_dot: Optional[Any] = Field(
        0.0,
        description="Size for expression levels at `dot_min`. Original type annotation: float",
        title="Smallest Dot",
    )
    kwds: Any = Field(
        ..., description="Parameters passed to matplotlib.pyplot.scatter", title="Kwds"
    )
    _api_name: str = PrivateAttr(default="sc.pl.dotplot")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlViolin(BaseAPI):
    """
    Violin plot. Wraps seaborn.violinplot for anndata.AnnData.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    keys: Any = Field(
        ...,
        description="Keys for accessing variables of `.var_names` or fields of `.obs`. Original type annotation: str | Sequence[str]",
        title="Keys",
    )
    groupby: Any = Field(
        None,
        description="The key of the observation grouping to consider. Original type annotation: str | None",
        title="Groupby",
    )
    log: Optional[Any] = Field(
        False,
        description="Plot on logarithmic axis. Original type annotation: bool",
        title="Log",
    )
    use_raw: Any = Field(
        None,
        description="Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present. Original type annotation: bool | None",
        title="Use Raw",
    )
    stripplot: Optional[Any] = Field(
        True,
        description="Add a stripplot on top of the violin plot. See :func:`~seaborn.stripplot`. Original type annotation: bool",
        title="Stripplot",
    )
    jitter: Optional[Any] = Field(
        True,
        description="Add jitter to the stripplot (only when stripplot is True). See :func:`~seaborn.stripplot`. Original type annotation: float | bool",
        title="Jitter",
    )
    size: Optional[Any] = Field(
        1,
        description="Size of the jitter points. Original type annotation: int",
        title="Size",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name, then the layer is plotted. `layer` takes precedence over `use_raw`. Original type annotation: str | None",
        title="Layer",
    )
    density_norm: Optional[Any] = Field(
        "width",
        description="The method used to scale the width of each violin. If 'width' (the default), each violin will have the same width. If 'area', each violin will have the same area. If 'count', a violin’s width corresponds to the number of observations. Original type annotation: DensityNorm",
        title="Density Norm",
    )
    order: Any = Field(
        None,
        description="Order in which to show the categories. Original type annotation: Sequence[str] | None",
        title="Order",
    )
    multi_panel: Any = Field(
        None,
        description="Display keys in multiple panels also when `groupby is not None`. Original type annotation: bool | None",
        title="Multi Panel",
    )
    xlabel: Optional[Any] = Field(
        "",
        description="Label of the x axis. Defaults to `groupby` if `rotation` is `None`, otherwise, no label is shown. Original type annotation: str",
        title="Xlabel",
    )
    ylabel: Any = Field(
        None,
        description="Label of the y axis. If `None` and `groupby` is `None`, defaults to `'value'`. If `None` and `groubpy` is not `None`, defaults to `keys`. Original type annotation: str | Sequence[str] | None",
        title="Ylabel",
    )
    rotation: Any = Field(
        None,
        description="Rotation of xtick labels. Original type annotation: float | None",
        title="Rotation",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: bool | str | None",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None",
        title="Ax",
    )
    scale: Optional[Any] = Field(
        0,
        description="No description available. Original type annotation: DensityNorm | Empty",
        title="Scale",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.violin")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlDendrogram(BaseAPI):
    """
    Plot a dendrogram of the categories defined in `groupby`. See :func:`~scanpy.tl.dendrogram`.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    groupby: Any = Field(
        ...,
        description="Categorical data column used to create the dendrogram. Original type annotation: str",
        title="Groupby",
    )
    dendrogram_key: Any = Field(
        None,
        description="Key under which the dendrogram information was stored. By default, the dendrogram information is stored under `.uns[f'dendrogram_{groupby}']`. Original type annotation: str | None",
        title="Dendrogram Key",
    )
    orientation: Optional[Any] = Field(
        "top",
        description="Origin of the tree. Will grow into the opposite direction. Original type annotation: Literal['top', 'bottom', 'left', 'right']",
        title="Orientation",
    )
    remove_labels: Optional[Any] = Field(
        False,
        description="Don’t draw labels. Used by functions like scanpy.pl.matrixplot to annotate matrix columns/rows. Original type annotation: bool",
        title="Remove Labels",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If True or a string, save the figure. A string is appended to the default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}. Original type annotation: str | bool | None",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None",
        title="Ax",
    )
    _api_name: str = PrivateAttr(default="sc.pl.dendrogram")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlDiffmap(BaseAPI):
    """
    Scatter plot in Diffusion Map basis.
    """

    adata: Any = Field(..., description="Annotated data matrix.", title="Adata")
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes.",
        title="Color",
    )
    mask_obs: Any = Field(
        None, description="No description available.", title="Mask Obs"
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
        title="Gene Symbols",
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
        title="Use Raw",
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points with higher values on top of others.",
        title="Sort Order",
    )
    edges: Optional[Any] = Field(
        False, description="No description available.", title="Edges"
    )
    edges_width: Optional[Any] = Field(
        0.1, description="No description available.", title="Edges Width"
    )
    edges_color: Optional[Any] = Field(
        "grey", description="No description available.", title="Edges Color"
    )
    neighbors_key: Any = Field(
        None, description="No description available.", title="Neighbors Key"
    )
    arrows: Optional[Any] = Field(
        False, description="No description available.", title="Arrows"
    )
    arrows_kwds: Any = Field(
        None, description="No description available.", title="Arrows Kwds"
    )
    groups: Any = Field(
        None,
        description="Restrict to a few categories in categorical observation annotation.",
        title="Groups",
    )
    components: Any = Field(
        None,
        description="For plotting available components or a subset specified by the user.",
        title="Components",
    )
    dimensions: Any = Field(
        None,
        description="0-indexed dimensions of the embedding to plot as integers.",
        title="Dimensions",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer that wants to be plotted.",
        title="Layer",
    )
    projection: Optional[Any] = Field(
        "2d", description="Projection of plot.", title="Projection"
    )
    scale_factor: Any = Field(
        None, description="No description available.", title="Scale Factor"
    )
    color_map: Any = Field(
        None,
        description="Color map to use for continuous variables.",
        title="Color Map",
    )
    cmap: Any = Field(None, description="No description available.", title="Cmap")
    palette: Any = Field(
        None,
        description="Colors to use for plotting categorical annotation groups.",
        title="Palette",
    )
    na_color: Optional[Any] = Field(
        "lightgray",
        description="Color to use for null or masked values.",
        title="Na Color",
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="Whether missing values get an entry in the legend.",
        title="Na In Legend",
    )
    size: Any = Field(None, description="Point size for each cell.", title="Size")
    frameon: Any = Field(
        None, description="Draw a frame around the scatter plot.", title="Frameon"
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string describing the size.",
        title="Legend Fontsize",
    )
    legend_fontweight: Optional[Any] = Field(
        "bold", description="Legend font weight.", title="Legend Fontweight"
    )
    legend_loc: Optional[Any] = Field(
        "right margin", description="Location of legend.", title="Legend Loc"
    )
    legend_fontoutline: Any = Field(
        None,
        description="Line width of the legend font outline.",
        title="Legend Fontoutline",
    )
    colorbar_loc: Optional[Any] = Field(
        "right",
        description="Location to place the colorbar for continuous variables.",
        title="Colorbar Loc",
    )
    vmax: Any = Field(
        None,
        description="The value representing the upper limit of the color scale.",
        title="Vmax",
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale.",
        title="Vmin",
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale.",
        title="Vcenter",
    )
    norm: Any = Field(None, description="No description available.", title="Norm")
    add_outline: Optional[Any] = Field(
        False,
        description="If set to True, this will add a thin border around groups of dots.",
        title="Add Outline",
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description="Width numbers used to adjust the outline.",
        title="Outline Width",
    )
    outline_color: Optional[Any] = Field(
        ["black", "white"],
        description="Color names used to adjust the add_outline.",
        title="Outline Color",
    )
    ncols: Optional[Any] = Field(
        4, description="Number of panels per row.", title="Ncols"
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.",
        title="Hspace",
    )
    wspace: Any = Field(
        None,
        description="Adjust the width of the space between multiple panels.",
        title="Wspace",
    )
    title: Any = Field(None, description="Provide title for panels.", title="Title")
    show: Any = Field(
        None, description="Show the plot, do not return axis.", title="Show"
    )
    save: Any = Field(
        None, description="If `True` or a `str`, save the figure.", title="Save"
    )
    ax: Any = Field(None, description="A matplotlib axes object.", title="Ax")
    return_fig: Any = Field(
        None, description="Return the matplotlib figure.", title="Return Fig"
    )
    marker: Optional[Any] = Field(
        ".", description="No description available.", title="Marker"
    )
    _api_name: str = PrivateAttr(default="sc.pl.diffmap")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlHighlyVariableGenes(BaseAPI):
    """
    Plot dispersions or normalized variance versus means for genes. Produces Supp. Fig. 5c of Zheng et al. (2017) and MeanVarPlot() and VariableFeaturePlot() of Seurat.
    """

    adata_or_result: Any = Field(
        ...,
        description="Adata Or Result can be either AnnData, pd.DataFrame, or np.recarray.",
        title="Adata Or Result",
    )
    log: Optional[Any] = Field(
        False, description="Plot on logarithmic axes if set to true.", title="Log"
    )
    show: Any = Field(
        None,
        description="If set to true, show the plot without returning axis; can also be set to None.",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If set to true or a string, save the figure with the specified filename; the filetype is inferred from the extension ('.pdf', '.png', '.svg').",
        title="Save",
    )
    highly_variable_genes: Optional[Any] = Field(
        True,
        description="Boolean value indicating whether the genes are highly variable or not.",
        title="Highly Variable Genes",
    )
    _api_name: str = PrivateAttr(default="sc.pl.highly_variable_genes")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlPca(BaseAPI):
    """
    Scatter plot in PCA coordinates. Use the parameter `annotate_var_explained` to annotate the explained variance.
    """

    adata: Any = Field(..., description="Annotated data matrix.", title="Adata")
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes.",
        title="Color",
    )
    mask_obs: Any = Field(
        None, description="No description available.", title="Mask Obs"
    )
    gene_symbols: Any = Field(
        None,
        description="Column name that stores gene symbols in `.var` DataFrame.",
        title="Gene Symbols",
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
        title="Use Raw",
    )
    sort_order: Optional[Any] = Field(
        True,
        description="Plot data points with higher values on top for continuous annotations.",
        title="Sort Order",
    )
    edges: Optional[Any] = Field(
        False, description="No description available.", title="Edges"
    )
    edges_width: Optional[Any] = Field(
        0.1, description="No description available.", title="Edges Width"
    )
    edges_color: Optional[Any] = Field(
        "grey", description="No description available.", title="Edges Color"
    )
    neighbors_key: Any = Field(
        None, description="No description available.", title="Neighbors Key"
    )
    arrows: Optional[Any] = Field(
        False, description="No description available.", title="Arrows"
    )
    arrows_kwds: Any = Field(
        None, description="No description available.", title="Arrows Kwds"
    )
    groups: Any = Field(
        None,
        description="Restrict to a few categories in categorical observation annotation.",
        title="Groups",
    )
    components: Any = Field(
        None, description="Plot specific components.", title="Components"
    )
    dimensions: Any = Field(
        None,
        description="Dimensions of the embedding to plot as integers.",
        title="Dimensions",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer to be plotted.",
        title="Layer",
    )
    projection: Optional[Any] = Field(
        "2d", description="Projection of plot (default: '2d').", title="Projection"
    )
    scale_factor: Any = Field(
        None, description="No description available.", title="Scale Factor"
    )
    color_map: Any = Field(
        None,
        description="Color map to use for continuous variables.",
        title="Color Map",
    )
    cmap: Any = Field(None, description="No description available.", title="Cmap")
    palette: Any = Field(
        None,
        description="Colors to use for plotting categorical annotation groups.",
        title="Palette",
    )
    na_color: Optional[Any] = Field(
        "lightgray",
        description="Color to use for null or masked values.",
        title="Na Color",
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="Whether missing values get an entry in the legend.",
        title="Na In Legend",
    )
    size: Any = Field(None, description="Point size for plotting.", title="Size")
    frameon: Any = Field(
        None, description="Draw a frame around the scatter plot.", title="Frameon"
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size or string describing the size of the legend font.",
        title="Legend Fontsize",
    )
    legend_fontweight: Optional[Any] = Field(
        "bold", description="Legend font weight.", title="Legend Fontweight"
    )
    legend_loc: Optional[Any] = Field(
        "right margin", description="Location of legend.", title="Legend Loc"
    )
    legend_fontoutline: Any = Field(
        None,
        description="Line width of the legend font outline.",
        title="Legend Fontoutline",
    )
    colorbar_loc: Optional[Any] = Field(
        "right",
        description="Location to place the colorbar for continuous variables.",
        title="Colorbar Loc",
    )
    vmax: Any = Field(
        None, description="The upper limit of the color scale.", title="Vmax"
    )
    vmin: Any = Field(
        None, description="The lower limit of the color scale.", title="Vmin"
    )
    vcenter: Any = Field(
        None, description="The center of the color scale.", title="Vcenter"
    )
    norm: Any = Field(None, description="No description available.", title="Norm")
    add_outline: Optional[Any] = Field(
        False, description="Add a border around groups of dots.", title="Add Outline"
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description="Width numbers used to adjust the outline.",
        title="Outline Width",
    )
    outline_color: Optional[Any] = Field(
        ["black", "white"],
        description="Color names used to adjust the outline.",
        title="Outline Color",
    )
    ncols: Optional[Any] = Field(
        4, description="Number of panels per row.", title="Ncols"
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.",
        title="Hspace",
    )
    wspace: Any = Field(
        None,
        description="Adjust the width of the space between multiple panels.",
        title="Wspace",
    )
    title: Any = Field(None, description="Title for panels.", title="Title")
    show: Any = Field(
        None, description="Show the plot, do not return axis.", title="Show"
    )
    save: Any = Field(None, description="Save the figure.", title="Save")
    ax: Any = Field(None, description="Matplotlib axes object.", title="Ax")
    return_fig: Any = Field(
        None, description="Return the matplotlib figure.", title="Return Fig"
    )
    marker: Optional[Any] = Field(
        ".", description="No description available.", title="Marker"
    )
    annotate_var_explained: Optional[Any] = Field(
        False, description="No description available.", title="Annotate Var Explained"
    )
    _api_name: str = PrivateAttr(default="sc.pl.pca")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlEmbeddingDensity(BaseAPI):
    """
    Plot the density of cells in an embedding (per condition). Plots the gaussian kernel density estimates (over condition) from the sc.tl.embedding_density() output. This function was written by Sophie Tritschler and implemented into Scanpy by Malte Luecken.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    basis: Optional[Any] = Field(
        "umap",
        description="The embedding over which the density was calculated. This embedded representation should be found in `adata.obsm['X_[basis]']`. Original type annotation: str",
        title="Basis",
    )
    key: Any = Field(
        None,
        description="Name of the `.obs` covariate that contains the density estimates. Alternatively, pass `groupby`. Original type annotation: str | None",
        title="Key",
    )
    groupby: Any = Field(
        None,
        description="Name of the condition used in `tl.embedding_density`. Alternatively, pass `key`. Original type annotation: str | None",
        title="Groupby",
    )
    group: Optional[Any] = Field(
        "all",
        description="The category in the categorical observation annotation to be plotted. For example, 'G1' in the cell cycle 'phase' covariate. If all categories are to be plotted use group='all' (default), If multiple categories want to be plotted use a list (e.g.: ['G1', 'S']. If the overall density wants to be ploted set group to 'None'. Original type annotation: str | Sequence[str] | None",
        title="Group",
    )
    color_map: Optional[Any] = Field(
        "YlOrRd",
        description="Matplolib color map to use for density plotting. Original type annotation: Colormap | str",
        title="Color Map",
    )
    bg_dotsize: Optional[Any] = Field(
        80,
        description="Dot size for background data points not in the `group`. Original type annotation: int | None",
        title="Bg Dotsize",
    )
    fg_dotsize: Optional[Any] = Field(
        180,
        description="Dot size for foreground data points in the `group`. Original type annotation: int | None",
        title="Fg Dotsize",
    )
    vmax: Optional[Any] = Field(
        1,
        description="The value representing the upper limit of the color scale. The format is the same as for `vmin`. Original type annotation: int | None",
        title="Vmax",
    )
    vmin: Optional[Any] = Field(
        0,
        description="The value representing the lower limit of the color scale. Values smaller than vmin are plotted with the same color as vmin. vmin can be a number, a string, a function or `None`. If vmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N). For example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then vmin is interpreted as the return value of the function over the list of values to plot. For example to set vmin to the mean of the values to plot, `def my_vmin(values): return np.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic minimum value is used as defined by matplotlib `scatter` function. When making multiple plots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`. Original type annotation: int | None",
        title="Vmin",
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Useful for diverging colormaps. The format is the same as for `vmin`. Example: `sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')`. Original type annotation: int | None",
        title="Vcenter",
    )
    norm: Any = Field(
        None,
        description="No description available. Original type annotation: Normalize | None",
        title="Norm",
    )
    ncols: Optional[Any] = Field(
        4,
        description="Number of panels per row. Original type annotation: int | None",
        title="Ncols",
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels. Original type annotation: float | None",
        title="Hspace",
    )
    wspace: Any = Field(
        None,
        description="Adjust the width of the space between multiple panels. Original type annotation: None",
        title="Wspace",
    )
    title: Any = Field(
        None,
        description="No description available. Original type annotation: str | None",
        title="Title",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}. Original type annotation: bool | str | None",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None",
        title="Ax",
    )
    return_fig: Any = Field(
        None,
        description="Return the matplotlib figure. Original type annotation: bool | None",
        title="Return Fig",
    )
    _api_name: str = PrivateAttr(default="sc.pl.embedding_density")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlRankGenesGroups(BaseAPI):
    """
    Plot ranking of genes.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix with original type annotation: AnnData.",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="The groups for which to show the gene ranking with original type annotation: str | Sequence[str] | None.",
        title="Groups",
    )
    n_genes: Optional[Any] = Field(
        20,
        description="Number of genes to show with original type annotation: int.",
        title="N Genes",
    )
    gene_symbols: Any = Field(
        None,
        description="Key for field in `.var` that stores gene symbols if you do not want to use `.var_names` with original type annotation: str | None.",
        title="Gene Symbols",
    )
    key: Optional[Any] = Field(
        "rank_genes_groups",
        description="No description available with original type annotation: str | None.",
        title="Key",
    )
    fontsize: Optional[Any] = Field(
        8,
        description="Fontsize for gene names with original type annotation: int.",
        title="Fontsize",
    )
    ncols: Optional[Any] = Field(
        4,
        description="Number of panels shown per row with original type annotation: int.",
        title="Ncols",
    )
    sharey: Optional[Any] = Field(
        True,
        description="Controls if the y-axis of each panel should be shared, with original type annotation: bool.",
        title="Sharey",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis with original type annotation: bool | None.",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`} with original type annotation: bool | None.",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Only works if plotting a single component with original type annotation: Axes | None.",
        title="Ax",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.rank_genes_groups")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlRankGenesGroupsDotplot(BaseAPI):
    """
    Plot ranking of genes using dotplot plot (see :func:`~scanpy.pl.dotplot`).
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="The groups for which to show the gene ranking. Original type annotation: str | Sequence[str] | None",
        title="Groups",
    )
    n_genes: Any = Field(
        None,
        description="Number of genes to show. This can be a negative number to show for example the down regulated genes. Original type annotation: int | None",
        title="N Genes",
    )
    groupby: Any = Field(
        None,
        description="The key of the observation grouping to consider. By default, the groupby is chosen from the rank genes groups parameter but other groupby options can be used. Original type annotation: str | None",
        title="Groupby",
    )
    values_to_plot: Any = Field(
        None,
        description="Instead of the mean gene value, plot the values computed by `sc.rank_genes_groups`. The options are: ['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj']. Original type annotation: Literal['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'] | None",
        title="Values To Plot",
    )
    var_names: Any = Field(
        None,
        description="Genes to plot. Sometimes is useful to pass a specific list of var names (e.g. genes) to check their fold changes or p-values, instead of the top/bottom genes. Original type annotation: Sequence[str] | Mapping[str, Sequence[str]] | None",
        title="Var Names",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols. By default `var_names` refer to the index column of the `.var` DataFrame. Original type annotation: str | None",
        title="Gene Symbols",
    )
    min_logfoldchange: Any = Field(
        None,
        description="Value to filter genes in groups if their logfoldchange is less than the min_logfoldchange. Original type annotation: float | None",
        title="Min Logfoldchange",
    )
    key: Any = Field(
        None,
        description="Key used to store the ranking results in `adata.uns`. Original type annotation: str | None",
        title="Key",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: bool | None",
        title="Save",
    )
    return_fig: Optional[Any] = Field(
        False,
        description="Returns :class:`DotPlot` object. Useful for fine-tuning the plot. Takes precedence over `show=False`. Original type annotation: bool",
        title="Return Fig",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.rank_genes_groups_dotplot")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlRankGenesGroupsViolin(BaseAPI):
    """
    Plot ranking of genes for all tested comparisons.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="List of group names. Original type annotation: Sequence[str] | None",
        title="Groups",
    )
    n_genes: Optional[Any] = Field(
        20,
        description="Number of genes to show. Is ignored if `gene_names` is passed. Original type annotation: int",
        title="N Genes",
    )
    gene_names: Any = Field(
        None,
        description="List of genes to plot. Is only useful if interested in a custom gene list, which is not the result of :func:`scanpy.tl.rank_genes_groups`. Original type annotation: Iterable[str] | None",
        title="Gene Names",
    )
    gene_symbols: Any = Field(
        None,
        description="Key for field in `.var` that stores gene symbols if you do not want to use `.var_names` displayed in the plot. Original type annotation: str | None",
        title="Gene Symbols",
    )
    use_raw: Any = Field(
        None,
        description="Use `raw` attribute of `adata` if present. Defaults to the value that was used in :func:`~scanpy.tl.rank_genes_groups`. Original type annotation: bool | None",
        title="Use Raw",
    )
    key: Any = Field(
        None,
        description="No description available. Original type annotation: str | None",
        title="Key",
    )
    split: Optional[Any] = Field(
        True,
        description="Whether to split the violins or not. Original type annotation: bool",
        title="Split",
    )
    density_norm: Optional[Any] = Field(
        "width",
        description="See :func:`~seaborn.violinplot`. Original type annotation: DensityNorm",
        title="Density Norm",
    )
    strip: Optional[Any] = Field(
        True,
        description="Show a strip plot on top of the violin plot. Original type annotation: bool",
        title="Strip",
    )
    jitter: Optional[Any] = Field(
        True,
        description="If set to 0, no points are drawn. See :func:`~seaborn.stripplot`. Original type annotation: float | bool",
        title="Jitter",
    )
    size: Optional[Any] = Field(
        1,
        description="Size of the jitter points. Original type annotation: int",
        title="Size",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None",
        title="Ax",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: bool | None",
        title="Save",
    )
    scale: Optional[Any] = Field(
        0,
        description="No description available. Original type annotation: DensityNorm | Empty",
        title="Scale",
    )
    _api_name: str = PrivateAttr(default="sc.pl.rank_genes_groups_violin")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlRankGenesGroupsHeatmap(BaseAPI):
    """
    Plot ranking of genes using heatmap plot (see :func:`~scanpy.pl.heatmap`).
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix with original type annotation as AnnData.",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="The groups for which to show the gene ranking; can be a string, a sequence of strings, or None.",
        title="Groups",
    )
    n_genes: Any = Field(
        None,
        description="Number of genes to display, which can be negative to show down-regulated genes; can be an integer or None.",
        title="N Genes",
    )
    groupby: Any = Field(
        None,
        description="The key of the observation grouping to consider, expected to be a categorical string or None.",
        title="Groupby",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in the .var DataFrame that contains gene symbols, allowing alternative names to be used.",
        title="Gene Symbols",
    )
    var_names: Any = Field(
        None,
        description="Sequence of strings or mapping in the .var DataFrame; no detailed description available.",
        title="Var Names",
    )
    min_logfoldchange: Any = Field(
        None,
        description="Value used to filter genes in groups based on logfoldchange; can be a float or None.",
        title="Min Logfoldchange",
    )
    key: Any = Field(
        None,
        description="Key to store the ranking results in adata.uns; can be a string or None.",
        title="Key",
    )
    show: Any = Field(
        None,
        description="Option to display the plot without returning the axis; can be a boolean or None.",
        title="Show",
    )
    save: Any = Field(
        None,
        description="Option to save the figure, where a string is appended to the default filename; can be a boolean or None.",
        title="Save",
    )
    kwds: Any = Field(
        ..., description="No detailed description available.", title="Kwds"
    )
    _api_name: str = PrivateAttr(default="sc.pl.rank_genes_groups_heatmap")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlRankGenesGroupsStackedViolin(BaseAPI):
    """
    Plot ranking of genes using stacked_violin plot.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="The groups for which to show the gene ranking. Original type annotation: str | Sequence[str] | None",
        title="Groups",
    )
    n_genes: Any = Field(
        None,
        description="Number of genes to show. This can be a negative number to show for example the down regulated genes. Is ignored if `gene_names` is passed. Original type annotation: int | None",
        title="N Genes",
    )
    groupby: Any = Field(
        None,
        description="The key of the observation grouping to consider. It is expected that groupby is a categorical. Original type annotation: str | None",
        title="Groupby",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols. Original type annotation: str | None",
        title="Gene Symbols",
    )
    var_names: Any = Field(
        None,
        description="No description available. Original type annotation: Sequence[str] | Mapping[str, Sequence[str]] | None",
        title="Var Names",
    )
    min_logfoldchange: Any = Field(
        None,
        description="Value to filter genes in groups if their logfoldchange is less than the min_logfoldchange. Original type annotation: float | None",
        title="Min Logfoldchange",
    )
    key: Any = Field(
        None,
        description="Key used to store the ranking results in `adata.uns`. Original type annotation: str | None",
        title="Key",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. Infer the filetype if ending on {'.pdf', '.png', '.svg'}. Original type annotation: bool | None",
        title="Save",
    )
    return_fig: Optional[Any] = Field(
        False,
        description="Returns :class:`StackedViolin` object. Useful for fine-tuning the plot. Original type annotation: bool",
        title="Return Fig",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.rank_genes_groups_stacked_violin")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlRankGenesGroupsMatrixplot(BaseAPI):
    """
    Plot ranking of genes using matrixplot plot (see :func:`~scanpy.pl.matrixplot`).
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="The groups for which to show the gene ranking. Original type annotation: str | Sequence[str] | None",
        title="Groups",
    )
    n_genes: Any = Field(
        None,
        description="Number of genes to show. This can be a negative number to show for example the down regulated genes. Original type annotation: int | None",
        title="N Genes",
    )
    groupby: Any = Field(
        None,
        description="The key of the observation grouping to consider. By default, the groupby is chosen from the rank genes groups parameter but other groupby options can be used. Original type annotation: str | None",
        title="Groupby",
    )
    values_to_plot: Any = Field(
        None,
        description="Instead of the mean gene value, plot the values computed by sc.rank_genes_groups. Original type annotation: Literal['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'] | None",
        title="Values To Plot",
    )
    var_names: Any = Field(
        None,
        description="Genes to plot. Sometimes is useful to pass a specific list of var names (e.g. genes) to check their fold changes or p-values. Original type annotation: Sequence[str] | Mapping[str, Sequence[str]] | None",
        title="Var Names",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in .var DataFrame that stores gene symbols. By default var_names refer to the index column of the .var DataFrame. Original type annotation: str | None",
        title="Gene Symbols",
    )
    min_logfoldchange: Any = Field(
        None,
        description="Value to filter genes in groups if their logfoldchange is less than the min_logfoldchange. Original type annotation: float | None",
        title="Min Logfoldchange",
    )
    key: Any = Field(
        None,
        description="Key used to store the ranking results in adata.uns. Original type annotation: str | None",
        title="Key",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If True or a str, save the figure. A string is appended to the default filename. Original type annotation: bool | None",
        title="Save",
    )
    return_fig: Optional[Any] = Field(
        False,
        description="Returns MatrixPlot object. Useful for fine-tuning the plot. Takes precedence over show=False. Original type annotation: bool",
        title="Return Fig",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.rank_genes_groups_matrixplot")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlRankGenesGroupsTracksplot(BaseAPI):
    """
    Plot ranking of genes using heatmap plot (see :func:`~scanpy.pl.heatmap`).
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="The groups for which to show the gene ranking. Original type annotation: str | Sequence[str] | None",
        title="Groups",
    )
    n_genes: Any = Field(
        None,
        description="Number of genes to show. This can be a negative number to show for example the down regulated genes. Is ignored if gene_names is passed. Original type annotation: int | None",
        title="N Genes",
    )
    groupby: Any = Field(
        None,
        description="The key of the observation grouping to consider. By default, the groupby is chosen from the rank genes groups parameter but other groupby options can be used. Expected to be a categorical. Original type annotation: str | None",
        title="Groupby",
    )
    var_names: Any = Field(
        None,
        description="No description available. Original type annotation: Sequence[str] | Mapping[str, Sequence[str]] | None",
        title="Var Names",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in .var DataFrame that stores gene symbols. By default var_names refer to the index column of the .var DataFrame. Original type annotation: str | None",
        title="Gene Symbols",
    )
    min_logfoldchange: Any = Field(
        None,
        description="Value to filter genes in groups if their logfoldchange is less than the min_logfoldchange. Original type annotation: float | None",
        title="Min Logfoldchange",
    )
    key: Any = Field(
        None,
        description="Key used to store the ranking results in adata.uns. Original type annotation: str | None",
        title="Key",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If True or a str, save the figure. A string is appended to the default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}. Original type annotation: bool | None",
        title="Save",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.rank_genes_groups_tracksplot")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlHighestExprGenes(BaseAPI):
    """
    Fraction of counts assigned to each gene over all cells. Computes, for each gene, the fraction of counts assigned to that gene within a cell. The `n_top` genes with the highest mean fraction over all cells are plotted as boxplots. This plot is similar to the `scater` package function `plotHighestExprs(type = 'highest-expression')`, see [here](https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette-qc.html). Quoting from there: 'We expect to see the “usual suspects”, i.e., mitochondrial genes, actin, ribosomal protein, MALAT1. A few spike-in transcripts may also be present here, though if all of the spike-ins are in the top 50, it suggests that too much spike-in RNA was added. A large number of pseudo-genes or predicted genes may indicate problems with alignment.' -- Davis McCarthy and Aaron Lun
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    n_top: Optional[Any] = Field(
        30, description="Number of top. Original type annotation: int", title="N Top"
    )
    layer: Any = Field(
        None,
        description="Layer from which to pull data. Original type annotation: str | None",
        title="Layer",
    )
    gene_symbols: Any = Field(
        None,
        description="Key for field in .var that stores gene symbols if you do not want to use .var_names. Original type annotation: str | None",
        title="Gene Symbols",
    )
    log: Optional[Any] = Field(
        False,
        description="Plot x-axis in log scale. Original type annotation: bool",
        title="Log",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: str | bool | None",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None",
        title="Ax",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.highest_expr_genes")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlTracksplot(BaseAPI):
    """
    Compact plot of expression of a list of genes. In this type of plot each var_name is plotted as a filled line plot where the y values correspond to the var_name values and x is each of the cells. Best results are obtained when using raw counts that are not log. 'groupby' is required to sort and order the values using the respective group and should be a categorical value.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`. If `var_names` is a mapping, then the key is used as a label to group the values. The mapping values should be sequences of valid `adata.var_names`. When `var_names` is a mapping, `var_group_labels` and `var_group_positions` are set. Original type annotation: _VarNames | Mapping[str, _VarNames]",
        title="Var Names",
    )
    groupby: Any = Field(
        ...,
        description="The key of the observation grouping to consider. Original type annotation: str",
        title="Groupby",
    )
    use_raw: Any = Field(
        None,
        description="Use `raw` attribute of `adata` if present. Original type annotation: bool | None",
        title="Use Raw",
    )
    log: Optional[Any] = Field(
        False,
        description="Plot on logarithmic axis. Original type annotation: bool",
        title="Log",
    )
    dendrogram: Optional[Any] = Field(
        False,
        description="If True or a valid dendrogram key, a dendrogram based on the hierarchical clustering between the `groupby` categories is added. The dendrogram information is computed using :func:`scanpy.tl.dendrogram`. Original type annotation: bool | str",
        title="Dendrogram",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols. By default `var_names` refer to the index column of the `.var` DataFrame. Setting this option allows alternative names to be used. Original type annotation: str | None",
        title="Gene Symbols",
    )
    var_group_positions: Any = Field(
        None,
        description="Use this parameter to highlight groups of `var_names`. This will draw a 'bracket' or a color block between the given start and end positions. If the parameter `var_group_labels` is set, the corresponding labels are added on top/left. Original type annotation: Sequence[tuple[int, int]] | None",
        title="Var Group Positions",
    )
    var_group_labels: Any = Field(
        None,
        description="Labels for each of the `var_group_positions` that want to be highlighted. Original type annotation: Sequence[str] | None",
        title="Var Group Labels",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name, then the layer is plotted. `layer` takes precedence over `use_raw`. Original type annotation: str | None",
        title="Layer",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: str | bool | None",
        title="Save",
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`. Otherwise the `rcParam['figure.figsize]` value is used. Format is (width, height). Original type annotation: tuple[float, float] | None",
        title="Figsize",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.tracksplot")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlClustermap(BaseAPI):
    """
    Hierarchically-clustered heatmap. Wraps seaborn.clustermap for anndata.AnnData.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    obs_keys: Any = Field(
        None,
        description="Categorical annotation to plot with a different color map. Currently, only a single key is supported. Original type annotation: str | None",
        title="Obs Keys",
    )
    use_raw: Any = Field(
        None,
        description="Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present. Original type annotation: bool | None",
        title="Use Raw",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: bool | str | None",
        title="Save",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.pl.clustermap")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlStackedViolin(BaseAPI):
    """
    Makes a compact image composed of individual violin plots stacked on top of each other. Useful for visualizing gene expression per cluster. Wraps seaborn.violinplot for AnnData. Provides a convenient interface to the StackedViolin class in scanpy.pl.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix representing the data to be visualized.",
        title="Adata",
    )
    var_names: Any = Field(
        ...,
        description="Subset of variable names to be used for grouping or coloring in the plot.",
        title="Var Names",
    )
    groupby: Any = Field(
        ...,
        description="Key specifying the observation grouping to consider for the plot.",
        title="Groupby",
    )
    log: Optional[Any] = Field(
        False,
        description="Boolean indicating whether to plot on a logarithmic axis.",
        title="Log",
    )
    use_raw: Any = Field(
        None,
        description="Boolean indicating whether to use the 'raw' attribute of the data if present.",
        title="Use Raw",
    )
    num_categories: Optional[Any] = Field(
        7,
        description="Number of groups to subdivide non-categorical groupby observations into.",
        title="Num Categories",
    )
    title: Any = Field(
        None, description="Title for the figure to be displayed.", title="Title"
    )
    colorbar_title: Optional[Any] = Field(
        "Median expression\nin group",
        description="Title for the color bar in the plot.",
        title="Colorbar Title",
    )
    figsize: Any = Field(
        None,
        description="Size of the figure when multiple panels are used.",
        title="Figsize",
    )
    dendrogram: Optional[Any] = Field(
        False,
        description="Option to add a dendrogram based on hierarchical clustering between groupby categories.",
        title="Dendrogram",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in the data frame that stores gene symbols.",
        title="Gene Symbols",
    )
    var_group_positions: Any = Field(
        None,
        description="Positions to highlight groups of variable names in the plot.",
        title="Var Group Positions",
    )
    var_group_labels: Any = Field(
        None,
        description="Labels for the var_group_positions to be highlighted.",
        title="Var Group Labels",
    )
    standard_scale: Any = Field(
        None,
        description="Option to standardize the given dimension between 0 and 1.",
        title="Standard Scale",
    )
    var_group_rotation: Any = Field(
        None,
        description="Degree of label rotation in the plot.",
        title="Var Group Rotation",
    )
    layer: Any = Field(
        None, description="Name of the data layer to be plotted.", title="Layer"
    )
    categories_order: Any = Field(
        None,
        description="Order in which to display the categories in the plot.",
        title="Categories Order",
    )
    swap_axes: Optional[Any] = Field(
        False,
        description="Boolean indicating whether to swap the x and y axes in the plot.",
        title="Swap Axes",
    )
    show: Any = Field(
        None,
        description="Boolean indicating whether to display the plot without returning the axis object.",
        title="Show",
    )
    save: Any = Field(
        None,
        description="Option to save the figure with a specified filename.",
        title="Save",
    )
    return_fig: Optional[Any] = Field(
        False,
        description="Boolean indicating whether to return the DotPlot object for fine-tuning.",
        title="Return Fig",
    )
    ax: Any = Field(
        None,
        description="Matplotlib axes object to be used for plotting a single component.",
        title="Ax",
    )
    vmin: Any = Field(
        None, description="Lower limit of the color scale in the plot.", title="Vmin"
    )
    vmax: Any = Field(
        None, description="Upper limit of the color scale in the plot.", title="Vmax"
    )
    vcenter: Any = Field(
        None,
        description="Value representing the center of the color scale for diverging colormaps.",
        title="Vcenter",
    )
    norm: Any = Field(
        None,
        description="Custom color normalization object from matplotlib.",
        title="Norm",
    )
    cmap: Optional[Any] = Field(
        "Blues",
        description="String denoting the matplotlib colormap for coloring.",
        title="Cmap",
    )
    stripplot: Optional[Any] = Field(
        False,
        description="Boolean indicating whether to add a stripplot on top of the violin plot.",
        title="Stripplot",
    )
    jitter: Optional[Any] = Field(
        False,
        description="Amount of jitter to be added to the stripplot points.",
        title="Jitter",
    )
    size: Optional[Any] = Field(
        1, description="Size of the jitter points in the stripplot.", title="Size"
    )
    row_palette: Any = Field(
        None,
        description="Palette to be used for coloring each violin plot row.",
        title="Row Palette",
    )
    density_norm: Optional[Any] = Field(
        0,
        description="Method used to scale the width of each violin in the plot.",
        title="Density Norm",
    )
    yticklabels: Optional[Any] = Field(
        False,
        description="Boolean indicating whether to display the y-axis tick labels.",
        title="Yticklabels",
    )
    order: Optional[Any] = Field(
        0, description="Order in which to display the data categories.", title="Order"
    )
    scale: Optional[Any] = Field(
        0, description="Scaling method for the plot.", title="Scale"
    )
    kwds: Any = Field(
        ...,
        description="Additional parameters for customization not explicitly described.",
        title="Kwds",
    )
    _api_name: str = PrivateAttr(default="sc.pl.stacked_violin")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlMatrixplot(BaseAPI):
    """
    Create a heatmap of the mean expression values per group of each var_names. This function provides a convenient interface to the scanpy.pl.MatrixPlot class. If you need more flexibility, you should use the scanpy.pl.MatrixPlot directly.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`. If `var_names` is a mapping, then the key is used as label to group the values. The mapping values should be sequences of valid `adata.var_names`. Original type annotation: _VarNames | Mapping[str, _VarNames]",
        title="Var Names",
    )
    groupby: Any = Field(
        ...,
        description="The key of the observation grouping to consider. Original type annotation: str | Sequence[str]",
        title="Groupby",
    )
    use_raw: Any = Field(
        None,
        description="Use `raw` attribute of `adata` if present. Original type annotation: bool | None",
        title="Use Raw",
    )
    log: Optional[Any] = Field(
        False,
        description="Plot on a logarithmic axis. Original type annotation: bool",
        title="Log",
    )
    num_categories: Optional[Any] = Field(
        7,
        description="Determines the number of groups into which the groupby observation should be subdivided. Original type annotation: int",
        title="Num Categories",
    )
    categories_order: Any = Field(
        None,
        description="Order in which to show the categories. Note: add_dendrogram or add_totals can change the categories order. Original type annotation: Sequence[str] | None",
        title="Categories Order",
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`. Original type annotation: tuple[float, float] | None",
        title="Figsize",
    )
    dendrogram: Optional[Any] = Field(
        False,
        description="If True or a valid dendrogram key, a dendrogram based on the hierarchical clustering between the `groupby` categories is added. Original type annotation: bool | str",
        title="Dendrogram",
    )
    title: Any = Field(
        None,
        description="Title for the figure. Original type annotation: str | None",
        title="Title",
    )
    cmap: Optional[Any] = Field(
        "viridis",
        description="String denoting matplotlib color map. Original type annotation: Colormap | str | None",
        title="Cmap",
    )
    colorbar_title: Optional[Any] = Field(
        "Mean expression\nin group",
        description="Title for the color bar. New line character (\\n) can be used. Original type annotation: str | None",
        title="Colorbar Title",
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols. Original type annotation: str | None",
        title="Gene Symbols",
    )
    var_group_positions: Any = Field(
        None,
        description="Use this parameter to highlight groups of `var_names`. Original type annotation: Sequence[tuple[int, int]] | None",
        title="Var Group Positions",
    )
    var_group_labels: Any = Field(
        None,
        description="Labels for each of the `var_group_positions` that want to be highlighted. Original type annotation: Sequence[str] | None",
        title="Var Group Labels",
    )
    var_group_rotation: Any = Field(
        None,
        description="Label rotation degrees. Original type annotation: float | None",
        title="Var Group Rotation",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer that wants to be plotted. Original type annotation: str | None",
        title="Layer",
    )
    standard_scale: Any = Field(
        None,
        description="Whether or not to standardize the given dimension between 0 and 1. Original type annotation: Literal['var', 'group'] | None",
        title="Standard Scale",
    )
    values_df: Any = Field(
        None,
        description="No description available. Original type annotation: pd.DataFrame | None",
        title="Values Df",
    )
    swap_axes: Optional[Any] = Field(
        False,
        description="By setting `swap_axes` then x are the `groupby` categories and y the `var_names`. Original type annotation: bool",
        title="Swap Axes",
    )
    show: Any = Field(
        None,
        description="Show the plot, do not return axis. Original type annotation: bool | None",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. Original type annotation: str | bool | None",
        title="Save",
    )
    ax: Any = Field(
        None,
        description="A matplotlib axes object. Original type annotation: _AxesSubplot | None",
        title="Ax",
    )
    return_fig: Optional[Any] = Field(
        False,
        description="Returns :class:`DotPlot` object. Original type annotation: bool | None",
        title="Return Fig",
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale. Original type annotation: float | None",
        title="Vmin",
    )
    vmax: Any = Field(
        None,
        description="The value representing the upper limit of the color scale. Original type annotation: float | None",
        title="Vmax",
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Original type annotation: float | None",
        title="Vcenter",
    )
    norm: Any = Field(
        None,
        description="Custom color normalization object from matplotlib. Original type annotation: Normalize | None",
        title="Norm",
    )
    kwds: Any = Field(
        ..., description="Are passed to :func:`matplotlib.pyplot.pcolor`.", title="Kwds"
    )
    _api_name: str = PrivateAttr(default="sc.pl.matrixplot")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlScrubletScoreDistribution(BaseAPI):
    """
    Plot histogram of doublet scores for observed transcriptomes and simulated doublets. The histogram for simulated doublets is useful for determining the correct doublet score threshold. Scrublet must have been run previously with the input object.
    """

    adata: Any = Field(
        ...,
        description="An AnnData object resulting from :func:`~scanpy.pp.scrublet`. Original type annotation: AnnData",
        title="Adata",
    )
    scale_hist_obs: Optional[Any] = Field(
        "log",
        description="Set y axis scale transformation in matplotlib for the plot of observed transcriptomes. Original type annotation: Scale",
        title="Scale Hist Obs",
    )
    scale_hist_sim: Optional[Any] = Field(
        "linear",
        description="Set y axis scale transformation in matplotlib for the plot of simulated doublets. Original type annotation: Scale",
        title="Scale Hist Sim",
    )
    figsize: Optional[Any] = Field(
        [8, 3],
        description="Width, height. Original type annotation: tuple[float | int, float | int]",
        title="Figsize",
    )
    return_fig: Optional[Any] = Field(
        False,
        description="No description available. Original type annotation: bool",
        title="Return Fig",
    )
    show: Optional[Any] = Field(
        True,
        description="Show the plot, do not return axis. Original type annotation: bool",
        title="Show",
    )
    save: Any = Field(
        None,
        description="If :data:`True` or a :class:`str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: str | bool | None",
        title="Save",
    )
    _api_name: str = PrivateAttr(default="sc.pl.scrublet_score_distribution")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlPcaLoadings(BaseAPI):
    """
    Rank genes according to contributions to PCs.
    """

    adata: Any = Field(..., description="Annotated data matrix.", title="Adata")
    components: Any = Field(
        None,
        description="For example, '1,2,3' means [1, 2, 3], first, second, third principal component.",
        title="Components",
    )
    include_lowest: Optional[Any] = Field(
        True,
        description="Whether to show the variables with both highest and lowest loadings.",
        title="Include Lowest",
    )
    n_points: Any = Field(
        None,
        description="Number of variables to plot for each component.",
        title="N Points",
    )
    show: Any = Field(
        None, description="Show the plot, do not return axis.", title="Show"
    )
    save: Any = Field(
        None,
        description="If True or a str, save the figure. A string is appended to the default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}.",
        title="Save",
    )
    _api_name: str = PrivateAttr(default="sc.pl.pca_loadings")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlDrawGraph(BaseAPI):
    """
    Scatter plot in graph-drawing basis.
    """

    adata: Any = Field(..., description="Annotated data matrix.", title="Adata")
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes.",
        title="Color",
    )
    mask_obs: Any = Field(
        None, description="No description available.", title="Mask Obs"
    )
    gene_symbols: Any = Field(
        None,
        description="Column name in `.var` DataFrame that stores gene symbols.",
        title="Gene Symbols",
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression.",
        title="Use Raw",
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points with higher values on top of others.",
        title="Sort Order",
    )
    edges: Optional[Any] = Field(False, description="Show edges.", title="Edges")
    edges_width: Optional[Any] = Field(
        0.1, description="Width of edges.", title="Edges Width"
    )
    edges_color: Optional[Any] = Field(
        "grey", description="Color of edges.", title="Edges Color"
    )
    neighbors_key: Any = Field(
        None,
        description="Where to look for neighbors connectivities.",
        title="Neighbors Key",
    )
    arrows: Optional[Any] = Field(False, description="Show arrows.", title="Arrows")
    arrows_kwds: Any = Field(
        None,
        description="Passed to :meth:`~matplotlib.axes.Axes.quiver`.",
        title="Arrows Kwds",
    )
    groups: Any = Field(
        None,
        description="Restrict to a few categories in categorical observation annotation.",
        title="Groups",
    )
    components: Any = Field(None, description="Components to plot.", title="Components")
    dimensions: Any = Field(
        None,
        description="0-indexed dimensions of the embedding to plot as integers.",
        title="Dimensions",
    )
    layer: Any = Field(
        None,
        description="Name of the AnnData object layer that wants to be plotted.",
        title="Layer",
    )
    projection: Optional[Any] = Field(
        "2d", description="Projection of plot.", title="Projection"
    )
    scale_factor: Any = Field(
        None, description="No description available.", title="Scale Factor"
    )
    color_map: Any = Field(
        None, description="Color map to use for continous variables.", title="Color Map"
    )
    cmap: Any = Field(None, description="No description available.", title="Cmap")
    palette: Any = Field(
        None,
        description="Colors to use for plotting categorical annotation groups.",
        title="Palette",
    )
    na_color: Optional[Any] = Field(
        "lightgray",
        description="Color to use for null or masked values.",
        title="Na Color",
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="If there are missing values, whether they get an entry in the legend.",
        title="Na In Legend",
    )
    size: Any = Field(None, description="Point size.", title="Size")
    frameon: Any = Field(
        None, description="Draw a frame around the scatter plot.", title="Frameon"
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string describing the size.",
        title="Legend Fontsize",
    )
    legend_fontweight: Optional[Any] = Field(
        "bold", description="Legend font weight.", title="Legend Fontweight"
    )
    legend_loc: Optional[Any] = Field(
        "right margin", description="Location of legend.", title="Legend Loc"
    )
    legend_fontoutline: Any = Field(
        None,
        description="Line width of the legend font outline in pt.",
        title="Legend Fontoutline",
    )
    colorbar_loc: Optional[Any] = Field(
        "right",
        description="Where to place the colorbar for continous variables.",
        title="Colorbar Loc",
    )
    vmax: Any = Field(
        None,
        description="The value representing the upper limit of the color scale.",
        title="Vmax",
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale.",
        title="Vmin",
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale.",
        title="Vcenter",
    )
    norm: Any = Field(None, description="No description available.", title="Norm")
    add_outline: Optional[Any] = Field(
        False,
        description="If set to True, this will add a thin border around groups of dots.",
        title="Add Outline",
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description="Width numbers used to adjust the outline.",
        title="Outline Width",
    )
    outline_color: Optional[Any] = Field(
        ["black", "white"],
        description="Valid color names used to adjust the add_outline.",
        title="Outline Color",
    )
    ncols: Optional[Any] = Field(
        4, description="Number of panels per row.", title="Ncols"
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.",
        title="Hspace",
    )
    wspace: Any = Field(
        None,
        description="Adjust the width of the space between multiple panels.",
        title="Wspace",
    )
    title: Any = Field(None, description="Provide title for panels.", title="Title")
    show: Any = Field(
        None, description="Show the plot, do not return axis.", title="Show"
    )
    save: Any = Field(
        None, description="If `True` or a `str`, save the figure.", title="Save"
    )
    ax: Any = Field(None, description="A matplotlib axes object.", title="Ax")
    return_fig: Any = Field(
        None, description="Return the matplotlib figure.", title="Return Fig"
    )
    marker: Optional[Any] = Field(
        ".", description="No description available.", title="Marker"
    )
    layout: Any = Field(
        None,
        description="One of the :func:`~scanpy.tl.draw_graph` layouts.",
        title="Layout",
    )
    _api_name: str = PrivateAttr(default="sc.pl.draw_graph")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


class ScPlPagaPath(BaseAPI):
    """
    Gene expression and annotation changes along paths in the abstracted graph.
    """

    adata: Any = Field(
        ..., description="An annotated data matrix representing AnnData.", title="Adata"
    )
    nodes: Any = Field(
        ...,
        description="A path through nodes of the abstracted graph, using names or indices of groups that were used in PAGA analysis.",
        title="Nodes",
    )
    keys: Any = Field(
        ...,
        description="Variables in adata.var_names or annotations in adata.obs to be plotted using the specified color_map.",
        title="Keys",
    )
    use_raw: Optional[Any] = Field(
        True,
        description="Boolean indicating whether to use adata.raw for retrieving gene expressions if set.",
        title="Use Raw",
    )
    annotations: Optional[Any] = Field(
        ["dpt_pseudotime"],
        description="Keys to plot with color_maps_annotations, must be keys for adata.obs.",
        title="Annotations",
    )
    color_map: Any = Field(
        None,
        description="Matplotlib colormap as a string or Colormap object.",
        title="Color Map",
    )
    color_maps_annotations: Optional[Any] = Field(
        {"dpt_pseudotime": "Greys"},
        description="Color maps for plotting annotations, keys must match annotations in a mapping structure.",
        title="Color Maps Annotations",
    )
    palette_groups: Any = Field(
        None,
        description="Palette groups to be used for coloring the abstracted graph.",
        title="Palette Groups",
    )
    n_avg: Optional[Any] = Field(
        1,
        description="Number of data points to include in the computation of running average.",
        title="N Avg",
    )
    groups_key: Any = Field(
        None,
        description="Key of the grouping used in PAGA analysis, defaults to adata.uns['paga']['groups'] if None.",
        title="Groups Key",
    )
    xlim: Optional[Any] = Field(
        [None, None], description="Tuple representing the x-axis limits.", title="Xlim"
    )
    title: Any = Field(None, description="Title of the plot.", title="Title")
    left_margin: Any = Field(
        None, description="Left margin of the plot.", title="Left Margin"
    )
    ytick_fontsize: Any = Field(
        None, description="Font size for y-axis ticks.", title="Ytick Fontsize"
    )
    title_fontsize: Any = Field(
        None, description="Font size for the plot title.", title="Title Fontsize"
    )
    show_node_names: Optional[Any] = Field(
        True,
        description="Boolean indicating whether to plot node names on the nodes bar.",
        title="Show Node Names",
    )
    show_yticks: Optional[Any] = Field(
        True,
        description="Boolean indicating whether to show the y-axis ticks.",
        title="Show Yticks",
    )
    show_colorbar: Optional[Any] = Field(
        True,
        description="Boolean indicating whether to show the colorbar.",
        title="Show Colorbar",
    )
    legend_fontsize: Any = Field(
        None, description="Font size for the legend.", title="Legend Fontsize"
    )
    legend_fontweight: Any = Field(
        None, description="Font weight for the legend.", title="Legend Fontweight"
    )
    normalize_to_zero_one: Optional[Any] = Field(
        False,
        description="Boolean indicating whether to shift and scale the running average to [0, 1] per gene.",
        title="Normalize To Zero One",
    )
    as_heatmap: Optional[Any] = Field(
        True,
        description="Boolean indicating whether to plot the timeseries as a heatmap.",
        title="As Heatmap",
    )
    return_data: Optional[Any] = Field(
        False,
        description="Boolean indicating whether to return the timeseries data in addition to the axes.",
        title="Return Data",
    )
    show: Any = Field(
        None,
        description="Boolean indicating whether to show the plot or not.",
        title="Show",
    )
    save: Any = Field(
        None,
        description="Boolean or string indicating whether to save the figure with an appended filename.",
        title="Save",
    )
    ax: Any = Field(None, description="A matplotlib axes object or None.", title="Ax")
    _api_name: str = PrivateAttr(default="sc.pl.paga_path")
    _products_original: list[str] = PrivateAttr(default=[])
    _data_name: str = PrivateAttr(default="adata")


TOOLS_DICT = {
    "sc.pl.paga": ScPlPaga,
    "sc.pl.scatter": ScPlScatter,
    "sc.pl.umap": ScPlUmap,
    "sc.pl.tsne": ScPlTsne,
    "sc.pl.heatmap": ScPlHeatmap,
    "sc.pl.dotplot": ScPlDotplot,
    "sc.pl.violin": ScPlViolin,
    "sc.pl.dendrogram": ScPlDendrogram,
    "sc.pl.diffmap": ScPlDiffmap,
    "sc.pl.highly_variable_genes": ScPlHighlyVariableGenes,
    "sc.pl.pca": ScPlPca,
    "sc.pl.embedding_density": ScPlEmbeddingDensity,
    "sc.pl.rank_genes_groups": ScPlRankGenesGroups,
    "sc.pl.rank_genes_groups_dotplot": ScPlRankGenesGroupsDotplot,
    "sc.pl.rank_genes_groups_violin": ScPlRankGenesGroupsViolin,
    "sc.pl.rank_genes_groups_heatmap": ScPlRankGenesGroupsHeatmap,
    "sc.pl.rank_genes_groups_stacked_violin": ScPlRankGenesGroupsStackedViolin,
    "sc.pl.rank_genes_groups_matrixplot": ScPlRankGenesGroupsMatrixplot,
    "sc.pl.rank_genes_groups_tracksplot": ScPlRankGenesGroupsTracksplot,
    "sc.pl.highest_expr_genes": ScPlHighestExprGenes,
    "sc.pl.tracksplot": ScPlTracksplot,
    "sc.pl.clustermap": ScPlClustermap,
    "sc.pl.stacked_violin": ScPlStackedViolin,
    "sc.pl.matrixplot": ScPlMatrixplot,
    "sc.pl.scrublet_score_distribution": ScPlScrubletScoreDistribution,
    "sc.pl.pca_loadings": ScPlPcaLoadings,
    "sc.pl.draw_graph": ScPlDrawGraph,
    "sc.pl.paga_path": ScPlPagaPath,
}
