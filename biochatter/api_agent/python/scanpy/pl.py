from __future__ import annotations

from typing import Any, Optional

from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScanpyPlottingPaga(BaseAPI):
    """
    Plot the PAGA graph through thresholding low-connectivity edges. Compute a coarse-grained layout of the data. Reuse this by passing `init_pos='paga'` to :func:`~scanpy.tl.umap` or :func:`~scanpy.tl.draw_graph` and obtain embeddings with more meaningful global topology :cite:p:`Wolf2019`. This uses ForceAtlas2 or igraph's layout algorithms for most layouts :cite:p:`Csardi2006`.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    threshold: Any = Field(
        None,
        description='Threshold for edge drawing based on weights. Can be a float or None.',
        title='Threshold',
    )
    color: Any = Field(
        None,
        description='Specifies node colors or degree visualization. Can also be used for pie chart visualization. Can be a string or a mapping.',
        title='Color',
    )
    layout: Any = Field(
        None,
        description="Plotting layout for computing node positions. Can be specific layouts like 'fa' or 'fr'. Can be _Layout type or None.",
        title='Layout',
    )
    layout_kwds: Optional[Any] = Field(
        {},
        description='Keywords for the layout. Original type annotation: Mapping[str, Any]',
        title='Layout Kwds',
    )
    init_pos: Any = Field(
        None,
        description='Initial x and y coordinates for layout. Can be an array or None.',
        title='Init Pos',
    )
    root: Optional[Any] = Field(
        0,
        description='Index of root node for tree layout or list of root node indices. Can be int, str, or Sequence[int].',
        title='Root',
    )
    labels: Any = Field(
        None,
        description='Node labels for visualization. Can be str, Sequence[str], or Mapping[str, str].',
        title='Labels',
    )
    single_component: Optional[Any] = Field(
        False,
        description='Boolean to restrict visualization to the largest connected component.',
        title='Single Component',
    )
    solid_edges: Optional[Any] = Field(
        'connectivities',
        description='Key that specifies the matrix storing solid black edges.',
        title='Solid Edges',
    )
    dashed_edges: Any = Field(
        None,
        description='Key that specifies the matrix storing dashed grey edges. Can be None.',
        title='Dashed Edges',
    )
    transitions: Any = Field(
        None,
        description='Key that specifies the matrix storing arrows. Can be None.',
        title='Transitions',
    )
    fontsize: Any = Field(
        None,
        description='Font size for node labels. Can be int or None.',
        title='Fontsize',
    )
    fontweight: Optional[Any] = Field(
        'bold',
        description='No description available. Original type annotation: str',
        title='Fontweight',
    )
    fontoutline: Any = Field(
        None,
        description='Width of white outline around fonts. Can be int or None.',
        title='Fontoutline',
    )
    text_kwds: Optional[Any] = Field(
        {},
        description='Keywords for text labels. Original type annotation: Mapping[str, Any]',
        title='Text Kwds',
    )
    node_size_scale: Optional[Any] = Field(
        1.0,
        description='Scale factor for node sizes. Must be a float.',
        title='Node Size Scale',
    )
    node_size_power: Optional[Any] = Field(
        0.5,
        description='Influence of group sizes on node radius. Must be a float.',
        title='Node Size Power',
    )
    edge_width_scale: Optional[Any] = Field(
        1.0,
        description='Scale factor for edge width. Must be a float.',
        title='Edge Width Scale',
    )
    min_edge_width: Any = Field(
        None,
        description='Minimum width of solid edges. Can be float or None.',
        title='Min Edge Width',
    )
    max_edge_width: Any = Field(
        None,
        description='Maximum width of solid and dashed edges. Can be float or None.',
        title='Max Edge Width',
    )
    arrowsize: Optional[Any] = Field(
        30,
        description='Size of arrow heads for directed graphs. Must be an int.',
        title='Arrowsize',
    )
    title: Any = Field(
        None,
        description='Title for the visualization. Can be str or None.',
        title='Title',
    )
    left_margin: Optional[Any] = Field(
        0.01,
        description='No description available. Must be a float.',
        title='Left Margin',
    )
    random_state: Optional[Any] = Field(
        0,
        description="Random initialization state for layouts like 'fr'. Must be an int or None.",
        title='Random State',
    )
    pos: Any = Field(
        None,
        description='Array or path to file for node coordinates. Can be np.ndarray, Path, str, or None.',
        title='Pos',
    )
    normalize_to_color: Optional[Any] = Field(
        False,
        description='Boolean for normalizing categorical plots to color or grouping.',
        title='Normalize To Color',
    )
    cmap: Any = Field(
        None,
        description='Color map for visualization. Can be str, Colormap, or None.',
        title='Cmap',
    )
    cax: Any = Field(
        None,
        description='Matplotlib axes object for colorbar. Can be Axes or None.',
        title='Cax',
    )
    colorbar: Any = Field(
        None, description='No description available.', title='Colorbar'
    )
    cb_kwds: Optional[Any] = Field(
        {},
        description='Keyword arguments for colorbar. Original type annotation: Mapping[str, Any]',
        title='Cb Kwds',
    )
    frameon: Any = Field(
        None,
        description='Boolean to draw a frame around the graph. Can be None.',
        title='Frameon',
    )
    add_pos: Optional[Any] = Field(
        True,
        description="Boolean to add positions to adata.uns['paga'].",
        title='Add Pos',
    )
    export_to_gexf: Optional[Any] = Field(
        False,
        description='Boolean to export to gexf format for external programs.',
        title='Export To Gexf',
    )
    use_raw: Optional[Any] = Field(
        True,
        description='No description available. Must be a boolean.',
        title='Use Raw',
    )
    colors: Any = Field(None, description='No description available.', title='Colors')
    groups: Any = Field(None, description='No description available.', title='Groups')
    plot: Optional[Any] = Field(
        True,
        description='Boolean to determine whether to create the figure.',
        title='Plot',
    )
    show: Any = Field(
        None, description='No description available. Can be bool or None.', title='Show'
    )
    save: Any = Field(
        None,
        description="Boolean or str to save the figure. Can infer filetype from string ending in '.pdf', '.png', or '.svg'.",
        title='Save',
    )
    ax: Any = Field(
        None, description='Matplotlib axes object. Can be Axes or None.', title='Ax'
    )

    _api_name = PrivateAttr(default='scanpy.plotting.paga')
    _products_original = PrivateAttr(default=['data.uns["paga"]["pos"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingScatter(BaseAPI):
    """
    Scatter plot along observations or variables axes. Color the plot using annotations of observations (.obs), variables (.var) or expression of genes (.var_names).
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    x: Any = Field(
        None,
        description='x coordinate. Original type annotation: str | None',
        title='X',
    )
    y: Any = Field(
        None,
        description='y coordinate. Original type annotation: str | None',
        title='Y',
    )
    color: Any = Field(
        None,
        description='Keys for annotations of observations/cells or variables/genes, or a hex color specification. Original type annotation: str | ColorLike | Collection[str | ColorLike] | None',
        title='Color',
    )
    use_raw: Any = Field(
        None,
        description="Whether to use 'raw' attribute of 'adata'. Defaults to 'True' if '.raw' is present. Original type annotation: bool | None",
        title='Use Raw',
    )
    layers: Any = Field(
        None,
        description="Use the 'layers' attribute of 'adata' if present, specify the layer for 'x', 'y' and 'color'. Original type annotation: str | Collection[str] | None",
        title='Layers',
    )
    sort_order: Optional[Any] = Field(
        True,
        description='For continuous annotations used as color parameter, plot data points with higher values on top of others. Original type annotation: bool',
        title='Sort Order',
    )
    alpha: Any = Field(
        None,
        description='No description available. Original type annotation: float | None',
        title='Alpha',
    )
    basis: Any = Field(
        None,
        description='String that denotes a plotting tool that computed coordinates. Original type annotation: _Basis | None',
        title='Basis',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation. Original type annotation: str | Iterable[str] | None',
        title='Groups',
    )
    components: Any = Field(
        None,
        description="For instance, ['1,2', '2,3']. To plot all available components use 'components='all''. Original type annotation: str | Collection[str] | None",
        title='Components',
    )
    projection: Optional[Any] = Field(
        '2d',
        description="Projection of plot (default: '2d'). Original type annotation: Literal['2d', '3d']",
        title='Projection',
    )
    legend_loc: Optional[Any] = Field(
        'right margin',
        description="Location of legend, either 'on data', 'right margin', 'None', or a valid keyword for the 'loc' parameter of :class:`~matplotlib.legend.Legend`. Original type annotation: _LegendLoc | None",
        title='Legend Loc',
    )
    legend_fontsize: Any = Field(
        None,
        description='Numeric size in pt or string describing the size. Original type annotation: float | _FontSize | None',
        title='Legend Fontsize',
    )
    legend_fontweight: Any = Field(
        None,
        description='Legend font weight. Original type annotation: int | _FontWeight | None',
        title='Legend Fontweight',
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline in pt. Original type annotation: float | None',
        title='Legend Fontoutline',
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continuous variables. Original type annotation: str | Colormap | None',
        title='Color Map',
    )
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups. Original type annotation: Cycler | ListedColormap | ColorLike | Sequence[ColorLike] | None',
        title='Palette',
    )
    frameon: Any = Field(
        None,
        description="Draw a frame around the scatter plot. Default to value set in :func:`~scanpy.set_figure_params`, defaults to 'True'. Original type annotation: bool | None",
        title='Frameon',
    )
    right_margin: Any = Field(
        None,
        description='No description available. Original type annotation: float | None',
        title='Right Margin',
    )
    left_margin: Any = Field(
        None,
        description='No description available. Original type annotation: float | None',
        title='Left Margin',
    )
    size: Any = Field(
        None,
        description='Point size. Original type annotation: float | None',
        title='Size',
    )
    marker: Optional[Any] = Field(
        '.', description='Original type annotation: str | Sequence[str]', title='Marker'
    )
    title: Any = Field(
        None,
        description='Provide title for panels either as string or list of strings. Original type annotation: str | Collection[str] | None',
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If 'True' or a 'str', save the figure. Original type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Original type annotation: Axes | None',
        title='Ax',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.scatter')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingUmap(BaseAPI):
    """
    Scatter plot in UMAP basis.
    """

    
    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: <class 'anndata._core.anndata.AnnData'>",
        title='Adata',
    )
    color: Any = Field(
        None,
        description='Keys for annotations of observations/cells or variables/genes.',
        title='Color',
    )
    mask_obs: Any = Field(
        None, description='No description available.', title='Mask Obs'
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description='Use `.raw` attribute of `adata` for coloring with gene expression.',
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description='For continuous annotations used as color parameter, plot data points with higher values on top of others.',
        title='Sort Order',
    )
    edges: Optional[Any] = Field(False, description='Show edges.', title='Edges')
    edges_width: Optional[Any] = Field(
        0.1, description='Width of edges.', title='Edges Width'
    )
    edges_color: Optional[Any] = Field(
        'grey', description='Color of edges.', title='Edges Color'
    )
    neighbors_key: Any = Field(
        None,
        description='Where to look for neighbors connectivities.',
        title='Neighbors Key',
    )
    arrows: Optional[Any] = Field(
        False,
        description='Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).',
        title='Arrows',
    )
    arrows_kwds: Any = Field(
        None,
        description='Passed to :meth:`~matplotlib.axes.Axes.quiver`.',
        title='Arrows Kwds',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.',
        title='Groups',
    )
    components: Any = Field(
        None, description='For plotting specific components.', title='Components'
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot as integers.',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted.',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d', description="Projection of plot (default: `'2d'`).", title='Projection'
    )
    scale_factor: Any = Field(
        None, description='No description available.', title='Scale Factor'
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continuous variables.',
        title='Color Map',
    )
    cmap: Any = Field(None, description='No description available.', title='Cmap')
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values.',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description='If there are missing values, whether they get an entry in the legend.',
        title='Na In Legend',
    )
    size: Any = Field(None, description='Point size.', title='Size')
    frameon: Any = Field(
        None, description='Draw a frame around the scatter plot.', title='Frameon'
    )
    legend_fontsize: Any = Field(
        None,
        description='Numeric size in pt or string describing the size.',
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold', description='Legend font weight.', title='Legend Fontweight'
    )
    legend_loc: Optional[Any] = Field(
        'right margin', description='Location of legend.', title='Legend Loc'
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline in pt.',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continuous variables.',
        title='Colorbar Loc',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale.',
        title='Vmax',
    )
    vmin: Any = Field(
        None,
        description='The value representing the lower limit of the color scale.',
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description='The value representing the center of the color scale.',
        title='Vcenter',
    )
    norm: Any = Field(None, description='No description available.', title='Norm')
    add_outline: Optional[Any] = Field(
        False,
        description='If set to True, this will add a thin border around groups of dots.',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Width numbers used to adjust the outline.',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Valid color names used to adjust the add_outline.',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4, description='Number of panels per row.', title='Ncols'
    )
    hspace: Optional[Any] = Field(
        0.25,
        description='Adjust the height of the space between multiple panels.',
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.',
        title='Wspace',
    )
    title: Any = Field(None, description='Provide title for panels.', title='Title')
    show: Any = Field(
        None, description='Show the plot, do not return axis.', title='Show'
    )
    save: Any = Field(
        None, description='If `True` or a `str`, save the figure.', title='Save'
    )
    ax: Any = Field(None, description='A matplotlib axes object.', title='Ax')
    return_fig: Any = Field(
        None, description='Return the matplotlib figure.', title='Return Fig'
    )
    marker: Optional[Any] = Field(
        '.', description='No description available.', title='Marker'
    )

    _api_name = PrivateAttr(default='scanpy.plotting.umap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingTsne(BaseAPI):
    """
    Scatter plot in tSNE basis.
    """

    
    adata: Any = Field(..., description='Annotated data matrix.', title='Adata')
    color: Any = Field(
        None,
        description='Keys for annotations of observations/cells or variables/genes.',
        title='Color',
    )
    mask_obs: Any = Field(
        None, description='No description available.', title='Mask Obs'
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description='Use `.raw` attribute of `adata` for coloring with gene expression.',
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description='For continuous annotations used as color parameter, plot data points with higher values on top of others.',
        title='Sort Order',
    )
    edges: Optional[Any] = Field(False, description='Show edges.', title='Edges')
    edges_width: Optional[Any] = Field(
        0.1, description='Width of edges.', title='Edges Width'
    )
    edges_color: Optional[Any] = Field(
        'grey', description='Color of edges.', title='Edges Color'
    )
    neighbors_key: Any = Field(
        None,
        description='Where to look for neighbors connectivities.',
        title='Neighbors Key',
    )
    arrows: Optional[Any] = Field(
        False,
        description='Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).',
        title='Arrows',
    )
    arrows_kwds: Any = Field(
        None,
        description='Passed to :meth:`~matplotlib.axes.Axes.quiver`.',
        title='Arrows Kwds',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.',
        title='Groups',
    )
    components: Any = Field(
        None,
        description='For plotting specific components or all available components.',
        title='Components',
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot.',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted.',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d', description='Projection of plot.', title='Projection'
    )
    scale_factor: Any = Field(
        None, description='No description available.', title='Scale Factor'
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continuous variables.',
        title='Color Map',
    )
    cmap: Any = Field(None, description='No description available.', title='Cmap')
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values.',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description='If missing values get an entry in the legend.',
        title='Na In Legend',
    )
    size: Any = Field(None, description='Point size for plotting.', title='Size')
    frameon: Any = Field(
        None, description='Draw a frame around the scatter plot.', title='Frameon'
    )
    legend_fontsize: Any = Field(
        None,
        description='Numeric size in pt or string describing the size.',
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold', description='Legend font weight.', title='Legend Fontweight'
    )
    legend_loc: Optional[Any] = Field(
        'right margin', description='Location of legend.', title='Legend Loc'
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline.',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continuous variables.',
        title='Colorbar Loc',
    )
    vmax: Any = Field(None, description='Upper limit of the color scale.', title='Vmax')
    vmin: Any = Field(None, description='Lower limit of the color scale.', title='Vmin')
    vcenter: Any = Field(
        None, description='Center of the color scale.', title='Vcenter'
    )
    norm: Any = Field(None, description='No description available.', title='Norm')
    add_outline: Optional[Any] = Field(
        False,
        description='If set to True, this will add a thin border around groups of dots.',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Width numbers used to adjust the outline.',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Valid color names used to adjust the outline.',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4, description='Number of panels per row.', title='Ncols'
    )
    hspace: Optional[Any] = Field(
        0.25,
        description='Adjust the height of the space between multiple panels.',
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.',
        title='Wspace',
    )
    title: Any = Field(None, description='Provide title for panels.', title='Title')
    show: Any = Field(None, description='Show the plot.', title='Show')
    save: Any = Field(None, description='Save the figure.', title='Save')
    ax: Any = Field(None, description='A matplotlib axes object.', title='Ax')
    return_fig: Any = Field(
        None, description='Return the matplotlib figure.', title='Return Fig'
    )
    marker: Optional[Any] = Field(
        '.', description='No description available.', title='Marker'
    )

    _api_name = PrivateAttr(default='scanpy.plotting.tsne')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingHeatmap(BaseAPI):
    """
    Heatmap of the expression values of genes. If `groupby` is given, the heatmap is ordered by the respective group. For example, a list of marker genes can be plotted, ordered by clustering. If the `groupby` observation annotation is not categorical the observation annotation is turned into a categorical by binning the data into the number specified in `num_categories`.
    """

    
    adata: Any = Field(..., description='Annotated data matrix.', title='Adata')
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`. If `var_names` is a mapping, the key is used to group the values, and it can be used for coloring or 'brackets' for grouping var names in the plot.",
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='The key of the observation grouping to consider.',
        title='Groupby',
    )
    use_raw: Any = Field(
        None, description='Use `raw` attribute of `adata` if present.', title='Use Raw'
    )
    log: Optional[Any] = Field(
        False, description='Plot on a logarithmic axis.', title='Log'
    )
    num_categories: Optional[Any] = Field(
        7,
        description="Determines the number of groups for groupby observation if it's not categorical.",
        title='Num Categories',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='Adds a dendrogram based on hierarchical clustering between the groupby categories.',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame storing gene symbols.',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description="Highlights groups of `var_names` using 'brackets' or color blocks between given start and end positions.",
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for each of the `var_group_positions` that are highlighted.',
        title='Var Group Labels',
    )
    var_group_rotation: Any = Field(
        None,
        description='Label rotation degrees for the plot.',
        title='Var Group Rotation',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer to be plotted.',
        title='Layer',
    )
    standard_scale: Any = Field(
        None,
        description='Standardizes the dimension between 0 and 1.',
        title='Standard Scale',
    )
    swap_axes: Optional[Any] = Field(
        False, description='Switches the x and y axes in the plot.', title='Swap Axes'
    )
    show_gene_labels: Any = Field(
        None,
        description='Controls the visibility of gene labels in the plot.',
        title='Show Gene Labels',
    )
    show: Any = Field(
        None, description='Shows the plot without returning the axis.', title='Show'
    )
    save: Any = Field(
        None,
        description='Saves the figure with an optional filename and filetype.',
        title='Save',
    )
    figsize: Any = Field(
        None, description='Figure size for multi-panel plots.', title='Figsize'
    )
    vmin: Any = Field(
        None, description='The lower limit of the color scale.', title='Vmin'
    )
    vmax: Any = Field(
        None, description='The upper limit of the color scale.', title='Vmax'
    )
    vcenter: Any = Field(
        None,
        description='The center of the color scale for diverging colormaps.',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib.',
        title='Norm',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.heatmap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingDotplot(BaseAPI):
    """
    Make a dot plot of the expression values of var_names. For each var_name and each groupby category, a dot is plotted representing the mean expression within each category (visualized by color) and the fraction of cells expressing the var_name in the category (visualized by the size of the dot). If groupby is not given, the dotplot assumes all data belongs to a single category. A gene is considered expressed if the expression value in the adata (or adata.raw) is above the specified threshold which is zero by default. This function provides a convenient interface to the scanpy.pl.DotPlot class for visualizing mean expression and percentage of cells expressing a gene across multiple clusters.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix with additional information about the data structure.',
        title='Adata',
    )
    var_names: Any = Field(
        ...,
        description='Subset of adata.var_names used for grouping and labeling in the plot.',
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='Key for grouping observations to consider in the plot.',
        title='Groupby',
    )
    use_raw: Any = Field(
        None,
        description='Boolean determining whether to use the raw attribute of adata if present.',
        title='Use Raw',
    )
    log: Optional[Any] = Field(
        False,
        description='Boolean indicating whether to plot data on a logarithmic axis.',
        title='Log',
    )
    num_categories: Optional[Any] = Field(
        7,
        description='Number of groups to subdivide the groupby observation into if it is not categorical.',
        title='Num Categories',
    )
    categories_order: Any = Field(
        None,
        description='Order in which to display the categories in the plot.',
        title='Categories Order',
    )
    expression_cutoff: Optional[Any] = Field(
        0.0,
        description='Threshold for gene expression values to determine if a gene is expressed.',
        title='Expression Cutoff',
    )
    mean_only_expressed: Optional[Any] = Field(
        False,
        description='Boolean to average gene expression only over cells expressing the genes.',
        title='Mean Only Expressed',
    )
    standard_scale: Any = Field(
        None,
        description='Option to standardize values between 0 and 1 for each variable or group.',
        title='Standard Scale',
    )
    title: Any = Field(None, description='Title for the figure.', title='Title')
    colorbar_title: Optional[Any] = Field(
        'Mean expression\nin group',
        description='Title for the color bar in the plot.',
        title='Colorbar Title',
    )
    size_title: Optional[Any] = Field(
        'Fraction of cells\nin group (%)',
        description='Title for the size legend in the plot.',
        title='Size Title',
    )
    figsize: Any = Field(
        None,
        description='Size of the figure when multi_panel=True, otherwise uses default rcParams value.',
        title='Figsize',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='Option to add a dendrogram based on hierarchical clustering to the plot.',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in .var DataFrame that stores gene symbols.',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description='Parameter to highlight groups of var_names with brackets or color blocks.',
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for var_group_positions to be highlighted in the plot.',
        title='Var Group Labels',
    )
    var_group_rotation: Any = Field(
        None,
        description='Rotation degrees for labels in the plot.',
        title='Var Group Rotation',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer to be plotted, with options for raw or specific layers.',
        title='Layer',
    )
    swap_axes: Optional[Any] = Field(
        False, description='Option to swap x and y axes in the plot.', title='Swap Axes'
    )
    dot_color_df: Any = Field(
        None, description='DataFrame for dot colors in the plot.', title='Dot Color Df'
    )
    show: Any = Field(
        None,
        description='Boolean to show the plot without returning the axis.',
        title='Show',
    )
    save: Any = Field(
        None,
        description='Option to save the figure with specified filename and filetype.',
        title='Save',
    )
    ax: Any = Field(
        None,
        description='Matplotlib axes object for plotting a single component.',
        title='Ax',
    )
    return_fig: Optional[Any] = Field(
        False,
        description='Option to return a DotPlot object for fine-tuning the plot.',
        title='Return Fig',
    )
    vmin: Any = Field(
        None, description='Lower limit of the color scale in the plot.', title='Vmin'
    )
    vmax: Any = Field(
        None, description='Upper limit of the color scale in the plot.', title='Vmax'
    )
    vcenter: Any = Field(
        None,
        description='Center of the color scale, useful for diverging colormaps.',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib.',
        title='Norm',
    )
    cmap: Optional[Any] = Field(
        'Reds', description='Matplotlib color map for the plot.', title='Cmap'
    )
    dot_max: Any = Field(
        None, description='Maximum size for dots in the plot.', title='Dot Max'
    )
    dot_min: Any = Field(
        None, description='Minimum size for dots in the plot.', title='Dot Min'
    )
    smallest_dot: Optional[Any] = Field(
        0.0,
        description='Size for expression levels with the minimum dot size.',
        title='Smallest Dot',
    )
    kwds: Any = Field(
        ...,
        description='Additional keyword arguments passed to matplotlib.pyplot.scatter.',
        title='Kwds',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.dotplot')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingViolin(BaseAPI):
    """
    Violin plot. Wraps seaborn.violinplot for AnnData.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    keys: Any = Field(
        ...,
        description='Keys for accessing variables of `.var_names` or fields of `.obs`. Original type annotation: str | Sequence[str]',
        title='Keys',
    )
    groupby: Any = Field(
        None,
        description='The key of the observation grouping to consider. Original type annotation: str | None',
        title='Groupby',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axis. Original type annotation: bool',
        title='Log',
    )
    use_raw: Any = Field(
        None,
        description='Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present. Original type annotation: bool | None',
        title='Use Raw',
    )
    stripplot: Optional[Any] = Field(
        True,
        description='Add a stripplot on top of the violin plot. See :func:`~seaborn.stripplot`. Original type annotation: bool',
        title='Stripplot',
    )
    jitter: Optional[Any] = Field(
        True,
        description='Add jitter to the stripplot (only when stripplot is True). See :func:`~seaborn.stripplot`. Original type annotation: float | bool',
        title='Jitter',
    )
    size: Optional[Any] = Field(
        1,
        description='Size of the jitter points. Original type annotation: int',
        title='Size',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name, then the layer is plotted. `layer` takes precedence over `use_raw`. Original type annotation: str | None',
        title='Layer',
    )
    density_norm: Optional[Any] = Field(
        'width',
        description="The method used to scale the width of each violin. If 'width' (the default), each violin will have the same width. If 'area', each violin will have the same area. If 'count', a violin’s width corresponds to the number of observations. Original type annotation: DensityNorm",
        title='Density Norm',
    )
    order: Any = Field(
        None,
        description='Order in which to show the categories. Original type annotation: Sequence[str] | None',
        title='Order',
    )
    multi_panel: Any = Field(
        None,
        description='Display keys in multiple panels also when `groupby is not None`. Original type annotation: bool | None',
        title='Multi Panel',
    )
    xlabel: Optional[Any] = Field(
        '',
        description='Label of the x axis. Defaults to `groupby` if `rotation` is `None`, otherwise, no label is shown. Original type annotation: str',
        title='Xlabel',
    )
    ylabel: Any = Field(
        None,
        description="Label of the y axis. If `None` and `groupby` is `None`, defaults to `'value'`. If `None` and `groubpy` is not `None`, defaults to `keys`. Original type annotation: str | Sequence[str] | None",
        title='Ylabel',
    )
    rotation: Any = Field(
        None,
        description='Rotation of xtick labels. Original type annotation: float | None',
        title='Rotation',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None',
        title='Ax',
    )
    scale: Optional[Any] = Field(
        0,
        description='No description available. Original type annotation: DensityNorm | Empty',
        title='Scale',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.violin')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingDendrogram(BaseAPI):
    """
    Plot a dendrogram of the categories defined in `groupby`. See :func:`~scanpy.tl.dendrogram`.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    groupby: Any = Field(
        ...,
        description='Categorical data column used to create the dendrogram. Original type annotation: str',
        title='Groupby',
    )
    dendrogram_key: Any = Field(
        None,
        description="Key under which the dendrogram information was stored. By default, the dendrogram information is stored under `.uns[f'dendrogram_{groupby}']`. Original type annotation: str | None",
        title='Dendrogram Key',
    )
    orientation: Optional[Any] = Field(
        'top',
        description="Origin of the tree. Will grow into the opposite direction. Original type annotation: Literal['top', 'bottom', 'left', 'right']",
        title='Orientation',
    )
    remove_labels: Optional[Any] = Field(
        False,
        description='Don’t draw labels. Used, for example, by scanpy.pl.matrixplot to annotate matrix columns/rows. Original type annotation: bool',
        title='Remove Labels',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If True or a str, save the figure. A string is appended to the default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}. Original type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None',
        title='Ax',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.dendrogram')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingDiffmap(BaseAPI):
    """
    Scatter plot in Diffusion Map basis.
    """

    
    adata: Any = Field(..., description='Annotated data matrix.', title='Adata')
    color: Any = Field(
        None,
        description='Keys for annotations of observations/cells or variables/genes.',
        title='Color',
    )
    mask_obs: Any = Field(
        None, description='No description available.', title='Mask Obs'
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description='Use `.raw` attribute of `adata` for coloring with gene expression.',
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description='For continuous annotations used as color parameter, plot data points with higher values on top of others.',
        title='Sort Order',
    )
    edges: Optional[Any] = Field(
        False, description='No description available.', title='Edges'
    )
    edges_width: Optional[Any] = Field(
        0.1, description='No description available.', title='Edges Width'
    )
    edges_color: Optional[Any] = Field(
        'grey', description='No description available.', title='Edges Color'
    )
    neighbors_key: Any = Field(
        None, description='No description available.', title='Neighbors Key'
    )
    arrows: Optional[Any] = Field(
        False, description='No description available.', title='Arrows'
    )
    arrows_kwds: Any = Field(
        None, description='No description available.', title='Arrows Kwds'
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.',
        title='Groups',
    )
    components: Any = Field(
        None,
        description='For plotting specific components or all available components.',
        title='Components',
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot as integers.',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted.',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d', description="Projection of plot (default: '2d').", title='Projection'
    )
    scale_factor: Any = Field(
        None, description='No description available.', title='Scale Factor'
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continuous variables.',
        title='Color Map',
    )
    cmap: Any = Field(None, description='No description available.', title='Cmap')
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values.',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description='Whether missing values get an entry in the legend.',
        title='Na In Legend',
    )
    size: Any = Field(
        None,
        description='Point size or a sequence containing the size for each cell.',
        title='Size',
    )
    frameon: Any = Field(
        None, description='Draw a frame around the scatter plot.', title='Frameon'
    )
    legend_fontsize: Any = Field(
        None,
        description='Numeric size or string describing the size of the legend font.',
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold', description='Legend font weight.', title='Legend Fontweight'
    )
    legend_loc: Optional[Any] = Field(
        'right margin', description='Location of the legend.', title='Legend Loc'
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline.',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continuous variables.',
        title='Colorbar Loc',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale.',
        title='Vmax',
    )
    vmin: Any = Field(
        None,
        description='The value representing the lower limit of the color scale.',
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description='The value representing the center of the color scale.',
        title='Vcenter',
    )
    norm: Any = Field(None, description='No description available.', title='Norm')
    add_outline: Optional[Any] = Field(
        False,
        description='Add a thin border around groups of dots.',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Adjust the width of the border color and gap color.',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Adjust the border and gap color for the outline.',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4, description='Number of panels per row.', title='Ncols'
    )
    hspace: Optional[Any] = Field(
        0.25,
        description='Adjust the height of the space between multiple panels.',
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.',
        title='Wspace',
    )
    title: Any = Field(None, description='Provide title for panels.', title='Title')
    show: Any = Field(
        None, description='Show the plot, do not return axis.', title='Show'
    )
    save: Any = Field(
        None, description='If `True` or a `str`, save the figure.', title='Save'
    )
    ax: Any = Field(None, description='A matplotlib axes object.', title='Ax')
    return_fig: Any = Field(
        None, description='Return the matplotlib figure.', title='Return Fig'
    )
    marker: Optional[Any] = Field(
        '.', description='No description available.', title='Marker'
    )

    _api_name = PrivateAttr(default='scanpy.plotting.diffmap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingHighlyVariableGenes(BaseAPI):
    """
    Plot dispersions or normalized variance versus means for genes. Produces Supp. Fig. 5c of Zheng et al. (2017) and MeanVarPlot() and VariableFeaturePlot() of Seurat.
    """

    
    adata_or_result: Any = Field(
        ...,
        description='No description available. Original type annotation: AnnData | pd.DataFrame | np.recarray',
        title='Adata Or Result',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axes. Original type annotation: bool',
        title='Log',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {'pdf', 'png', 'svg'}. Original type annotation: bool | str | None",
        title='Save',
    )
    highly_variable_genes: Optional[Any] = Field(
        True,
        description='No description available. Original type annotation: bool',
        title='Highly Variable Genes',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.highly_variable_genes')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingPca(BaseAPI):
    """
    Scatter plot in PCA coordinates. Use the parameter `annotate_var_explained` to annotate the explained variance.
    """

    
    adata: Any = Field(..., description='Annotated data matrix.', title='Adata')
    color: Any = Field(
        None,
        description='Keys for annotations of observations/cells or variables/genes.',
        title='Color',
    )
    mask_obs: Any = Field(
        None, description='No description available.', title='Mask Obs'
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description='Use `.raw` attribute of `adata` for coloring with gene expression.',
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description='For continuous annotations used as color parameter, plot data points with higher values on top of others.',
        title='Sort Order',
    )
    edges: Optional[Any] = Field(
        False, description='No description available.', title='Edges'
    )
    edges_width: Optional[Any] = Field(
        0.1, description='No description available.', title='Edges Width'
    )
    edges_color: Optional[Any] = Field(
        'grey', description='No description available.', title='Edges Color'
    )
    neighbors_key: Any = Field(
        None, description='No description available.', title='Neighbors Key'
    )
    arrows: Optional[Any] = Field(
        False, description='No description available.', title='Arrows'
    )
    arrows_kwds: Any = Field(
        None, description='No description available.', title='Arrows Kwds'
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.',
        title='Groups',
    )
    components: Any = Field(
        None,
        description="For plotting components, e.g., ['1,2', '2,3'].",
        title='Components',
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot as integers.',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted.',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d', description="Projection of plot (default: '2d').", title='Projection'
    )
    scale_factor: Any = Field(
        None, description='No description available.', title='Scale Factor'
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continuous variables.',
        title='Color Map',
    )
    cmap: Any = Field(None, description='No description available.', title='Cmap')
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values.',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description='Whether missing values get an entry in the legend.',
        title='Na In Legend',
    )
    size: Any = Field(None, description='Point size for plotting.', title='Size')
    frameon: Any = Field(
        None, description='Draw a frame around the scatter plot.', title='Frameon'
    )
    legend_fontsize: Any = Field(
        None,
        description='Numeric size or string describing the size of the legend font.',
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold', description='Legend font weight.', title='Legend Fontweight'
    )
    legend_loc: Optional[Any] = Field(
        'right margin', description='Location of the legend.', title='Legend Loc'
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline.',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continuous variables.',
        title='Colorbar Loc',
    )
    vmax: Any = Field(
        None, description='The upper limit of the color scale.', title='Vmax'
    )
    vmin: Any = Field(
        None, description='The lower limit of the color scale.', title='Vmin'
    )
    vcenter: Any = Field(
        None, description='The center of the color scale.', title='Vcenter'
    )
    norm: Any = Field(None, description='No description available.', title='Norm')
    add_outline: Optional[Any] = Field(
        False,
        description='If set to true, this will add a thin border around groups of dots.',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Width of the border and gap color for outline.',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Border and gap color for the outline.',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4, description='Number of panels per row.', title='Ncols'
    )
    hspace: Optional[Any] = Field(
        0.25,
        description='Adjust the height of the space between multiple panels.',
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.',
        title='Wspace',
    )
    title: Any = Field(None, description='Provide title for panels.', title='Title')
    show: Any = Field(
        None, description='Show the plot without returning axis.', title='Show'
    )
    save: Any = Field(
        None,
        description='Save the figure with options for filename and type.',
        title='Save',
    )
    ax: Any = Field(None, description='A matplotlib axes object.', title='Ax')
    return_fig: Any = Field(
        None, description='Return the matplotlib figure.', title='Return Fig'
    )
    marker: Optional[Any] = Field(
        '.', description='No description available.', title='Marker'
    )
    annotate_var_explained: Optional[Any] = Field(
        False, description='No description available.', title='Annotate Var Explained'
    )

    _api_name = PrivateAttr(default='scanpy.plotting.pca')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingEmbeddingDensity(BaseAPI):
    """
    Plot the density of cells in an embedding (per condition). Plots the gaussian kernel density estimates (over condition) from the sc.tl.embedding_density() output. This function was written by Sophie Tritschler and implemented into Scanpy by Malte Luecken.
    """

    
    adata: Any = Field(
        ...,
        description='The annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    basis: Optional[Any] = Field(
        'umap',
        description="The embedding over which the density was calculated. This embedded representation should be found in `adata.obsm['X_[basis]']`. Original type annotation: str",
        title='Basis',
    )
    key: Any = Field(
        None,
        description='Name of the `.obs` covariate that contains the density estimates. Alternatively, pass `groupby`. Original type annotation: str | None',
        title='Key',
    )
    groupby: Any = Field(
        None,
        description='Name of the condition used in `tl.embedding_density`. Alternatively, pass `key`. Original type annotation: str | None',
        title='Groupby',
    )
    group: Optional[Any] = Field(
        'all',
        description='The category in the categorical observation annotation to be plotted. Original type annotation: str | Sequence[str] | None',
        title='Group',
    )
    color_map: Optional[Any] = Field(
        'YlOrRd',
        description='Matplolib color map to use for density plotting. Original type annotation: Colormap | str',
        title='Color Map',
    )
    bg_dotsize: Optional[Any] = Field(
        80,
        description='Dot size for background data points not in the `group`. Original type annotation: int | None',
        title='Bg Dotsize',
    )
    fg_dotsize: Optional[Any] = Field(
        180,
        description='Dot size for foreground data points in the `group`. Original type annotation: int | None',
        title='Fg Dotsize',
    )
    vmax: Optional[Any] = Field(
        1,
        description='The value representing the upper limit of the color scale. Original type annotation: int | None',
        title='Vmax',
    )
    vmin: Optional[Any] = Field(
        0,
        description='The value representing the lower limit of the color scale. Values smaller than vmin are plotted with the same color as vmin. Original type annotation: int | None',
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description='The value representing the center of the color scale. Useful for diverging colormaps. Original type annotation: int | None',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='No description available. Original type annotation: Normalize | None',
        title='Norm',
    )
    ncols: Optional[Any] = Field(
        4,
        description='Number of panels per row. Original type annotation: int | None',
        title='Ncols',
    )
    hspace: Optional[Any] = Field(
        0.25,
        description='Adjust the height of the space between multiple panels. Original type annotation: float | None',
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels. Original type annotation: None',
        title='Wspace',
    )
    title: Any = Field(
        None,
        description='No description available. Original type annotation: str | None',
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}. Original type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None',
        title='Ax',
    )
    return_fig: Any = Field(
        None,
        description='Return the matplotlib figure. Original type annotation: bool | None',
        title='Return Fig',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.embedding_density')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingRankGenesGroupsDotplot(BaseAPI):
    """
    Plot ranking of genes using dotplot plot (see :func:`~scanpy.pl.dotplot`).
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    groups: Any = Field(
        None,
        description='The groups for which to show the gene ranking. Original type annotation: str | Sequence[str] | None',
        title='Groups',
    )
    n_genes: Any = Field(
        None,
        description='Number of genes to show. This can be a negative number to show for example the down regulated genes. Original type annotation: int | None',
        title='N Genes',
    )
    groupby: Any = Field(
        None,
        description='The key of the observation grouping to consider. By default, the groupby is chosen from the rank genes groups parameter but other groupby options can be used. Original type annotation: str | None',
        title='Groupby',
    )
    values_to_plot: Any = Field(
        None,
        description="Instead of the mean gene value, plot the values computed by `sc.rank_genes_groups`. Original type annotation: Literal['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'] | None",
        title='Values To Plot',
    )
    var_names: Any = Field(
        None,
        description='Genes to plot. Sometimes is useful to pass a specific list of var names (e.g. genes) to check their fold changes or p-values. Original type annotation: Sequence[str] | Mapping[str, Sequence[str]] | None',
        title='Var Names',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols. Original type annotation: str | None',
        title='Gene Symbols',
    )
    min_logfoldchange: Any = Field(
        None,
        description='Value to filter genes in groups if their logfoldchange is less than the min_logfoldchange. Original type annotation: float | None',
        title='Min Logfoldchange',
    )
    key: Any = Field(
        None,
        description='Key used to store the ranking results in `adata.uns`. Original type annotation: str | None',
        title='Key',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description='If `True` or a `str`, save the figure. Original type annotation: bool | None',
        title='Save',
    )
    return_fig: Optional[Any] = Field(
        False,
        description='Returns :class:`DotPlot` object. Useful for fine-tuning the plot. Original type annotation: bool',
        title='Return Fig',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.rank_genes_groups_dotplot')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingCorrelationMatrix(BaseAPI):
    """
    Plot the correlation matrix computed as part of `sc.tl.dendrogram`.
    """

    
    adata: Any = Field(
        ...,
        description='AnnData type object representing data to be plotted.',
        title='Adata',
    )
    groupby: Any = Field(
        ...,
        description='Categorical data column used to create the dendrogram.',
        title='Groupby',
    )
    show_correlation_numbers: Optional[Any] = Field(
        False,
        description='Option to plot the correlation on top of each cell.',
        title='Show Correlation Numbers',
    )
    dendrogram: Any = Field(
        None,
        description='Option to add a dendrogram based on hierarchical clustering between groupby categories.',
        title='Dendrogram',
    )
    figsize: Any = Field(
        None,
        description='Figure size for the correlation matrix plot.',
        title='Figsize',
    )
    show: Any = Field(
        None,
        description='Option to display the plot without returning the axis.',
        title='Show',
    )
    save: Any = Field(
        None,
        description='Option to save the figure with a specified filename and inferred filetype.',
        title='Save',
    )
    ax: Any = Field(
        None,
        description='Matplotlib axes object for plotting a single component.',
        title='Ax',
    )
    vmin: Any = Field(
        None, description='Lower limit of the color scale for plotting.', title='Vmin'
    )
    vmax: Any = Field(
        None, description='Upper limit of the color scale for plotting.', title='Vmax'
    )
    vcenter: Any = Field(
        None,
        description='Center value of the color scale, useful for diverging colormaps.',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib.',
        title='Norm',
    )
    kwds: Any = Field(
        ..., description='Additional keyword arguments for customization.', title='Kwds'
    )

    _api_name = PrivateAttr(default='scanpy.plotting.correlation_matrix')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingHighestExprGenes(BaseAPI):
    """
    Fraction of counts assigned to each gene over all cells. Computes, for each gene, the fraction of counts assigned to that gene within a cell. The `n_top` genes with the highest mean fraction over all cells are plotted as boxplots. This plot is similar to the `scater` package function `plotHighestExprs(type = "highest-expression")`, see `here <https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette-qc.html>`. Quoting from there: *We expect to see the “usual suspects”, i.e., mitochondrial genes, actin, ribosomal protein, MALAT1. A few spike-in transcripts may also be present here, though if all of the spike-ins are in the top 50, it suggests that too much spike-in RNA was added. A large number of pseudo-genes or predicted genes may indicate problems with alignment.* -- Davis McCarthy and Aaron Lun
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    n_top: Optional[Any] = Field(
        30, description='Number of top. Original type annotation: int', title='N Top'
    )
    layer: Any = Field(
        None,
        description='Layer from which to pull data. Original type annotation: str | None',
        title='Layer',
    )
    gene_symbols: Any = Field(
        None,
        description='Key for field in .var that stores gene symbols if you do not want to use .var_names. Original type annotation: str | None',
        title='Gene Symbols',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot x-axis in log scale. Original type annotation: bool',
        title='Log',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component. Original type annotation: Axes | None',
        title='Ax',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.highest_expr_genes')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingTracksplot(BaseAPI):
    """
    Compact plot of expression of a list of genes. In this type of plot each var_name is plotted as a filled line plot where the y values correspond to the var_name values and x is each of the cells. Best results are obtained when using raw counts that are not log. `groupby` is required to sort and order the values using the respective group and should be a categorical value.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    var_names: Any = Field(
        ...,
        description='`var_names` should be a valid subset of `adata.var_names`. If `var_names` is a mapping, the key is used as a label to group the values. The mapping values should be sequences of valid `adata.var_names`. Original type annotation: _VarNames | Mapping[str, _VarNames]',
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='The key of the observation grouping to consider. Original type annotation: str',
        title='Groupby',
    )
    use_raw: Any = Field(
        None,
        description='Use `raw` attribute of `adata` if present. Original type annotation: bool | None',
        title='Use Raw',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on a logarithmic axis. Original type annotation: bool',
        title='Log',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='If True or a valid dendrogram key, a dendrogram based on hierarchical clustering between the `groupby` categories is added. The dendrogram information is computed using scanpy.tl.dendrogram. Original type annotation: bool | str',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols. By default, `var_names` refer to the index column of the `.var` DataFrame. Setting this option allows alternative names to be used. Original type annotation: str | None',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description="Use this parameter to highlight groups of `var_names` by drawing 'brackets' or color blocks between given start and end positions. Original type annotation: Sequence[tuple[int, int]] | None",
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for each of the `var_group_positions` that want to be highlighted. Original type annotation: Sequence[str] | None',
        title='Var Group Labels',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default, adata.raw.X is plotted. Original type annotation: str | None',
        title='Layer',
    )
    show: Any = Field(
        None,
        description='Show the plot without returning the axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure with a string appended to the default filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}. Original type annotation: str | bool | None",
        title='Save',
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`. Otherwise, the `rcParam['figure.figsize]` value is used. Format is (width, height). Original type annotation: tuple[float, float] | None",
        title='Figsize',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.tracksplot')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingClustermap(BaseAPI):
    """
    Hierarchically-clustered heatmap. Wraps seaborn.clustermap for anndata.AnnData.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    obs_keys: Any = Field(
        None,
        description='Categorical annotation to plot with a different color map. Currently, only a single key is supported. Original type annotation: str | None',
        title='Obs Keys',
    )
    use_raw: Any = Field(
        None,
        description='Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present. Original type annotation: bool | None',
        title='Use Raw',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis. Original type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}. Original type annotation: bool | str | None",
        title='Save',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.clustermap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingStackedViolin(BaseAPI):
    """
    Stacked violin plots. Makes a compact image composed of individual violin plots stacked on top of each other. Useful for visualizing gene expression per cluster. Wraps seaborn.violinplot for anndata.AnnData. Provides a convenient interface to the scanpy.pl.StackedViolin class for more flexibility.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix with original type annotation as AnnData.',
        title='Adata',
    )
    var_names: Any = Field(
        ...,
        description='`var_names` should be a valid subset of `adata.var_names`, allowing for grouping and coloring based on a mapping if provided.',
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='Key of the observation grouping to consider, with original type annotation as string or sequence of strings.',
        title='Groupby',
    )
    log: Optional[Any] = Field(
        False,
        description='Boolean flag indicating whether to plot on a logarithmic axis.',
        title='Log',
    )
    use_raw: Any = Field(
        None,
        description='Boolean flag or None to specify whether to use `raw` attribute of `adata` if present.',
        title='Use Raw',
    )
    num_categories: Optional[Any] = Field(
        7,
        description='Integer value determining the number of groups to divide the groupby observation into if it is not categorical.',
        title='Num Categories',
    )
    title: Any = Field(
        None,
        description='Title for the figure with original type annotation as string or None.',
        title='Title',
    )
    colorbar_title: Optional[Any] = Field(
        'Median expression\nin group',
        description='Title for the color bar with support for newline characters and original type annotation as string or None.',
        title='Colorbar Title',
    )
    figsize: Any = Field(
        None,
        description='Figure size when `multi_panel=True`, otherwise using default values with original type annotation as tuple of floats or None.',
        title='Figsize',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='Option to include a dendrogram based on hierarchical clustering between groupby categories with original type annotation as bool or str.',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols with original type annotation as string or None.',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description='Parameter to highlight groups of `var_names` with original type annotation as sequence of tuples or None.',
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for highlighted `var_group_positions` with original type annotation as sequence of strings or None.',
        title='Var Group Labels',
    )
    standard_scale: Any = Field(
        None,
        description='Option to standardize dimension between 0 and 1 with original type annotation as specific literals or None.',
        title='Standard Scale',
    )
    var_group_rotation: Any = Field(
        None,
        description='Degree of label rotation for visualization with original type annotation as float or None.',
        title='Var Group Rotation',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer to plot with original type annotation as string or None.',
        title='Layer',
    )
    categories_order: Any = Field(
        None,
        description='Order in which to display the categories with original type annotation as sequence of strings or None.',
        title='Categories Order',
    )
    swap_axes: Optional[Any] = Field(
        False,
        description='Boolean flag to swap x and y axes for visualization.',
        title='Swap Axes',
    )
    show: Any = Field(
        None,
        description='Boolean flag or None to indicate whether to show the plot without returning axis.',
        title='Show',
    )
    save: Any = Field(
        None,
        description='Option to save the figure with a given filename extension and original type annotation as bool, string, or None.',
        title='Save',
    )
    return_fig: Optional[Any] = Field(
        False,
        description='Boolean flag or None to return a DotPlot object instead of showing the plot.',
        title='Return Fig',
    )
    ax: Any = Field(
        None,
        description='Matplotlib axes object to use for plotting a single component or None.',
        title='Ax',
    )
    vmin: Any = Field(
        None,
        description='Lower limit of the color scale with original type annotation as float or None.',
        title='Vmin',
    )
    vmax: Any = Field(
        None,
        description='Upper limit of the color scale with original type annotation as float or None.',
        title='Vmax',
    )
    vcenter: Any = Field(
        None,
        description='Center value of the color scale, particularly useful for diverging colormaps with original type annotation as float or None.',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib with original type annotation as specific object or None.',
        title='Norm',
    )
    cmap: Optional[Any] = Field(
        'Blues',
        description='String denoting matplotlib color map with original type annotation as colormap, string, or None.',
        title='Cmap',
    )
    stripplot: Optional[Any] = Field(
        False,
        description='Boolean flag to add a stripplot on top of the violin plot.',
        title='Stripplot',
    )
    jitter: Optional[Any] = Field(
        False,
        description='Degree of jitter to add to the stripplot if enabled with original type annotation as float, bool, or None.',
        title='Jitter',
    )
    size: Optional[Any] = Field(
        1,
        description='Size of the jitter points on the plot with original type annotation as float.',
        title='Size',
    )
    row_palette: Any = Field(
        None,
        description='Palette to assign colors to violin plot rows with original type annotation as string or None.',
        title='Row Palette',
    )
    density_norm: Optional[Any] = Field(
        0,
        description='Method used to scale the width of each violin plot with original type annotation as specific enum or empty.',
        title='Density Norm',
    )
    yticklabels: Optional[Any] = Field(
        False,
        description='Boolean flag to show or hide y tick labels on the plot.',
        title='Yticklabels',
    )
    order: Optional[Any] = Field(
        0,
        description='Sequence of strings or None/Empty type with no specific description available.',
        title='Order',
    )
    scale: Optional[Any] = Field(
        0,
        description='Density normalization method or empty with no specific description available.',
        title='Scale',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.stacked_violin')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')

TOOLS_DICT = {
    "scanpy.plotting.paga": ScanpyPlottingPaga,
    "scanpy.plotting.scatter": ScanpyPlottingScatter,
    "scanpy.plotting.umap": ScanpyPlottingUmap,
    "scanpy.plotting.tsne": ScanpyPlottingTsne,
    "scanpy.plotting.heatmap": ScanpyPlottingHeatmap,
    "scanpy.plotting.dotplot": ScanpyPlottingDotplot,
    "scanpy.plotting.violin": ScanpyPlottingViolin,
    "scanpy.plotting.dendrogram": ScanpyPlottingDendrogram,
    "scanpy.plotting.diffmap": ScanpyPlottingDiffmap,
    "scanpy.plotting.highly_variable_genes": ScanpyPlottingHighlyVariableGenes,
    "scanpy.plotting.pca": ScanpyPlottingPca,
    "scanpy.plotting.embedding_density": ScanpyPlottingEmbeddingDensity,
    "scanpy.plotting.rank_genes_groups_dotplot": ScanpyPlottingRankGenesGroupsDotplot,
    "scanpy.plotting.correlation_matrix": ScanpyPlottingCorrelationMatrix,
    "scanpy.plotting.highest_expr_genes": ScanpyPlottingHighestExprGenes,
    "scanpy.plotting.tracksplot": ScanpyPlottingTracksplot,
    "scanpy.plotting.clustermap": ScanpyPlottingClustermap,
    "scanpy.plotting.stacked_violin": ScanpyPlottingStackedViolin,
}