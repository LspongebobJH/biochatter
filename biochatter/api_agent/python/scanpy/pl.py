from __future__ import annotations

from typing import Any, Optional

from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScanpyPlottingPaga(BaseAPI):
    """
    Plot the PAGA graph through thresholding low-connectivity edges.
    
    Compute a coarse-grained layout of the data. Reuse this by passing
    `init_pos='paga'` to :func:`~scanpy.tl.umap` or
    :func:`~scanpy.tl.draw_graph` and obtain embeddings with more meaningful
    global topology :cite:p:`Wolf2019`.
    
    This uses ForceAtlas2 or igraph's layout algorithms for most layouts :cite:p:`Csardi2006`.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    threshold
        Do not draw edges for weights below this threshold. Set to 0 if you want
        all edges. Discarding low-connectivity edges helps in getting a much
        clearer picture of the graph.
    color
        Gene name or `obs` annotation defining the node colors.
        Also plots the degree of the abstracted graph when
        passing {`'degree_dashed'`, `'degree_solid'`}.
    
        Can be also used to visualize pie chart at each node in the following form:
        `{<group name or index>: {<color>: <fraction>, ...}, ...}`. If the fractions
        do not sum to 1, a new category called `'rest'` colored grey will be created.
    labels
        The node labels. If `None`, this defaults to the group labels stored in
        the categorical for which :func:`~scanpy.tl.paga` has been computed.
    pos
        Two-column array-like storing the x and y coordinates for drawing.
        Otherwise, path to a `.gdf` file that has been exported from Gephi or
        a similar graph visualization software.
    layout
        Plotting layout that computes positions.
        `'fa'` stands for “ForceAtlas2”,
        `'fr'` stands for “Fruchterman-Reingold”,
        `'rt'` stands for “Reingold-Tilford”,
        `'eq_tree'` stands for “eqally spaced tree”.
        All but `'fa'` and `'eq_tree'` are igraph layouts.
        All other igraph layouts are also permitted.
        See also parameter `pos` and :func:`~scanpy.tl.draw_graph`.
    layout_kwds
        Keywords for the layout.
    init_pos
        Two-column array storing the x and y coordinates for initializing the
        layout.
    random_state
        For layouts with random initialization like `'fr'`, change this to use
        different intial states for the optimization. If `None`, the initial
        state is not reproducible.
    root
        If choosing a tree layout, this is the index of the root node or a list
        of root node indices. If this is a non-empty vector then the supplied
        node IDs are used as the roots of the trees (or a single tree if the
        graph is connected). If this is `None` or an empty list, the root
        vertices are automatically calculated based on topological sorting.
    transitions
        Key for `.uns['paga']` that specifies the matrix that stores the
        arrows, for instance `'transitions_confidence'`.
    solid_edges
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn solid black.
    dashed_edges
        Key for `.uns['paga']` that specifies the matrix that stores the edges
        to be drawn dashed grey. If `None`, no dashed edges are drawn.
    single_component
        Restrict to largest connected component.
    fontsize
        Font size for node labels.
    fontoutline
        Width of the white outline around fonts.
    text_kwds
        Keywords for :meth:`~matplotlib.axes.Axes.text`.
    node_size_scale
        Increase or decrease the size of the nodes.
    node_size_power
        The power with which groups sizes influence the radius of the nodes.
    edge_width_scale
        Edge with scale in units of `rcParams['lines.linewidth']`.
    min_edge_width
        Min width of solid edges.
    max_edge_width
        Max width of solid and dashed edges.
    arrowsize
       For directed graphs, choose the size of the arrow head head's length and
       width. See :py:class: `matplotlib.patches.FancyArrowPatch` for attribute
       `mutation_scale` for more info.
    export_to_gexf
        Export to gexf format to be read by graph visualization programs such as
        Gephi.
    normalize_to_color
        Whether to normalize categorical plots to `color` or the underlying
        grouping.
    cmap
        The color map.
    cax
        A matplotlib axes object for a potential colorbar.
    cb_kwds
        Keyword arguments for :class:`~matplotlib.colorbar.Colorbar`,
        for instance, `ticks`.
    add_pos
        Add the positions to `adata.uns['paga']`.
    title
        Provide a title.
    frameon
        Draw a frame around the PAGA graph.
    plot
        If `False`, do not create the figure, simply compute the layout.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on \{`'.pdf'`, `'.png'`, `'.svg'`\}.
    ax
        A matplotlib axes object.
    
    Returns
    -------
    If `show==False`, one or more :class:`~matplotlib.axes.Axes` objects.
    Adds `'pos'` to `adata.uns['paga']` if `add_pos` is `True`.
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
    
        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()
        sc.tl.paga(adata, groups='louvain')
        sc.pl.paga(adata)
    
    You can increase node and edge sizes by specifying additional arguments.
    
    .. plot::
        :context: close-figs
    
        sc.pl.paga(adata, node_size_scale=10, edge_width_scale=2)
    
    Notes
    -----
    When initializing the positions, note that – for some reason – igraph
    mirrors coordinates along the x axis... that is, you should increase the
    `maxiter` parameter by 1 if the layout is flipped.
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    tl.paga
    pl.paga_compare
    pl.paga_path
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    threshold: Any = Field(
        None,
        description='Do not draw edges for weights below this threshold. Set to 0 if you want\nall edges. Discarding low-connectivity edges helps in getting a much\nclearer picture of the graph.\nOriginal type annotation: float | None',
        title='Threshold',
    )
    color: Any = Field(
        None,
        description="Gene name or `obs` annotation defining the node colors.\nAlso plots the degree of the abstracted graph when\npassing {`'degree_dashed'`, `'degree_solid'`}.\n\nCan be also used to visualize pie chart at each node in the following form:\n`{<group name or index>: {<color>: <fraction>, ...}, ...}`. If the fractions\ndo not sum to 1, a new category called `'rest'` colored grey will be created.\nOriginal type annotation: str | Mapping[str | int, Mapping[Any, float]] | None",
        title='Color',
    )
    layout: Any = Field(
        None,
        description="Plotting layout that computes positions.\n`'fa'` stands for “ForceAtlas2”,\n`'fr'` stands for “Fruchterman-Reingold”,\n`'rt'` stands for “Reingold-Tilford”,\n`'eq_tree'` stands for “eqally spaced tree”.\nAll but `'fa'` and `'eq_tree'` are igraph layouts.\nAll other igraph layouts are also permitted.\nSee also parameter `pos` and :func:`~scanpy.tl.draw_graph`.\nOriginal type annotation: _Layout | None",
        title='Layout',
    )
    layout_kwds: Optional[Any] = Field(
        {},
        description='Keywords for the layout.\nOriginal type annotation: Mapping[str, Any]',
        title='Layout Kwds',
    )
    init_pos: Any = Field(
        None,
        description='Two-column array storing the x and y coordinates for initializing the\nlayout.\nOriginal type annotation: np.ndarray | None',
        title='Init Pos',
    )
    root: Optional[Any] = Field(
        0,
        description='If choosing a tree layout, this is the index of the root node or a list\nof root node indices. If this is a non-empty vector then the supplied\nnode IDs are used as the roots of the trees (or a single tree if the\ngraph is connected). If this is `None` or an empty list, the root\nvertices are automatically calculated based on topological sorting.\nOriginal type annotation: int | str | Sequence[int] | None',
        title='Root',
    )
    labels: Any = Field(
        None,
        description='The node labels. If `None`, this defaults to the group labels stored in\nthe categorical for which :func:`~scanpy.tl.paga` has been computed.\nOriginal type annotation: str | Sequence[str] | Mapping[str, str] | None',
        title='Labels',
    )
    single_component: Optional[Any] = Field(
        False,
        description='Restrict to largest connected component.\nOriginal type annotation: bool',
        title='Single Component',
    )
    solid_edges: Optional[Any] = Field(
        'connectivities',
        description="Key for `.uns['paga']` that specifies the matrix that stores the edges\nto be drawn solid black.\nOriginal type annotation: str",
        title='Solid Edges',
    )
    dashed_edges: Any = Field(
        None,
        description="Key for `.uns['paga']` that specifies the matrix that stores the edges\nto be drawn dashed grey. If `None`, no dashed edges are drawn.\nOriginal type annotation: str | None",
        title='Dashed Edges',
    )
    transitions: Any = Field(
        None,
        description="Key for `.uns['paga']` that specifies the matrix that stores the\narrows, for instance `'transitions_confidence'`.\nOriginal type annotation: str | None",
        title='Transitions',
    )
    fontsize: Any = Field(
        None,
        description='Font size for node labels.\nOriginal type annotation: int | None',
        title='Fontsize',
    )
    fontweight: Optional[Any] = Field(
        'bold',
        description='No description available.\nOriginal type annotation: str',
        title='Fontweight',
    )
    fontoutline: Any = Field(
        None,
        description='Width of the white outline around fonts.\nOriginal type annotation: int | None',
        title='Fontoutline',
    )
    text_kwds: Optional[Any] = Field(
        {},
        description='Keywords for :meth:`~matplotlib.axes.Axes.text`.\nOriginal type annotation: Mapping[str, Any]',
        title='Text Kwds',
    )
    node_size_scale: Optional[Any] = Field(
        1.0,
        description='Increase or decrease the size of the nodes.\nOriginal type annotation: float',
        title='Node Size Scale',
    )
    node_size_power: Optional[Any] = Field(
        0.5,
        description='The power with which groups sizes influence the radius of the nodes.\nOriginal type annotation: float',
        title='Node Size Power',
    )
    edge_width_scale: Optional[Any] = Field(
        1.0,
        description="Edge with scale in units of `rcParams['lines.linewidth']`.\nOriginal type annotation: float",
        title='Edge Width Scale',
    )
    min_edge_width: Any = Field(
        None,
        description='Min width of solid edges.\nOriginal type annotation: float | None',
        title='Min Edge Width',
    )
    max_edge_width: Any = Field(
        None,
        description='Max width of solid and dashed edges.\nOriginal type annotation: float | None',
        title='Max Edge Width',
    )
    arrowsize: Optional[Any] = Field(
        30,
        description="For directed graphs, choose the size of the arrow head head's length and\nwidth. See :py:class: `matplotlib.patches.FancyArrowPatch` for attribute\n`mutation_scale` for more info.\nOriginal type annotation: int",
        title='Arrowsize',
    )
    title: Any = Field(
        None,
        description='Provide a title.\nOriginal type annotation: str | None',
        title='Title',
    )
    left_margin: Optional[Any] = Field(
        0.01,
        description='No description available.\nOriginal type annotation: float',
        title='Left Margin',
    )
    random_state: Optional[Any] = Field(
        0,
        description="For layouts with random initialization like `'fr'`, change this to use\ndifferent intial states for the optimization. If `None`, the initial\nstate is not reproducible.\nOriginal type annotation: int | None",
        title='Random State',
    )
    pos: Any = Field(
        None,
        description='Two-column array-like storing the x and y coordinates for drawing.\nOtherwise, path to a `.gdf` file that has been exported from Gephi or\na similar graph visualization software.\nOriginal type annotation: np.ndarray | Path | str | None',
        title='Pos',
    )
    normalize_to_color: Optional[Any] = Field(
        False,
        description='Whether to normalize categorical plots to `color` or the underlying\ngrouping.\nOriginal type annotation: bool',
        title='Normalize To Color',
    )
    cmap: Any = Field(
        None,
        description='The color map.\nOriginal type annotation: str | Colormap | None',
        title='Cmap',
    )
    cax: Any = Field(
        None,
        description='A matplotlib axes object for a potential colorbar.\nOriginal type annotation: Axes | None',
        title='Cax',
    )
    colorbar: Any = Field(
        None, description='No description available.', title='Colorbar'
    )
    cb_kwds: Optional[Any] = Field(
        {},
        description='Keyword arguments for :class:`~matplotlib.colorbar.Colorbar`,\nfor instance, `ticks`.\nOriginal type annotation: Mapping[str, Any]',
        title='Cb Kwds',
    )
    frameon: Any = Field(
        None,
        description='Draw a frame around the PAGA graph.\nOriginal type annotation: bool | None',
        title='Frameon',
    )
    add_pos: Optional[Any] = Field(
        True,
        description="Add the positions to `adata.uns['paga']`.\nOriginal type annotation: bool",
        title='Add Pos',
    )
    export_to_gexf: Optional[Any] = Field(
        False,
        description='Export to gexf format to be read by graph visualization programs such as\nGephi.\nOriginal type annotation: bool',
        title='Export To Gexf',
    )
    use_raw: Optional[Any] = Field(
        True,
        description='No description available.\nOriginal type annotation: bool',
        title='Use Raw',
    )
    colors: Any = Field(None, description='No description available.', title='Colors')
    groups: Any = Field(None, description='No description available.', title='Groups')
    plot: Optional[Any] = Field(
        True,
        description='If `False`, do not create the figure, simply compute the layout.\nOriginal type annotation: bool',
        title='Plot',
    )
    show: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on \\{`'.pdf'`, `'.png'`, `'.svg'`\\}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object.\nOriginal type annotation: Axes | None',
        title='Ax',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.paga')
    _products_original = PrivateAttr(default=['data.uns["paga"]["pos"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingScatter(BaseAPI):
    """
    Scatter plot along observations or variables axes.
    
    Color the plot using annotations of observations (`.obs`), variables
    (`.var`) or expression of genes (`.var_names`).
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    x
        x coordinate.
    y
        y coordinate.
    color
        Keys for annotations of observations/cells or variables/genes,
        or a hex color specification, e.g.,
        `'ann1'`, `'#fe57a1'`, or `['ann1', 'ann2']`.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
    layers
        Use the `layers` attribute of `adata` if present: specify the layer for
        `x`, `y` and `color`. If `layers` is a string, then it is expanded to
        `(layers, layers, layers)`.
    basis
        String that denotes a plotting tool that computed coordinates.
    sort_order
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups
        Restrict to a few categories in categorical observation annotation.
        The default is not to restrict to any groups.
    dimensions
        0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].
        Unlike `components`, this argument is used in the same way as `colors`, e.g. is
        used to specify a single plot at a time. Will eventually replace the components
        argument.
    components
        For instance, `['1,2', '2,3']`. To plot all available components use
        `components='all'`.
    projection
        Projection of plot (default: `'2d'`).
    legend_loc
        Location of legend, either `'on data'`, `'right margin'`, `None`,
        or a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.
    legend_fontsize
        Numeric size in pt or string describing the size.
        See :meth:`~matplotlib.text.Text.set_fontsize`.
    legend_fontweight
        Legend font weight. A numeric value in range 0-1000 or a string.
        Defaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.
        See :meth:`~matplotlib.text.Text.set_fontweight`.
    legend_fontoutline
        Line width of the legend font outline in pt. Draws a white outline using
        the path effect :class:`~matplotlib.patheffects.withStroke`.
    colorbar_loc
        Where to place the colorbar for continous variables. If `None`, no colorbar
        is added.
    size
        Point size. If `None`, is automatically computed as 120000 / n_cells.
        Can be a sequence containing the size for each cell. The order should be
        the same as in adata.obs.
    color_map
        Color map to use for continous variables. Can be a name or a
        :class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`
        or `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.
        If `None`, the value of `mpl.rcParams["image.cmap"]` is used.
        The default `color_map` can be set using :func:`~scanpy.set_figure_params`.
    palette
        Colors to use for plotting categorical annotation groups.
        The palette can be a valid :class:`~matplotlib.colors.ListedColormap` name
        (`'Set2'`, `'tab20'`, …), a :class:`~cycler.Cycler` object, a dict mapping
        categories to colors, or a sequence of colors. Colors must be valid to
        matplotlib. (see :func:`~matplotlib.colors.is_color_like`).
        If `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
        variable already has colors stored in `adata.uns["{var}_colors"]`.
        If provided, values of `adata.uns["{var}_colors"]` will be set.
    na_color
        Color to use for null or masked values. Can be anything matplotlib accepts as a
        color. Used for all points if `color=None`.
    na_in_legend
        If there are missing values, whether they get an entry in the legend. Currently
        only implemented for categorical legends.
    frameon
        Draw a frame around the scatter plot. Defaults to value set in
        :func:`~scanpy.set_figure_params`, defaults to `True`.
    title
        Provide title for panels either as string or list of strings,
        e.g. `['title1', 'title2', ...]`.
    
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    
    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    x: Any = Field(
        None,
        description='x coordinate.\nOriginal type annotation: str | None',
        title='X',
    )
    y: Any = Field(
        None,
        description='y coordinate.\nOriginal type annotation: str | None',
        title='Y',
    )
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes,\nor a hex color specification, e.g.,\n`'ann1'`, `'#fe57a1'`, or `['ann1', 'ann2']`.\nOriginal type annotation: str | ColorLike | Collection[str | ColorLike] | None",
        title='Color',
    )
    use_raw: Any = Field(
        None,
        description='Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    layers: Any = Field(
        None,
        description='Use the `layers` attribute of `adata` if present: specify the layer for\n`x`, `y` and `color`. If `layers` is a string, then it is expanded to\n`(layers, layers, layers)`.\nOriginal type annotation: str | Collection[str] | None',
        title='Layers',
    )
    sort_order: Optional[Any] = Field(
        True,
        description='For continuous annotations used as color parameter, plot data points\nwith higher values on top of others.\nOriginal type annotation: bool',
        title='Sort Order',
    )
    alpha: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: float | None',
        title='Alpha',
    )
    basis: Any = Field(
        None,
        description='String that denotes a plotting tool that computed coordinates.\nOriginal type annotation: _Basis | None',
        title='Basis',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.\nThe default is not to restrict to any groups.\nOriginal type annotation: str | Iterable[str] | None',
        title='Groups',
    )
    components: Any = Field(
        None,
        description="For instance, `['1,2', '2,3']`. To plot all available components use\n`components='all'`.\nOriginal type annotation: str | Collection[str] | None",
        title='Components',
    )
    projection: Optional[Any] = Field(
        '2d',
        description="Projection of plot (default: `'2d'`).\nOriginal type annotation: Literal['2d', '3d']",
        title='Projection',
    )
    legend_loc: Optional[Any] = Field(
        'right margin',
        description="Location of legend, either `'on data'`, `'right margin'`, `None`,\nor a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.\nOriginal type annotation: _LegendLoc | None",
        title='Legend Loc',
    )
    legend_fontsize: Any = Field(
        None,
        description='Numeric size in pt or string describing the size.\nSee :meth:`~matplotlib.text.Text.set_fontsize`.\nOriginal type annotation: float | _FontSize | None',
        title='Legend Fontsize',
    )
    legend_fontweight: Any = Field(
        None,
        description="Legend font weight. A numeric value in range 0-1000 or a string.\nDefaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.\nSee :meth:`~matplotlib.text.Text.set_fontweight`.\nOriginal type annotation: int | _FontWeight | None",
        title='Legend Fontweight',
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline in pt. Draws a white outline using\nthe path effect :class:`~matplotlib.patheffects.withStroke`.\nOriginal type annotation: float | None',
        title='Legend Fontoutline',
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continous variables. Can be a name or a\n:class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`\nor `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.\nIf `None`, the value of `mpl.rcParams["image.cmap"]` is used.\nThe default `color_map` can be set using :func:`~scanpy.set_figure_params`.\nOriginal type annotation: str | Colormap | None',
        title='Color Map',
    )
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.\nThe palette can be a valid :class:`~matplotlib.colors.ListedColormap` name\n(`\'Set2\'`, `\'tab20\'`, …), a :class:`~cycler.Cycler` object, a dict mapping\ncategories to colors, or a sequence of colors. Colors must be valid to\nmatplotlib. (see :func:`~matplotlib.colors.is_color_like`).\nIf `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical\nvariable already has colors stored in `adata.uns["{var}_colors"]`.\nIf provided, values of `adata.uns["{var}_colors"]` will be set.\nOriginal type annotation: Cycler | ListedColormap | ColorLike | Sequence[ColorLike] | None',
        title='Palette',
    )
    frameon: Any = Field(
        None,
        description='Draw a frame around the scatter plot. Defaults to value set in\n:func:`~scanpy.set_figure_params`, defaults to `True`.\nOriginal type annotation: bool | None',
        title='Frameon',
    )
    right_margin: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: float | None',
        title='Right Margin',
    )
    left_margin: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: float | None',
        title='Left Margin',
    )
    size: Any = Field(
        None,
        description='Point size. If `None`, is automatically computed as 120000 / n_cells.\nCan be a sequence containing the size for each cell. The order should be\nthe same as in adata.obs.\nOriginal type annotation: float | None',
        title='Size',
    )
    marker: Optional[Any] = Field(
        '.',
        description='No description available.\nOriginal type annotation: str | Sequence[str]',
        title='Marker',
    )
    title: Any = Field(
        None,
        description="Provide title for panels either as string or list of strings,\ne.g. `['title1', 'title2', ...]`.\nOriginal type annotation: str | Collection[str] | None",
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: Axes | None',
        title='Ax',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.scatter')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingUmap(BaseAPI):
    """
    Scatter plot in UMAP basis.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    color
        Keys for annotations of observations/cells or variables/genes, e.g.,
        `'ann1'` or `['ann1', 'ann2']`.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols. By default `var_names`
        refer to the index column of the `.var` DataFrame. Setting this option allows
        alternative names to be used.
    use_raw
        Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,
        defaults to `True` if `layer` isn't provided and `adata.raw` is present.
    layer
        Name of the AnnData object layer that wants to be plotted. By default
        adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.
        If `layer` is set to a valid layer name, then the layer is plotted. `layer`
        takes precedence over `use_raw`.
    edges
        Show edges.
    edges_width
        Width of edges.
    edges_color
        Color of edges. See :func:`~networkx.drawing.nx_pylab.draw_networkx_edges`.
    neighbors_key
        Where to look for neighbors connectivities.
        If not specified, this looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, this looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    arrows
        Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).
    arrows_kwds
        Passed to :meth:`~matplotlib.axes.Axes.quiver`
    sort_order
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups
        Restrict to a few categories in categorical observation annotation.
        The default is not to restrict to any groups.
    dimensions
        0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].
        Unlike `components`, this argument is used in the same way as `colors`, e.g. is
        used to specify a single plot at a time. Will eventually replace the components
        argument.
    components
        For instance, `['1,2', '2,3']`. To plot all available components use
        `components='all'`.
    projection
        Projection of plot (default: `'2d'`).
    legend_loc
        Location of legend, either `'on data'`, `'right margin'`, `None`,
        or a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.
    legend_fontsize
        Numeric size in pt or string describing the size.
        See :meth:`~matplotlib.text.Text.set_fontsize`.
    legend_fontweight
        Legend font weight. A numeric value in range 0-1000 or a string.
        Defaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.
        See :meth:`~matplotlib.text.Text.set_fontweight`.
    legend_fontoutline
        Line width of the legend font outline in pt. Draws a white outline using
        the path effect :class:`~matplotlib.patheffects.withStroke`.
    colorbar_loc
        Where to place the colorbar for continous variables. If `None`, no colorbar
        is added.
    size
        Point size. If `None`, is automatically computed as 120000 / n_cells.
        Can be a sequence containing the size for each cell. The order should be
        the same as in adata.obs.
    color_map
        Color map to use for continous variables. Can be a name or a
        :class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`
        or `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.
        If `None`, the value of `mpl.rcParams["image.cmap"]` is used.
        The default `color_map` can be set using :func:`~scanpy.set_figure_params`.
    palette
        Colors to use for plotting categorical annotation groups.
        The palette can be a valid :class:`~matplotlib.colors.ListedColormap` name
        (`'Set2'`, `'tab20'`, …), a :class:`~cycler.Cycler` object, a dict mapping
        categories to colors, or a sequence of colors. Colors must be valid to
        matplotlib. (see :func:`~matplotlib.colors.is_color_like`).
        If `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
        variable already has colors stored in `adata.uns["{var}_colors"]`.
        If provided, values of `adata.uns["{var}_colors"]` will be set.
    na_color
        Color to use for null or masked values. Can be anything matplotlib accepts as a
        color. Used for all points if `color=None`.
    na_in_legend
        If there are missing values, whether they get an entry in the legend. Currently
        only implemented for categorical legends.
    frameon
        Draw a frame around the scatter plot. Defaults to value set in
        :func:`~scanpy.set_figure_params`, defaults to `True`.
    title
        Provide title for panels either as string or list of strings,
        e.g. `['title1', 'title2', ...]`.
    
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin. vmin can be a number, a string, a function or `None`. If
        vmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).
        For example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then
        vmin is interpreted as the return value of the function over the list of values to plot.
        For example to set vmin tp the mean of the values to plot, `def my_vmin(values): return
        np.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic
        minimum value is used as defined by matplotlib `scatter` function. When making multiple
        plots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`
    vmax
        The value representing the upper limit of the color scale. The format is the same as for `vmin`.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
        The format is the same as for `vmin`.
        Example: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``
    add_outline
        If set to True, this will add a thin border around groups of dots. In some situations
        this can enhance the aesthetics of the resulting image
    outline_color
        Tuple with two valid color names used to adjust the add_outline. The first color is the
        border color (default: black), while the second color is a gap color between the
        border color and the scatter dot (default: white).
    outline_width
        Tuple with two width numbers used to adjust the outline. The first value is the width
        of the border color as a fraction of the scatter dot size (default: 0.3). The second value is
        width of the gap color (default: 0.05).
    ncols
        Number of panels per row.
    wspace
        Adjust the width of the space between multiple panels.
    hspace
        Adjust the height of the space between multiple panels.
    return_fig
        Return the matplotlib figure.
    kwargs
        Arguments to pass to :func:`matplotlib.pyplot.scatter`,
        for instance: the maximum and minimum values (e.g. `vmin=-2, vmax=5`).
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    
    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.pl.umap(adata)
    
    Colour points by discrete variable (Louvain clusters).
    
    .. plot::
        :context: close-figs
    
        sc.pl.umap(adata, color="louvain")
    
    Colour points by gene expression.
    
    .. plot::
        :context: close-figs
    
        sc.pl.umap(adata, color="HES4")
    
    Plot muliple umaps for different gene expressions.
    
    .. plot::
        :context: close-figs
    
        sc.pl.umap(adata, color=["HES4", "TNFRSF4"])
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    tl.umap
    """

    
    adata: Any = Field(
        ...,
        description="Annotated data matrix.\nOriginal type annotation: <class 'anndata._core.anndata.AnnData'>",
        title='Adata',
    )
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes, e.g.,\n`'ann1'` or `['ann1', 'ann2']`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Color',
    )
    mask_obs: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: numpy.ndarray[typing.Any, numpy.dtype[numpy.bool_]] | str | None',
        title='Mask Obs',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols. By default `var_names`\nrefer to the index column of the `.var` DataFrame. Setting this option allows\nalternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,\ndefaults to `True` if `layer` isn't provided and `adata.raw` is present.\nOriginal type annotation: bool | None",
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points\nwith higher values on top of others.\nOriginal type annotation: <class 'bool'>",
        title='Sort Order',
    )
    edges: Optional[Any] = Field(
        False,
        description="Show edges.\nOriginal type annotation: <class 'bool'>",
        title='Edges',
    )
    edges_width: Optional[Any] = Field(
        0.1,
        description="Width of edges.\nOriginal type annotation: <class 'float'>",
        title='Edges Width',
    )
    edges_color: Optional[Any] = Field(
        'grey',
        description='Color of edges. See :func:`~networkx.drawing.nx_pylab.draw_networkx_edges`.\nOriginal type annotation: str | collections.abc.Sequence[float] | collections.abc.Sequence[str]',
        title='Edges Color',
    )
    neighbors_key: Any = Field(
        None,
        description="Where to look for neighbors connectivities.\nIf not specified, this looks .obsp['connectivities'] for connectivities\n(default storage place for pp.neighbors).\nIf specified, this looks\n.obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.\nOriginal type annotation: str | None",
        title='Neighbors Key',
    )
    arrows: Optional[Any] = Field(
        False,
        description="Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).\nOriginal type annotation: <class 'bool'>",
        title='Arrows',
    )
    arrows_kwds: Any = Field(
        None,
        description='Passed to :meth:`~matplotlib.axes.Axes.quiver`\nOriginal type annotation: collections.abc.Mapping[str, typing.Any] | None',
        title='Arrows Kwds',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.\nThe default is not to restrict to any groups.\nOriginal type annotation: str | collections.abc.Sequence[str] | None',
        title='Groups',
    )
    components: Any = Field(
        None,
        description="For instance, `['1,2', '2,3']`. To plot all available components use\n`components='all'`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Components',
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].\nUnlike `components`, this argument is used in the same way as `colors`, e.g. is\nused to specify a single plot at a time. Will eventually replace the components\nargument.\nOriginal type annotation: tuple[int, int] | collections.abc.Sequence[tuple[int, int]] | None',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default\nadata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.\nIf `layer` is set to a valid layer name, then the layer is plotted. `layer`\ntakes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d',
        description="Projection of plot (default: `'2d'`).\nOriginal type annotation: typing.Literal['2d', '3d']",
        title='Projection',
    )
    scale_factor: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: float | None',
        title='Scale Factor',
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continous variables. Can be a name or a\n:class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`\nor `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.\nIf `None`, the value of `mpl.rcParams["image.cmap"]` is used.\nThe default `color_map` can be set using :func:`~scanpy.set_figure_params`.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Color Map',
    )
    cmap: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Cmap',
    )
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.\nThe palette can be a valid :class:`~matplotlib.colors.ListedColormap` name\n(`\'Set2\'`, `\'tab20\'`, …), a :class:`~cycler.Cycler` object, a dict mapping\ncategories to colors, or a sequence of colors. Colors must be valid to\nmatplotlib. (see :func:`~matplotlib.colors.is_color_like`).\nIf `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical\nvariable already has colors stored in `adata.uns["{var}_colors"]`.\nIf provided, values of `adata.uns["{var}_colors"]` will be set.\nOriginal type annotation: str | collections.abc.Sequence[str] | cycler.Cycler | None',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values. Can be anything matplotlib accepts as a\ncolor. Used for all points if `color=None`.\nOriginal type annotation: str | tuple[float, ...]',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="If there are missing values, whether they get an entry in the legend. Currently\nonly implemented for categorical legends.\nOriginal type annotation: <class 'bool'>",
        title='Na In Legend',
    )
    size: Any = Field(
        None,
        description='Point size. If `None`, is automatically computed as 120000 / n_cells.\nCan be a sequence containing the size for each cell. The order should be\nthe same as in adata.obs.\nOriginal type annotation: float | collections.abc.Sequence[float] | None',
        title='Size',
    )
    frameon: Any = Field(
        None,
        description='Draw a frame around the scatter plot. Defaults to value set in\n:func:`~scanpy.set_figure_params`, defaults to `True`.\nOriginal type annotation: bool | None',
        title='Frameon',
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string describing the size.\nSee :meth:`~matplotlib.text.Text.set_fontsize`.\nOriginal type annotation: typing.Union[float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], NoneType]",
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold',
        description="Legend font weight. A numeric value in range 0-1000 or a string.\nDefaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.\nSee :meth:`~matplotlib.text.Text.set_fontweight`.\nOriginal type annotation: typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]",
        title='Legend Fontweight',
    )
    legend_loc: Optional[Any] = Field(
        'right margin',
        description="Location of legend, either `'on data'`, `'right margin'`, `None`,\nor a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.\nOriginal type annotation: typing.Optional[typing.Literal['none', 'right margin', 'on data', 'on data export', 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center']]",
        title='Legend Loc',
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline in pt. Draws a white outline using\nthe path effect :class:`~matplotlib.patheffects.withStroke`.\nOriginal type annotation: int | None',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continous variables. If `None`, no colorbar\nis added.\nOriginal type annotation: str | None',
        title='Colorbar Loc',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. The format is the same as for `vmin`.\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None',
        title='Vmax',
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin. vmin can be a number, a string, a function or `None`. If\nvmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).\nFor example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then\nvmin is interpreted as the return value of the function over the list of values to plot.\nFor example to set vmin tp the mean of the values to plot, `def my_vmin(values): return\nnp.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic\nminimum value is used as defined by matplotlib `scatter` function. When making multiple\nplots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Useful for diverging colormaps.\nThe format is the same as for `vmin`.\nExample: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Normalize | collections.abc.Sequence[matplotlib.colors.Normalize] | None',
        title='Norm',
    )
    add_outline: Optional[Any] = Field(
        False,
        description='If set to True, this will add a thin border around groups of dots. In some situations\nthis can enhance the aesthetics of the resulting image\nOriginal type annotation: bool | None',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Tuple with two width numbers used to adjust the outline. The first value is the width\nof the border color as a fraction of the scatter dot size (default: 0.3). The second value is\nwidth of the gap color (default: 0.05).\nOriginal type annotation: tuple[float, float]',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Tuple with two valid color names used to adjust the add_outline. The first color is the\nborder color (default: black), while the second color is a gap color between the\nborder color and the scatter dot (default: white).\nOriginal type annotation: tuple[str, str]',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4,
        description="Number of panels per row.\nOriginal type annotation: <class 'int'>",
        title='Ncols',
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.\nOriginal type annotation: <class 'float'>",
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.\nOriginal type annotation: float | None',
        title='Wspace',
    )
    title: Any = Field(
        None,
        description="Provide title for panels either as string or list of strings,\ne.g. `['title1', 'title2', ...]`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: matplotlib.axes._axes.Axes | None',
        title='Ax',
    )
    return_fig: Any = Field(
        None,
        description='Return the matplotlib figure.\nOriginal type annotation: bool | None',
        title='Return Fig',
    )
    marker: Optional[Any] = Field(
        '.',
        description='No description available.\nOriginal type annotation: str | collections.abc.Sequence[str]',
        title='Marker',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.umap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingTsne(BaseAPI):
    """
    Scatter plot in tSNE basis.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    color
        Keys for annotations of observations/cells or variables/genes, e.g.,
        `'ann1'` or `['ann1', 'ann2']`.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols. By default `var_names`
        refer to the index column of the `.var` DataFrame. Setting this option allows
        alternative names to be used.
    use_raw
        Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,
        defaults to `True` if `layer` isn't provided and `adata.raw` is present.
    layer
        Name of the AnnData object layer that wants to be plotted. By default
        adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.
        If `layer` is set to a valid layer name, then the layer is plotted. `layer`
        takes precedence over `use_raw`.
    edges
        Show edges.
    edges_width
        Width of edges.
    edges_color
        Color of edges. See :func:`~networkx.drawing.nx_pylab.draw_networkx_edges`.
    neighbors_key
        Where to look for neighbors connectivities.
        If not specified, this looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, this looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    arrows
        Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).
    arrows_kwds
        Passed to :meth:`~matplotlib.axes.Axes.quiver`
    sort_order
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups
        Restrict to a few categories in categorical observation annotation.
        The default is not to restrict to any groups.
    dimensions
        0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].
        Unlike `components`, this argument is used in the same way as `colors`, e.g. is
        used to specify a single plot at a time. Will eventually replace the components
        argument.
    components
        For instance, `['1,2', '2,3']`. To plot all available components use
        `components='all'`.
    projection
        Projection of plot (default: `'2d'`).
    legend_loc
        Location of legend, either `'on data'`, `'right margin'`, `None`,
        or a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.
    legend_fontsize
        Numeric size in pt or string describing the size.
        See :meth:`~matplotlib.text.Text.set_fontsize`.
    legend_fontweight
        Legend font weight. A numeric value in range 0-1000 or a string.
        Defaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.
        See :meth:`~matplotlib.text.Text.set_fontweight`.
    legend_fontoutline
        Line width of the legend font outline in pt. Draws a white outline using
        the path effect :class:`~matplotlib.patheffects.withStroke`.
    colorbar_loc
        Where to place the colorbar for continous variables. If `None`, no colorbar
        is added.
    size
        Point size. If `None`, is automatically computed as 120000 / n_cells.
        Can be a sequence containing the size for each cell. The order should be
        the same as in adata.obs.
    color_map
        Color map to use for continous variables. Can be a name or a
        :class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`
        or `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.
        If `None`, the value of `mpl.rcParams["image.cmap"]` is used.
        The default `color_map` can be set using :func:`~scanpy.set_figure_params`.
    palette
        Colors to use for plotting categorical annotation groups.
        The palette can be a valid :class:`~matplotlib.colors.ListedColormap` name
        (`'Set2'`, `'tab20'`, …), a :class:`~cycler.Cycler` object, a dict mapping
        categories to colors, or a sequence of colors. Colors must be valid to
        matplotlib. (see :func:`~matplotlib.colors.is_color_like`).
        If `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
        variable already has colors stored in `adata.uns["{var}_colors"]`.
        If provided, values of `adata.uns["{var}_colors"]` will be set.
    na_color
        Color to use for null or masked values. Can be anything matplotlib accepts as a
        color. Used for all points if `color=None`.
    na_in_legend
        If there are missing values, whether they get an entry in the legend. Currently
        only implemented for categorical legends.
    frameon
        Draw a frame around the scatter plot. Defaults to value set in
        :func:`~scanpy.set_figure_params`, defaults to `True`.
    title
        Provide title for panels either as string or list of strings,
        e.g. `['title1', 'title2', ...]`.
    
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin. vmin can be a number, a string, a function or `None`. If
        vmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).
        For example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then
        vmin is interpreted as the return value of the function over the list of values to plot.
        For example to set vmin tp the mean of the values to plot, `def my_vmin(values): return
        np.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic
        minimum value is used as defined by matplotlib `scatter` function. When making multiple
        plots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`
    vmax
        The value representing the upper limit of the color scale. The format is the same as for `vmin`.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
        The format is the same as for `vmin`.
        Example: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``
    add_outline
        If set to True, this will add a thin border around groups of dots. In some situations
        this can enhance the aesthetics of the resulting image
    outline_color
        Tuple with two valid color names used to adjust the add_outline. The first color is the
        border color (default: black), while the second color is a gap color between the
        border color and the scatter dot (default: white).
    outline_width
        Tuple with two width numbers used to adjust the outline. The first value is the width
        of the border color as a fraction of the scatter dot size (default: 0.3). The second value is
        width of the gap color (default: 0.05).
    ncols
        Number of panels per row.
    wspace
        Adjust the width of the space between multiple panels.
    hspace
        Adjust the height of the space between multiple panels.
    return_fig
        Return the matplotlib figure.
    kwargs
        Arguments to pass to :func:`matplotlib.pyplot.scatter`,
        for instance: the maximum and minimum values (e.g. `vmin=-2, vmax=5`).
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    
    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    
    Examples
    --------
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.tsne(adata)
        sc.pl.tsne(adata, color='bulk_labels')
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    tl.tsne
    """

    
    adata: Any = Field(
        ...,
        description="Annotated data matrix.\nOriginal type annotation: <class 'anndata._core.anndata.AnnData'>",
        title='Adata',
    )
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes, e.g.,\n`'ann1'` or `['ann1', 'ann2']`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Color',
    )
    mask_obs: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: numpy.ndarray[typing.Any, numpy.dtype[numpy.bool_]] | str | None',
        title='Mask Obs',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols. By default `var_names`\nrefer to the index column of the `.var` DataFrame. Setting this option allows\nalternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,\ndefaults to `True` if `layer` isn't provided and `adata.raw` is present.\nOriginal type annotation: bool | None",
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points\nwith higher values on top of others.\nOriginal type annotation: <class 'bool'>",
        title='Sort Order',
    )
    edges: Optional[Any] = Field(
        False,
        description="Show edges.\nOriginal type annotation: <class 'bool'>",
        title='Edges',
    )
    edges_width: Optional[Any] = Field(
        0.1,
        description="Width of edges.\nOriginal type annotation: <class 'float'>",
        title='Edges Width',
    )
    edges_color: Optional[Any] = Field(
        'grey',
        description='Color of edges. See :func:`~networkx.drawing.nx_pylab.draw_networkx_edges`.\nOriginal type annotation: str | collections.abc.Sequence[float] | collections.abc.Sequence[str]',
        title='Edges Color',
    )
    neighbors_key: Any = Field(
        None,
        description="Where to look for neighbors connectivities.\nIf not specified, this looks .obsp['connectivities'] for connectivities\n(default storage place for pp.neighbors).\nIf specified, this looks\n.obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.\nOriginal type annotation: str | None",
        title='Neighbors Key',
    )
    arrows: Optional[Any] = Field(
        False,
        description="Show arrows (deprecated in favour of `scvelo.pl.velocity_embedding`).\nOriginal type annotation: <class 'bool'>",
        title='Arrows',
    )
    arrows_kwds: Any = Field(
        None,
        description='Passed to :meth:`~matplotlib.axes.Axes.quiver`\nOriginal type annotation: collections.abc.Mapping[str, typing.Any] | None',
        title='Arrows Kwds',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.\nThe default is not to restrict to any groups.\nOriginal type annotation: str | collections.abc.Sequence[str] | None',
        title='Groups',
    )
    components: Any = Field(
        None,
        description="For instance, `['1,2', '2,3']`. To plot all available components use\n`components='all'`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Components',
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].\nUnlike `components`, this argument is used in the same way as `colors`, e.g. is\nused to specify a single plot at a time. Will eventually replace the components\nargument.\nOriginal type annotation: tuple[int, int] | collections.abc.Sequence[tuple[int, int]] | None',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default\nadata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.\nIf `layer` is set to a valid layer name, then the layer is plotted. `layer`\ntakes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d',
        description="Projection of plot (default: `'2d'`).\nOriginal type annotation: typing.Literal['2d', '3d']",
        title='Projection',
    )
    scale_factor: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: float | None',
        title='Scale Factor',
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continous variables. Can be a name or a\n:class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`\nor `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.\nIf `None`, the value of `mpl.rcParams["image.cmap"]` is used.\nThe default `color_map` can be set using :func:`~scanpy.set_figure_params`.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Color Map',
    )
    cmap: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Cmap',
    )
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.\nThe palette can be a valid :class:`~matplotlib.colors.ListedColormap` name\n(`\'Set2\'`, `\'tab20\'`, …), a :class:`~cycler.Cycler` object, a dict mapping\ncategories to colors, or a sequence of colors. Colors must be valid to\nmatplotlib. (see :func:`~matplotlib.colors.is_color_like`).\nIf `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical\nvariable already has colors stored in `adata.uns["{var}_colors"]`.\nIf provided, values of `adata.uns["{var}_colors"]` will be set.\nOriginal type annotation: str | collections.abc.Sequence[str] | cycler.Cycler | None',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values. Can be anything matplotlib accepts as a\ncolor. Used for all points if `color=None`.\nOriginal type annotation: str | tuple[float, ...]',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="If there are missing values, whether they get an entry in the legend. Currently\nonly implemented for categorical legends.\nOriginal type annotation: <class 'bool'>",
        title='Na In Legend',
    )
    size: Any = Field(
        None,
        description='Point size. If `None`, is automatically computed as 120000 / n_cells.\nCan be a sequence containing the size for each cell. The order should be\nthe same as in adata.obs.\nOriginal type annotation: float | collections.abc.Sequence[float] | None',
        title='Size',
    )
    frameon: Any = Field(
        None,
        description='Draw a frame around the scatter plot. Defaults to value set in\n:func:`~scanpy.set_figure_params`, defaults to `True`.\nOriginal type annotation: bool | None',
        title='Frameon',
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string describing the size.\nSee :meth:`~matplotlib.text.Text.set_fontsize`.\nOriginal type annotation: typing.Union[float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], NoneType]",
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold',
        description="Legend font weight. A numeric value in range 0-1000 or a string.\nDefaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.\nSee :meth:`~matplotlib.text.Text.set_fontweight`.\nOriginal type annotation: typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]",
        title='Legend Fontweight',
    )
    legend_loc: Optional[Any] = Field(
        'right margin',
        description="Location of legend, either `'on data'`, `'right margin'`, `None`,\nor a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.\nOriginal type annotation: typing.Optional[typing.Literal['none', 'right margin', 'on data', 'on data export', 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center']]",
        title='Legend Loc',
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline in pt. Draws a white outline using\nthe path effect :class:`~matplotlib.patheffects.withStroke`.\nOriginal type annotation: int | None',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continous variables. If `None`, no colorbar\nis added.\nOriginal type annotation: str | None',
        title='Colorbar Loc',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. The format is the same as for `vmin`.\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None',
        title='Vmax',
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin. vmin can be a number, a string, a function or `None`. If\nvmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).\nFor example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then\nvmin is interpreted as the return value of the function over the list of values to plot.\nFor example to set vmin tp the mean of the values to plot, `def my_vmin(values): return\nnp.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic\nminimum value is used as defined by matplotlib `scatter` function. When making multiple\nplots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Useful for diverging colormaps.\nThe format is the same as for `vmin`.\nExample: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Normalize | collections.abc.Sequence[matplotlib.colors.Normalize] | None',
        title='Norm',
    )
    add_outline: Optional[Any] = Field(
        False,
        description='If set to True, this will add a thin border around groups of dots. In some situations\nthis can enhance the aesthetics of the resulting image\nOriginal type annotation: bool | None',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Tuple with two width numbers used to adjust the outline. The first value is the width\nof the border color as a fraction of the scatter dot size (default: 0.3). The second value is\nwidth of the gap color (default: 0.05).\nOriginal type annotation: tuple[float, float]',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Tuple with two valid color names used to adjust the add_outline. The first color is the\nborder color (default: black), while the second color is a gap color between the\nborder color and the scatter dot (default: white).\nOriginal type annotation: tuple[str, str]',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4,
        description="Number of panels per row.\nOriginal type annotation: <class 'int'>",
        title='Ncols',
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.\nOriginal type annotation: <class 'float'>",
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.\nOriginal type annotation: float | None',
        title='Wspace',
    )
    title: Any = Field(
        None,
        description="Provide title for panels either as string or list of strings,\ne.g. `['title1', 'title2', ...]`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: matplotlib.axes._axes.Axes | None',
        title='Ax',
    )
    return_fig: Any = Field(
        None,
        description='Return the matplotlib figure.\nOriginal type annotation: bool | None',
        title='Return Fig',
    )
    marker: Optional[Any] = Field(
        '.',
        description='No description available.\nOriginal type annotation: str | collections.abc.Sequence[str]',
        title='Marker',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.tsne')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingHeatmap(BaseAPI):
    """
    Heatmap of the expression values of genes.
    
    If `groupby` is given, the heatmap is ordered by the respective group. For
    example, a list of marker genes can be plotted, ordered by clustering. If
    the `groupby` observation annotation is not categorical the observation
    annotation is turned into a categorical by binning the data into the number
    specified in `num_categories`.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    var_names
        `var_names` should be a valid subset of `adata.var_names`.
        If `var_names` is a mapping, then the key is used as label
        to group the values (see `var_group_labels`). The mapping values
        should be sequences of valid `adata.var_names`. In this
        case either coloring or 'brackets' are used for the grouping
        of var names depending on the plot. When `var_names` is a mapping,
        then the `var_group_labels` and `var_group_positions` are set.
    groupby
        The key of the observation grouping to consider.
    use_raw
        Use `raw` attribute of `adata` if present.
    log
        Plot on logarithmic axis.
    num_categories
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    categories_order
        Order in which to show the categories. Note: add_dendrogram or add_totals
        can change the categories order.
    figsize
        Figure size when `multi_panel=True`.
        Otherwise the `rcParam['figure.figsize]` value is used.
        Format is (width, height)
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the hierarchical
        clustering between the `groupby` categories is added.
        The dendrogram information is computed using :func:`scanpy.tl.dendrogram`.
        If `tl.dendrogram` has not been called previously the function is called
        with default parameters.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols.
        By default `var_names` refer to the index column of the `.var` DataFrame.
        Setting this option allows alternative names to be used.
    var_group_positions
        Use this parameter to highlight groups of `var_names`.
        This will draw a 'bracket' or a color block between the given start and end
        positions. If the parameter `var_group_labels` is set, the corresponding
        labels are added on top/left. E.g. `var_group_positions=[(4,10)]`
        will add a bracket between the fourth `var_name` and the tenth `var_name`.
        By giving more positions, more brackets/color blocks are drawn.
    var_group_labels
        Labels for each of the `var_group_positions` that want to be highlighted.
    var_group_rotation
        Label rotation degrees.
        By default, labels larger than 4 characters are rotated 90 degrees.
    layer
        Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.
        If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,
        then the layer is plotted. `layer` takes precedence over `use_raw`.
    standard_scale
        Whether or not to standardize that dimension between 0 and 1, meaning for each variable or observation,
        subtract the minimum and divide each by its maximum.
    swap_axes
         By default, the x axis contains `var_names` (e.g. genes) and the y axis the `groupby`
         categories (if any). By setting `swap_axes` then x are the `groupby` categories and y the `var_names`.
    show_gene_labels
         By default gene labels are shown when there are 50 or less genes. Otherwise the labels are removed.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin.
    vmax
        The value representing the upper limit of the color scale. Values larger than vmax are plotted
        with the same color as vmax.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
    norm
        Custom color normalization object from matplotlib. See
        `https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.
    **kwds
        Are passed to :func:`matplotlib.pyplot.imshow`.
    
    Returns
    -------
    Dict of :class:`~matplotlib.axes.Axes`
    
    Examples
    --------
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.heatmap(adata, markers, groupby='bulk_labels', swap_axes=True)
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    pl.rank_genes_groups_heatmap
    tl.rank_genes_groups
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`.\nIf `var_names` is a mapping, then the key is used as label\nto group the values (see `var_group_labels`). The mapping values\nshould be sequences of valid `adata.var_names`. In this\ncase either coloring or 'brackets' are used for the grouping\nof var names depending on the plot. When `var_names` is a mapping,\nthen the `var_group_labels` and `var_group_positions` are set.\nOriginal type annotation: _VarNames | Mapping[str, _VarNames]",
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='The key of the observation grouping to consider.\nOriginal type annotation: str | Sequence[str]',
        title='Groupby',
    )
    use_raw: Any = Field(
        None,
        description='Use `raw` attribute of `adata` if present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axis.\nOriginal type annotation: bool',
        title='Log',
    )
    num_categories: Optional[Any] = Field(
        7,
        description='Only used if groupby observation is not categorical. This value\ndetermines the number of groups into which the groupby observation\nshould be subdivided.\nOriginal type annotation: int',
        title='Num Categories',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='If True or a valid dendrogram key, a dendrogram based on the hierarchical\nclustering between the `groupby` categories is added.\nThe dendrogram information is computed using :func:`scanpy.tl.dendrogram`.\nIf `tl.dendrogram` has not been called previously the function is called\nwith default parameters.\nOriginal type annotation: bool | str',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.\nBy default `var_names` refer to the index column of the `.var` DataFrame.\nSetting this option allows alternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description="Use this parameter to highlight groups of `var_names`.\nThis will draw a 'bracket' or a color block between the given start and end\npositions. If the parameter `var_group_labels` is set, the corresponding\nlabels are added on top/left. E.g. `var_group_positions=[(4,10)]`\nwill add a bracket between the fourth `var_name` and the tenth `var_name`.\nBy giving more positions, more brackets/color blocks are drawn.\nOriginal type annotation: Sequence[tuple[int, int]] | None",
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for each of the `var_group_positions` that want to be highlighted.\nOriginal type annotation: Sequence[str] | None',
        title='Var Group Labels',
    )
    var_group_rotation: Any = Field(
        None,
        description='Label rotation degrees.\nBy default, labels larger than 4 characters are rotated 90 degrees.\nOriginal type annotation: float | None',
        title='Var Group Rotation',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.\nIf `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,\nthen the layer is plotted. `layer` takes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    standard_scale: Any = Field(
        None,
        description="Whether or not to standardize that dimension between 0 and 1, meaning for each variable or observation,\nsubtract the minimum and divide each by its maximum.\nOriginal type annotation: Literal['var', 'obs'] | None",
        title='Standard Scale',
    )
    swap_axes: Optional[Any] = Field(
        False,
        description='By default, the x axis contains `var_names` (e.g. genes) and the y axis the `groupby`\ncategories (if any). By setting `swap_axes` then x are the `groupby` categories and y the `var_names`.\nOriginal type annotation: bool',
        title='Swap Axes',
    )
    show_gene_labels: Any = Field(
        None,
        description='By default gene labels are shown when there are 50 or less genes. Otherwise the labels are removed.\nOriginal type annotation: bool | None',
        title='Show Gene Labels',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: str | bool | None",
        title='Save',
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`.\nOtherwise the `rcParam['figure.figsize]` value is used.\nFormat is (width, height)\nOriginal type annotation: tuple[float, float] | None",
        title='Figsize',
    )
    vmin: Any = Field(
        None,
        description='The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin.\nOriginal type annotation: float | None',
        title='Vmin',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. Values larger than vmax are plotted\nwith the same color as vmax.\nOriginal type annotation: float | None',
        title='Vmax',
    )
    vcenter: Any = Field(
        None,
        description='The value representing the center of the color scale. Useful for diverging colormaps.\nOriginal type annotation: float | None',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib. See\n`https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.\nOriginal type annotation: Normalize | None',
        title='Norm',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.heatmap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingDotplot(BaseAPI):
    """
    Make a *dot plot* of the expression values of `var_names`.
    
    For each var_name and each `groupby` category a dot is plotted.
    Each dot represents two values: mean expression within each category
    (visualized by color) and fraction of cells expressing the `var_name` in the
    category (visualized by the size of the dot). If `groupby` is not given,
    the dotplot assumes that all data belongs to a single category.
    
    .. note::
       A gene is considered expressed if the expression value in the `adata` (or
       `adata.raw`) is above the specified threshold which is zero by default.
    
    An example of dotplot usage is to visualize, for multiple marker genes,
    the mean value and the percentage of cells expressing the gene
    across  multiple clusters.
    
    This function provides a convenient interface to the :class:`~scanpy.pl.DotPlot`
    class. If you need more flexibility, you should use :class:`~scanpy.pl.DotPlot`
    directly.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    var_names
        `var_names` should be a valid subset of `adata.var_names`.
        If `var_names` is a mapping, then the key is used as label
        to group the values (see `var_group_labels`). The mapping values
        should be sequences of valid `adata.var_names`. In this
        case either coloring or 'brackets' are used for the grouping
        of var names depending on the plot. When `var_names` is a mapping,
        then the `var_group_labels` and `var_group_positions` are set.
    groupby
        The key of the observation grouping to consider.
    use_raw
        Use `raw` attribute of `adata` if present.
    log
        Plot on logarithmic axis.
    num_categories
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    categories_order
        Order in which to show the categories. Note: add_dendrogram or add_totals
        can change the categories order.
    figsize
        Figure size when `multi_panel=True`.
        Otherwise the `rcParam['figure.figsize]` value is used.
        Format is (width, height)
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the hierarchical
        clustering between the `groupby` categories is added.
        The dendrogram information is computed using :func:`scanpy.tl.dendrogram`.
        If `tl.dendrogram` has not been called previously the function is called
        with default parameters.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols.
        By default `var_names` refer to the index column of the `.var` DataFrame.
        Setting this option allows alternative names to be used.
    var_group_positions
        Use this parameter to highlight groups of `var_names`.
        This will draw a 'bracket' or a color block between the given start and end
        positions. If the parameter `var_group_labels` is set, the corresponding
        labels are added on top/left. E.g. `var_group_positions=[(4,10)]`
        will add a bracket between the fourth `var_name` and the tenth `var_name`.
        By giving more positions, more brackets/color blocks are drawn.
    var_group_labels
        Labels for each of the `var_group_positions` that want to be highlighted.
    var_group_rotation
        Label rotation degrees.
        By default, labels larger than 4 characters are rotated 90 degrees.
    layer
        Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.
        If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,
        then the layer is plotted. `layer` takes precedence over `use_raw`.
    title
        Title for the figure
    colorbar_title
        Title for the color bar. New line character (\n) can be used.
    cmap
        String denoting matplotlib color map.
    standard_scale
        Whether or not to standardize the given dimension between 0 and 1, meaning for
        each variable or group, subtract the minimum and divide each by its maximum.
    swap_axes
         By default, the x axis contains `var_names` (e.g. genes) and the y axis
         the `groupby` categories. By setting `swap_axes` then x are the
         `groupby` categories and y the `var_names`.
    return_fig
        Returns :class:`DotPlot` object. Useful for fine-tuning
        the plot. Takes precedence over `show=False`.
    
    size_title
        Title for the size legend. New line character (\n) can be used.
    expression_cutoff
        Expression cutoff that is used for binarizing the gene expression and
        determining the fraction of cells expressing given genes. A gene is
        expressed only if the expression value is greater than this threshold.
    mean_only_expressed
        If True, gene expression is averaged only over the cells
        expressing the given genes.
    dot_max
        If ``None``, the maximum dot size is set to the maximum fraction value found
        (e.g. 0.6). If given, the value should be a number between 0 and 1.
        All fractions larger than dot_max are clipped to this value.
    dot_min
        If ``None``, the minimum dot size is set to 0. If given,
        the value should be a number between 0 and 1.
        All fractions smaller than dot_min are clipped to this value.
    smallest_dot
        All expression levels with `dot_min` are plotted with this size.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin.
    vmax
        The value representing the upper limit of the color scale. Values larger than vmax are plotted
        with the same color as vmax.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
    norm
        Custom color normalization object from matplotlib. See
        `https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.
    kwds
        Are passed to :func:`matplotlib.pyplot.scatter`.
    
    Returns
    -------
    If `return_fig` is `True`, returns a :class:`~scanpy.pl.DotPlot` object,
    else if `show` is false, return axes dict
    
    See Also
    --------
    :class:`~scanpy.pl.DotPlot`: The DotPlot class can be used to to control
        several visual parameters not available in this function.
    :func:`~scanpy.pl.rank_genes_groups_dotplot`: to plot marker genes
        identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    
    Examples
    --------
    Create a dot plot using the given markers and the PBMC example dataset grouped by
    the category 'bulk_labels'.
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)
    
    Using var_names as dict:
    
    .. plot::
        :context: close-figs
    
        markers = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
        sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)
    
    Get DotPlot object for fine tuning
    
    .. plot::
        :context: close-figs
    
        dp = sc.pl.dotplot(adata, markers, 'bulk_labels', return_fig=True)
        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).show()
    
    The axes used can be obtained using the get_axes() method
    
    .. code-block:: python
    
        axes_dict = dp.get_axes()
        print(axes_dict)
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`.\nIf `var_names` is a mapping, then the key is used as label\nto group the values (see `var_group_labels`). The mapping values\nshould be sequences of valid `adata.var_names`. In this\ncase either coloring or 'brackets' are used for the grouping\nof var names depending on the plot. When `var_names` is a mapping,\nthen the `var_group_labels` and `var_group_positions` are set.\nOriginal type annotation: _VarNames | Mapping[str, _VarNames]",
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='The key of the observation grouping to consider.\nOriginal type annotation: str | Sequence[str]',
        title='Groupby',
    )
    use_raw: Any = Field(
        None,
        description='Use `raw` attribute of `adata` if present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axis.\nOriginal type annotation: bool',
        title='Log',
    )
    num_categories: Optional[Any] = Field(
        7,
        description='Only used if groupby observation is not categorical. This value\ndetermines the number of groups into which the groupby observation\nshould be subdivided.\nOriginal type annotation: int',
        title='Num Categories',
    )
    categories_order: Any = Field(
        None,
        description='Order in which to show the categories. Note: add_dendrogram or add_totals\ncan change the categories order.\nOriginal type annotation: Sequence[str] | None',
        title='Categories Order',
    )
    expression_cutoff: Optional[Any] = Field(
        0.0,
        description='Expression cutoff that is used for binarizing the gene expression and\ndetermining the fraction of cells expressing given genes. A gene is\nexpressed only if the expression value is greater than this threshold.\nOriginal type annotation: float',
        title='Expression Cutoff',
    )
    mean_only_expressed: Optional[Any] = Field(
        False,
        description='If True, gene expression is averaged only over the cells\nexpressing the given genes.\nOriginal type annotation: bool',
        title='Mean Only Expressed',
    )
    standard_scale: Any = Field(
        None,
        description="Whether or not to standardize the given dimension between 0 and 1, meaning for\neach variable or group, subtract the minimum and divide each by its maximum.\nOriginal type annotation: Literal['var', 'group'] | None",
        title='Standard Scale',
    )
    title: Any = Field(
        None,
        description='Title for the figure\nOriginal type annotation: str | None',
        title='Title',
    )
    colorbar_title: Optional[Any] = Field(
        'Mean expression\nin group',
        description='Title for the color bar. New line character (\\n) can be used.\nOriginal type annotation: str | None',
        title='Colorbar Title',
    )
    size_title: Optional[Any] = Field(
        'Fraction of cells\nin group (%)',
        description='Title for the size legend. New line character (\\n) can be used.\nOriginal type annotation: str | None',
        title='Size Title',
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`.\nOtherwise the `rcParam['figure.figsize]` value is used.\nFormat is (width, height)\nOriginal type annotation: tuple[float, float] | None",
        title='Figsize',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='If True or a valid dendrogram key, a dendrogram based on the hierarchical\nclustering between the `groupby` categories is added.\nThe dendrogram information is computed using :func:`scanpy.tl.dendrogram`.\nIf `tl.dendrogram` has not been called previously the function is called\nwith default parameters.\nOriginal type annotation: bool | str',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.\nBy default `var_names` refer to the index column of the `.var` DataFrame.\nSetting this option allows alternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description="Use this parameter to highlight groups of `var_names`.\nThis will draw a 'bracket' or a color block between the given start and end\npositions. If the parameter `var_group_labels` is set, the corresponding\nlabels are added on top/left. E.g. `var_group_positions=[(4,10)]`\nwill add a bracket between the fourth `var_name` and the tenth `var_name`.\nBy giving more positions, more brackets/color blocks are drawn.\nOriginal type annotation: Sequence[tuple[int, int]] | None",
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for each of the `var_group_positions` that want to be highlighted.\nOriginal type annotation: Sequence[str] | None',
        title='Var Group Labels',
    )
    var_group_rotation: Any = Field(
        None,
        description='Label rotation degrees.\nBy default, labels larger than 4 characters are rotated 90 degrees.\nOriginal type annotation: float | None',
        title='Var Group Rotation',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.\nIf `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,\nthen the layer is plotted. `layer` takes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    swap_axes: Optional[Any] = Field(
        False,
        description='By default, the x axis contains `var_names` (e.g. genes) and the y axis\nthe `groupby` categories. By setting `swap_axes` then x are the\n`groupby` categories and y the `var_names`.\nOriginal type annotation: bool | None',
        title='Swap Axes',
    )
    dot_color_df: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: pd.DataFrame | None',
        title='Dot Color Df',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: _AxesSubplot | None',
        title='Ax',
    )
    return_fig: Optional[Any] = Field(
        False,
        description='Returns :class:`DotPlot` object. Useful for fine-tuning\nthe plot. Takes precedence over `show=False`.\nOriginal type annotation: bool | None',
        title='Return Fig',
    )
    vmin: Any = Field(
        None,
        description='The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin.\nOriginal type annotation: float | None',
        title='Vmin',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. Values larger than vmax are plotted\nwith the same color as vmax.\nOriginal type annotation: float | None',
        title='Vmax',
    )
    vcenter: Any = Field(
        None,
        description='The value representing the center of the color scale. Useful for diverging colormaps.\nOriginal type annotation: float | None',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib. See\n`https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.\nOriginal type annotation: Normalize | None',
        title='Norm',
    )
    cmap: Optional[Any] = Field(
        'Reds',
        description='String denoting matplotlib color map.\nOriginal type annotation: Colormap | str | None',
        title='Cmap',
    )
    dot_max: Any = Field(
        None,
        description='If ``None``, the maximum dot size is set to the maximum fraction value found\n(e.g. 0.6). If given, the value should be a number between 0 and 1.\nAll fractions larger than dot_max are clipped to this value.\nOriginal type annotation: float | None',
        title='Dot Max',
    )
    dot_min: Any = Field(
        None,
        description='If ``None``, the minimum dot size is set to 0. If given,\nthe value should be a number between 0 and 1.\nAll fractions smaller than dot_min are clipped to this value.\nOriginal type annotation: float | None',
        title='Dot Min',
    )
    smallest_dot: Optional[Any] = Field(
        0.0,
        description='All expression levels with `dot_min` are plotted with this size.\nOriginal type annotation: float',
        title='Smallest Dot',
    )
    kwds: Any = Field(
        ...,
        description='Are passed to :func:`matplotlib.pyplot.scatter`.',
        title='Kwds',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.dotplot')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingViolin(BaseAPI):
    """
    Violin plot.
    
    Wraps :func:`seaborn.violinplot` for :class:`~anndata.AnnData`.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    keys
        Keys for accessing variables of `.var_names` or fields of `.obs`.
    groupby
        The key of the observation grouping to consider.
    log
        Plot on logarithmic axis.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
    stripplot
        Add a stripplot on top of the violin plot.
        See :func:`~seaborn.stripplot`.
    jitter
        Add jitter to the stripplot (only when stripplot is True)
        See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    layer
        Name of the AnnData object layer that wants to be plotted. By
        default adata.raw.X is plotted. If `use_raw=False` is set,
        then `adata.X` is plotted. If `layer` is set to a valid layer name,
        then the layer is plotted. `layer` takes precedence over `use_raw`.
    density_norm
        The method used to scale the width of each violin.
        If 'width' (the default), each violin will have the same width.
        If 'area', each violin will have the same area.
        If 'count', a violin’s width corresponds to the number of observations.
    order
        Order in which to show the categories.
    multi_panel
        Display keys in multiple panels also when `groupby is not None`.
    xlabel
        Label of the x axis. Defaults to `groupby` if `rotation` is `None`,
        otherwise, no label is shown.
    ylabel
        Label of the y axis. If `None` and `groupby` is `None`, defaults
        to `'value'`. If `None` and `groubpy` is not `None`, defaults to `keys`.
    rotation
        Rotation of xtick labels.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    **kwds
        Are passed to :func:`~seaborn.violinplot`.
    
    Returns
    -------
    A :class:`~matplotlib.axes.Axes` object if `ax` is `None` else `None`.
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.pl.violin(adata, keys='S_score')
    
    Plot by category. Rotate x-axis labels so that they do not overlap.
    
    .. plot::
        :context: close-figs
    
        sc.pl.violin(adata, keys='S_score', groupby='bulk_labels', rotation=90)
    
    Set order of categories to be plotted or select specific categories to be plotted.
    
    .. plot::
        :context: close-figs
    
        groupby_order = ['CD34+', 'CD19+ B']
        sc.pl.violin(adata, keys='S_score', groupby='bulk_labels', rotation=90,
            order=groupby_order)
    
    Plot multiple keys.
    
    .. plot::
        :context: close-figs
    
        sc.pl.violin(adata, keys=['S_score', 'G2M_score'], groupby='bulk_labels',
            rotation=90)
    
    For large datasets consider omitting the overlaid scatter plot.
    
    .. plot::
        :context: close-figs
    
        sc.pl.violin(adata, keys='S_score', stripplot=False)
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    pl.stacked_violin
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    keys: Any = Field(
        ...,
        description='Keys for accessing variables of `.var_names` or fields of `.obs`.\nOriginal type annotation: str | Sequence[str]',
        title='Keys',
    )
    groupby: Any = Field(
        None,
        description='The key of the observation grouping to consider.\nOriginal type annotation: str | None',
        title='Groupby',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axis.\nOriginal type annotation: bool',
        title='Log',
    )
    use_raw: Any = Field(
        None,
        description='Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    stripplot: Optional[Any] = Field(
        True,
        description='Add a stripplot on top of the violin plot.\nSee :func:`~seaborn.stripplot`.\nOriginal type annotation: bool',
        title='Stripplot',
    )
    jitter: Optional[Any] = Field(
        True,
        description='Add jitter to the stripplot (only when stripplot is True)\nSee :func:`~seaborn.stripplot`.\nOriginal type annotation: float | bool',
        title='Jitter',
    )
    size: Optional[Any] = Field(
        1,
        description='Size of the jitter points.\nOriginal type annotation: int',
        title='Size',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By\ndefault adata.raw.X is plotted. If `use_raw=False` is set,\nthen `adata.X` is plotted. If `layer` is set to a valid layer name,\nthen the layer is plotted. `layer` takes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    density_norm: Optional[Any] = Field(
        'width',
        description="The method used to scale the width of each violin.\nIf 'width' (the default), each violin will have the same width.\nIf 'area', each violin will have the same area.\nIf 'count', a violin’s width corresponds to the number of observations.\nOriginal type annotation: DensityNorm",
        title='Density Norm',
    )
    order: Any = Field(
        None,
        description='Order in which to show the categories.\nOriginal type annotation: Sequence[str] | None',
        title='Order',
    )
    multi_panel: Any = Field(
        None,
        description='Display keys in multiple panels also when `groupby is not None`.\nOriginal type annotation: bool | None',
        title='Multi Panel',
    )
    xlabel: Optional[Any] = Field(
        '',
        description='Label of the x axis. Defaults to `groupby` if `rotation` is `None`,\notherwise, no label is shown.\nOriginal type annotation: str',
        title='Xlabel',
    )
    ylabel: Any = Field(
        None,
        description="Label of the y axis. If `None` and `groupby` is `None`, defaults\nto `'value'`. If `None` and `groubpy` is not `None`, defaults to `keys`.\nOriginal type annotation: str | Sequence[str] | None",
        title='Ylabel',
    )
    rotation: Any = Field(
        None,
        description='Rotation of xtick labels.\nOriginal type annotation: float | None',
        title='Rotation',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: Axes | None',
        title='Ax',
    )
    scale: Optional[Any] = Field(
        0,
        description='No description available.\nOriginal type annotation: DensityNorm | Empty',
        title='Scale',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.violin')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingDendrogram(BaseAPI):
    """
    Plot a dendrogram of the categories defined in `groupby`.
    
    See :func:`~scanpy.tl.dendrogram`.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        Categorical data column used to create the dendrogram
    dendrogram_key
        Key under with the dendrogram information was stored.
        By default the dendrogram information is stored under
        `.uns[f'dendrogram_{groupby}']`.
    orientation
        Origin of the tree. Will grow into the opposite direction.
    remove_labels
        Don’t draw labels. Used e.g. by :func:`scanpy.pl.matrixplot`
        to annotate matrix columns/rows.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    
    Returns
    -------
    :class:`matplotlib.axes.Axes`
    
    Examples
    --------
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.dendrogram(adata, 'bulk_labels')
        sc.pl.dendrogram(adata, 'bulk_labels')
    
    .. currentmodule:: scanpy
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    groupby: Any = Field(
        ...,
        description='Categorical data column used to create the dendrogram\nOriginal type annotation: str',
        title='Groupby',
    )
    dendrogram_key: Any = Field(
        None,
        description="Key under with the dendrogram information was stored.\nBy default the dendrogram information is stored under\n`.uns[f'dendrogram_{groupby}']`.\nOriginal type annotation: str | None",
        title='Dendrogram Key',
    )
    orientation: Optional[Any] = Field(
        'top',
        description="Origin of the tree. Will grow into the opposite direction.\nOriginal type annotation: Literal['top', 'bottom', 'left', 'right']",
        title='Orientation',
    )
    remove_labels: Optional[Any] = Field(
        False,
        description='Don’t draw labels. Used e.g. by :func:`scanpy.pl.matrixplot`\nto annotate matrix columns/rows.\nOriginal type annotation: bool',
        title='Remove Labels',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: Axes | None',
        title='Ax',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.dendrogram')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingDiffmap(BaseAPI):
    """
    Scatter plot in Diffusion Map basis.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    color
        Keys for annotations of observations/cells or variables/genes, e.g.,
        `'ann1'` or `['ann1', 'ann2']`.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols. By default `var_names`
        refer to the index column of the `.var` DataFrame. Setting this option allows
        alternative names to be used.
    use_raw
        Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,
        defaults to `True` if `layer` isn't provided and `adata.raw` is present.
    layer
        Name of the AnnData object layer that wants to be plotted. By default
        adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.
        If `layer` is set to a valid layer name, then the layer is plotted. `layer`
        takes precedence over `use_raw`.
    sort_order
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups
        Restrict to a few categories in categorical observation annotation.
        The default is not to restrict to any groups.
    dimensions
        0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].
        Unlike `components`, this argument is used in the same way as `colors`, e.g. is
        used to specify a single plot at a time. Will eventually replace the components
        argument.
    components
        For instance, `['1,2', '2,3']`. To plot all available components use
        `components='all'`.
    projection
        Projection of plot (default: `'2d'`).
    legend_loc
        Location of legend, either `'on data'`, `'right margin'`, `None`,
        or a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.
    legend_fontsize
        Numeric size in pt or string describing the size.
        See :meth:`~matplotlib.text.Text.set_fontsize`.
    legend_fontweight
        Legend font weight. A numeric value in range 0-1000 or a string.
        Defaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.
        See :meth:`~matplotlib.text.Text.set_fontweight`.
    legend_fontoutline
        Line width of the legend font outline in pt. Draws a white outline using
        the path effect :class:`~matplotlib.patheffects.withStroke`.
    colorbar_loc
        Where to place the colorbar for continous variables. If `None`, no colorbar
        is added.
    size
        Point size. If `None`, is automatically computed as 120000 / n_cells.
        Can be a sequence containing the size for each cell. The order should be
        the same as in adata.obs.
    color_map
        Color map to use for continous variables. Can be a name or a
        :class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`
        or `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.
        If `None`, the value of `mpl.rcParams["image.cmap"]` is used.
        The default `color_map` can be set using :func:`~scanpy.set_figure_params`.
    palette
        Colors to use for plotting categorical annotation groups.
        The palette can be a valid :class:`~matplotlib.colors.ListedColormap` name
        (`'Set2'`, `'tab20'`, …), a :class:`~cycler.Cycler` object, a dict mapping
        categories to colors, or a sequence of colors. Colors must be valid to
        matplotlib. (see :func:`~matplotlib.colors.is_color_like`).
        If `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
        variable already has colors stored in `adata.uns["{var}_colors"]`.
        If provided, values of `adata.uns["{var}_colors"]` will be set.
    na_color
        Color to use for null or masked values. Can be anything matplotlib accepts as a
        color. Used for all points if `color=None`.
    na_in_legend
        If there are missing values, whether they get an entry in the legend. Currently
        only implemented for categorical legends.
    frameon
        Draw a frame around the scatter plot. Defaults to value set in
        :func:`~scanpy.set_figure_params`, defaults to `True`.
    title
        Provide title for panels either as string or list of strings,
        e.g. `['title1', 'title2', ...]`.
    
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin. vmin can be a number, a string, a function or `None`. If
        vmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).
        For example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then
        vmin is interpreted as the return value of the function over the list of values to plot.
        For example to set vmin tp the mean of the values to plot, `def my_vmin(values): return
        np.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic
        minimum value is used as defined by matplotlib `scatter` function. When making multiple
        plots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`
    vmax
        The value representing the upper limit of the color scale. The format is the same as for `vmin`.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
        The format is the same as for `vmin`.
        Example: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``
    add_outline
        If set to True, this will add a thin border around groups of dots. In some situations
        this can enhance the aesthetics of the resulting image
    outline_color
        Tuple with two valid color names used to adjust the add_outline. The first color is the
        border color (default: black), while the second color is a gap color between the
        border color and the scatter dot (default: white).
    outline_width
        Tuple with two width numbers used to adjust the outline. The first value is the width
        of the border color as a fraction of the scatter dot size (default: 0.3). The second value is
        width of the gap color (default: 0.05).
    ncols
        Number of panels per row.
    wspace
        Adjust the width of the space between multiple panels.
    hspace
        Adjust the height of the space between multiple panels.
    return_fig
        Return the matplotlib figure.
    kwargs
        Arguments to pass to :func:`matplotlib.pyplot.scatter`,
        for instance: the maximum and minimum values (e.g. `vmin=-2, vmax=5`).
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    
    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    
    Examples
    --------
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.diffmap(adata)
        sc.pl.diffmap(adata, color='bulk_labels')
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    tl.diffmap
    """

    
    adata: Any = Field(
        ...,
        description="Annotated data matrix.\nOriginal type annotation: <class 'anndata._core.anndata.AnnData'>",
        title='Adata',
    )
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes, e.g.,\n`'ann1'` or `['ann1', 'ann2']`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Color',
    )
    mask_obs: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: numpy.ndarray[typing.Any, numpy.dtype[numpy.bool_]] | str | None',
        title='Mask Obs',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols. By default `var_names`\nrefer to the index column of the `.var` DataFrame. Setting this option allows\nalternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,\ndefaults to `True` if `layer` isn't provided and `adata.raw` is present.\nOriginal type annotation: bool | None",
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points\nwith higher values on top of others.\nOriginal type annotation: <class 'bool'>",
        title='Sort Order',
    )
    edges: Optional[Any] = Field(
        False,
        description="No description available.\nOriginal type annotation: <class 'bool'>",
        title='Edges',
    )
    edges_width: Optional[Any] = Field(
        0.1,
        description="No description available.\nOriginal type annotation: <class 'float'>",
        title='Edges Width',
    )
    edges_color: Optional[Any] = Field(
        'grey',
        description='No description available.\nOriginal type annotation: str | collections.abc.Sequence[float] | collections.abc.Sequence[str]',
        title='Edges Color',
    )
    neighbors_key: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: str | None',
        title='Neighbors Key',
    )
    arrows: Optional[Any] = Field(
        False,
        description="No description available.\nOriginal type annotation: <class 'bool'>",
        title='Arrows',
    )
    arrows_kwds: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: collections.abc.Mapping[str, typing.Any] | None',
        title='Arrows Kwds',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.\nThe default is not to restrict to any groups.\nOriginal type annotation: str | collections.abc.Sequence[str] | None',
        title='Groups',
    )
    components: Any = Field(
        None,
        description="For instance, `['1,2', '2,3']`. To plot all available components use\n`components='all'`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Components',
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].\nUnlike `components`, this argument is used in the same way as `colors`, e.g. is\nused to specify a single plot at a time. Will eventually replace the components\nargument.\nOriginal type annotation: tuple[int, int] | collections.abc.Sequence[tuple[int, int]] | None',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default\nadata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.\nIf `layer` is set to a valid layer name, then the layer is plotted. `layer`\ntakes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d',
        description="Projection of plot (default: `'2d'`).\nOriginal type annotation: typing.Literal['2d', '3d']",
        title='Projection',
    )
    scale_factor: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: float | None',
        title='Scale Factor',
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continous variables. Can be a name or a\n:class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`\nor `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.\nIf `None`, the value of `mpl.rcParams["image.cmap"]` is used.\nThe default `color_map` can be set using :func:`~scanpy.set_figure_params`.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Color Map',
    )
    cmap: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Cmap',
    )
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.\nThe palette can be a valid :class:`~matplotlib.colors.ListedColormap` name\n(`\'Set2\'`, `\'tab20\'`, …), a :class:`~cycler.Cycler` object, a dict mapping\ncategories to colors, or a sequence of colors. Colors must be valid to\nmatplotlib. (see :func:`~matplotlib.colors.is_color_like`).\nIf `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical\nvariable already has colors stored in `adata.uns["{var}_colors"]`.\nIf provided, values of `adata.uns["{var}_colors"]` will be set.\nOriginal type annotation: str | collections.abc.Sequence[str] | cycler.Cycler | None',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values. Can be anything matplotlib accepts as a\ncolor. Used for all points if `color=None`.\nOriginal type annotation: str | tuple[float, ...]',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="If there are missing values, whether they get an entry in the legend. Currently\nonly implemented for categorical legends.\nOriginal type annotation: <class 'bool'>",
        title='Na In Legend',
    )
    size: Any = Field(
        None,
        description='Point size. If `None`, is automatically computed as 120000 / n_cells.\nCan be a sequence containing the size for each cell. The order should be\nthe same as in adata.obs.\nOriginal type annotation: float | collections.abc.Sequence[float] | None',
        title='Size',
    )
    frameon: Any = Field(
        None,
        description='Draw a frame around the scatter plot. Defaults to value set in\n:func:`~scanpy.set_figure_params`, defaults to `True`.\nOriginal type annotation: bool | None',
        title='Frameon',
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string describing the size.\nSee :meth:`~matplotlib.text.Text.set_fontsize`.\nOriginal type annotation: typing.Union[float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], NoneType]",
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold',
        description="Legend font weight. A numeric value in range 0-1000 or a string.\nDefaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.\nSee :meth:`~matplotlib.text.Text.set_fontweight`.\nOriginal type annotation: typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]",
        title='Legend Fontweight',
    )
    legend_loc: Optional[Any] = Field(
        'right margin',
        description="Location of legend, either `'on data'`, `'right margin'`, `None`,\nor a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.\nOriginal type annotation: typing.Optional[typing.Literal['none', 'right margin', 'on data', 'on data export', 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center']]",
        title='Legend Loc',
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline in pt. Draws a white outline using\nthe path effect :class:`~matplotlib.patheffects.withStroke`.\nOriginal type annotation: int | None',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continous variables. If `None`, no colorbar\nis added.\nOriginal type annotation: str | None',
        title='Colorbar Loc',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. The format is the same as for `vmin`.\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None',
        title='Vmax',
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin. vmin can be a number, a string, a function or `None`. If\nvmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).\nFor example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then\nvmin is interpreted as the return value of the function over the list of values to plot.\nFor example to set vmin tp the mean of the values to plot, `def my_vmin(values): return\nnp.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic\nminimum value is used as defined by matplotlib `scatter` function. When making multiple\nplots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Useful for diverging colormaps.\nThe format is the same as for `vmin`.\nExample: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Normalize | collections.abc.Sequence[matplotlib.colors.Normalize] | None',
        title='Norm',
    )
    add_outline: Optional[Any] = Field(
        False,
        description='If set to True, this will add a thin border around groups of dots. In some situations\nthis can enhance the aesthetics of the resulting image\nOriginal type annotation: bool | None',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Tuple with two width numbers used to adjust the outline. The first value is the width\nof the border color as a fraction of the scatter dot size (default: 0.3). The second value is\nwidth of the gap color (default: 0.05).\nOriginal type annotation: tuple[float, float]',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Tuple with two valid color names used to adjust the add_outline. The first color is the\nborder color (default: black), while the second color is a gap color between the\nborder color and the scatter dot (default: white).\nOriginal type annotation: tuple[str, str]',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4,
        description="Number of panels per row.\nOriginal type annotation: <class 'int'>",
        title='Ncols',
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.\nOriginal type annotation: <class 'float'>",
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.\nOriginal type annotation: float | None',
        title='Wspace',
    )
    title: Any = Field(
        None,
        description="Provide title for panels either as string or list of strings,\ne.g. `['title1', 'title2', ...]`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: matplotlib.axes._axes.Axes | None',
        title='Ax',
    )
    return_fig: Any = Field(
        None,
        description='Return the matplotlib figure.\nOriginal type annotation: bool | None',
        title='Return Fig',
    )
    marker: Optional[Any] = Field(
        '.',
        description='No description available.\nOriginal type annotation: str | collections.abc.Sequence[str]',
        title='Marker',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.diffmap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingHighlyVariableGenes(BaseAPI):
    """
    Plot dispersions or normalized variance versus means for genes.
    
    Produces Supp. Fig. 5c of Zheng et al. (2017) and MeanVarPlot() and
    VariableFeaturePlot() of Seurat.
    
    Parameters
    ----------
    adata
        Result of :func:`~scanpy.pp.highly_variable_genes`.
    log
        Plot on logarithmic axes.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.
    """

    
    adata_or_result: Any = Field(
        ...,
        description='No description available.\nOriginal type annotation: AnnData | pd.DataFrame | np.recarray',
        title='Adata Or Result',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axes.\nOriginal type annotation: bool',
        title='Log',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    highly_variable_genes: Optional[Any] = Field(
        True,
        description='No description available.\nOriginal type annotation: bool',
        title='Highly Variable Genes',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.highly_variable_genes')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingPca(BaseAPI):
    """
    Scatter plot in PCA coordinates.
    
    Use the parameter `annotate_var_explained` to annotate the explained variance.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    color
        Keys for annotations of observations/cells or variables/genes, e.g.,
        `'ann1'` or `['ann1', 'ann2']`.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols. By default `var_names`
        refer to the index column of the `.var` DataFrame. Setting this option allows
        alternative names to be used.
    use_raw
        Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,
        defaults to `True` if `layer` isn't provided and `adata.raw` is present.
    layer
        Name of the AnnData object layer that wants to be plotted. By default
        adata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.
        If `layer` is set to a valid layer name, then the layer is plotted. `layer`
        takes precedence over `use_raw`.
    annotate_var_explained
    sort_order
        For continuous annotations used as color parameter, plot data points
        with higher values on top of others.
    groups
        Restrict to a few categories in categorical observation annotation.
        The default is not to restrict to any groups.
    dimensions
        0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].
        Unlike `components`, this argument is used in the same way as `colors`, e.g. is
        used to specify a single plot at a time. Will eventually replace the components
        argument.
    components
        For instance, `['1,2', '2,3']`. To plot all available components use
        `components='all'`.
    projection
        Projection of plot (default: `'2d'`).
    legend_loc
        Location of legend, either `'on data'`, `'right margin'`, `None`,
        or a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.
    legend_fontsize
        Numeric size in pt or string describing the size.
        See :meth:`~matplotlib.text.Text.set_fontsize`.
    legend_fontweight
        Legend font weight. A numeric value in range 0-1000 or a string.
        Defaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.
        See :meth:`~matplotlib.text.Text.set_fontweight`.
    legend_fontoutline
        Line width of the legend font outline in pt. Draws a white outline using
        the path effect :class:`~matplotlib.patheffects.withStroke`.
    colorbar_loc
        Where to place the colorbar for continous variables. If `None`, no colorbar
        is added.
    size
        Point size. If `None`, is automatically computed as 120000 / n_cells.
        Can be a sequence containing the size for each cell. The order should be
        the same as in adata.obs.
    color_map
        Color map to use for continous variables. Can be a name or a
        :class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`
        or `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.
        If `None`, the value of `mpl.rcParams["image.cmap"]` is used.
        The default `color_map` can be set using :func:`~scanpy.set_figure_params`.
    palette
        Colors to use for plotting categorical annotation groups.
        The palette can be a valid :class:`~matplotlib.colors.ListedColormap` name
        (`'Set2'`, `'tab20'`, …), a :class:`~cycler.Cycler` object, a dict mapping
        categories to colors, or a sequence of colors. Colors must be valid to
        matplotlib. (see :func:`~matplotlib.colors.is_color_like`).
        If `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical
        variable already has colors stored in `adata.uns["{var}_colors"]`.
        If provided, values of `adata.uns["{var}_colors"]` will be set.
    na_color
        Color to use for null or masked values. Can be anything matplotlib accepts as a
        color. Used for all points if `color=None`.
    na_in_legend
        If there are missing values, whether they get an entry in the legend. Currently
        only implemented for categorical legends.
    frameon
        Draw a frame around the scatter plot. Defaults to value set in
        :func:`~scanpy.set_figure_params`, defaults to `True`.
    title
        Provide title for panels either as string or list of strings,
        e.g. `['title1', 'title2', ...]`.
    
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin. vmin can be a number, a string, a function or `None`. If
        vmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).
        For example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then
        vmin is interpreted as the return value of the function over the list of values to plot.
        For example to set vmin tp the mean of the values to plot, `def my_vmin(values): return
        np.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic
        minimum value is used as defined by matplotlib `scatter` function. When making multiple
        plots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`
    vmax
        The value representing the upper limit of the color scale. The format is the same as for `vmin`.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
        The format is the same as for `vmin`.
        Example: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``
    add_outline
        If set to True, this will add a thin border around groups of dots. In some situations
        this can enhance the aesthetics of the resulting image
    outline_color
        Tuple with two valid color names used to adjust the add_outline. The first color is the
        border color (default: black), while the second color is a gap color between the
        border color and the scatter dot (default: white).
    outline_width
        Tuple with two width numbers used to adjust the outline. The first value is the width
        of the border color as a fraction of the scatter dot size (default: 0.3). The second value is
        width of the gap color (default: 0.05).
    ncols
        Number of panels per row.
    wspace
        Adjust the width of the space between multiple panels.
    hspace
        Adjust the height of the space between multiple panels.
    return_fig
        Return the matplotlib figure.
    kwargs
        Arguments to pass to :func:`matplotlib.pyplot.scatter`,
        for instance: the maximum and minimum values (e.g. `vmin=-2, vmax=5`).
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    
    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc3k_processed()
        sc.pl.pca(adata)
    
    Colour points by discrete variable (Louvain clusters).
    
    .. plot::
        :context: close-figs
    
        sc.pl.pca(adata, color="louvain")
    
    Colour points by gene expression.
    
    .. plot::
        :context: close-figs
    
        sc.pl.pca(adata, color="CST3")
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    pp.pca
    """

    
    adata: Any = Field(
        ...,
        description="Annotated data matrix.\nOriginal type annotation: <class 'anndata._core.anndata.AnnData'>",
        title='Adata',
    )
    color: Any = Field(
        None,
        description="Keys for annotations of observations/cells or variables/genes, e.g.,\n`'ann1'` or `['ann1', 'ann2']`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Color',
    )
    mask_obs: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: numpy.ndarray[typing.Any, numpy.dtype[numpy.bool_]] | str | None',
        title='Mask Obs',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols. By default `var_names`\nrefer to the index column of the `.var` DataFrame. Setting this option allows\nalternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    use_raw: Any = Field(
        None,
        description="Use `.raw` attribute of `adata` for coloring with gene expression. If `None`,\ndefaults to `True` if `layer` isn't provided and `adata.raw` is present.\nOriginal type annotation: bool | None",
        title='Use Raw',
    )
    sort_order: Optional[Any] = Field(
        True,
        description="For continuous annotations used as color parameter, plot data points\nwith higher values on top of others.\nOriginal type annotation: <class 'bool'>",
        title='Sort Order',
    )
    edges: Optional[Any] = Field(
        False,
        description="No description available.\nOriginal type annotation: <class 'bool'>",
        title='Edges',
    )
    edges_width: Optional[Any] = Field(
        0.1,
        description="No description available.\nOriginal type annotation: <class 'float'>",
        title='Edges Width',
    )
    edges_color: Optional[Any] = Field(
        'grey',
        description='No description available.\nOriginal type annotation: str | collections.abc.Sequence[float] | collections.abc.Sequence[str]',
        title='Edges Color',
    )
    neighbors_key: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: str | None',
        title='Neighbors Key',
    )
    arrows: Optional[Any] = Field(
        False,
        description="No description available.\nOriginal type annotation: <class 'bool'>",
        title='Arrows',
    )
    arrows_kwds: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: collections.abc.Mapping[str, typing.Any] | None',
        title='Arrows Kwds',
    )
    groups: Any = Field(
        None,
        description='Restrict to a few categories in categorical observation annotation.\nThe default is not to restrict to any groups.\nOriginal type annotation: str | collections.abc.Sequence[str] | None',
        title='Groups',
    )
    components: Any = Field(
        None,
        description="For instance, `['1,2', '2,3']`. To plot all available components use\n`components='all'`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Components',
    )
    dimensions: Any = Field(
        None,
        description='0-indexed dimensions of the embedding to plot as integers. E.g. [(0, 1), (1, 2)].\nUnlike `components`, this argument is used in the same way as `colors`, e.g. is\nused to specify a single plot at a time. Will eventually replace the components\nargument.\nOriginal type annotation: tuple[int, int] | collections.abc.Sequence[tuple[int, int]] | None',
        title='Dimensions',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default\nadata.raw.X is plotted. If `use_raw=False` is set, then `adata.X` is plotted.\nIf `layer` is set to a valid layer name, then the layer is plotted. `layer`\ntakes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    projection: Optional[Any] = Field(
        '2d',
        description="Projection of plot (default: `'2d'`).\nOriginal type annotation: typing.Literal['2d', '3d']",
        title='Projection',
    )
    scale_factor: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: float | None',
        title='Scale Factor',
    )
    color_map: Any = Field(
        None,
        description='Color map to use for continous variables. Can be a name or a\n:class:`~matplotlib.colors.Colormap` instance (e.g. `"magma`", `"viridis"`\nor `mpl.cm.cividis`), see :meth:`~matplotlib.cm.ColormapRegistry.get_cmap`.\nIf `None`, the value of `mpl.rcParams["image.cmap"]` is used.\nThe default `color_map` can be set using :func:`~scanpy.set_figure_params`.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Color Map',
    )
    cmap: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Colormap | str | None',
        title='Cmap',
    )
    palette: Any = Field(
        None,
        description='Colors to use for plotting categorical annotation groups.\nThe palette can be a valid :class:`~matplotlib.colors.ListedColormap` name\n(`\'Set2\'`, `\'tab20\'`, …), a :class:`~cycler.Cycler` object, a dict mapping\ncategories to colors, or a sequence of colors. Colors must be valid to\nmatplotlib. (see :func:`~matplotlib.colors.is_color_like`).\nIf `None`, `mpl.rcParams["axes.prop_cycle"]` is used unless the categorical\nvariable already has colors stored in `adata.uns["{var}_colors"]`.\nIf provided, values of `adata.uns["{var}_colors"]` will be set.\nOriginal type annotation: str | collections.abc.Sequence[str] | cycler.Cycler | None',
        title='Palette',
    )
    na_color: Optional[Any] = Field(
        'lightgray',
        description='Color to use for null or masked values. Can be anything matplotlib accepts as a\ncolor. Used for all points if `color=None`.\nOriginal type annotation: str | tuple[float, ...]',
        title='Na Color',
    )
    na_in_legend: Optional[Any] = Field(
        True,
        description="If there are missing values, whether they get an entry in the legend. Currently\nonly implemented for categorical legends.\nOriginal type annotation: <class 'bool'>",
        title='Na In Legend',
    )
    size: Any = Field(
        None,
        description='Point size. If `None`, is automatically computed as 120000 / n_cells.\nCan be a sequence containing the size for each cell. The order should be\nthe same as in adata.obs.\nOriginal type annotation: float | collections.abc.Sequence[float] | None',
        title='Size',
    )
    frameon: Any = Field(
        None,
        description='Draw a frame around the scatter plot. Defaults to value set in\n:func:`~scanpy.set_figure_params`, defaults to `True`.\nOriginal type annotation: bool | None',
        title='Frameon',
    )
    legend_fontsize: Any = Field(
        None,
        description="Numeric size in pt or string describing the size.\nSee :meth:`~matplotlib.text.Text.set_fontsize`.\nOriginal type annotation: typing.Union[float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], NoneType]",
        title='Legend Fontsize',
    )
    legend_fontweight: Optional[Any] = Field(
        'bold',
        description="Legend font weight. A numeric value in range 0-1000 or a string.\nDefaults to `'bold'` if `legend_loc == 'on data'`, otherwise to `'normal'`.\nSee :meth:`~matplotlib.text.Text.set_fontweight`.\nOriginal type annotation: typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]",
        title='Legend Fontweight',
    )
    legend_loc: Optional[Any] = Field(
        'right margin',
        description="Location of legend, either `'on data'`, `'right margin'`, `None`,\nor a valid keyword for the `loc` parameter of :class:`~matplotlib.legend.Legend`.\nOriginal type annotation: typing.Optional[typing.Literal['none', 'right margin', 'on data', 'on data export', 'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center']]",
        title='Legend Loc',
    )
    legend_fontoutline: Any = Field(
        None,
        description='Line width of the legend font outline in pt. Draws a white outline using\nthe path effect :class:`~matplotlib.patheffects.withStroke`.\nOriginal type annotation: int | None',
        title='Legend Fontoutline',
    )
    colorbar_loc: Optional[Any] = Field(
        'right',
        description='Where to place the colorbar for continous variables. If `None`, no colorbar\nis added.\nOriginal type annotation: str | None',
        title='Colorbar Loc',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. The format is the same as for `vmin`.\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None',
        title='Vmax',
    )
    vmin: Any = Field(
        None,
        description="The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin. vmin can be a number, a string, a function or `None`. If\nvmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).\nFor example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then\nvmin is interpreted as the return value of the function over the list of values to plot.\nFor example to set vmin tp the mean of the values to plot, `def my_vmin(values): return\nnp.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic\nminimum value is used as defined by matplotlib `scatter` function. When making multiple\nplots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Useful for diverging colormaps.\nThe format is the same as for `vmin`.\nExample: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``\nOriginal type annotation: str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float] | collections.abc.Sequence[str | float | collections.abc.Callable[[collections.abc.Sequence[float]], float]] | None",
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: matplotlib.colors.Normalize | collections.abc.Sequence[matplotlib.colors.Normalize] | None',
        title='Norm',
    )
    add_outline: Optional[Any] = Field(
        False,
        description='If set to True, this will add a thin border around groups of dots. In some situations\nthis can enhance the aesthetics of the resulting image\nOriginal type annotation: bool | None',
        title='Add Outline',
    )
    outline_width: Optional[Any] = Field(
        [0.3, 0.05],
        description='Tuple with two width numbers used to adjust the outline. The first value is the width\nof the border color as a fraction of the scatter dot size (default: 0.3). The second value is\nwidth of the gap color (default: 0.05).\nOriginal type annotation: tuple[float, float]',
        title='Outline Width',
    )
    outline_color: Optional[Any] = Field(
        ['black', 'white'],
        description='Tuple with two valid color names used to adjust the add_outline. The first color is the\nborder color (default: black), while the second color is a gap color between the\nborder color and the scatter dot (default: white).\nOriginal type annotation: tuple[str, str]',
        title='Outline Color',
    )
    ncols: Optional[Any] = Field(
        4,
        description="Number of panels per row.\nOriginal type annotation: <class 'int'>",
        title='Ncols',
    )
    hspace: Optional[Any] = Field(
        0.25,
        description="Adjust the height of the space between multiple panels.\nOriginal type annotation: <class 'float'>",
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.\nOriginal type annotation: float | None',
        title='Wspace',
    )
    title: Any = Field(
        None,
        description="Provide title for panels either as string or list of strings,\ne.g. `['title1', 'title2', ...]`.\nOriginal type annotation: str | collections.abc.Sequence[str] | None",
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: matplotlib.axes._axes.Axes | None',
        title='Ax',
    )
    return_fig: Any = Field(
        None,
        description='Return the matplotlib figure.\nOriginal type annotation: bool | None',
        title='Return Fig',
    )
    marker: Optional[Any] = Field(
        '.',
        description='No description available.\nOriginal type annotation: str | collections.abc.Sequence[str]',
        title='Marker',
    )
    annotate_var_explained: Optional[Any] = Field(
        False,
        description="No description available.\nOriginal type annotation: <class 'bool'>",
        title='Annotate Var Explained',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.pca')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingEmbeddingDensity(BaseAPI):
    """
    Plot the density of cells in an embedding (per condition).
    
    Plots the gaussian kernel density estimates (over condition) from the
    `sc.tl.embedding_density()` output.
    
    This function was written by Sophie Tritschler and implemented into
    Scanpy by Malte Luecken.
    
    Parameters
    ----------
    adata
        The annotated data matrix.
    basis
        The embedding over which the density was calculated. This embedded
        representation should be found in `adata.obsm['X_[basis]']``.
    key
        Name of the `.obs` covariate that contains the density estimates. Alternatively, pass `groupby`.
    groupby
        Name of the condition used in `tl.embedding_density`. Alternatively, pass `key`.
    group
        The category in the categorical observation annotation to be plotted.
        For example, 'G1' in the cell cycle 'phase' covariate. If all categories
        are to be plotted use group='all' (default), If multiple categories
        want to be plotted use a list (e.g.: ['G1', 'S']. If the overall density
        wants to be ploted set group to 'None'.
    color_map
        Matplolib color map to use for density plotting.
    bg_dotsize
        Dot size for background data points not in the `group`.
    fg_dotsize
        Dot size for foreground data points in the `group`.
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin. vmin can be a number, a string, a function or `None`. If
        vmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).
        For example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then
        vmin is interpreted as the return value of the function over the list of values to plot.
        For example to set vmin tp the mean of the values to plot, `def my_vmin(values): return
        np.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic
        minimum value is used as defined by matplotlib `scatter` function. When making multiple
        plots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`
    vmax
        The value representing the upper limit of the color scale. The format is the same as for `vmin`.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
        The format is the same as for `vmin`.
        Example: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``
    ncols
        Number of panels per row.
    wspace
        Adjust the width of the space between multiple panels.
    hspace
        Adjust the height of the space between multiple panels.
    return_fig
        Return the matplotlib figure.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.umap(adata)
        sc.tl.embedding_density(adata, basis='umap', groupby='phase')
    
    Plot all categories be default
    
    .. plot::
        :context: close-figs
    
        sc.pl.embedding_density(adata, basis='umap', key='umap_density_phase')
    
    Plot selected categories
    
    .. plot::
        :context: close-figs
    
        sc.pl.embedding_density(
            adata,
            basis='umap',
            key='umap_density_phase',
            group=['G1', 'S'],
        )
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    tl.embedding_density
    """

    
    adata: Any = Field(
        ...,
        description='The annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    basis: Optional[Any] = Field(
        'umap',
        description="The embedding over which the density was calculated. This embedded\nrepresentation should be found in `adata.obsm['X_[basis]']``.\nOriginal type annotation: str",
        title='Basis',
    )
    key: Any = Field(
        None,
        description='Name of the `.obs` covariate that contains the density estimates. Alternatively, pass `groupby`.\nOriginal type annotation: str | None',
        title='Key',
    )
    groupby: Any = Field(
        None,
        description='Name of the condition used in `tl.embedding_density`. Alternatively, pass `key`.\nOriginal type annotation: str | None',
        title='Groupby',
    )
    group: Optional[Any] = Field(
        'all',
        description="The category in the categorical observation annotation to be plotted.\nFor example, 'G1' in the cell cycle 'phase' covariate. If all categories\nare to be plotted use group='all' (default), If multiple categories\nwant to be plotted use a list (e.g.: ['G1', 'S']. If the overall density\nwants to be ploted set group to 'None'.\nOriginal type annotation: str | Sequence[str] | None",
        title='Group',
    )
    color_map: Optional[Any] = Field(
        'YlOrRd',
        description='Matplolib color map to use for density plotting.\nOriginal type annotation: Colormap | str',
        title='Color Map',
    )
    bg_dotsize: Optional[Any] = Field(
        80,
        description='Dot size for background data points not in the `group`.\nOriginal type annotation: int | None',
        title='Bg Dotsize',
    )
    fg_dotsize: Optional[Any] = Field(
        180,
        description='Dot size for foreground data points in the `group`.\nOriginal type annotation: int | None',
        title='Fg Dotsize',
    )
    vmax: Optional[Any] = Field(
        1,
        description='The value representing the upper limit of the color scale. The format is the same as for `vmin`.\nOriginal type annotation: int | None',
        title='Vmax',
    )
    vmin: Optional[Any] = Field(
        0,
        description="The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin. vmin can be a number, a string, a function or `None`. If\nvmin is a string and has the format `pN`, this is interpreted as a vmin=percentile(N).\nFor example vmin='p1.5' is interpreted as the 1.5 percentile. If vmin is function, then\nvmin is interpreted as the return value of the function over the list of values to plot.\nFor example to set vmin tp the mean of the values to plot, `def my_vmin(values): return\nnp.mean(values)` and then set `vmin=my_vmin`. If vmin is None (default) an automatic\nminimum value is used as defined by matplotlib `scatter` function. When making multiple\nplots, vmin can be a list of values, one for each plot. For example `vmin=[0.1, 'p1', None, my_vmin]`\nOriginal type annotation: int | None",
        title='Vmin',
    )
    vcenter: Any = Field(
        None,
        description="The value representing the center of the color scale. Useful for diverging colormaps.\nThe format is the same as for `vmin`.\nExample: ``sc.pl.umap(adata, color='TREM2', vcenter='p50', cmap='RdBu_r')``\nOriginal type annotation: int | None",
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: Normalize | None',
        title='Norm',
    )
    ncols: Optional[Any] = Field(
        4,
        description='Number of panels per row.\nOriginal type annotation: int | None',
        title='Ncols',
    )
    hspace: Optional[Any] = Field(
        0.25,
        description='Adjust the height of the space between multiple panels.\nOriginal type annotation: float | None',
        title='Hspace',
    )
    wspace: Any = Field(
        None,
        description='Adjust the width of the space between multiple panels.\nOriginal type annotation: None',
        title='Wspace',
    )
    title: Any = Field(
        None,
        description='No description available.\nOriginal type annotation: str | None',
        title='Title',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: Axes | None',
        title='Ax',
    )
    return_fig: Any = Field(
        None,
        description='Return the matplotlib figure.\nOriginal type annotation: bool | None',
        title='Return Fig',
    )

    _api_name = PrivateAttr(default='scanpy.plotting.embedding_density')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingRankGenesGroupsDotplot(BaseAPI):
    """
    Plot ranking of genes using dotplot plot (see :func:`~scanpy.pl.dotplot`).
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        The groups for which to show the gene ranking.
    n_genes
        Number of genes to show. This can be a negative number to show for
        example the down regulated genes. eg: num_genes=-10. Is ignored if
        `gene_names` is passed.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols. By default `var_names`
        refer to the index column of the `.var` DataFrame. Setting this option allows
        alternative names to be used.
    groupby
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank genes groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see :func:`~scanpy.pl.dotplot`).
    min_logfoldchange
        Value to filter genes in groups if their logfoldchange is less than the
        min_logfoldchange
    key
        Key used to store the ranking results in `adata.uns`.
    values_to_plot
        Instead of the mean gene value, plot the values computed by `sc.rank_genes_groups`.
        The options are: ['scores', 'logfoldchanges', 'pvals', 'pvals_adj',
        'log10_pvals', 'log10_pvals_adj']. When plotting logfoldchanges a divergent
        colormap is recommended. See examples below.
    var_names
        Genes to plot. Sometimes is useful to pass a specific list of var names (e.g. genes)
        to check their fold changes or p-values, instead of the top/bottom genes. The
        var_names could be a dictionary or a list as in :func:`~scanpy.pl.dotplot` or
        :func:`~scanpy.pl.matrixplot`. See examples below.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    return_fig
        Returns :class:`DotPlot` object. Useful for fine-tuning
        the plot. Takes precedence over `show=False`.
    **kwds
        Are passed to :func:`~scanpy.pl.dotplot`.
    
    Returns
    -------
    If `return_fig` is `True`, returns a :class:`DotPlot` object,
    else if `show` is false, return axes dict
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.rank_genes_groups(adata, 'bulk_labels', n_genes=adata.raw.shape[1])
    
    Plot top 2 genes per group.
    
    .. plot::
        :context: close-figs
    
        sc.pl.rank_genes_groups_dotplot(adata,n_genes=2)
    
    Plot with scaled expressions for easier identification of differences.
    
    .. plot::
        :context: close-figs
    
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=2, standard_scale='var')
    
    Plot `logfoldchanges` instead of gene expression. In this case a diverging colormap
    like `bwr` or `seismic` works better. To center the colormap in zero, the minimum
    and maximum values to plot are set to -4 and 4 respectively.
    Also, only genes with a log fold change of 3 or more are shown.
    
    .. plot::
        :context: close-figs
    
        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=4,
            values_to_plot="logfoldchanges", cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change'
        )
    
    Also, the last genes can be plotted. This can be useful to identify genes
    that are lowly expressed in a group. For this `n_genes=-4` is used
    
    .. plot::
        :context: close-figs
    
        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=-4,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
        )
    
    A list specific genes can be given to check their log fold change. If a
    dictionary, the dictionary keys will be added as labels in the plot.
    
    .. plot::
        :context: close-figs
    
        var_names = {'T-cell': ['CD3D', 'CD3E', 'IL32'],
                      'B-cell': ['CD79A', 'CD79B', 'MS4A1'],
                      'myeloid': ['CST3', 'LYZ'] }
        sc.pl.rank_genes_groups_dotplot(
            adata,
            var_names=var_names,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
        )
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    tl.rank_genes_groups
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    groups: Any = Field(
        None,
        description='The groups for which to show the gene ranking.\nOriginal type annotation: str | Sequence[str] | None',
        title='Groups',
    )
    n_genes: Any = Field(
        None,
        description='Number of genes to show. This can be a negative number to show for\nexample the down regulated genes. eg: num_genes=-10. Is ignored if\n`gene_names` is passed.\nOriginal type annotation: int | None',
        title='N Genes',
    )
    groupby: Any = Field(
        None,
        description='The key of the observation grouping to consider. By default,\nthe groupby is chosen from the rank genes groups parameter but\nother groupby options can be used.  It is expected that\ngroupby is a categorical. If groupby is not a categorical observation,\nit would be subdivided into `num_categories` (see :func:`~scanpy.pl.dotplot`).\nOriginal type annotation: str | None',
        title='Groupby',
    )
    values_to_plot: Any = Field(
        None,
        description="Instead of the mean gene value, plot the values computed by `sc.rank_genes_groups`.\nThe options are: ['scores', 'logfoldchanges', 'pvals', 'pvals_adj',\n'log10_pvals', 'log10_pvals_adj']. When plotting logfoldchanges a divergent\ncolormap is recommended. See examples below.\nOriginal type annotation: Literal['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'] | None",
        title='Values To Plot',
    )
    var_names: Any = Field(
        None,
        description='Genes to plot. Sometimes is useful to pass a specific list of var names (e.g. genes)\nto check their fold changes or p-values, instead of the top/bottom genes. The\nvar_names could be a dictionary or a list as in :func:`~scanpy.pl.dotplot` or\n:func:`~scanpy.pl.matrixplot`. See examples below.\nOriginal type annotation: Sequence[str] | Mapping[str, Sequence[str]] | None',
        title='Var Names',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols. By default `var_names`\nrefer to the index column of the `.var` DataFrame. Setting this option allows\nalternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    min_logfoldchange: Any = Field(
        None,
        description='Value to filter genes in groups if their logfoldchange is less than the\nmin_logfoldchange\nOriginal type annotation: float | None',
        title='Min Logfoldchange',
    )
    key: Any = Field(
        None,
        description='Key used to store the ranking results in `adata.uns`.\nOriginal type annotation: str | None',
        title='Key',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | None",
        title='Save',
    )
    return_fig: Optional[Any] = Field(
        False,
        description='Returns :class:`DotPlot` object. Useful for fine-tuning\nthe plot. Takes precedence over `show=False`.\nOriginal type annotation: bool',
        title='Return Fig',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.rank_genes_groups_dotplot')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingCorrelationMatrix(BaseAPI):
    """
    Plot the correlation matrix computed as part of `sc.tl.dendrogram`.
    
    Parameters
    ----------
    adata
    groupby
        Categorical data column used to create the dendrogram
    show_correlation_numbers
        If `show_correlation=True`, plot the correlation on top of each cell.
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the
        hierarchical clustering between the `groupby` categories is added.
        The dendrogram is computed using :func:`scanpy.tl.dendrogram`.
        If `tl.dendrogram` has not been called previously,
        the function is called with default parameters.
    figsize
        By default a figure size that aims to produce a squared correlation
        matrix plot is used. Format is (width, height)
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin.
    vmax
        The value representing the upper limit of the color scale. Values larger than vmax are plotted
        with the same color as vmax.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
    norm
        Custom color normalization object from matplotlib. See
        `https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.
    **kwds
        Only if `show_correlation` is True:
        Are passed to :func:`matplotlib.pyplot.pcolormesh` when plotting the
        correlation heatmap. `cmap` can be used to change the color palette.
    
    Returns
    -------
    If `show=False`, returns a list of :class:`matplotlib.axes.Axes` objects.
    
    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, "bulk_labels")
    >>> sc.pl.correlation_matrix(adata, "bulk_labels")
    """

    
    adata: Any = Field(
        ...,
        description='No description available.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    groupby: Any = Field(
        ...,
        description='Categorical data column used to create the dendrogram\nOriginal type annotation: str',
        title='Groupby',
    )
    show_correlation_numbers: Optional[Any] = Field(
        False,
        description='If `show_correlation=True`, plot the correlation on top of each cell.\nOriginal type annotation: bool',
        title='Show Correlation Numbers',
    )
    dendrogram: Any = Field(
        None,
        description='If True or a valid dendrogram key, a dendrogram based on the\nhierarchical clustering between the `groupby` categories is added.\nThe dendrogram is computed using :func:`scanpy.tl.dendrogram`.\nIf `tl.dendrogram` has not been called previously,\nthe function is called with default parameters.\nOriginal type annotation: bool | str | None',
        title='Dendrogram',
    )
    figsize: Any = Field(
        None,
        description='By default a figure size that aims to produce a squared correlation\nmatrix plot is used. Format is (width, height)\nOriginal type annotation: tuple[float, float] | None',
        title='Figsize',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: Axes | None',
        title='Ax',
    )
    vmin: Any = Field(
        None,
        description='The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin.\nOriginal type annotation: float | None',
        title='Vmin',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. Values larger than vmax are plotted\nwith the same color as vmax.\nOriginal type annotation: float | None',
        title='Vmax',
    )
    vcenter: Any = Field(
        None,
        description='The value representing the center of the color scale. Useful for diverging colormaps.\nOriginal type annotation: float | None',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib. See\n`https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.\nOriginal type annotation: Normalize | None',
        title='Norm',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.correlation_matrix')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingHighestExprGenes(BaseAPI):
    """
    Fraction of counts assigned to each gene over all cells.
    
    Computes, for each gene, the fraction of counts assigned to that gene within
    a cell. The `n_top` genes with the highest mean fraction over all cells are
    plotted as boxplots.
    
    This plot is similar to the `scater` package function `plotHighestExprs(type
    = "highest-expression")`, see `here
    <https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette-qc.html>`__. Quoting
    
        *We expect to see the “usual suspects”, i.e., mitochondrial genes, actin,
        ribosomal protein, MALAT1. A few spike-in transcripts may also be
        present here, though if all of the spike-ins are in the top 50, it
        suggests that too much spike-in RNA was added. A large number of
        pseudo-genes or predicted genes may indicate problems with alignment.*
        -- Davis McCarthy and Aaron Lun
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_top
        Number of top
    layer
        Layer from which to pull data.
    gene_symbols
        Key for field in .var that stores gene symbols if you do not want to use .var_names.
    log
        Plot x-axis in log scale
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    **kwds
        Are passed to :func:`~seaborn.boxplot`.
    
    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes`.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    n_top: Optional[Any] = Field(
        30, description='Number of top\nOriginal type annotation: int', title='N Top'
    )
    layer: Any = Field(
        None,
        description='Layer from which to pull data.\nOriginal type annotation: str | None',
        title='Layer',
    )
    gene_symbols: Any = Field(
        None,
        description='Key for field in .var that stores gene symbols if you do not want to use .var_names.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot x-axis in log scale\nOriginal type annotation: bool',
        title='Log',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: str | bool | None",
        title='Save',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: Axes | None',
        title='Ax',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.highest_expr_genes')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingTracksplot(BaseAPI):
    """
    Compact plot of expression of a list of genes.
    
    In this type of plot each var_name is plotted as a filled line plot where the
    y values correspond to the var_name values and x is each of the cells. Best results
    are obtained when using raw counts that are not log.
    
    `groupby` is required to sort and order the values using the respective group
    and should be a categorical value.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    var_names
        `var_names` should be a valid subset of `adata.var_names`.
        If `var_names` is a mapping, then the key is used as label
        to group the values (see `var_group_labels`). The mapping values
        should be sequences of valid `adata.var_names`. In this
        case either coloring or 'brackets' are used for the grouping
        of var names depending on the plot. When `var_names` is a mapping,
        then the `var_group_labels` and `var_group_positions` are set.
    groupby
        The key of the observation grouping to consider.
    use_raw
        Use `raw` attribute of `adata` if present.
    log
        Plot on logarithmic axis.
    num_categories
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    categories_order
        Order in which to show the categories. Note: add_dendrogram or add_totals
        can change the categories order.
    figsize
        Figure size when `multi_panel=True`.
        Otherwise the `rcParam['figure.figsize]` value is used.
        Format is (width, height)
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the hierarchical
        clustering between the `groupby` categories is added.
        The dendrogram information is computed using :func:`scanpy.tl.dendrogram`.
        If `tl.dendrogram` has not been called previously the function is called
        with default parameters.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols.
        By default `var_names` refer to the index column of the `.var` DataFrame.
        Setting this option allows alternative names to be used.
    var_group_positions
        Use this parameter to highlight groups of `var_names`.
        This will draw a 'bracket' or a color block between the given start and end
        positions. If the parameter `var_group_labels` is set, the corresponding
        labels are added on top/left. E.g. `var_group_positions=[(4,10)]`
        will add a bracket between the fourth `var_name` and the tenth `var_name`.
        By giving more positions, more brackets/color blocks are drawn.
    var_group_labels
        Labels for each of the `var_group_positions` that want to be highlighted.
    var_group_rotation
        Label rotation degrees.
        By default, labels larger than 4 characters are rotated 90 degrees.
    layer
        Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.
        If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,
        then the layer is plotted. `layer` takes precedence over `use_raw`.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    **kwds
        Are passed to :func:`~seaborn.heatmap`.
    
    Returns
    -------
    A list of :class:`~matplotlib.axes.Axes`.
    
    Examples
    --------
    Using var_names as list:
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.tracksplot(adata, markers, groupby='bulk_labels', dendrogram=True)
    
    Using var_names as dict:
    
    .. plot::
        :context: close-figs
    
        markers = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
        sc.pl.tracksplot(adata, markers, groupby='bulk_labels', dendrogram=True)
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    pl.rank_genes_groups_tracksplot: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`.\nIf `var_names` is a mapping, then the key is used as label\nto group the values (see `var_group_labels`). The mapping values\nshould be sequences of valid `adata.var_names`. In this\ncase either coloring or 'brackets' are used for the grouping\nof var names depending on the plot. When `var_names` is a mapping,\nthen the `var_group_labels` and `var_group_positions` are set.\nOriginal type annotation: _VarNames | Mapping[str, _VarNames]",
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='The key of the observation grouping to consider.\nOriginal type annotation: str',
        title='Groupby',
    )
    use_raw: Any = Field(
        None,
        description='Use `raw` attribute of `adata` if present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axis.\nOriginal type annotation: bool',
        title='Log',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='If True or a valid dendrogram key, a dendrogram based on the hierarchical\nclustering between the `groupby` categories is added.\nThe dendrogram information is computed using :func:`scanpy.tl.dendrogram`.\nIf `tl.dendrogram` has not been called previously the function is called\nwith default parameters.\nOriginal type annotation: bool | str',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.\nBy default `var_names` refer to the index column of the `.var` DataFrame.\nSetting this option allows alternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description="Use this parameter to highlight groups of `var_names`.\nThis will draw a 'bracket' or a color block between the given start and end\npositions. If the parameter `var_group_labels` is set, the corresponding\nlabels are added on top/left. E.g. `var_group_positions=[(4,10)]`\nwill add a bracket between the fourth `var_name` and the tenth `var_name`.\nBy giving more positions, more brackets/color blocks are drawn.\nOriginal type annotation: Sequence[tuple[int, int]] | None",
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for each of the `var_group_positions` that want to be highlighted.\nOriginal type annotation: Sequence[str] | None',
        title='Var Group Labels',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.\nIf `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,\nthen the layer is plotted. `layer` takes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: str | bool | None",
        title='Save',
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`.\nOtherwise the `rcParam['figure.figsize]` value is used.\nFormat is (width, height)\nOriginal type annotation: tuple[float, float] | None",
        title='Figsize',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.tracksplot')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingClustermap(BaseAPI):
    """
    Hierarchically-clustered heatmap.
    
    Wraps :func:`seaborn.clustermap` for :class:`~anndata.AnnData`.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    obs_keys
        Categorical annotation to plot with a different color map.
        Currently, only a single key is supported.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    **kwds
        Keyword arguments passed to :func:`~seaborn.clustermap`.
    
    Returns
    -------
    If `show` is `False`, a :class:`~seaborn.matrix.ClusterGrid` object
    (see :func:`~seaborn.clustermap`).
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.krumsiek11()
        sc.pl.clustermap(adata)
    
    .. plot::
        :context: close-figs
    
        sc.pl.clustermap(adata, obs_keys='cell_type')
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    obs_keys: Any = Field(
        None,
        description='Categorical annotation to plot with a different color map.\nCurrently, only a single key is supported.\nOriginal type annotation: str | None',
        title='Obs Keys',
    )
    use_raw: Any = Field(
        None,
        description='Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    kwds: Any = Field(..., description='No description available.', title='Kwds')

    _api_name = PrivateAttr(default='scanpy.plotting.clustermap')
    _products_original = PrivateAttr(default=[])
    _data_name = PrivateAttr(default='adata')



class ScanpyPlottingStackedViolin(BaseAPI):
    """
    Stacked violin plots.
    
    Makes a compact image composed of individual violin plots
    (from :func:`~seaborn.violinplot`) stacked on top of each other.
    Useful to visualize gene expression per cluster.
    
    Wraps :func:`seaborn.violinplot` for :class:`~anndata.AnnData`.
    
    This function provides a convenient interface to the
    :class:`~scanpy.pl.StackedViolin` class. If you need more flexibility,
    you should use :class:`~scanpy.pl.StackedViolin` directly.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    var_names
        `var_names` should be a valid subset of `adata.var_names`.
        If `var_names` is a mapping, then the key is used as label
        to group the values (see `var_group_labels`). The mapping values
        should be sequences of valid `adata.var_names`. In this
        case either coloring or 'brackets' are used for the grouping
        of var names depending on the plot. When `var_names` is a mapping,
        then the `var_group_labels` and `var_group_positions` are set.
    groupby
        The key of the observation grouping to consider.
    use_raw
        Use `raw` attribute of `adata` if present.
    log
        Plot on logarithmic axis.
    num_categories
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    categories_order
        Order in which to show the categories. Note: add_dendrogram or add_totals
        can change the categories order.
    figsize
        Figure size when `multi_panel=True`.
        Otherwise the `rcParam['figure.figsize]` value is used.
        Format is (width, height)
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the hierarchical
        clustering between the `groupby` categories is added.
        The dendrogram information is computed using :func:`scanpy.tl.dendrogram`.
        If `tl.dendrogram` has not been called previously the function is called
        with default parameters.
    gene_symbols
        Column name in `.var` DataFrame that stores gene symbols.
        By default `var_names` refer to the index column of the `.var` DataFrame.
        Setting this option allows alternative names to be used.
    var_group_positions
        Use this parameter to highlight groups of `var_names`.
        This will draw a 'bracket' or a color block between the given start and end
        positions. If the parameter `var_group_labels` is set, the corresponding
        labels are added on top/left. E.g. `var_group_positions=[(4,10)]`
        will add a bracket between the fourth `var_name` and the tenth `var_name`.
        By giving more positions, more brackets/color blocks are drawn.
    var_group_labels
        Labels for each of the `var_group_positions` that want to be highlighted.
    var_group_rotation
        Label rotation degrees.
        By default, labels larger than 4 characters are rotated 90 degrees.
    layer
        Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.
        If `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,
        then the layer is plotted. `layer` takes precedence over `use_raw`.
    title
        Title for the figure
    colorbar_title
        Title for the color bar. New line character (\n) can be used.
    cmap
        String denoting matplotlib color map.
    standard_scale
        Whether or not to standardize the given dimension between 0 and 1, meaning for
        each variable or group, subtract the minimum and divide each by its maximum.
    swap_axes
         By default, the x axis contains `var_names` (e.g. genes) and the y axis
         the `groupby` categories. By setting `swap_axes` then x are the
         `groupby` categories and y the `var_names`.
    return_fig
        Returns :class:`DotPlot` object. Useful for fine-tuning
        the plot. Takes precedence over `show=False`.
    
    stripplot
        Add a stripplot on top of the violin plot.
        See :func:`~seaborn.stripplot`.
    jitter
        Add jitter to the stripplot (only when stripplot is True)
        See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    density_norm
        The method used to scale the width of each violin.
        If 'width' (the default), each violin will have the same width.
        If 'area', each violin will have the same area.
        If 'count', a violin’s width corresponds to the number of observations.
    yticklabels
        Set to true to view the y tick labels.
    row_palette
        Be default, median values are mapped to the violin color using a
        color map (see `cmap` argument). Alternatively, a 'row_palette` can
        be given to color each violin plot row using a different colors.
        The value should be a valid seaborn or matplotlib palette name
        (see :func:`~seaborn.color_palette`).
        Alternatively, a single color name or hex value can be passed,
        e.g. `'red'` or `'#cc33ff'`.
    show
         Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    ax
        A matplotlib axes object. Only works if plotting a single component.
    vmin
        The value representing the lower limit of the color scale. Values smaller than vmin are plotted
        with the same color as vmin.
    vmax
        The value representing the upper limit of the color scale. Values larger than vmax are plotted
        with the same color as vmax.
    vcenter
        The value representing the center of the color scale. Useful for diverging colormaps.
    norm
        Custom color normalization object from matplotlib. See
        `https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.
    **kwds
        Are passed to :func:`~seaborn.violinplot`.
    
    Returns
    -------
    If `return_fig` is `True`, returns a :class:`~scanpy.pl.StackedViolin` object,
    else if `show` is false, return axes dict
    
    See Also
    --------
    :class:`~scanpy.pl.StackedViolin`: The StackedViolin class can be used to to control
        several visual parameters not available in this function.
    :func:`~scanpy.pl.rank_genes_groups_stacked_violin` to plot marker genes identified
        using the :func:`~scanpy.tl.rank_genes_groups` function.
    
    Examples
    --------
    Visualization of violin plots of a few genes grouped by the category `bulk_labels`:
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)
    
    Same visualization but passing var_names as dict, which adds a grouping of
    the genes on top of the image:
    
    .. plot::
        :context: close-figs
    
        markers = {'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}
        sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)
    
    Get StackedViolin object for fine tuning
    
    .. plot::
        :context: close-figs
    
        vp = sc.pl.stacked_violin(adata, markers, 'bulk_labels', return_fig=True)
        vp.add_totals().style(ylim=(0,5)).show()
    
    The axes used can be obtained using the get_axes() method:
    
    .. code-block:: python
    
        axes_dict = vp.get_axes()
        print(axes_dict)
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    var_names: Any = Field(
        ...,
        description="`var_names` should be a valid subset of `adata.var_names`.\nIf `var_names` is a mapping, then the key is used as label\nto group the values (see `var_group_labels`). The mapping values\nshould be sequences of valid `adata.var_names`. In this\ncase either coloring or 'brackets' are used for the grouping\nof var names depending on the plot. When `var_names` is a mapping,\nthen the `var_group_labels` and `var_group_positions` are set.\nOriginal type annotation: _VarNames | Mapping[str, _VarNames]",
        title='Var Names',
    )
    groupby: Any = Field(
        ...,
        description='The key of the observation grouping to consider.\nOriginal type annotation: str | Sequence[str]',
        title='Groupby',
    )
    log: Optional[Any] = Field(
        False,
        description='Plot on logarithmic axis.\nOriginal type annotation: bool',
        title='Log',
    )
    use_raw: Any = Field(
        None,
        description='Use `raw` attribute of `adata` if present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    num_categories: Optional[Any] = Field(
        7,
        description='Only used if groupby observation is not categorical. This value\ndetermines the number of groups into which the groupby observation\nshould be subdivided.\nOriginal type annotation: int',
        title='Num Categories',
    )
    title: Any = Field(
        None,
        description='Title for the figure\nOriginal type annotation: str | None',
        title='Title',
    )
    colorbar_title: Optional[Any] = Field(
        'Median expression\nin group',
        description='Title for the color bar. New line character (\\n) can be used.\nOriginal type annotation: str | None',
        title='Colorbar Title',
    )
    figsize: Any = Field(
        None,
        description="Figure size when `multi_panel=True`.\nOtherwise the `rcParam['figure.figsize]` value is used.\nFormat is (width, height)\nOriginal type annotation: tuple[float, float] | None",
        title='Figsize',
    )
    dendrogram: Optional[Any] = Field(
        False,
        description='If True or a valid dendrogram key, a dendrogram based on the hierarchical\nclustering between the `groupby` categories is added.\nThe dendrogram information is computed using :func:`scanpy.tl.dendrogram`.\nIf `tl.dendrogram` has not been called previously the function is called\nwith default parameters.\nOriginal type annotation: bool | str',
        title='Dendrogram',
    )
    gene_symbols: Any = Field(
        None,
        description='Column name in `.var` DataFrame that stores gene symbols.\nBy default `var_names` refer to the index column of the `.var` DataFrame.\nSetting this option allows alternative names to be used.\nOriginal type annotation: str | None',
        title='Gene Symbols',
    )
    var_group_positions: Any = Field(
        None,
        description="Use this parameter to highlight groups of `var_names`.\nThis will draw a 'bracket' or a color block between the given start and end\npositions. If the parameter `var_group_labels` is set, the corresponding\nlabels are added on top/left. E.g. `var_group_positions=[(4,10)]`\nwill add a bracket between the fourth `var_name` and the tenth `var_name`.\nBy giving more positions, more brackets/color blocks are drawn.\nOriginal type annotation: Sequence[tuple[int, int]] | None",
        title='Var Group Positions',
    )
    var_group_labels: Any = Field(
        None,
        description='Labels for each of the `var_group_positions` that want to be highlighted.\nOriginal type annotation: Sequence[str] | None',
        title='Var Group Labels',
    )
    standard_scale: Any = Field(
        None,
        description="Whether or not to standardize the given dimension between 0 and 1, meaning for\neach variable or group, subtract the minimum and divide each by its maximum.\nOriginal type annotation: Literal['var', 'group'] | None",
        title='Standard Scale',
    )
    var_group_rotation: Any = Field(
        None,
        description='Label rotation degrees.\nBy default, labels larger than 4 characters are rotated 90 degrees.\nOriginal type annotation: float | None',
        title='Var Group Rotation',
    )
    layer: Any = Field(
        None,
        description='Name of the AnnData object layer that wants to be plotted. By default adata.raw.X is plotted.\nIf `use_raw=False` is set, then `adata.X` is plotted. If `layer` is set to a valid layer name,\nthen the layer is plotted. `layer` takes precedence over `use_raw`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    categories_order: Any = Field(
        None,
        description='Order in which to show the categories. Note: add_dendrogram or add_totals\ncan change the categories order.\nOriginal type annotation: Sequence[str] | None',
        title='Categories Order',
    )
    swap_axes: Optional[Any] = Field(
        False,
        description='By default, the x axis contains `var_names` (e.g. genes) and the y axis\nthe `groupby` categories. By setting `swap_axes` then x are the\n`groupby` categories and y the `var_names`.\nOriginal type annotation: bool',
        title='Swap Axes',
    )
    show: Any = Field(
        None,
        description='Show the plot, do not return axis.\nOriginal type annotation: bool | None',
        title='Show',
    )
    save: Any = Field(
        None,
        description="If `True` or a `str`, save the figure.\nA string is appended to the default filename.\nInfer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.\nOriginal type annotation: bool | str | None",
        title='Save',
    )
    return_fig: Optional[Any] = Field(
        False,
        description='Returns :class:`DotPlot` object. Useful for fine-tuning\nthe plot. Takes precedence over `show=False`.\nOriginal type annotation: bool | None',
        title='Return Fig',
    )
    ax: Any = Field(
        None,
        description='A matplotlib axes object. Only works if plotting a single component.\nOriginal type annotation: _AxesSubplot | None',
        title='Ax',
    )
    vmin: Any = Field(
        None,
        description='The value representing the lower limit of the color scale. Values smaller than vmin are plotted\nwith the same color as vmin.\nOriginal type annotation: float | None',
        title='Vmin',
    )
    vmax: Any = Field(
        None,
        description='The value representing the upper limit of the color scale. Values larger than vmax are plotted\nwith the same color as vmax.\nOriginal type annotation: float | None',
        title='Vmax',
    )
    vcenter: Any = Field(
        None,
        description='The value representing the center of the color scale. Useful for diverging colormaps.\nOriginal type annotation: float | None',
        title='Vcenter',
    )
    norm: Any = Field(
        None,
        description='Custom color normalization object from matplotlib. See\n`https://matplotlib.org/stable/tutorials/colors/colormapnorms.html` for details.\nOriginal type annotation: Normalize | None',
        title='Norm',
    )
    cmap: Optional[Any] = Field(
        'Blues',
        description='String denoting matplotlib color map.\nOriginal type annotation: Colormap | str | None',
        title='Cmap',
    )
    stripplot: Optional[Any] = Field(
        False,
        description='Add a stripplot on top of the violin plot.\nSee :func:`~seaborn.stripplot`.\nOriginal type annotation: bool',
        title='Stripplot',
    )
    jitter: Optional[Any] = Field(
        False,
        description='Add jitter to the stripplot (only when stripplot is True)\nSee :func:`~seaborn.stripplot`.\nOriginal type annotation: float | bool',
        title='Jitter',
    )
    size: Optional[Any] = Field(
        1,
        description='Size of the jitter points.\nOriginal type annotation: float',
        title='Size',
    )
    row_palette: Any = Field(
        None,
        description="Be default, median values are mapped to the violin color using a\ncolor map (see `cmap` argument). Alternatively, a 'row_palette` can\nbe given to color each violin plot row using a different colors.\nThe value should be a valid seaborn or matplotlib palette name\n(see :func:`~seaborn.color_palette`).\nAlternatively, a single color name or hex value can be passed,\ne.g. `'red'` or `'#cc33ff'`.\nOriginal type annotation: str | None",
        title='Row Palette',
    )
    density_norm: Optional[Any] = Field(
        0,
        description="The method used to scale the width of each violin.\nIf 'width' (the default), each violin will have the same width.\nIf 'area', each violin will have the same area.\nIf 'count', a violin’s width corresponds to the number of observations.\nOriginal type annotation: DensityNorm | Empty",
        title='Density Norm',
    )
    yticklabels: Optional[Any] = Field(
        False,
        description='Set to true to view the y tick labels.\nOriginal type annotation: bool',
        title='Yticklabels',
    )
    order: Optional[Any] = Field(
        0,
        description='No description available.\nOriginal type annotation: Sequence[str] | None | Empty',
        title='Order',
    )
    scale: Optional[Any] = Field(
        0,
        description='No description available.\nOriginal type annotation: DensityNorm | Empty',
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