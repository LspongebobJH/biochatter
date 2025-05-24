from __future__ import annotations

from typing import Any, Optional

from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScanpyToolsPaga(BaseAPI):
    """
    Map out the coarse-grained connectivity structures of complex manifolds :cite:p:`Wolf2019`.
    
    By quantifying the connectivity of partitions (groups, clusters) of the
    single-cell graph, partition-based graph abstraction (PAGA) generates a much
    simpler abstracted graph (*PAGA graph*) of partitions, in which edge weights
    represent confidence in the presence of connections. By thresholding this
    confidence in :func:`~scanpy.pl.paga`, a much simpler representation of the
    manifold data is obtained, which is nonetheless faithful to the topology of
    the manifold.
    
    The confidence should be interpreted as the ratio of the actual versus the
    expected value of connections under the null model of randomly connecting
    partitions. We do not provide a p-value as this null model does not
    precisely capture what one would consider "connected" in real data, hence it
    strongly overestimates the expected value. See an extensive discussion of
    this in :cite:t:`Wolf2019`.
    
    .. note::
        Note that you can use the result of :func:`~scanpy.pl.paga` in
        :func:`~scanpy.tl.umap` and :func:`~scanpy.tl.draw_graph` via
        `init_pos='paga'` to get single-cell embeddings that are typically more
        faithful to the global topology.
    
    Parameters
    ----------
    adata
        An annotated data matrix.
    groups
        Key for categorical in `adata.obs`. You can pass your predefined groups
        by choosing any categorical annotation of observations. Default:
        The first present key of `'leiden'` or `'louvain'`.
    use_rna_velocity
        Use RNA velocity to orient edges in the abstracted graph and estimate
        transitions. Requires that `adata.uns` contains a directed single-cell
        graph with key `['velocity_graph']`. This feature might be subject
        to change in the future.
    model
        The PAGA connectivity model.
    neighbors_key
        If not specified, paga looks `.uns['neighbors']` for neighbors settings
        and `.obsp['connectivities']`, `.obsp['distances']` for connectivities and
        distances respectively (default storage places for `pp.neighbors`).
        If specified, paga looks `.uns[neighbors_key]` for neighbors settings and
        `.obsp[.uns[neighbors_key]['connectivities_key']]`,
        `.obsp[.uns[neighbors_key]['distances_key']]` for connectivities and distances
        respectively.
    copy
        Copy `adata` before computation and return a copy. Otherwise, perform
        computation inplace and return `None`.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.uns['connectivities']` : :class:`numpy.ndarray` (dtype `float`)
        The full adjacency matrix of the abstracted graph, weights correspond to
        confidence in the connectivities of partitions.
    `adata.uns['connectivities_tree']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        The adjacency matrix of the tree-like subgraph that best explains
        the topology.
    
    Notes
    -----
    Together with a random walk-based distance measure
    (e.g. :func:`scanpy.tl.dpt`) this generates a partial coordinatization of
    data useful for exploring and explaining its variation.
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    pl.paga
    pl.paga_path
    pl.paga_compare
    """

    
    adata: Any = Field(
        ...,
        description='An annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    groups: Any = Field(
        None,
        description="Key for categorical in `adata.obs`. You can pass your predefined groups\nby choosing any categorical annotation of observations. Default:\nThe first present key of `'leiden'` or `'louvain'`.\nOriginal type annotation: str | None",
        title='Groups',
    )
    use_rna_velocity: Optional[Any] = Field(
        False,
        description="Use RNA velocity to orient edges in the abstracted graph and estimate\ntransitions. Requires that `adata.uns` contains a directed single-cell\ngraph with key `['velocity_graph']`. This feature might be subject\nto change in the future.\nOriginal type annotation: bool",
        title='Use Rna Velocity',
    )
    model: Optional[Any] = Field(
        'v1.2',
        description="The PAGA connectivity model.\nOriginal type annotation: Literal['v1.2', 'v1.0']",
        title='Model',
    )
    neighbors_key: Any = Field(
        None,
        description="If not specified, paga looks `.uns['neighbors']` for neighbors settings\nand `.obsp['connectivities']`, `.obsp['distances']` for connectivities and\ndistances respectively (default storage places for `pp.neighbors`).\nIf specified, paga looks `.uns[neighbors_key]` for neighbors settings and\n`.obsp[.uns[neighbors_key]['connectivities_key']]`,\n`.obsp[.uns[neighbors_key]['distances_key']]` for connectivities and distances\nrespectively.\nOriginal type annotation: str | None",
        title='Neighbors Key',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Copy `adata` before computation and return a copy. Otherwise, perform\ncomputation inplace and return `None`.\nOriginal type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.tools.paga')
    _products_original = PrivateAttr(default=['data.uns["paga"]["connectivities"]', 'data.uns["paga"]["connectivities_tree"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyToolsLeiden(BaseAPI):
    """
    Cluster cells into subgroups :cite:p:`Traag2019`.
    
    Cluster cells using the Leiden algorithm :cite:p:`Traag2019`,
    an improved version of the Louvain algorithm :cite:p:`Blondel2008`.
    It was proposed for single-cell analysis by :cite:t:`Levine2015`.
    
    This requires having run :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.
    
    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters.
        Set to `None` if overriding `partition_type`
        to one that doesn’t accept a `resolution_parameter`.
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain `(obs_key, list_of_categories)`.
    key_added
        `adata.obs` key under which to add the cluster labels.
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    n_iterations
        How many iterations of the Leiden clustering algorithm to perform.
        Positive values above 2 define the total number of iterations to perform,
        -1 has the algorithm run until it reaches its optimal clustering.
        2 is faster and the default for underlying packages.
    partition_type
        Type of partition to use.
        Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
        For the available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    neighbors_key
        Use neighbors connectivities as adjacency.
        If not specified, leiden looks at .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, leiden looks at
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    obsp
        Use .obsp[obsp] as adjacency. You can't specify both
        `obsp` and `neighbors_key` at the same time.
    copy
        Whether to copy `adata` or modify it inplace.
    flavor
        Which package's implementation to use.
    **clustering_args
        Any further arguments to pass to :func:`~leidenalg.find_partition` (which in turn passes arguments to the `partition_type`)
        or :meth:`igraph.Graph.community_leiden` from `igraph`.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obs['leiden' | key_added]` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.
    
    `adata.uns['leiden' | key_added]['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    """

    
    adata: Any = Field(
        ...,
        description='The annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    resolution: Optional[Any] = Field(
        1,
        description='A parameter value controlling the coarseness of the clustering.\nHigher values lead to more clusters.\nSet to `None` if overriding `partition_type`\nto one that doesn’t accept a `resolution_parameter`.\nOriginal type annotation: float',
        title='Resolution',
    )
    restrict_to: Any = Field(
        None,
        description='Restrict the clustering to the categories within the key for sample\nannotation, tuple needs to contain `(obs_key, list_of_categories)`.\nOriginal type annotation: tuple[str, Sequence[str]] | None',
        title='Restrict To',
    )
    random_state: Optional[Any] = Field(
        0,
        description='Change the initialization of the optimization.\nOriginal type annotation: _LegacyRandom',
        title='Random State',
    )
    key_added: Optional[Any] = Field(
        'leiden',
        description='`adata.obs` key under which to add the cluster labels.\nOriginal type annotation: str',
        title='Key Added',
    )
    adjacency: Any = Field(
        None,
        description='Sparse adjacency matrix of the graph, defaults to neighbors connectivities.\nOriginal type annotation: _CSMatrix | None',
        title='Adjacency',
    )
    directed: Any = Field(
        None,
        description='Whether to treat the graph as directed or undirected.\nOriginal type annotation: bool | None',
        title='Directed',
    )
    use_weights: Optional[Any] = Field(
        True,
        description='If `True`, edge weights from the graph are used in the computation\n(placing more emphasis on stronger edges).\nOriginal type annotation: bool',
        title='Use Weights',
    )
    n_iterations: Optional[Any] = Field(
        -1,
        description='How many iterations of the Leiden clustering algorithm to perform.\nPositive values above 2 define the total number of iterations to perform,\n-1 has the algorithm run until it reaches its optimal clustering.\n2 is faster and the default for underlying packages.\nOriginal type annotation: int',
        title='N Iterations',
    )
    partition_type: Any = Field(
        None,
        description='Type of partition to use.\nDefaults to :class:`~leidenalg.RBConfigurationVertexPartition`.\nFor the available options, consult the documentation for\n:func:`~leidenalg.find_partition`.\nOriginal type annotation: type[MutableVertexPartition] | None',
        title='Partition Type',
    )
    neighbors_key: Any = Field(
        None,
        description="Use neighbors connectivities as adjacency.\nIf not specified, leiden looks at .obsp['connectivities'] for connectivities\n(default storage place for pp.neighbors).\nIf specified, leiden looks at\n.obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.\nOriginal type annotation: str | None",
        title='Neighbors Key',
    )
    obsp: Any = Field(
        None,
        description="Use .obsp[obsp] as adjacency. You can't specify both\n`obsp` and `neighbors_key` at the same time.\nOriginal type annotation: str | None",
        title='Obsp',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Whether to copy `adata` or modify it inplace.\nOriginal type annotation: bool',
        title='Copy',
    )
    flavor: Optional[Any] = Field(
        'leidenalg',
        description="Which package's implementation to use.\nOriginal type annotation: Literal['leidenalg', 'igraph']",
        title='Flavor',
    )
    clustering_args: Any = Field(
        ..., description='No description available.', title='Clustering Args'
    )

    _api_name = PrivateAttr(default='scanpy.tools.leiden')
    _products_original = PrivateAttr(default=['data.obs["leiden"]', 'data.uns["leiden"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyToolsLouvain(BaseAPI):
    """
    Cluster cells into subgroups :cite:p:`Blondel2008,Levine2015,Traag2017`.
    
    Cluster cells using the Louvain algorithm :cite:p:`Blondel2008` in the implementation
    of :cite:t:`Traag2017`. The Louvain algorithm was proposed for single-cell
    analysis by :cite:t:`Levine2015`.
    
    This requires having run :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first,
    or explicitly passing a ``adjacency`` matrix.
    
    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        For the default flavor (``'vtraag'``) or for ```RAPIDS```, you can provide a
        resolution (higher resolution means finding more and smaller clusters),
        which defaults to 1.0.
        See “Time as a resolution parameter” in :cite:t:`Lambiotte2014`.
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain ``(obs_key, list_of_categories)``.
    key_added
        Key under which to add the cluster labels. (default: ``'louvain'``)
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    flavor
        Choose between to packages for computing the clustering.
    
        ``'vtraag'``
            Much more powerful than ``'igraph'``, and the default.
        ``'igraph'``
            Built in ``igraph`` method.
        ``'rapids'``
            GPU accelerated implementation.
    
            .. deprecated:: 1.10.0
                Use :func:`rapids_singlecell.tl.louvain` instead.
    directed
        Interpret the ``adjacency`` matrix as directed graph?
    use_weights
        Use weights from knn graph.
    partition_type
        Type of partition to use.
        Only a valid argument if ``flavor`` is ``'vtraag'``.
    partition_kwargs
        Key word arguments to pass to partitioning,
        if ``vtraag`` method is being used.
    neighbors_key
        Use neighbors connectivities as adjacency.
        If not specified, louvain looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, louvain looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    obsp
        Use .obsp[obsp] as adjacency. You can't specify both
        `obsp` and `neighbors_key` at the same time.
    copy
        Copy adata or modify it inplace.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obs['louvain' | key_added]` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.
    
    `adata.uns['louvain' | key_added]['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    """

    
    adata: Any = Field(
        ...,
        description='The annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    resolution: Any = Field(
        None,
        description="For the default flavor (``'vtraag'``) or for ```RAPIDS```, you can provide a\nresolution (higher resolution means finding more and smaller clusters),\nwhich defaults to 1.0.\nSee “Time as a resolution parameter” in :cite:t:`Lambiotte2014`.\nOriginal type annotation: float | None",
        title='Resolution',
    )
    random_state: Optional[Any] = Field(
        0,
        description='Change the initialization of the optimization.\nOriginal type annotation: _LegacyRandom',
        title='Random State',
    )
    restrict_to: Any = Field(
        None,
        description='Restrict the clustering to the categories within the key for sample\nannotation, tuple needs to contain ``(obs_key, list_of_categories)``.\nOriginal type annotation: tuple[str, Sequence[str]] | None',
        title='Restrict To',
    )
    key_added: Optional[Any] = Field(
        'louvain',
        description="Key under which to add the cluster labels. (default: ``'louvain'``)\nOriginal type annotation: str",
        title='Key Added',
    )
    adjacency: Any = Field(
        None,
        description='Sparse adjacency matrix of the graph, defaults to neighbors connectivities.\nOriginal type annotation: _CSMatrix | None',
        title='Adjacency',
    )
    flavor: Optional[Any] = Field(
        'vtraag',
        description="Choose between to packages for computing the clustering.\n\n``'vtraag'``\n    Much more powerful than ``'igraph'``, and the default.\n``'igraph'``\n    Built in ``igraph`` method.\n``'rapids'``\n    GPU accelerated implementation.\n\n    .. deprecated:: 1.10.0\n        Use :func:`rapids_singlecell.tl.louvain` instead.\nOriginal type annotation: Literal['vtraag', 'igraph', 'rapids']",
        title='Flavor',
    )
    directed: Optional[Any] = Field(
        True,
        description='Interpret the ``adjacency`` matrix as directed graph?\nOriginal type annotation: bool',
        title='Directed',
    )
    use_weights: Optional[Any] = Field(
        False,
        description='Use weights from knn graph.\nOriginal type annotation: bool',
        title='Use Weights',
    )
    partition_type: Any = Field(
        None,
        description="Type of partition to use.\nOnly a valid argument if ``flavor`` is ``'vtraag'``.\nOriginal type annotation: type[MutableVertexPartition] | None",
        title='Partition Type',
    )
    partition_kwargs: Optional[Any] = Field(
        {},
        description='Key word arguments to pass to partitioning,\nif ``vtraag`` method is being used.\nOriginal type annotation: Mapping[str, Any]',
        title='Partition Kwargs',
    )
    neighbors_key: Any = Field(
        None,
        description="Use neighbors connectivities as adjacency.\nIf not specified, louvain looks .obsp['connectivities'] for connectivities\n(default storage place for pp.neighbors).\nIf specified, louvain looks\n.obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.\nOriginal type annotation: str | None",
        title='Neighbors Key',
    )
    obsp: Any = Field(
        None,
        description="Use .obsp[obsp] as adjacency. You can't specify both\n`obsp` and `neighbors_key` at the same time.\nOriginal type annotation: str | None",
        title='Obsp',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Copy adata or modify it inplace.\nOriginal type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.tools.louvain')
    _products_original = PrivateAttr(default=['data.obs["louvain"]', 'data.uns["louvain"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyToolsUmap(BaseAPI):
    """
    Embed the neighborhood graph using UMAP :cite:p:`McInnes2018`.
    
    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph. tSNE, by contrast, optimizes the distribution of
    nearest-neighbor distances in the embedding such that these best match the
    distribution of distances in the high-dimensional space.
    We use the implementation of umap-learn_ :cite:p:`McInnes2018`.
    For a few comparisons of UMAP with tSNE, see :cite:t:`Becht2018`.
    
    .. _umap-learn: https://github.com/lmcinnes/umap
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    min_dist
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
    spread
        The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are.
    n_components
        The number of dimensions of the embedding.
    maxiter
        The number of iterations (epochs) of the optimization. Called `n_epochs`
        in the original UMAP.
    alpha
        The initial learning rate for the embedding optimization.
    gamma
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
    init_pos
        How to initialize the low dimensional embedding. Called `init` in the
        original UMAP. Options are:
    
        * Any key for `adata.obsm`.
        * 'paga': positions from :func:`~scanpy.pl.paga`.
        * 'spectral': use a spectral embedding of the graph.
        * 'random': assign initial embedding positions at random.
        * A numpy array of initial embedding positions.
    random_state
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState` or `Generator`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    a
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    b
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    method
        Chosen implementation.
    
        ``'umap'``
            Umap’s simplical set embedding.
        ``'rapids'``
            GPU accelerated implementation.
    
            .. deprecated:: 1.10.0
                Use :func:`rapids_singlecell.tl.umap` instead.
    key_added
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ `['X_umap']` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ `['umap']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ ``[key_added]`` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ ``[key_added]``.
    neighbors_key
        Umap looks in
        :attr:`~anndata.AnnData.uns`\ ``[neighbors_key]`` for neighbors settings and
        :attr:`~anndata.AnnData.obsp`\ ``[.uns[neighbors_key]['connectivities_key']]`` for connectivities.
    copy
        Return a copy instead of writing to adata.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obsm['X_umap' | key_added]` : :class:`numpy.ndarray` (dtype `float`)
        UMAP coordinates of data.
    `adata.uns['umap' | key_added]` : :class:`dict`
        UMAP parameters.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    min_dist: Optional[Any] = Field(
        0.5,
        description='The effective minimum distance between embedded points. Smaller values\nwill result in a more clustered/clumped embedding where nearby points on\nthe manifold are drawn closer together, while larger values will result\non a more even dispersal of points. The value should be set relative to\nthe ``spread`` value, which determines the scale at which embedded\npoints will be spread out. The default of in the `umap-learn` package is\n0.1.\nOriginal type annotation: float',
        title='Min Dist',
    )
    spread: Optional[Any] = Field(
        1.0,
        description='The effective scale of embedded points. In combination with `min_dist`\nthis determines how clustered/clumped the embedded points are.\nOriginal type annotation: float',
        title='Spread',
    )
    n_components: Optional[Any] = Field(
        2,
        description='The number of dimensions of the embedding.\nOriginal type annotation: int',
        title='N Components',
    )
    maxiter: Any = Field(
        None,
        description='The number of iterations (epochs) of the optimization. Called `n_epochs`\nin the original UMAP.\nOriginal type annotation: int | None',
        title='Maxiter',
    )
    alpha: Optional[Any] = Field(
        1.0,
        description='The initial learning rate for the embedding optimization.\nOriginal type annotation: float',
        title='Alpha',
    )
    gamma: Optional[Any] = Field(
        1.0,
        description='Weighting applied to negative samples in low dimensional embedding\noptimization. Values higher than one will result in greater weight\nbeing given to negative samples.\nOriginal type annotation: float',
        title='Gamma',
    )
    negative_sample_rate: Optional[Any] = Field(
        5,
        description='The number of negative edge/1-simplex samples to use per positive\nedge/1-simplex sample in optimizing the low dimensional embedding.\nOriginal type annotation: int',
        title='Negative Sample Rate',
    )
    init_pos: Optional[Any] = Field(
        'spectral',
        description="How to initialize the low dimensional embedding. Called `init` in the\noriginal UMAP. Options are:\n\n* Any key for `adata.obsm`.\n* 'paga': positions from :func:`~scanpy.pl.paga`.\n* 'spectral': use a spectral embedding of the graph.\n* 'random': assign initial embedding positions at random.\n* A numpy array of initial embedding positions.\nOriginal type annotation: _InitPos | np.ndarray | None",
        title='Init Pos',
    )
    random_state: Optional[Any] = Field(
        0,
        description='If `int`, `random_state` is the seed used by the random number generator;\nIf `RandomState` or `Generator`, `random_state` is the random number generator;\nIf `None`, the random number generator is the `RandomState` instance used\nby `np.random`.\nOriginal type annotation: _LegacyRandom',
        title='Random State',
    )
    a: Any = Field(
        None,
        description='More specific parameters controlling the embedding. If `None` these\nvalues are set automatically as determined by `min_dist` and\n`spread`.\nOriginal type annotation: float | None',
        title='A',
    )
    b: Any = Field(
        None,
        description='More specific parameters controlling the embedding. If `None` these\nvalues are set automatically as determined by `min_dist` and\n`spread`.\nOriginal type annotation: float | None',
        title='B',
    )
    method: Optional[Any] = Field(
        'umap',
        description="Chosen implementation.\n\n``'umap'``\n    Umap’s simplical set embedding.\n``'rapids'``\n    GPU accelerated implementation.\n\n    .. deprecated:: 1.10.0\n        Use :func:`rapids_singlecell.tl.umap` instead.\nOriginal type annotation: Literal['umap', 'rapids']",
        title='Method',
    )
    key_added: Any = Field(
        None,
        description="If not specified, the embedding is stored as\n:attr:`~anndata.AnnData.obsm`\\ `['X_umap']` and the the parameters in\n:attr:`~anndata.AnnData.uns`\\ `['umap']`.\nIf specified, the embedding is stored as\n:attr:`~anndata.AnnData.obsm`\\ ``[key_added]`` and the the parameters in\n:attr:`~anndata.AnnData.uns`\\ ``[key_added]``.\nOriginal type annotation: str | None",
        title='Key Added',
    )
    neighbors_key: Optional[Any] = Field(
        'neighbors',
        description="Umap looks in\n:attr:`~anndata.AnnData.uns`\\ ``[neighbors_key]`` for neighbors settings and\n:attr:`~anndata.AnnData.obsp`\\ ``[.uns[neighbors_key]['connectivities_key']]`` for connectivities.\nOriginal type annotation: str",
        title='Neighbors Key',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Return a copy instead of writing to adata.\nOriginal type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.tools.umap')
    _products_original = PrivateAttr(default=['data.obsm["X_umap"]', 'data.uns["umap"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyToolsTsne(BaseAPI):
    """
    t-SNE :cite:p:`vanDerMaaten2008,Amir2013,Pedregosa2011`.
    
    t-distributed stochastic neighborhood embedding (tSNE, :cite:t:`vanDerMaaten2008`) was
    proposed for visualizating single-cell data by :cite:t:`Amir2013`. Here, by default,
    we use the implementation of *scikit-learn* :cite:p:`Pedregosa2011`. You can achieve
    a huge speedup and better convergence if you install Multicore-tSNE_
    by :cite:t:`Ulyanov2016`, which will be automatically detected by Scanpy.
    
    .. _multicore-tsne: https://github.com/DmitryUlyanov/Multicore-TSNE
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_pcs
        Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.
    use_rep
        Use the indicated representation. `'X'` or any key for `.obsm` is valid.
        If `None`, the representation is chosen automatically:
        For `.n_vars` < :attr:`~scanpy._settings.ScanpyConfig.N_PCS` (default: 50), `.X` is used, otherwise 'X_pca' is used.
        If 'X_pca' is not present, it’s computed with default parameters or `n_pcs` if present.
    perplexity
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    metric
        Distance metric calculate neighbors on.
    early_exaggeration
        Controls how tight natural clusters in the original space are in the
        embedded space and how much space will be between them. For larger
        values, the space between natural clusters will be larger in the
        embedded space. Again, the choice of this parameter is not very
        critical. If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
    learning_rate
        Note that the R-package "Rtsne" uses a default of 200.
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.
    random_state
        Change this to use different intial states for the optimization.
        If `None`, the initial state is not reproducible.
    n_jobs
        Number of jobs for parallel computation.
        `None` means using :attr:`scanpy._settings.ScanpyConfig.n_jobs`.
    key_added
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ `['X_tsne']` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ `['tsne']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ ``[key_added]`` and the the parameters in
        :attr:`~anndata.AnnData.uns`\ ``[key_added]``.
    copy
        Return a copy instead of writing to `adata`.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obsm['X_tsne' | key_added]` : :class:`numpy.ndarray` (dtype `float`)
        tSNE coordinates of data.
    `adata.uns['tsne' | key_added]` : :class:`dict`
        tSNE parameters.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    n_pcs: Any = Field(
        None,
        description='Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.\nOriginal type annotation: int | None',
        title='N Pcs',
    )
    use_rep: Any = Field(
        None,
        description="Use the indicated representation. `'X'` or any key for `.obsm` is valid.\nIf `None`, the representation is chosen automatically:\nFor `.n_vars` < :attr:`~scanpy._settings.ScanpyConfig.N_PCS` (default: 50), `.X` is used, otherwise 'X_pca' is used.\nIf 'X_pca' is not present, it’s computed with default parameters or `n_pcs` if present.\nOriginal type annotation: str | None",
        title='Use Rep',
    )
    perplexity: Optional[Any] = Field(
        30,
        description='The perplexity is related to the number of nearest neighbors that\nis used in other manifold learning algorithms. Larger datasets\nusually require a larger perplexity. Consider selecting a value\nbetween 5 and 50. The choice is not extremely critical since t-SNE\nis quite insensitive to this parameter.\nOriginal type annotation: float',
        title='Perplexity',
    )
    metric: Optional[Any] = Field(
        'euclidean',
        description='Distance metric calculate neighbors on.\nOriginal type annotation: str',
        title='Metric',
    )
    early_exaggeration: Optional[Any] = Field(
        12,
        description='Controls how tight natural clusters in the original space are in the\nembedded space and how much space will be between them. For larger\nvalues, the space between natural clusters will be larger in the\nembedded space. Again, the choice of this parameter is not very\ncritical. If the cost function increases during initial optimization,\nthe early exaggeration factor or the learning rate might be too high.\nOriginal type annotation: float',
        title='Early Exaggeration',
    )
    learning_rate: Optional[Any] = Field(
        1000,
        description='Note that the R-package "Rtsne" uses a default of 200.\nThe learning rate can be a critical parameter. It should be\nbetween 100 and 1000. If the cost function increases during initial\noptimization, the early exaggeration factor or the learning rate\nmight be too high. If the cost function gets stuck in a bad local\nminimum increasing the learning rate helps sometimes.\nOriginal type annotation: float',
        title='Learning Rate',
    )
    random_state: Optional[Any] = Field(
        0,
        description='Change this to use different intial states for the optimization.\nIf `None`, the initial state is not reproducible.\nOriginal type annotation: _LegacyRandom',
        title='Random State',
    )
    use_fast_tsne: Optional[Any] = Field(
        False,
        description='No description available.\nOriginal type annotation: bool',
        title='Use Fast Tsne',
    )
    n_jobs: Any = Field(
        None,
        description='Number of jobs for parallel computation.\n`None` means using :attr:`scanpy._settings.ScanpyConfig.n_jobs`.\nOriginal type annotation: int | None',
        title='N Jobs',
    )
    key_added: Any = Field(
        None,
        description="If not specified, the embedding is stored as\n:attr:`~anndata.AnnData.obsm`\\ `['X_tsne']` and the the parameters in\n:attr:`~anndata.AnnData.uns`\\ `['tsne']`.\nIf specified, the embedding is stored as\n:attr:`~anndata.AnnData.obsm`\\ ``[key_added]`` and the the parameters in\n:attr:`~anndata.AnnData.uns`\\ ``[key_added]``.\nOriginal type annotation: str | None",
        title='Key Added',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Return a copy instead of writing to `adata`.\nOriginal type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.tools.tsne')
    _products_original = PrivateAttr(default=['data.obsm["X_tsne"]', 'data.uns["tsne"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyToolsDiffmap(BaseAPI):
    """
    Diffusion Maps :cite:p:`Coifman2005,Haghverdi2015,Wolf2018`.
    
    Diffusion maps :cite:p:`Coifman2005` have been proposed for visualizing single-cell
    data by :cite:t:`Haghverdi2015`. This tool uses the adapted Gaussian kernel suggested
    by :cite:t:`Haghverdi2016` with the implementation of :cite:t:`Wolf2018`.
    
    The width ("sigma") of the connectivity kernel is implicitly determined by
    the number of neighbors used to compute the single-cell graph in
    :func:`~scanpy.pp.neighbors`. To reproduce the original implementation
    using a Gaussian kernel, use `method=='gauss'` in
    :func:`~scanpy.pp.neighbors`. To use an exponential kernel, use the default
    `method=='umap'`. Differences between these options shouldn't usually be
    dramatic.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_comps
        The number of dimensions of the representation.
    neighbors_key
        If not specified, diffmap looks in .uns['neighbors'] for neighbors settings
        and .obsp['connectivities'] and .obsp['distances'] for connectivities and
        distances, respectively (default storage places for pp.neighbors).
        If specified, diffmap looks in .uns[neighbors_key] for neighbors settings and
        .obsp[.uns[neighbors_key]['connectivities_key']] and
        .obsp[.uns[neighbors_key]['distances_key']] for connectivities and distances,
        respectively.
    random_state
        A numpy random seed
    copy
        Return a copy instead of writing to adata.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obsm['X_diffmap']` : :class:`numpy.ndarray` (dtype `float`)
        Diffusion map representation of data, which is the right eigen basis of
        the transition matrix with eigenvectors as columns.
    
    `adata.uns['diffmap_evals']` : :class:`numpy.ndarray` (dtype `float`)
        Array of size (number of eigen vectors).
        Eigenvalues of transition matrix.
    
    Notes
    -----
    The 0-th column in `adata.obsm["X_diffmap"]` is the steady-state solution,
    which is non-informative in diffusion maps.
    Therefore, the first diffusion component is at index 1,
    e.g. `adata.obsm["X_diffmap"][:,1]`
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    n_comps: Optional[Any] = Field(
        15,
        description='The number of dimensions of the representation.\nOriginal type annotation: int',
        title='N Comps',
    )
    neighbors_key: Any = Field(
        None,
        description="If not specified, diffmap looks in .uns['neighbors'] for neighbors settings\nand .obsp['connectivities'] and .obsp['distances'] for connectivities and\ndistances, respectively (default storage places for pp.neighbors).\nIf specified, diffmap looks in .uns[neighbors_key] for neighbors settings and\n.obsp[.uns[neighbors_key]['connectivities_key']] and\n.obsp[.uns[neighbors_key]['distances_key']] for connectivities and distances,\nrespectively.\nOriginal type annotation: str | None",
        title='Neighbors Key',
    )
    random_state: Optional[Any] = Field(
        0,
        description='A numpy random seed\nOriginal type annotation: _LegacyRandom',
        title='Random State',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Return a copy instead of writing to adata.\nOriginal type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.tools.diffmap')
    _products_original = PrivateAttr(default=['data.obsm["X_diffmap"]', 'data.uns["diffmap_evals"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyToolsEmbeddingDensity(BaseAPI):
    """
    Calculate the density of cells in an embedding (per condition).
    
    Gaussian kernel density estimation is used to calculate the density of
    cells in an embedded space. This can be performed per category over a
    categorical cell annotation. The cell density can be plotted using the
    `pl.embedding_density` function.
    
    Note that density values are scaled to be between 0 and 1. Thus, the
    density value at each cell is only comparable to densities in
    the same category.
    
    Beware that the KDE estimate used (`scipy.stats.gaussian_kde`) becomes
    unreliable if you don't have enough cells in a category.
    
    This function was written by Sophie Tritschler and implemented into
    Scanpy by Malte Luecken.
    
    Parameters
    ----------
    adata
        The annotated data matrix.
    basis
        The embedding over which the density will be calculated. This embedded
        representation is found in `adata.obsm['X_[basis]']``.
    groupby
        Key for categorical observation/cell annotation for which densities
        are calculated per category.
    key_added
        Name of the `.obs` covariate that will be added with the density
        estimates.
    components
        The embedding dimensions over which the density should be calculated.
        This is limited to two components.
    
    Returns
    -------
    Sets the following fields (`key_added` defaults to `[basis]_density_[groupby]`, where `[basis]` is one of `umap`, `diffmap`, `pca`, `tsne`, or `draw_graph_fa` and `[groupby]` denotes the parameter input):
    
    `adata.obs[key_added]` : :class:`numpy.ndarray` (dtype `float`)
        Embedding density values for each cell.
    `adata.uns['[key_added]_params']` : :class:`dict`
        A dict with the values for the parameters `covariate` (for the `groupby` parameter) and `components`.
    
    Examples
    --------
    
    .. plot::
        :context: close-figs
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.umap(adata)
        sc.tl.embedding_density(adata, basis='umap', groupby='phase')
        sc.pl.embedding_density(
            adata, basis='umap', key='umap_density_phase', group='G1'
        )
    
    .. plot::
        :context: close-figs
    
        sc.pl.embedding_density(
            adata, basis='umap', key='umap_density_phase', group='S'
        )
    
    .. currentmodule:: scanpy
    
    See Also
    --------
    pl.embedding_density
    """

    
    adata: Any = Field(
        ...,
        description='The annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    basis: Optional[Any] = Field(
        'umap',
        description="The embedding over which the density will be calculated. This embedded\nrepresentation is found in `adata.obsm['X_[basis]']``.\nOriginal type annotation: str",
        title='Basis',
    )
    groupby: Any = Field(
        None,
        description='Key for categorical observation/cell annotation for which densities\nare calculated per category.\nOriginal type annotation: str | None',
        title='Groupby',
    )
    key_added: Any = Field(
        None,
        description='Name of the `.obs` covariate that will be added with the density\nestimates.\nOriginal type annotation: str | None',
        title='Key Added',
    )
    components: Any = Field(
        None,
        description='The embedding dimensions over which the density should be calculated.\nThis is limited to two components.\nOriginal type annotation: str | Sequence[str] | None',
        title='Components',
    )

    _api_name = PrivateAttr(default='scanpy.tools.embedding_density')
    _products_original = PrivateAttr(default=['data.obs["umap_density"]', 'data.uns["umap_density_params"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyToolsRankGenesGroups(BaseAPI):
    """
    Rank genes for characterizing groups.
    
    Expects logarithmized data.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The key of the observations grouping to consider.
    mask_var
        Select subset of genes to use in statistical tests.
    use_raw
        Use `raw` attribute of `adata` if present. The default behavior is to use `raw` if present.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.
    groups
        Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
        shall be restricted, or `'all'` (default), for all groups. Note that if
        `reference='rest'` all groups will still be used as the reference, not
        just those specified in `groups`.
    reference
        If `'rest'`, compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.
    n_genes
        The number of genes that appear in the returned tables.
        Defaults to all genes.
    method
        The default method is `'t-test'`,
        `'t-test_overestim_var'` overestimates variance of each group,
        `'wilcoxon'` uses Wilcoxon rank-sum,
        `'logreg'` uses logistic regression. See :cite:t:`Ntranos2019`,
        `here <https://github.com/scverse/scanpy/issues/95>`__ and `here
        <https://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters>`__,
        for why this is meaningful.
    corr_method
        p-value correction method.
        Used only for `'t-test'`, `'t-test_overestim_var'`, and `'wilcoxon'`.
    tie_correct
        Use tie correction for `'wilcoxon'` scores.
        Used only for `'wilcoxon'`.
    rankby_abs
        Rank genes by the absolute value of the score, not by the
        score. The returned scores are never the absolute values.
    pts
        Compute the fraction of cells expressing the genes.
    key_added
        The key in `adata.uns` information is saved to.
    copy
        Whether to copy `adata` or modify it inplace.
    kwds
        Are passed to test methods. Currently this affects only parameters that
        are passed to :class:`sklearn.linear_model.LogisticRegression`.
        For instance, you can pass `penalty='l1'` to try to come up with a
        minimal set of genes that are good predictors (sparse solution meaning
        few non-zero fitted coefficients).
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.uns['rank_genes_groups' | key_added]['names']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    `adata.uns['rank_genes_groups' | key_added]['scores']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the z-score
        underlying the computation of a p-value for each gene for each
        group. Ordered according to scores.
    `adata.uns['rank_genes_groups' | key_added]['logfoldchanges']` : structured :class:`numpy.ndarray` (dtype `object`)
        Structured array to be indexed by group id storing the log2
        fold change for each gene for each group. Ordered according to
        scores. Only provided if method is 't-test' like.
        Note: this is an approximation calculated from mean-log values.
    `adata.uns['rank_genes_groups' | key_added]['pvals']` : structured :class:`numpy.ndarray` (dtype `float`)
        p-values.
    `adata.uns['rank_genes_groups' | key_added]['pvals_adj']` : structured :class:`numpy.ndarray` (dtype `float`)
        Corrected p-values.
    `adata.uns['rank_genes_groups' | key_added]['pts']` : :class:`pandas.DataFrame` (dtype `float`)
        Fraction of cells expressing the genes for each group.
    `adata.uns['rank_genes_groups' | key_added]['pts_rest']` : :class:`pandas.DataFrame` (dtype `float`)
        Only if `reference` is set to `'rest'`.
        Fraction of cells from the union of the rest of each group
        expressing the genes.
    
    Notes
    -----
    There are slight inconsistencies depending on whether sparse
    or dense data are passed. See `here <https://github.com/scverse/scanpy/blob/main/tests/test_rank_genes_groups.py>`__.
    
    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")
    >>> # to visualize the results
    >>> sc.pl.rank_genes_groups(adata)
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    groupby: Any = Field(
        ...,
        description='The key of the observations grouping to consider.\nOriginal type annotation: str',
        title='Groupby',
    )
    mask_var: Any = Field(
        None,
        description='Select subset of genes to use in statistical tests.\nOriginal type annotation: NDArray[np.bool_] | str | None',
        title='Mask Var',
    )
    use_raw: Any = Field(
        None,
        description='Use `raw` attribute of `adata` if present. The default behavior is to use `raw` if present.\nOriginal type annotation: bool | None',
        title='Use Raw',
    )
    groups: Optional[Any] = Field(
        'all',
        description="Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison\nshall be restricted, or `'all'` (default), for all groups. Note that if\n`reference='rest'` all groups will still be used as the reference, not\njust those specified in `groups`.\nOriginal type annotation: Literal['all'] | Iterable[str]",
        title='Groups',
    )
    reference: Optional[Any] = Field(
        'rest',
        description="If `'rest'`, compare each group to the union of the rest of the group.\nIf a group identifier, compare with respect to this group.\nOriginal type annotation: str",
        title='Reference',
    )
    n_genes: Any = Field(
        None,
        description='The number of genes that appear in the returned tables.\nDefaults to all genes.\nOriginal type annotation: int | None',
        title='N Genes',
    )
    rankby_abs: Optional[Any] = Field(
        False,
        description='Rank genes by the absolute value of the score, not by the\nscore. The returned scores are never the absolute values.\nOriginal type annotation: bool',
        title='Rankby Abs',
    )
    pts: Optional[Any] = Field(
        False,
        description='Compute the fraction of cells expressing the genes.\nOriginal type annotation: bool',
        title='Pts',
    )
    key_added: Any = Field(
        None,
        description='The key in `adata.uns` information is saved to.\nOriginal type annotation: str | None',
        title='Key Added',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Whether to copy `adata` or modify it inplace.\nOriginal type annotation: bool',
        title='Copy',
    )
    method: Any = Field(
        None,
        description="The default method is `'t-test'`,\n`'t-test_overestim_var'` overestimates variance of each group,\n`'wilcoxon'` uses Wilcoxon rank-sum,\n`'logreg'` uses logistic regression. See :cite:t:`Ntranos2019`,\n`here <https://github.com/scverse/scanpy/issues/95>`__ and `here\n<https://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters>`__,\nfor why this is meaningful.\nOriginal type annotation: _Method | None",
        title='Method',
    )
    corr_method: Optional[Any] = Field(
        'benjamini-hochberg',
        description="p-value correction method.\nUsed only for `'t-test'`, `'t-test_overestim_var'`, and `'wilcoxon'`.\nOriginal type annotation: _CorrMethod",
        title='Corr Method',
    )
    tie_correct: Optional[Any] = Field(
        False,
        description="Use tie correction for `'wilcoxon'` scores.\nUsed only for `'wilcoxon'`.\nOriginal type annotation: bool",
        title='Tie Correct',
    )
    layer: Any = Field(
        None,
        description='Key from `adata.layers` whose value will be used to perform tests on.\nOriginal type annotation: str | None',
        title='Layer',
    )
    kwds: Any = Field(
        ...,
        description="Are passed to test methods. Currently this affects only parameters that\nare passed to :class:`sklearn.linear_model.LogisticRegression`.\nFor instance, you can pass `penalty='l1'` to try to come up with a\nminimal set of genes that are good predictors (sparse solution meaning\nfew non-zero fitted coefficients).",
        title='Kwds',
    )

    _api_name = PrivateAttr(default='scanpy.tools.rank_genes_groups')
    _products_original = PrivateAttr(default=['data.uns["rank_genes_groups"]'])
    _data_name = PrivateAttr(default='adata')

TOOLS_DICT = {
    "scanpy.tools.paga": ScanpyToolsPaga,
    "scanpy.tools.leiden": ScanpyToolsLeiden,
    "scanpy.tools.louvain": ScanpyToolsLouvain,
    "scanpy.tools.umap": ScanpyToolsUmap,
    "scanpy.tools.tsne": ScanpyToolsTsne,
    "scanpy.tools.diffmap": ScanpyToolsDiffmap,
    "scanpy.tools.embedding_density": ScanpyToolsEmbeddingDensity,
    "scanpy.tools.rank_genes_groups": ScanpyToolsRankGenesGroups,
}