from __future__ import annotations

from typing import Any, Optional

from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScanpyPreprocessingNeighbors(BaseAPI):
    """
    Compute the nearest neighbors distance matrix and a neighborhood graph of observations :cite:p:`McInnes2018`.
    
    The neighbor search efficiency of this heavily relies on UMAP :cite:p:`McInnes2018`,
    which also provides a method for estimating connectivities of data points -
    the connectivity of the manifold (`method=='umap'`). If `method=='gauss'`,
    connectivities are computed according to :cite:t:`Coifman2005`, in the adaption of
    :cite:t:`Haghverdi2016`.
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_neighbors
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.
    
        *ignored if ``transformer`` is an instance.*
    n_pcs
        Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.
    use_rep
        Use the indicated representation. `'X'` or any key for `.obsm` is valid.
        If `None`, the representation is chosen automatically:
        For `.n_vars` < :attr:`~scanpy._settings.ScanpyConfig.N_PCS` (default: 50), `.X` is used, otherwise 'X_pca' is used.
        If 'X_pca' is not present, it’s computed with default parameters or `n_pcs` if present.
    knn
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    method
        Use 'umap' :cite:p:`McInnes2018` or 'gauss' (Gauss kernel following :cite:t:`Coifman2005`
        with adaptive width :cite:t:`Haghverdi2016`) for computing connectivities.
    transformer
        Approximate kNN search implementation following the API of
        :class:`~sklearn.neighbors.KNeighborsTransformer`.
        See :doc:`/how-to/knn-transformers` for more details.
        Also accepts the following known options:
    
        `None` (the default)
            Behavior depends on data size.
            For small data, we will calculate exact kNN, otherwise we use
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
        `'pynndescent'`
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
        `'rapids'`
            A transformer based on :class:`cuml.neighbors.NearestNeighbors`.
    
            .. deprecated:: 1.10.0
               Use :func:`rapids_singlecell.pp.neighbors` instead.
    metric
        A known metric’s name or a callable that returns a distance.
    
        *ignored if ``transformer`` is an instance.*
    metric_kwds
        Options for the metric.
    
        *ignored if ``transformer`` is an instance.*
    random_state
        A numpy random seed.
    
        *ignored if ``transformer`` is an instance.*
    key_added
        If not specified, the neighbors data is stored in `.uns['neighbors']`,
        distances and connectivities are stored in `.obsp['distances']` and
        `.obsp['connectivities']` respectively.
        If specified, the neighbors data is added to .uns[key_added],
        distances are stored in `.obsp[key_added+'_distances']` and
        connectivities in `.obsp[key_added+'_connectivities']`.
    copy
        Return a copy instead of writing to adata.
    
    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:
    
    `adata.obsp['distances' | key_added+'_distances']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Distance matrix of the nearest neighbors search. Each row (cell) has `n_neighbors`-1 non-zero entries. These are the distances to their `n_neighbors`-1 nearest neighbors (excluding the cell itself).
    `adata.obsp['connectivities' | key_added+'_connectivities']` : :class:`scipy.sparse._csr.csr_matrix` (dtype `float`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    `adata.uns['neighbors' | key_added]` : :class:`dict`
        neighbors parameters.
    
    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> # Basic usage
    >>> sc.pp.neighbors(adata, 20, metric="cosine")
    >>> # Provide your own transformer for more control and flexibility
    >>> from sklearn.neighbors import KNeighborsTransformer
    >>> transformer = KNeighborsTransformer(
    ...     n_neighbors=10, metric="manhattan", algorithm="kd_tree"
    ... )
    >>> sc.pp.neighbors(adata, transformer=transformer)
    >>> # now you can e.g. access the index: `transformer._tree`
    
    See Also
    --------
    :doc:`/how-to/knn-transformers`
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    n_neighbors: Optional[Any] = Field(
        15,
        description='The size of local neighborhood (in terms of number of neighboring data\npoints) used for manifold approximation. Larger values result in more\nglobal views of the manifold, while smaller values result in more local\ndata being preserved. In general values should be in the range 2 to 100.\nIf `knn` is `True`, number of nearest neighbors to be searched. If `knn`\nis `False`, a Gaussian kernel width is set to the distance of the\n`n_neighbors` neighbor.\n\n*ignored if ``transformer`` is an instance.*\nOriginal type annotation: int',
        title='N Neighbors',
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
    knn: Optional[Any] = Field(
        True,
        description='If `True`, use a hard threshold to restrict the number of neighbors to\n`n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian\nKernel to assign low weights to neighbors more distant than the\n`n_neighbors` nearest neighbor.\nOriginal type annotation: bool',
        title='Knn',
    )
    method: Optional[Any] = Field(
        'umap',
        description="Use 'umap' :cite:p:`McInnes2018` or 'gauss' (Gauss kernel following :cite:t:`Coifman2005`\nwith adaptive width :cite:t:`Haghverdi2016`) for computing connectivities.\nOriginal type annotation: _Method",
        title='Method',
    )
    transformer: Any = Field(
        None,
        description="Approximate kNN search implementation following the API of\n:class:`~sklearn.neighbors.KNeighborsTransformer`.\nSee :doc:`/how-to/knn-transformers` for more details.\nAlso accepts the following known options:\n\n`None` (the default)\n    Behavior depends on data size.\n    For small data, we will calculate exact kNN, otherwise we use\n    :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`\n`'pynndescent'`\n    :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`\n`'rapids'`\n    A transformer based on :class:`cuml.neighbors.NearestNeighbors`.\n\n    .. deprecated:: 1.10.0\n       Use :func:`rapids_singlecell.pp.neighbors` instead.\nOriginal type annotation: KnnTransformerLike | _KnownTransformer | None",
        title='Transformer',
    )
    metric: Optional[Any] = Field(
        'euclidean',
        description='A known metric’s name or a callable that returns a distance.\n\n*ignored if ``transformer`` is an instance.*\nOriginal type annotation: _Metric | _MetricFn',
        title='Metric',
    )
    metric_kwds: Optional[Any] = Field(
        {},
        description='Options for the metric.\n\n*ignored if ``transformer`` is an instance.*\nOriginal type annotation: Mapping[str, Any]',
        title='Metric Kwds',
    )
    random_state: Optional[Any] = Field(
        0,
        description='A numpy random seed.\n\n*ignored if ``transformer`` is an instance.*\nOriginal type annotation: _LegacyRandom',
        title='Random State',
    )
    key_added: Any = Field(
        None,
        description="If not specified, the neighbors data is stored in `.uns['neighbors']`,\ndistances and connectivities are stored in `.obsp['distances']` and\n`.obsp['connectivities']` respectively.\nIf specified, the neighbors data is added to .uns[key_added],\ndistances are stored in `.obsp[key_added+'_distances']` and\nconnectivities in `.obsp[key_added+'_connectivities']`.\nOriginal type annotation: str | None",
        title='Key Added',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Return a copy instead of writing to adata.\nOriginal type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.preprocessing.neighbors')
    _products_original = PrivateAttr(default=['data.uns["neighbors"]', 'data.obsp["distances"]', 'data.obsp["connectivities"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyPreprocessingLogP(BaseAPI):
    """
    Logarithmize the data matrix.
    
    Computes :math:`X = \log(X + 1)`,
    where :math:`log` denotes the natural logarithm unless a different base is given.
    
    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    base
        Base of the logarithm. Natural logarithm is used by default.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    chunked
        Process the data matrix in chunks, which will save memory.
        Applies only to :class:`~anndata.AnnData`.
    chunk_size
        `n_obs` of the chunks to process the data in.
    layer
        Entry of layers to transform.
    obsm
        Entry of obsm to transform.
    
    Returns
    -------
    Returns or updates `data`, depending on `copy`.
    """

    
    data: Any = Field(
        ...,
        description='The (annotated) data matrix of shape `n_obs` × `n_vars`.\nRows correspond to cells and columns to genes.\nOriginal type annotation: AnnData | np.ndarray | _CSMatrix',
        title='Data',
    )
    base: Any = Field(
        None,
        description='Base of the logarithm. Natural logarithm is used by default.\nOriginal type annotation: Number | None',
        title='Base',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='If an :class:`~anndata.AnnData` is passed, determines whether a copy\nis returned.\nOriginal type annotation: bool',
        title='Copy',
    )
    chunked: Any = Field(
        None,
        description='Process the data matrix in chunks, which will save memory.\nApplies only to :class:`~anndata.AnnData`.\nOriginal type annotation: bool | None',
        title='Chunked',
    )
    chunk_size: Any = Field(
        None,
        description='`n_obs` of the chunks to process the data in.\nOriginal type annotation: int | None',
        title='Chunk Size',
    )
    layer: Any = Field(
        None,
        description='Entry of layers to transform.\nOriginal type annotation: str | None',
        title='Layer',
    )
    obsm: Any = Field(
        None,
        description='Entry of obsm to transform.\nOriginal type annotation: str | None',
        title='Obsm',
    )

    _api_name = PrivateAttr(default='scanpy.preprocessing.log1p')
    _products_original = PrivateAttr(default=['data.X'])
    _data_name = PrivateAttr(default='data')



class ScanpyPreprocessingHighlyVariableGenes(BaseAPI):
    """
    Annotate highly variable genes :cite:p:`Satija2015,Zheng2017,Stuart2019`.
    
    Expects logarithmized data, except when `flavor='seurat_v3'`/`'seurat_v3_paper'`, in which count
    data is expected.
    
    Depending on `flavor`, this reproduces the R-implementations of Seurat
    :cite:p:`Satija2015`, Cell Ranger :cite:p:`Zheng2017`, and Seurat v3 :cite:p:`Stuart2019`.
    
    `'seurat_v3'`/`'seurat_v3_paper'` requires `scikit-misc` package. If you plan to use this flavor, consider
    installing `scanpy` with this optional dependency: `scanpy[skmisc]`.
    
    For the dispersion-based methods (`flavor='seurat'` :cite:t:`Satija2015` and
    `flavor='cell_ranger'` :cite:t:`Zheng2017`), the normalized dispersion is obtained
    by scaling with the mean and standard deviation of the dispersions for genes
    falling into a given bin for mean expression of genes. This means that for each
    bin of mean expression, highly variable genes are selected.
    
    For `flavor='seurat_v3'`/`'seurat_v3_paper'` :cite:p:`Stuart2019`, a normalized variance for each gene
    is computed. First, the data are standardized (i.e., z-score normalization
    per feature) with a regularized standard deviation. Next, the normalized variance
    is computed as the variance of each gene after the transformation. Genes are ranked
    by the normalized variance.
    Only if `batch_key` is not `None`, the two flavors differ: For `flavor='seurat_v3'`, genes are first sorted by the median (across batches) rank, with ties broken by the number of batches a gene is a HVG.
    For `flavor='seurat_v3_paper'`, genes are first sorted by the number of batches a gene is a HVG, with ties broken by the median (across batches) rank.
    
    The following may help when comparing to Seurat's naming:
    If `batch_key=None` and `flavor='seurat'`, this mimics Seurat's `FindVariableFeatures(…, method='mean.var.plot')`.
    If `batch_key=None` and `flavor='seurat_v3'`/`flavor='seurat_v3_paper'`, this mimics Seurat's `FindVariableFeatures(..., method='vst')`.
    If `batch_key` is not `None` and `flavor='seurat_v3_paper'`, this mimics Seurat's `SelectIntegrationFeatures`.
    
    See also `scanpy.experimental.pp._highly_variable_genes` for additional flavors
    (e.g. Pearson residuals).
    
    Parameters
    ----------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    layer
        If provided, use `adata.layers[layer]` for expression values instead of `adata.X`.
    n_top_genes
        Number of highly-variable genes to keep. Mandatory if `flavor='seurat_v3'`.
    min_mean
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.
    max_mean
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.
    min_disp
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.
    max_disp
        If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the
        normalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.
    span
        The fraction of the data (cells) used when estimating the variance in the loess
        model fit if `flavor='seurat_v3'`.
    n_bins
        Number of bins for binning the mean gene expression. Normalization is
        done with respect to each bin. If just a single gene falls into a bin,
        the normalized dispersion is artificially set to 1. You'll be informed
        about this if you set `settings.verbosity = 4`.
    flavor
        Choose the flavor for identifying highly variable genes. For the dispersion
        based methods in their default workflows, Seurat passes the cutoffs whereas
        Cell Ranger passes `n_top_genes`.
    subset
        Inplace subset to highly-variable genes if `True` otherwise merely indicate
        highly variable genes.
    inplace
        Whether to place calculated metrics in `.var` or return them.
    batch_key
        If specified, highly-variable genes are selected within each batch separately and merged.
        This simple process avoids the selection of batch-specific genes and acts as a
        lightweight batch correction method. For all flavors, except `seurat_v3`, genes are first sorted
        by how many batches they are a HVG. For dispersion-based flavors ties are broken
        by normalized dispersion. For `flavor = 'seurat_v3_paper'`, ties are broken by the median
        (across batches) rank based on within-batch normalized variance.
    check_values
        Check if counts in selected layer are integers. A Warning is returned if set to True.
        Only used if `flavor='seurat_v3'`/`'seurat_v3_paper'`.
    
    Returns
    -------
    Returns a :class:`pandas.DataFrame` with calculated metrics if `inplace=False`, else returns an `AnnData` object where it sets the following field:
    
    `adata.var['highly_variable']` : :class:`pandas.Series` (dtype `bool`)
        boolean indicator of highly-variable genes
    `adata.var['means']` : :class:`pandas.Series` (dtype `float`)
        means per gene
    `adata.var['dispersions']` : :class:`pandas.Series` (dtype `float`)
        For dispersion-based flavors, dispersions per gene
    `adata.var['dispersions_norm']` : :class:`pandas.Series` (dtype `float`)
        For dispersion-based flavors, normalized dispersions per gene
    `adata.var['variances']` : :class:`pandas.Series` (dtype `float`)
        For `flavor='seurat_v3'`/`'seurat_v3_paper'`, variance per gene
    `adata.var['variances_norm']`/`'seurat_v3_paper'` : :class:`pandas.Series` (dtype `float`)
        For `flavor='seurat_v3'`/`'seurat_v3_paper'`, normalized variance per gene, averaged in
        the case of multiple batches
    `adata.var['highly_variable_rank']` : :class:`pandas.Series` (dtype `float`)
        For `flavor='seurat_v3'`/`'seurat_v3_paper'`, rank of the gene according to normalized
        variance, in case of multiple batches description above
    `adata.var['highly_variable_nbatches']` : :class:`pandas.Series` (dtype `int`)
        If `batch_key` is given, this denotes in how many batches genes are detected as HVG
    `adata.var['highly_variable_intersection']` : :class:`pandas.Series` (dtype `bool`)
        If `batch_key` is given, this denotes the genes that are highly variable in all batches
    
    Notes
    -----
    This function replaces :func:`~scanpy.pp.filter_genes_dispersion`.
    """

    
    adata: Any = Field(
        ...,
        description='The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond\nto cells and columns to genes.\nOriginal type annotation: AnnData',
        title='Adata',
    )
    layer: Any = Field(
        None,
        description='If provided, use `adata.layers[layer]` for expression values instead of `adata.X`.\nOriginal type annotation: str | None',
        title='Layer',
    )
    n_top_genes: Any = Field(
        None,
        description="Number of highly-variable genes to keep. Mandatory if `flavor='seurat_v3'`.\nOriginal type annotation: int | None",
        title='N Top Genes',
    )
    min_disp: Optional[Any] = Field(
        0.5,
        description="If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the\nnormalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.\nOriginal type annotation: float",
        title='Min Disp',
    )
    max_disp: Optional[Any] = Field(
        'inf',
        description="If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the\nnormalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.\nOriginal type annotation: float",
        title='Max Disp',
    )
    min_mean: Optional[Any] = Field(
        0.0125,
        description="If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the\nnormalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.\nOriginal type annotation: float",
        title='Min Mean',
    )
    max_mean: Optional[Any] = Field(
        3,
        description="If `n_top_genes` unequals `None`, this and all other cutoffs for the means and the\nnormalized dispersions are ignored. Ignored if `flavor='seurat_v3'`.\nOriginal type annotation: float",
        title='Max Mean',
    )
    span: Optional[Any] = Field(
        0.3,
        description="The fraction of the data (cells) used when estimating the variance in the loess\nmodel fit if `flavor='seurat_v3'`.\nOriginal type annotation: float",
        title='Span',
    )
    n_bins: Optional[Any] = Field(
        20,
        description="Number of bins for binning the mean gene expression. Normalization is\ndone with respect to each bin. If just a single gene falls into a bin,\nthe normalized dispersion is artificially set to 1. You'll be informed\nabout this if you set `settings.verbosity = 4`.\nOriginal type annotation: int",
        title='N Bins',
    )
    flavor: Optional[Any] = Field(
        'seurat',
        description="Choose the flavor for identifying highly variable genes. For the dispersion\nbased methods in their default workflows, Seurat passes the cutoffs whereas\nCell Ranger passes `n_top_genes`.\nOriginal type annotation: Literal['seurat', 'cell_ranger', 'seurat_v3', 'seurat_v3_paper']",
        title='Flavor',
    )
    subset: Optional[Any] = Field(
        False,
        description='Inplace subset to highly-variable genes if `True` otherwise merely indicate\nhighly variable genes.\nOriginal type annotation: bool',
        title='Subset',
    )
    inplace: Optional[Any] = Field(
        True,
        description='Whether to place calculated metrics in `.var` or return them.\nOriginal type annotation: bool',
        title='Inplace',
    )
    batch_key: Any = Field(
        None,
        description="If specified, highly-variable genes are selected within each batch separately and merged.\nThis simple process avoids the selection of batch-specific genes and acts as a\nlightweight batch correction method. For all flavors, except `seurat_v3`, genes are first sorted\nby how many batches they are a HVG. For dispersion-based flavors ties are broken\nby normalized dispersion. For `flavor = 'seurat_v3_paper'`, ties are broken by the median\n(across batches) rank based on within-batch normalized variance.\nOriginal type annotation: str | None",
        title='Batch Key',
    )
    check_values: Optional[Any] = Field(
        True,
        description="Check if counts in selected layer are integers. A Warning is returned if set to True.\nOnly used if `flavor='seurat_v3'`/`'seurat_v3_paper'`.\nOriginal type annotation: bool",
        title='Check Values',
    )

    _api_name = PrivateAttr(default='scanpy.preprocessing.highly_variable_genes')
    _products_original = PrivateAttr(default=['data.var["highly_variable"]', 'data.var["means"]', 'data.var["dispersions"]', 'data.var["dispersions_norm"]', 'data.var["variances"]', 'data.var["variances_norm"]', 'data.var["highly_variable_rank"]', 'data.var["highly_variable_nbatches"]', 'data.var["highly_variable_intersection"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyPreprocessingPca(BaseAPI):
    """
    Principal component analysis :cite:p:`Pedregosa2011`.
    
    Computes PCA coordinates, loadings and variance decomposition.
    Uses the implementation of *scikit-learn* :cite:p:`Pedregosa2011`.
    
    .. versionchanged:: 1.5.0
    
        In previous versions, computing a PCA on a sparse matrix would make
        a dense copy of the array for mean centering.
        As of scanpy 1.5.0, mean centering is implicit.
        While results are extremely similar, they are not exactly the same.
        If you would like to reproduce the old results, pass a dense array.
    
    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    n_comps
        Number of principal components to compute. Defaults to 50, or 1 - minimum
        dimension size of selected representation.
    layer
        If provided, which element of layers to use for PCA.
    zero_center
        If `True`, compute standard PCA from covariance matrix.
        If `False`, omit zero-centering variables
        (uses *scikit-learn* :class:`~sklearn.decomposition.TruncatedSVD` or
        *dask-ml* :class:`~dask_ml.decomposition.TruncatedSVD`),
        which allows to handle sparse input efficiently.
        Passing `None` decides automatically based on sparseness of the data.
    svd_solver
        SVD solver to use:
    
        `None`
            See `chunked` and `zero_center` descriptions to determine which class will be used.
            Depending on the class and the type of X different values for default will be set.
            For sparse *dask* arrays, will use `'covariance_eigh'`.
            If *scikit-learn* :class:`~sklearn.decomposition.PCA` is used, will give `'arpack'`,
            if *scikit-learn* :class:`~sklearn.decomposition.TruncatedSVD` is used, will give `'randomized'`,
            if *dask-ml* :class:`~dask_ml.decomposition.PCA` or :class:`~dask_ml.decomposition.IncrementalPCA` is used, will give `'auto'`,
            if *dask-ml* :class:`~dask_ml.decomposition.TruncatedSVD` is used, will give `'tsqr'`
        `'arpack'`
            for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)
            Not available with *dask* arrays.
        `'covariance_eigh'`
            Classic eigendecomposition of the covariance matrix, suited for tall-and-skinny matrices.
            With dask, array must be CSR or dense and chunked as (N, adata.shape[1]).
        `'randomized'`
            for the randomized algorithm due to Halko (2009). For *dask* arrays,
            this will use :func:`~dask.array.linalg.svd_compressed`.
        `'auto'`
            chooses automatically depending on the size of the problem.
        `'tsqr'`
            Only available with dense *dask* arrays. "tsqr"
            algorithm from Benson et. al. (2013).
    
        .. versionchanged:: 1.9.3
           Default value changed from `'arpack'` to None.
        .. versionchanged:: 1.4.5
           Default value changed from `'auto'` to `'arpack'`.
    
        Efficient computation of the principal components of a sparse matrix
        currently only works with the `'arpack`' or `'covariance_eigh`' solver.
    
        If X is a sparse *dask* array, a custom `'covariance_eigh'` solver will be used.
        If X is a dense *dask* array, *dask-ml* classes :class:`~dask_ml.decomposition.PCA`,
        :class:`~dask_ml.decomposition.IncrementalPCA`, or
        :class:`~dask_ml.decomposition.TruncatedSVD` will be used.
        Otherwise their *scikit-learn* counterparts :class:`~sklearn.decomposition.PCA`,
        :class:`~sklearn.decomposition.IncrementalPCA`, or
        :class:`~sklearn.decomposition.TruncatedSVD` will be used.
    random_state
        Change to use different initial states for the optimization.
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “Returns”.
    mask_var
        To run only on a certain set of genes given by a boolean array
        or a string referring to an array in :attr:`~anndata.AnnData.var`.
        By default, uses `.var['highly_variable']` if available, else everything.
    use_highly_variable
        Whether to use highly variable genes only, stored in
        `.var['highly_variable']`.
        By default uses them if they have been determined beforehand.
    
        .. deprecated:: 1.10.0
           Use `mask_var` instead
    
    layer
        Layer of `adata` to use as expression values.
    dtype
        Numpy data type string to which to convert the result.
    chunked
        If `True`, perform an incremental PCA on segments of `chunk_size`.
        The incremental PCA automatically zero centers and ignores settings of
        `random_seed` and `svd_solver`. Uses sklearn :class:`~sklearn.decomposition.IncrementalPCA` or
        *dask-ml* :class:`~dask_ml.decomposition.IncrementalPCA`. If `False`, perform a full PCA and
        use sklearn :class:`~sklearn.decomposition.PCA` or
        *dask-ml* :class:`~dask_ml.decomposition.PCA`
    chunk_size
        Number of observations to include in each chunk.
        Required if `chunked=True` was passed.
    key_added
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ `['X_pca']`, the loadings as
        :attr:`~anndata.AnnData.varm`\ `['PCs']`, and the the parameters in
        :attr:`~anndata.AnnData.uns`\ `['pca']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ ``[key_added]``, the loadings as
        :attr:`~anndata.AnnData.varm`\ ``[key_added]``, and the the parameters in
        :attr:`~anndata.AnnData.uns`\ ``[key_added]``.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.
    
    Returns
    -------
    If `data` is array-like and `return_info=False` was passed,
    this function returns the PCA representation of `data` as an
    array of the same type as the input array.
    
    Otherwise, it returns `None` if `copy=False`, else an updated `AnnData` object.
    Sets the following fields:
    
    `.obsm['X_pca' | key_added]` : :class:`~scipy.sparse.csr_matrix` | :class:`~scipy.sparse.csc_matrix` | :class:`~numpy.ndarray` (shape `(adata.n_obs, n_comps)`)
        PCA representation of data.
    `.varm['PCs' | key_added]` : :class:`~numpy.ndarray` (shape `(adata.n_vars, n_comps)`)
        The principal components containing the loadings.
    `.uns['pca' | key_added]['variance_ratio']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Ratio of explained variance.
    `.uns['pca' | key_added]['variance']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Explained variance, equivalent to the eigenvalues of the
        covariance matrix.
    """

    
    data: Any = Field(
        ...,
        description='The (annotated) data matrix of shape `n_obs` × `n_vars`.\nRows correspond to cells and columns to genes.\nOriginal type annotation: AnnData | np.ndarray | _CSMatrix',
        title='Data',
    )
    n_comps: Any = Field(
        None,
        description='Number of principal components to compute. Defaults to 50, or 1 - minimum\ndimension size of selected representation.\nOriginal type annotation: int | None',
        title='N Comps',
    )
    layer: Any = Field(
        None,
        description='Layer of `adata` to use as expression values.\nOriginal type annotation: str | None',
        title='Layer',
    )
    zero_center: Optional[Any] = Field(
        True,
        description='If `True`, compute standard PCA from covariance matrix.\nIf `False`, omit zero-centering variables\n(uses *scikit-learn* :class:`~sklearn.decomposition.TruncatedSVD` or\n*dask-ml* :class:`~dask_ml.decomposition.TruncatedSVD`),\nwhich allows to handle sparse input efficiently.\nPassing `None` decides automatically based on sparseness of the data.\nOriginal type annotation: bool | None',
        title='Zero Center',
    )
    svd_solver: Any = Field(
        None,
        description='SVD solver to use:\n\n`None`\n    See `chunked` and `zero_center` descriptions to determine which class will be used.\n    Depending on the class and the type of X different values for default will be set.\n    For sparse *dask* arrays, will use `\'covariance_eigh\'`.\n    If *scikit-learn* :class:`~sklearn.decomposition.PCA` is used, will give `\'arpack\'`,\n    if *scikit-learn* :class:`~sklearn.decomposition.TruncatedSVD` is used, will give `\'randomized\'`,\n    if *dask-ml* :class:`~dask_ml.decomposition.PCA` or :class:`~dask_ml.decomposition.IncrementalPCA` is used, will give `\'auto\'`,\n    if *dask-ml* :class:`~dask_ml.decomposition.TruncatedSVD` is used, will give `\'tsqr\'`\n`\'arpack\'`\n    for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)\n    Not available with *dask* arrays.\n`\'covariance_eigh\'`\n    Classic eigendecomposition of the covariance matrix, suited for tall-and-skinny matrices.\n    With dask, array must be CSR or dense and chunked as (N, adata.shape[1]).\n`\'randomized\'`\n    for the randomized algorithm due to Halko (2009). For *dask* arrays,\n    this will use :func:`~dask.array.linalg.svd_compressed`.\n`\'auto\'`\n    chooses automatically depending on the size of the problem.\n`\'tsqr\'`\n    Only available with dense *dask* arrays. "tsqr"\n    algorithm from Benson et. al. (2013).\n\n.. versionchanged:: 1.9.3\n   Default value changed from `\'arpack\'` to None.\n.. versionchanged:: 1.4.5\n   Default value changed from `\'auto\'` to `\'arpack\'`.\n\nEfficient computation of the principal components of a sparse matrix\ncurrently only works with the `\'arpack`\' or `\'covariance_eigh`\' solver.\n\nIf X is a sparse *dask* array, a custom `\'covariance_eigh\'` solver will be used.\nIf X is a dense *dask* array, *dask-ml* classes :class:`~dask_ml.decomposition.PCA`,\n:class:`~dask_ml.decomposition.IncrementalPCA`, or\n:class:`~dask_ml.decomposition.TruncatedSVD` will be used.\nOtherwise their *scikit-learn* counterparts :class:`~sklearn.decomposition.PCA`,\n:class:`~sklearn.decomposition.IncrementalPCA`, or\n:class:`~sklearn.decomposition.TruncatedSVD` will be used.\nOriginal type annotation: SvdSolver | None',
        title='Svd Solver',
    )
    random_state: Optional[Any] = Field(
        0,
        description='Change to use different initial states for the optimization.\nOriginal type annotation: _LegacyRandom',
        title='Random State',
    )
    return_info: Optional[Any] = Field(
        False,
        description='Only relevant when not passing an :class:`~anndata.AnnData`:\nsee “Returns”.\nOriginal type annotation: bool',
        title='Return Info',
    )
    mask_var: Optional[Any] = Field(
        0,
        description="To run only on a certain set of genes given by a boolean array\nor a string referring to an array in :attr:`~anndata.AnnData.var`.\nBy default, uses `.var['highly_variable']` if available, else everything.\nOriginal type annotation: NDArray[np.bool_] | str | None | Empty",
        title='Mask Var',
    )
    use_highly_variable: Any = Field(
        None,
        description="Whether to use highly variable genes only, stored in\n`.var['highly_variable']`.\nBy default uses them if they have been determined beforehand.\n\n.. deprecated:: 1.10.0\n   Use `mask_var` instead\nOriginal type annotation: bool | None",
        title='Use Highly Variable',
    )
    dtype: Optional[Any] = Field(
        'float32',
        description='Numpy data type string to which to convert the result.\nOriginal type annotation: DTypeLike',
        title='Dtype',
    )
    chunked: Optional[Any] = Field(
        False,
        description='If `True`, perform an incremental PCA on segments of `chunk_size`.\nThe incremental PCA automatically zero centers and ignores settings of\n`random_seed` and `svd_solver`. Uses sklearn :class:`~sklearn.decomposition.IncrementalPCA` or\n*dask-ml* :class:`~dask_ml.decomposition.IncrementalPCA`. If `False`, perform a full PCA and\nuse sklearn :class:`~sklearn.decomposition.PCA` or\n*dask-ml* :class:`~dask_ml.decomposition.PCA`\nOriginal type annotation: bool',
        title='Chunked',
    )
    chunk_size: Any = Field(
        None,
        description='Number of observations to include in each chunk.\nRequired if `chunked=True` was passed.\nOriginal type annotation: int | None',
        title='Chunk Size',
    )
    key_added: Any = Field(
        None,
        description="If not specified, the embedding is stored as\n:attr:`~anndata.AnnData.obsm`\\ `['X_pca']`, the loadings as\n:attr:`~anndata.AnnData.varm`\\ `['PCs']`, and the the parameters in\n:attr:`~anndata.AnnData.uns`\\ `['pca']`.\nIf specified, the embedding is stored as\n:attr:`~anndata.AnnData.obsm`\\ ``[key_added]``, the loadings as\n:attr:`~anndata.AnnData.varm`\\ ``[key_added]``, and the the parameters in\n:attr:`~anndata.AnnData.uns`\\ ``[key_added]``.\nOriginal type annotation: str | None",
        title='Key Added',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='If an :class:`~anndata.AnnData` is passed, determines whether a copy\nis returned. Is ignored otherwise.\nOriginal type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.preprocessing.pca')
    _products_original = PrivateAttr(default=['data.obsm["X_pca"]', 'data.varm["PCs"]', 'data.uns["pca"]["variance_ratio"]', 'data.uns["pca"]["variance"]'])
    _data_name = PrivateAttr(default='data')

TOOLS_DICT = {
    "scanpy.preprocessing.neighbors": ScanpyPreprocessingNeighbors,
    "scanpy.preprocessing.log1p": ScanpyPreprocessingLogP,
    "scanpy.preprocessing.highly_variable_genes": ScanpyPreprocessingHighlyVariableGenes,
    "scanpy.preprocessing.pca": ScanpyPreprocessingPca,
}