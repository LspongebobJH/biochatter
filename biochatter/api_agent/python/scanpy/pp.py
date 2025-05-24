from __future__ import annotations

from typing import Any, Optional

from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScanpyPreprocessingNeighbors(BaseAPI):
    
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