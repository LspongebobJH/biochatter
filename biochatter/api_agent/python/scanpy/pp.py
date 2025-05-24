from __future__ import annotations

from typing import Any, Optional

from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScanpyPreprocessingNeighbors(BaseAPI):
    """
    Compute the nearest neighbors distance matrix and a neighborhood graph of observations. The efficiency of neighbor search heavily relies on UMAP, which also estimates connectivities of data points. Connectivities are computed differently based on the chosen method.
    """

    
    adata: Any = Field(
        ...,
        description='Annotated data matrix. Original type annotation: AnnData',
        title='Adata',
    )
    n_neighbors: Optional[Any] = Field(
        15,
        description='The size of local neighborhood used for manifold approximation. Larger values provide more global views, while smaller values preserve more local data. Should be in the range 2 to 100. Ignored if transformer is an instance. Original type annotation: int',
        title='N Neighbors',
    )
    n_pcs: Any = Field(
        None,
        description='Specifies the number of principal components to use. If n_pcs==0, uses .X if use_rep is None. Original type annotation: int | None',
        title='N Pcs',
    )
    use_rep: Any = Field(
        None,
        description='Indicates the representation to use. If None, representation is chosen automatically based on data size. Original type annotation: str | None',
        title='Use Rep',
    )
    knn: Optional[Any] = Field(
        True,
        description='Determines whether to use a hard threshold for neighbors or a Gaussian Kernel. Original type annotation: bool',
        title='Knn',
    )
    method: Optional[Any] = Field(
        'umap',
        description="Specifies the method for computing connectivities, either 'umap' or 'gauss'. Original type annotation: _Method",
        title='Method',
    )
    transformer: Any = Field(
        None,
        description="Implementation for approximate kNN search. Accepts various options like 'pynndescent' or 'rapids'. Original type annotation: KnnTransformerLike | _KnownTransformer | None",
        title='Transformer',
    )
    metric: Optional[Any] = Field(
        'euclidean',
        description='Name of a known metric or a callable that returns a distance. Ignored if transformer is an instance. Original type annotation: _Metric | _MetricFn',
        title='Metric',
    )
    metric_kwds: Optional[Any] = Field(
        {},
        description='Options for the metric. Ignored if transformer is an instance. Original type annotation: Mapping[str, Any]',
        title='Metric Kwds',
    )
    random_state: Optional[Any] = Field(
        0,
        description='Numpy random seed. Ignored if transformer is an instance. Original type annotation: _LegacyRandom',
        title='Random State',
    )
    key_added: Any = Field(
        None,
        description='Specifies where the neighbors data is stored. If not specified, stored in default locations. Original type annotation: str | None',
        title='Key Added',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='Determines whether to return a copy instead of writing to adata. Original type annotation: bool',
        title='Copy',
    )

    _api_name = PrivateAttr(default='scanpy.preprocessing.neighbors')
    _products_original = PrivateAttr(default=['data.uns["neighbors"]', 'data.obsp["distances"]', 'data.obsp["connectivities"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyPreprocessingLogP(BaseAPI):
    """
    Logarithmize the data matrix. Computes :math:`X = \log(X + 1)`, where :math:`log` denotes the natural logarithm unless a different base is given.
    """

    
    data: Any = Field(
        ...,
        description='The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.',
        title='Data',
    )
    base: Any = Field(
        None,
        description='Base of the logarithm. Natural logarithm is used by default.',
        title='Base',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='If an :class:`~anndata.AnnData` is passed, determines whether a copy is returned.',
        title='Copy',
    )
    chunked: Any = Field(
        None,
        description='Process the data matrix in chunks, which will save memory. Applies only to :class:`~anndata.AnnData`.',
        title='Chunked',
    )
    chunk_size: Any = Field(
        None,
        description='`n_obs` of the chunks to process the data in.',
        title='Chunk Size',
    )
    layer: Any = Field(None, description='Entry of layers to transform.', title='Layer')
    obsm: Any = Field(None, description='Entry of obsm to transform.', title='Obsm')

    _api_name = PrivateAttr(default='scanpy.preprocessing.log1p')
    _products_original = PrivateAttr(default=['data.X'])
    _data_name = PrivateAttr(default='data')



class ScanpyPreprocessingHighlyVariableGenes(BaseAPI):
    """
    Annotate highly variable genes based on different flavors such as 'seurat', 'cell_ranger', 'seurat_v3', or 'seurat_v3_paper'. The method involves selecting highly variable genes based on dispersion or normalized variance, with variations depending on the flavor and batch_key parameter. Additional flavors and comparisons to Seurat's naming are also provided.
    """

    
    adata: Any = Field(
        ...,
        description='The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.',
        title='Adata',
    )
    layer: Any = Field(
        None,
        description='If provided, use `adata.layers[layer]` for expression values instead of `adata.X`.',
        title='Layer',
    )
    n_top_genes: Any = Field(
        None,
        description="Number of highly-variable genes to keep. Mandatory if `flavor='seurat_v3'.",
        title='N Top Genes',
    )
    min_disp: Optional[Any] = Field(
        0.5,
        description="Cutoff for the minimum dispersion, ignored if `flavor='seurat_v3'`.",
        title='Min Disp',
    )
    max_disp: Optional[Any] = Field(
        'inf',
        description="Cutoff for the maximum dispersion, ignored if `flavor='seurat_v3'`.",
        title='Max Disp',
    )
    min_mean: Optional[Any] = Field(
        0.0125,
        description="Cutoff for the minimum mean, ignored if `flavor='seurat_v3'`.",
        title='Min Mean',
    )
    max_mean: Optional[Any] = Field(
        3,
        description="Cutoff for the maximum mean, ignored if `flavor='seurat_v3'`.",
        title='Max Mean',
    )
    span: Optional[Any] = Field(
        0.3,
        description="Fraction of the data used when estimating variance in the loess model fit if `flavor='seurat_v3'`.",
        title='Span',
    )
    n_bins: Optional[Any] = Field(
        20,
        description='Number of bins for binning the mean gene expression for normalization.',
        title='N Bins',
    )
    flavor: Optional[Any] = Field(
        'seurat',
        description='Choose the method for identifying highly variable genes.',
        title='Flavor',
    )
    subset: Optional[Any] = Field(
        False,
        description='Subset to highly-variable genes if `True`, otherwise merely indicate highly variable genes.',
        title='Subset',
    )
    inplace: Optional[Any] = Field(
        True,
        description='Whether to place calculated metrics in `.var` or return them.',
        title='Inplace',
    )
    batch_key: Any = Field(
        None,
        description='Key to select highly-variable genes within each batch separately and merge them.',
        title='Batch Key',
    )
    check_values: Optional[Any] = Field(
        True,
        description='Check if counts in selected layer are integers, used for specific flavors.',
        title='Check Values',
    )

    _api_name = PrivateAttr(default='scanpy.preprocessing.highly_variable_genes')
    _products_original = PrivateAttr(default=['data.var["highly_variable"]', 'data.var["means"]', 'data.var["dispersions"]', 'data.var["dispersions_norm"]', 'data.var["variances"]', 'data.var["variances_norm"]', 'data.var["highly_variable_rank"]', 'data.var["highly_variable_nbatches"]', 'data.var["highly_variable_intersection"]'])
    _data_name = PrivateAttr(default='adata')



class ScanpyPreprocessingPca(BaseAPI):
    """
    Principal component analysis: Computes PCA coordinates, loadings, and variance decomposition using the implementation of scikit-learn. In previous versions, computing PCA on a sparse matrix would make a dense copy of the array for mean centering. As of scanpy 1.5.0, mean centering is implicit, resulting in extremely similar but not exactly the same results. To reproduce the old results, pass a dense array.
    """

    
    data: Any = Field(
        ...,
        description='The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.',
        title='Data',
    )
    n_comps: Any = Field(
        None,
        description='Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation.',
        title='N Comps',
    )
    layer: Any = Field(
        None, description='Layer of `adata` to use as expression values.', title='Layer'
    )
    zero_center: Optional[Any] = Field(
        True,
        description='If `True`, compute standard PCA from covariance matrix. If `False`, omit zero-centering variables.',
        title='Zero Center',
    )
    svd_solver: Any = Field(
        None,
        description="SVD solver to use based on the type of input data. Options include 'arpack', 'covariance_eigh', 'randomized', and 'auto'.",
        title='Svd Solver',
    )
    random_state: Optional[Any] = Field(
        0,
        description='Change to use different initial states for the optimization.',
        title='Random State',
    )
    return_info: Optional[Any] = Field(
        False,
        description="Only relevant when not passing an AnnData: see 'Returns'.",
        title='Return Info',
    )
    mask_var: Optional[Any] = Field(
        0,
        description='To run only on a certain set of genes given by a boolean array or a string referring to an array in AnnData.var.',
        title='Mask Var',
    )
    use_highly_variable: Any = Field(
        None,
        description="Whether to use highly variable genes only, stored in .var['highly_variable'].",
        title='Use Highly Variable',
    )
    dtype: Optional[Any] = Field(
        'float32',
        description='Numpy data type string to which to convert the result.',
        title='Dtype',
    )
    chunked: Optional[Any] = Field(
        False,
        description='If `True`, perform an incremental PCA on segments of `chunk_size`. If `False`, perform a full PCA.',
        title='Chunked',
    )
    chunk_size: Any = Field(
        None,
        description='Number of observations to include in each chunk. Required if `chunked=True` was passed.',
        title='Chunk Size',
    )
    key_added: Any = Field(
        None,
        description="If not specified, the embedding is stored in AnnData.obsm['X_pca'], the loadings in AnnData.varm['PCs'], and the parameters in AnnData.uns['pca'].",
        title='Key Added',
    )
    copy_: Optional[Any] = Field(
        False,
        alias='copy',
        description='If an AnnData is passed, determines whether a copy is returned.',
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