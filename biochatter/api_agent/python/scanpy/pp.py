from __future__ import annotations

from typing import Any, Optional
from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScPpNeighbors(BaseAPI):
    """
    Compute the nearest neighbors distance matrix and a neighborhood graph of observations. The efficiency heavily relies on UMAP, which also provides a method for estimating connectivities of data points. Connectivities are computed based on different methods such as UMAP, Gaussian, or adaptations of specific references.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    n_neighbors: Optional[Any] = Field(
        15,
        description="The size of local neighborhood used for manifold approximation, with options for knn or Gaussian kernel. Original type annotation: int",
        title="N Neighbors",
    )
    n_pcs: Any = Field(
        None,
        description="Specifies the number of PCs to use, with an option for automatic selection. Original type annotation: int | None",
        title="N Pcs",
    )
    use_rep: Any = Field(
        None,
        description="Indicates the type of representation to use, with options for automatic selection based on data size. Original type annotation: str | None",
        title="Use Rep",
    )
    knn: Optional[Any] = Field(
        True,
        description="Determines whether to use a hard threshold or Gaussian kernel for neighbor assignment. Original type annotation: bool",
        title="Knn",
    )
    method: Optional[Any] = Field(
        "umap",
        description="Specifies the method for computing connectivities, with options for 'umap' or 'gauss'. Original type annotation: _Method",
        title="Method",
    )
    transformer: Any = Field(
        None,
        description="Defines the kNN search implementation, with options for various transformers. Original type annotation: KnnTransformerLike | _KnownTransformer | None",
        title="Transformer",
    )
    metric: Optional[Any] = Field(
        "euclidean",
        description="Defines the metric used for distance calculation, with options for a known metric's name or a callable function. Original type annotation: _Metric | _MetricFn",
        title="Metric",
    )
    metric_kwds: Optional[Any] = Field(
        {},
        description="Options for the metric, to be ignored if a transformer is used. Original type annotation: Mapping[str, Any]",
        title="Metric Kwds",
    )
    random_state: Optional[Any] = Field(
        0,
        description="Specifies a numpy random seed, to be ignored if a transformer is used. Original type annotation: _LegacyRandom",
        title="Random State",
    )
    key_added: Any = Field(
        None,
        description="Specifies where the neighbors data is stored within the data structure. Original type annotation: str | None",
        title="Key Added",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Determines whether to return a copy of the data instead of modifying the original data. Original type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.pp.neighbors")
    _products_original: list[str] = PrivateAttr(
        default=[
            'data.uns["neighbors"]',
            'data.obsp["distances"]',
            'data.obsp["connectivities"]',
        ]
    )
    _data_name: str = PrivateAttr(default="adata")


class ScPpLogP(BaseAPI):
    """
    Logarithmize the data matrix. Computes :math:`X = \\log(X + 1)`, where :math:`log` denotes the natural logarithm unless a different base is given.
    """

    data: Any = Field(
        ...,
        description="The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.",
        title="Data",
    )
    base: Any = Field(
        None,
        description="Base of the logarithm. Natural logarithm is used by default.",
        title="Base",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="If an :class:`~anndata.AnnData` is passed, determines whether a copy is returned.",
        title="Copy",
    )
    chunked: Any = Field(
        None,
        description="Process the data matrix in chunks, which will save memory. Applies only to :class:`~anndata.AnnData`.",
        title="Chunked",
    )
    chunk_size: Any = Field(
        None,
        description="`n_obs` of the chunks to process the data in.",
        title="Chunk Size",
    )
    layer: Any = Field(None, description="Entry of layers to transform.", title="Layer")
    obsm: Any = Field(None, description="Entry of obsm to transform.", title="Obsm")
    _api_name: str = PrivateAttr(default="sc.pp.log1p")
    _products_original: list[str] = PrivateAttr(default=["data.X"])
    _data_name: str = PrivateAttr(default="data")


class ScPpHighlyVariableGenes(BaseAPI):
    """
    Annotate highly variable genes by reproducing R-implementations of Seurat, Cell Ranger, and Seurat v3, depending on the chosen flavor. Various normalization and selection methods are applied to identify highly variable genes, with differences in handling batch effects and ranking depending on the flavor chosen. Additional flavors and comparisons to Seurat's naming conventions are also provided.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.",
        title="Adata",
    )
    layer: Any = Field(
        None,
        description="If provided, use `adata.layers[layer]` for expression values instead of `adata.X`.",
        title="Layer",
    )
    n_top_genes: Any = Field(
        None,
        description="Number of highly-variable genes to keep. Mandatory if `flavor='seurat_v3'.",
        title="N Top Genes",
    )
    min_disp: Optional[Any] = Field(
        0.5,
        description="Cutoff for the minimum dispersion, ignored if `flavor='seurat_v3'`.",
        title="Min Disp",
    )
    max_disp: Optional[Any] = Field(
        "inf",
        description="Cutoff for the maximum dispersion, ignored if `flavor='seurat_v3'`.",
        title="Max Disp",
    )
    min_mean: Optional[Any] = Field(
        0.0125,
        description="Cutoff for the minimum mean, ignored if `flavor='seurat_v3'`.",
        title="Min Mean",
    )
    max_mean: Optional[Any] = Field(
        3,
        description="Cutoff for the maximum mean, ignored if `flavor='seurat_v3'`.",
        title="Max Mean",
    )
    span: Optional[Any] = Field(
        0.3,
        description="The fraction of data used when estimating variance in the loess model if `flavor='seurat_v3'`.",
        title="Span",
    )
    n_bins: Optional[Any] = Field(
        20,
        description="Number of bins for binning the mean gene expression for normalization.",
        title="N Bins",
    )
    flavor: Optional[Any] = Field(
        "seurat",
        description="Choose the method for identifying highly variable genes.",
        title="Flavor",
    )
    subset: Optional[Any] = Field(
        False,
        description="Subset to highly-variable genes if `True`, otherwise indicate highly variable genes.",
        title="Subset",
    )
    inplace: Optional[Any] = Field(
        True,
        description="Whether to place calculated metrics in `.var` or return them.",
        title="Inplace",
    )
    batch_key: Any = Field(
        None,
        description="Key for selecting highly-variable genes within each batch separately and merging them.",
        title="Batch Key",
    )
    check_values: Optional[Any] = Field(
        True,
        description="Check if counts in selected layer are integers. Only used for specific flavors.",
        title="Check Values",
    )
    _api_name: str = PrivateAttr(default="sc.pp.highly_variable_genes")
    _products_original: list[str] = PrivateAttr(
        default=[
            'data.var["highly_variable"]',
            'data.var["means"]',
            'data.var["dispersions"]',
            'data.var["dispersions_norm"]',
            'data.var["variances"]',
            'data.var["variances_norm"]',
            'data.var["highly_variable_rank"]',
            'data.var["highly_variable_nbatches"]',
            'data.var["highly_variable_intersection"]',
        ]
    )
    _data_name: str = PrivateAttr(default="adata")


class ScPpPca(BaseAPI):
    """
    Principal component analysis: Computes PCA coordinates, loadings, and variance decomposition using the implementation of scikit-learn. Note the version change in 1.5.0 regarding mean centering in sparse matrices.
    """

    data: Any = Field(
        ...,
        description="The (annotated) data matrix of shape `n_obs` × `n_vars`, where rows are cells and columns are genes.",
        title="Data",
    )
    n_comps: Any = Field(
        None,
        description="Number of principal components to compute, with a default of 50 or 1 minus the minimum dimension size of the selected representation.",
        title="N Comps",
    )
    layer: Any = Field(
        None, description="Layer of `adata` to use as expression values.", title="Layer"
    )
    zero_center: Optional[Any] = Field(
        True,
        description="Determines whether to compute standard PCA from a covariance matrix. If `True`, standard PCA is computed; if `False`, zero-centering variables are omitted.",
        title="Zero Center",
    )
    svd_solver: Any = Field(
        None,
        description="Specifies the SVD solver to use based on different conditions and types of data input.",
        title="Svd Solver",
    )
    random_state: Optional[Any] = Field(
        0,
        description="Parameter used to change initial states for optimization.",
        title="Random State",
    )
    return_info: Optional[Any] = Field(
        False,
        description="Relevant when not passing an `AnnData` object, indicating details in 'Returns'.",
        title="Return Info",
    )
    mask_var: Optional[Any] = Field(
        0,
        description="Specifies a certain set of genes to run the computation on, with options for boolean array, string reference, or default behavior.",
        title="Mask Var",
    )
    use_highly_variable: Any = Field(
        None,
        description="Indicates whether to use highly variable genes only, with a deprecation notice to use `mask_var` instead.",
        title="Use Highly Variable",
    )
    dtype: Optional[Any] = Field(
        "float32",
        description="Specifies the Numpy data type string to convert the result into.",
        title="Dtype",
    )
    chunked: Optional[Any] = Field(
        False,
        description="Determines whether to perform incremental PCA on segments of a specified chunk size, affecting zero-centering and solver settings.",
        title="Chunked",
    )
    chunk_size: Any = Field(
        None,
        description="Number of observations to include in each chunk when `chunked=True`.",
        title="Chunk Size",
    )
    key_added: Any = Field(
        None,
        description="Specifies where the computed embedding, loadings, and parameters are stored within an `AnnData` object.",
        title="Key Added",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Determines whether a copy of the result is returned when an `AnnData` object is passed.",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.pp.pca")
    _products_original: list[str] = PrivateAttr(
        default=[
            'data.obsm["X_pca"]',
            'data.varm["PCs"]',
            'data.uns["pca"]["variance_ratio"]',
            'data.uns["pca"]["variance"]',
        ]
    )
    _data_name: str = PrivateAttr(default="data")


class ScPpCalculateQcMetrics(BaseAPI):
    """
    Calculate quality control metrics. Calculates a number of qc metrics for an AnnData object, largely based on calculateQCMetrics from scater. Most efficient on a sparse CSR or dense matrix. Note that this method can take a while to compile on the first call and is cached to disk for later use.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    expr_type: Optional[Any] = Field(
        "counts",
        description="Name of kind of values in X. Original type annotation: str",
        title="Expr Type",
    )
    var_type: Optional[Any] = Field(
        "genes",
        description="The kind of thing the variables are. Original type annotation: str",
        title="Var Type",
    )
    qc_vars: Optional[Any] = Field(
        [],
        description="Keys for boolean columns of .var which identify variables you could want to control for (e.g. 'ERCC' or 'mito'). Original type annotation: Collection[str] | str",
        title="Qc Vars",
    )
    percent_top: Optional[Any] = Field(
        [50, 100, 200, 500],
        description="List of ranks at which the cumulative proportion of expression will be reported as a percentage. Original type annotation: Collection[int] | None",
        title="Percent Top",
    )
    layer: Any = Field(
        None,
        description="If provided, use adata.layers[layer] for expression values instead of adata.X. Original type annotation: str | None",
        title="Layer",
    )
    use_raw: Optional[Any] = Field(
        False,
        description="If True, use adata.raw.X for expression values instead of adata.X. Original type annotation: bool",
        title="Use Raw",
    )
    inplace: Optional[Any] = Field(
        False,
        description="Whether to place calculated metrics in adata's .obs and .var. Original type annotation: bool",
        title="Inplace",
    )
    log1p: Optional[Any] = Field(
        True,
        description="Set to False to skip computing log1p transformed annotations. Original type annotation: bool",
        title="Log1P",
    )
    parallel: Any = Field(
        None,
        description="No description available. Original type annotation: bool | None",
        title="Parallel",
    )
    _api_name: str = PrivateAttr(default="sc.pp.calculate_qc_metrics")
    _products_original: list[str] = PrivateAttr(
        default=[
            'data.obs["total_genes_by_counts"]',
            'data.obs["total_counts"]',
            'data.obs["pct_counts_in_top_50_genes"]',
            'data.obs["pct_counts_in_top_100_genes"]',
            'data.obs["pct_counts_in_top_200_genes"]',
        ]
    )
    _data_name: str = PrivateAttr(default="adata")


class ScPpFilterCells(BaseAPI):
    """
    Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_genes` genes expressed. This is to filter measurement outliers,
    i.e. “unreliable” observations.

    Only provide one of the optional parameters `min_counts`, `min_genes`,
    `max_counts`, `max_genes` per call.
    """

    data: Any = Field(
        ...,
        description="The (annotated) data matrix of shape `n_obs` × `n_vars`.\nRows correspond to cells and columns to genes.\nOriginal type annotation: AnnData | _CSMatrix | np.ndarray | DaskArray",
        title="Data",
    )
    min_counts: Any = Field(
        None,
        description="Minimum number of counts required for a cell to pass filtering.\nOriginal type annotation: int | None",
        title="Min Counts",
    )
    min_genes: Any = Field(
        None,
        description="Minimum number of genes expressed required for a cell to pass filtering.\nOriginal type annotation: int | None",
        title="Min Genes",
    )
    max_counts: Any = Field(
        None,
        description="Maximum number of counts required for a cell to pass filtering.\nOriginal type annotation: int | None",
        title="Max Counts",
    )
    max_genes: Any = Field(
        None,
        description="Maximum number of genes expressed required for a cell to pass filtering.\nOriginal type annotation: int | None",
        title="Max Genes",
    )
    inplace: Optional[Any] = Field(
        True,
        description="Perform computation inplace or return result.\nOriginal type annotation: bool",
        title="Inplace",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="No description available.\nOriginal type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.pp.filter_cells")
    _products_original: list[str] = PrivateAttr(
        default=["data.X", 'data.obs["n_counts"]', 'data.obs["n_genes"]']
    )
    _data_name: str = PrivateAttr(default="data")


class ScPpFilterGenes(BaseAPI):
    """
    Filter genes based on number of cells or counts. Keep genes that have at least `min_counts` counts or are expressed in at least `min_cells` cells or have at most `max_counts` counts or are expressed in at most `max_cells` cells. Only provide one of the optional parameters `min_counts`, `min_cells`, `max_counts`, `max_cells` per call.
    """

    data: Any = Field(
        ...,
        description="An annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes. Original type annotation: AnnData | _CSMatrix | np.ndarray | DaskArray",
        title="Data",
    )
    min_counts: Any = Field(
        None,
        description="Minimum number of counts required for a gene to pass filtering. Original type annotation: int | None",
        title="Min Counts",
    )
    min_cells: Any = Field(
        None,
        description="Minimum number of cells expressed required for a gene to pass filtering. Original type annotation: int | None",
        title="Min Cells",
    )
    max_counts: Any = Field(
        None,
        description="Maximum number of counts required for a gene to pass filtering. Original type annotation: int | None",
        title="Max Counts",
    )
    max_cells: Any = Field(
        None,
        description="Maximum number of cells expressed required for a gene to pass filtering. Original type annotation: int | None",
        title="Max Cells",
    )
    inplace: Optional[Any] = Field(
        True,
        description="Perform computation inplace or return result. Original type annotation: bool",
        title="Inplace",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="No description available. Original type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.pp.filter_genes")
    _products_original: list[str] = PrivateAttr(
        default=["data.X", 'data.var["n_counts"]', 'data.var["n_genes"]']
    )
    _data_name: str = PrivateAttr(default="data")


class ScPpNormalizeTotal(BaseAPI):
    """
    Normalize counts per cell. Normalize each cell by total counts over all genes, so that every cell has the same total count after normalization. If choosing `target_sum=1e6`, this is CPM normalization. If `exclude_highly_expressed=True`, very highly expressed genes are excluded from the computation of the normalization factor (size factor) for each cell. This is meaningful as these can strongly influence the resulting normalized values for all other genes. Similar functions are used by Seurat, Cell Ranger, or SPRING. Note: When used with a Dask array in adata.X, this function will have to call functions that trigger `.compute()` on the Dask array if `exclude_highly_expressed` is `True`, `layer_norm` is not `None`, or if `key_added` is not `None`.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.",
        title="Adata",
    )
    target_sum: Any = Field(
        None,
        description="If `None`, after normalization, each observation (cell) has a total count equal to the median of total counts for observations (cells) before normalization.",
        title="Target Sum",
    )
    exclude_highly_expressed: Optional[Any] = Field(
        False,
        description="Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell based on a specified maximum fraction of total counts in a cell.",
        title="Exclude Highly Expressed",
    )
    max_fraction: Optional[Any] = Field(
        0.05,
        description="Specifies the threshold fraction of counts in a cell above which genes are considered highly expressed.",
        title="Max Fraction",
    )
    key_added: Any = Field(
        None,
        description="Name of the field in `adata.obs` where the normalization factor is stored.",
        title="Key Added",
    )
    layer: Any = Field(
        None,
        description="Specifies the layer to normalize instead of `X`. If `None`, `X` is normalized.",
        title="Layer",
    )
    layers: Any = Field(
        None,
        description="Specifies the available layers for normalization.",
        title="Layers",
    )
    layer_norm: Any = Field(
        None,
        description="Specifies the type of normalization for the layers.",
        title="Layer Norm",
    )
    inplace: Optional[Any] = Field(
        True,
        description="Determines whether to update `adata` in place or return a dictionary with normalized copies of `adata.X` and `adata.layers`.",
        title="Inplace",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Determines whether to modify the copied input object. Not compatible with `inplace=False`.",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.pp.normalize_total")
    _products_original: list[str] = PrivateAttr(default=["data.X"])
    _data_name: str = PrivateAttr(default="adata")


class ScPpRegressOut(BaseAPI):
    """
    Regress out (mostly) unwanted sources of variation. Uses simple linear regression. This is inspired by Seurat's `regressOut` function in R (Satija et al., 2015). Note that this function tends to overcorrect in certain circumstances as described in issue #526.
    """

    adata: Any = Field(..., description="", title="Adata")
    keys: Any = Field(..., description="", title="Keys")
    layer: Any = Field(None, description="", title="Layer")
    n_jobs: Any = Field(None, description="", title="N Jobs")
    copy_: Optional[Any] = Field(False, alias="copy", description="", title="Copy")
    _api_name: str = PrivateAttr(default="sc.pp.regress_out")
    _products_original: list[str] = PrivateAttr(default=["data.X"])
    _data_name: str = PrivateAttr(default="adata")


class ScPpScale(BaseAPI):
    """
    Scale data to unit variance and zero mean. Variables that do not display any variation are retained and set to 0 during this operation, with the possibility of being set to NaNs in the future.
    """

    data: Any = Field(
        ...,
        description="The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.",
        title="Data",
    )
    zero_center: Optional[Any] = Field(
        True,
        description="If `False`, omit zero-centering variables, which allows to handle sparse input efficiently.",
        title="Zero Center",
    )
    max_value: Any = Field(
        None,
        description="Clip (truncate) to this value after scaling. If `None`, do not clip.",
        title="Max Value",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Whether this function should be performed inplace. If an AnnData object is passed, this also determines if a copy is returned.",
        title="Copy",
    )
    layer: Any = Field(
        None,
        description="If provided, specifies which element of layers to scale.",
        title="Layer",
    )
    obsm: Any = Field(
        None,
        description="If provided, specifies which element of obsm to scale.",
        title="Obsm",
    )
    mask_obs: Any = Field(
        None,
        description="Restrict both the derivation of scaling parameters and the scaling itself to a certain set of observations. The mask is specified as a boolean array or a string referring to an array in AnnData.obs. This will transform data from csc to csr format if `issparse(data)`.",
        title="Mask Obs",
    )
    _api_name: str = PrivateAttr(default="sc.pp.scale")
    _products_original: list[str] = PrivateAttr(
        default=["data.X", 'data.var["mean"]', 'data.var["std"]', 'data.var["var"]']
    )
    _data_name: str = PrivateAttr(default="data")


class ScPpSample(BaseAPI):
    """
    Sample observations or variables with or without replacement.
    """

    data: Any = Field(
        ...,
        description="The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.",
        title="Data",
    )
    fraction: Any = Field(
        None,
        description="Sample to this `fraction` of the number of observations or variables. This can be larger than 1.0 if `replace=True`. See `axis` and `replace`.",
        title="Fraction",
    )
    n: Any = Field(
        None,
        description="Sample to this number of observations or variables. See `axis`.",
        title="N",
    )
    rng: Any = Field(
        None, description="Random seed to change subsampling.", title="Rng"
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="If an :class:`~anndata.AnnData` is passed, determines whether a copy is returned.",
        title="Copy",
    )
    replace: Optional[Any] = Field(
        False,
        description="If True, samples are drawn with replacement.",
        title="Replace",
    )
    axis: Optional[Any] = Field(
        "obs",
        description="Sample `obs`ervations (axis 0) or `var`iables (axis 1).",
        title="Axis",
    )
    p: Any = Field(
        None,
        description="Drawing probabilities (floats) or mask (bools). Either an `axis`-sized array or the name of a column. If `p` is an array of probabilities, it must sum to 1.",
        title="P",
    )
    _api_name: str = PrivateAttr(default="sc.pp.sample")
    _products_original: list[str] = PrivateAttr(default=["data.X"])
    _data_name: str = PrivateAttr(default="data")


class ScPpDownsampleCounts(BaseAPI):
    """
    Downsample counts from count matrix. If `counts_per_cell` is specified, each cell will be downsampled. If `total_counts` is specified, the expression matrix will be downsampled to contain at most `total_counts`.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    counts_per_cell: Any = Field(
        None,
        description="Target total counts per cell. If a cell has more than 'counts_per_cell', it will be downsampled to this number. Resulting counts can be specified on a per cell basis by passing an array. Should be an integer or integer ndarray with the same length as the number of obs. Original type annotation: int | Collection[int] | None",
        title="Counts Per Cell",
    )
    total_counts: Any = Field(
        None,
        description="Target total counts. If the count matrix has more than `total_counts` it will be downsampled to have this number. Original type annotation: int | None",
        title="Total Counts",
    )
    random_state: Optional[Any] = Field(
        0,
        description="Random seed for subsampling. Original type annotation: _LegacyRandom",
        title="Random State",
    )
    replace: Optional[Any] = Field(
        False,
        description="Whether to sample the counts with replacement. Original type annotation: bool",
        title="Replace",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Determines whether a copy of `adata` is returned. Original type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.pp.downsample_counts")
    _products_original: list[str] = PrivateAttr(default=["data.X"])
    _data_name: str = PrivateAttr(default="adata")


class ScPpCombat(BaseAPI):
    """
    ComBat function for batch effect correction. Corrects for batch effects by fitting linear models, gains statistical power via an EB framework where information is borrowed across genes. Implementation available at: https://github.com/brentp/combat.py
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix with original type annotation as AnnData.",
        title="Adata",
    )
    key: Optional[Any] = Field(
        "batch",
        description="Categorical annotation key from AnnData.obs for batch effect removal with original type annotation as str.",
        title="Key",
    )
    covariates: Any = Field(
        None,
        description="Additional covariates for batch effect removal including adjustment variables or biological conditions, referred to as design matrix X in Equation 2.1 in Johnson2006 and mod argument in combat function in sva R package. Missing covariates may introduce bias or remove biological signal in unbalanced designs with original type annotation as Collection[str] or None.",
        title="Covariates",
    )
    inplace: Optional[Any] = Field(
        True,
        description="Indicator for replacing adata.X or returning corrected data with original type annotation as bool.",
        title="Inplace",
    )
    _api_name: str = PrivateAttr(default="sc.pp.combat")
    _products_original: list[str] = PrivateAttr(default=["data.X"])
    _data_name: str = PrivateAttr(default="adata")


class ScPpScrublet(BaseAPI):
    """
    Predict doublets using Scrublet. Predict cell doublets using a nearest-neighbor classifier of observed transcriptomes and simulated doublets. Works best with raw counts matrix from a single sample or similar samples from the same experiment. This function is a wrapper around functions that pre-process using Scanpy and directly call functions of Scrublet. You may also undertake your own preprocessing, simulate doublets with 'scanpy.pp.scrublet_simulate_doublets', and run the core scrublet function 'scanpy.pp.scrublet' with adata_sim set.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix of cells and genes, expected to be un-normalized unless adata_sim is provided for consistent processing.",
        title="Adata",
    )
    adata_sim: Any = Field(
        None,
        description="An optional annData object generated for advanced use, with the same number of genes as adata, typically built from adata_obs after filtering and selecting highly variable genes.",
        title="Adata Sim",
    )
    batch_key: Any = Field(
        None,
        description="An optional column name in AnnData.obs used to discriminate between batches.",
        title="Batch Key",
    )
    sim_doublet_ratio: Optional[Any] = Field(
        2.0,
        description="The number of doublets to simulate relative to the number of observed transcriptomes.",
        title="Sim Doublet Ratio",
    )
    expected_doublet_rate: Optional[Any] = Field(
        0.05,
        description="The estimated doublet rate for the experiment if adata_sim is not provided.",
        title="Expected Doublet Rate",
    )
    stdev_doublet_rate: Optional[Any] = Field(
        0.02,
        description="The uncertainty in the expected doublet rate if adata_sim is not provided.",
        title="Stdev Doublet Rate",
    )
    synthetic_doublet_umi_subsampling: Optional[Any] = Field(
        1.0,
        description="The rate for sampling UMIs when creating synthetic doublets if adata_sim is not provided.",
        title="Synthetic Doublet Umi Subsampling",
    )
    knn_dist_metric: Optional[Any] = Field(
        "euclidean",
        description="The distance metric used for finding nearest neighbors.",
        title="Knn Dist Metric",
    )
    normalize_variance: Optional[Any] = Field(
        True,
        description="A boolean indicating whether to normalize the data such that each gene has a variance of 1.",
        title="Normalize Variance",
    )
    log_transform: Optional[Any] = Field(
        False,
        description="A boolean indicating whether to log-transform the data prior to PCA.",
        title="Log Transform",
    )
    mean_center: Optional[Any] = Field(
        True,
        description="A boolean indicating whether to center the data such that each gene has a mean of 0.",
        title="Mean Center",
    )
    n_prin_comps: Optional[Any] = Field(
        30,
        description="The number of principal components used to embed the transcriptomes before k-nearest-neighbor graph construction.",
        title="N Prin Comps",
    )
    use_approx_neighbors: Any = Field(
        None,
        description="A boolean indicating whether to use approximate nearest neighbor method for the KNN classifier.",
        title="Use Approx Neighbors",
    )
    get_doublet_neighbor_parents: Optional[Any] = Field(
        False,
        description="A boolean indicating whether to return the parent transcriptomes that generated the doublet neighbors of each observed transcriptome.",
        title="Get Doublet Neighbor Parents",
    )
    n_neighbors: Any = Field(
        None,
        description="The number of neighbors used to construct the KNN graph of observed transcriptomes and simulated doublets.",
        title="N Neighbors",
    )
    threshold: Any = Field(
        None,
        description="The doublet score threshold for calling a transcriptome a doublet.",
        title="Threshold",
    )
    verbose: Optional[Any] = Field(
        True,
        description="A boolean indicating whether to log progress updates.",
        title="Verbose",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="A boolean indicating whether to return a copy of the input adata with Scrublet results added.",
        title="Copy",
    )
    random_state: Optional[Any] = Field(
        0,
        description="The initial state for doublet simulation and nearest neighbors.",
        title="Random State",
    )
    _api_name: str = PrivateAttr(default="sc.pp.scrublet")
    _products_original: list[str] = PrivateAttr(
        default=[
            'data.obs["doublet_score"]',
            'data.obs["predicted_doublet"]',
            'data.uns["scrublet"]["doublet_scores_sim"]',
            'data.uns["scrublet"]["doublet_parents"]',
            'data.uns["scrublet"]["parameters"]',
        ]
    )
    _data_name: str = PrivateAttr(default="adata")


class ScPpScrubletSimulateDoublets(BaseAPI):
    """
    Simulate doublets by adding the counts of random observed transcriptome pairs.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix of shape ``n_obs`` × ``n_vars``. Rows correspond to cells and columns to genes. Genes should have been filtered for expression and variability, and the object should contain raw expression of the same dimensions.",
        title="Adata",
    )
    layer: Any = Field(
        None,
        description="Layer of adata where raw values are stored, or 'X' if values are in .X.",
        title="Layer",
    )
    sim_doublet_ratio: Optional[Any] = Field(
        2.0,
        description="Number of doublets to simulate relative to the number of observed transcriptomes. If `None`, self.sim_doublet_ratio is used.",
        title="Sim Doublet Ratio",
    )
    synthetic_doublet_umi_subsampling: Optional[Any] = Field(
        1.0,
        description="Rate for sampling UMIs when creating synthetic doublets. If 1.0, each doublet is created by simply adding the UMIs from two randomly sampled observed transcriptomes. For values less than 1, the UMI counts are added and then randomly sampled at the specified rate.",
        title="Synthetic Doublet Umi Subsampling",
    )
    random_seed: Optional[Any] = Field(
        0, description="No description available.", title="Random Seed"
    )
    _api_name: str = PrivateAttr(default="sc.pp.scrublet_simulate_doublets")
    _products_original: list[str] = PrivateAttr(
        default=[
            'data.obsm["scrublet"]["doublet_parents"]',
            'data.uns["scrublet"]["parameters"]',
        ]
    )
    _data_name: str = PrivateAttr(default="adata")


TOOLS_DICT = {
    "sc.pp.neighbors": ScPpNeighbors,
    "sc.pp.log1p": ScPpLogP,
    "sc.pp.highly_variable_genes": ScPpHighlyVariableGenes,
    "sc.pp.pca": ScPpPca,
    "sc.pp.calculate_qc_metrics": ScPpCalculateQcMetrics,
    "sc.pp.filter_cells": ScPpFilterCells,
    "sc.pp.filter_genes": ScPpFilterGenes,
    "sc.pp.normalize_total": ScPpNormalizeTotal,
    "sc.pp.regress_out": ScPpRegressOut,
    "sc.pp.scale": ScPpScale,
    "sc.pp.sample": ScPpSample,
    "sc.pp.downsample_counts": ScPpDownsampleCounts,
    "sc.pp.combat": ScPpCombat,
    "sc.pp.scrublet": ScPpScrublet,
    "sc.pp.scrublet_simulate_doublets": ScPpScrubletSimulateDoublets,
}
