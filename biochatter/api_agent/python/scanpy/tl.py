from __future__ import annotations

from typing import Any, Optional
from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI


class ScTlPaga(BaseAPI):
    """
    Map out the coarse-grained connectivity structures of complex manifolds. Explains the concept of partition-based graph abstraction (PAGA) and how it simplifies the representation of single-cell graph connectivity. Discusses the interpretation of confidence in PAGA and its application for obtaining a simplified representation of manifold data.
    """

    adata: Any = Field(
        ...,
        description="An annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    groups: Any = Field(
        None,
        description="Key for categorical in `adata.obs`. You can pass your predefined groups by choosing any categorical annotation of observations. Default: The first present key of `'leiden'` or `'louvain'`. Original type annotation: str | None",
        title="Groups",
    )
    use_rna_velocity: Optional[Any] = Field(
        False,
        description="Use RNA velocity to orient edges in the abstracted graph and estimate transitions. Requires that `adata.uns` contains a directed single-cell graph with key `['velocity_graph']`. This feature might be subject to change in the future. Original type annotation: bool",
        title="Use Rna Velocity",
    )
    model: Optional[Any] = Field(
        "v1.2",
        description="The PAGA connectivity model. Original type annotation: Literal['v1.2', 'v1.0']",
        title="Model",
    )
    neighbors_key: Any = Field(
        None,
        description="If not specified, paga looks `.uns['neighbors']` for neighbors settings and `.obsp['connectivities']`, `.obsp['distances']` for connectivities and distances respectively (default storage places for `pp.neighbors`). If specified, paga looks `.uns[neighbors_key]` for neighbors settings and `.obsp[.uns[neighbors_key]['connectivities_key']]`, `.obsp[.uns[neighbors_key]['distances_key']]` for connectivities and distances respectively. Original type annotation: str | None",
        title="Neighbors Key",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Copy `adata` before computation and return a copy. Otherwise, perform computation inplace and return `None`. Original type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.tl.paga")
    _products_original: list[str] = PrivateAttr(
        default=[
            'data.uns["paga"]["connectivities"]',
            'data.uns["paga"]["connectivities_tree"]',
        ]
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlLeiden(BaseAPI):
    """
    Cluster cells into subgroups using the Leiden algorithm, an improved version of the Louvain algorithm, proposed for single-cell analysis. Prior execution of functions scanpy.pp.neighbors or scanpy.external.pp.bbknn is required.
    """

    adata: Any = Field(..., description="The annotated data matrix.", title="Adata")
    resolution: Optional[Any] = Field(
        1,
        description="A parameter value controlling the coarseness of the clustering.",
        title="Resolution",
    )
    restrict_to: Any = Field(
        None,
        description="Restrict the clustering to specific categories within the sample annotation.",
        title="Restrict To",
    )
    random_state: Optional[Any] = Field(
        0,
        description="Change the initialization of the optimization.",
        title="Random State",
    )
    key_added: Optional[Any] = Field(
        "leiden",
        description="Key under which to add the cluster labels.",
        title="Key Added",
    )
    adjacency: Any = Field(
        None,
        description="Sparse adjacency matrix of the graph, defaults to neighbors connectivities.",
        title="Adjacency",
    )
    directed: Any = Field(
        None,
        description="Indicates whether the graph is directed or undirected.",
        title="Directed",
    )
    use_weights: Optional[Any] = Field(
        True,
        description="Determines if edge weights from the graph are used in the computation.",
        title="Use Weights",
    )
    n_iterations: Optional[Any] = Field(
        -1,
        description="Number of iterations of the Leiden clustering algorithm to perform.",
        title="N Iterations",
    )
    partition_type: Any = Field(
        None,
        description="Type of partition to use in clustering.",
        title="Partition Type",
    )
    neighbors_key: Any = Field(
        None,
        description="Specifies the key to use for neighbors connectivities as adjacency.",
        title="Neighbors Key",
    )
    obsp: Any = Field(
        None, description="Specifies the adjacency to use in clustering.", title="Obsp"
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Whether to copy the data or modify it in place.",
        title="Copy",
    )
    flavor: Optional[Any] = Field(
        "leidenalg",
        description="Specifies which package's implementation to use.",
        title="Flavor",
    )
    clustering_args: Any = Field(
        ..., description="No description available.", title="Clustering Args"
    )
    _api_name: str = PrivateAttr(default="sc.tl.leiden")
    _products_original: list[str] = PrivateAttr(
        default=['data.obs["leiden"]', 'data.uns["leiden"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlLouvain(BaseAPI):
    """
    Cluster cells into subgroups using the Louvain algorithm with references to relevant papers and functions.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    resolution: Any = Field(
        None,
        description="For clustering, you can specify a resolution parameter that defaults to 1.0. Original type annotation: float | None",
        title="Resolution",
    )
    random_state: Optional[Any] = Field(
        0,
        description="Change the initialization of the optimization. Original type annotation: _LegacyRandom",
        title="Random State",
    )
    restrict_to: Any = Field(
        None,
        description="Restrict clustering to specific categories within the sample annotation. Original type annotation: tuple[str, Sequence[str]] | None",
        title="Restrict To",
    )
    key_added: Optional[Any] = Field(
        "louvain",
        description="Key to add the cluster labels under. Original type annotation: str",
        title="Key Added",
    )
    adjacency: Any = Field(
        None,
        description="Sparse adjacency matrix of the graph, defaults to neighbors connectivities. Original type annotation: _CSMatrix | None",
        title="Adjacency",
    )
    flavor: Optional[Any] = Field(
        "vtraag",
        description="Choose between different packages for computing clustering, including 'vtraag', 'igraph', and 'rapids'. Original type annotation: Literal['vtraag', 'igraph', 'rapids']",
        title="Flavor",
    )
    directed: Optional[Any] = Field(
        True,
        description="Indicate whether the adjacency matrix represents a directed graph. Original type annotation: bool",
        title="Directed",
    )
    use_weights: Optional[Any] = Field(
        False,
        description="Specify whether to use weights from knn graph. Original type annotation: bool",
        title="Use Weights",
    )
    partition_type: Any = Field(
        None,
        description="Type of partition to use, applicable only if the flavor is 'vtraag'. Original type annotation: type[MutableVertexPartition] | None",
        title="Partition Type",
    )
    partition_kwargs: Optional[Any] = Field(
        {},
        description="Keyword arguments to pass to partitioning, particularly for the 'vtraag' method. Original type annotation: Mapping[str, Any]",
        title="Partition Kwargs",
    )
    neighbors_key: Any = Field(
        None,
        description="Specify the key to use for neighbors connectivities as adjacency. Original type annotation: str | None",
        title="Neighbors Key",
    )
    obsp: Any = Field(
        None,
        description="Specify the adjacency using .obsp[obsp], cannot be used simultaneously with 'neighbors_key'. Original type annotation: str | None",
        title="Obsp",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Indicate whether to copy the data or modify it in place. Original type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.tl.louvain")
    _products_original: list[str] = PrivateAttr(
        default=['data.obs["louvain"]', 'data.uns["louvain"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlUmap(BaseAPI):
    """
    Embed the neighborhood graph using UMAP. UMAP is a technique for visualizing high-dimensional data that is faster than tSNE and optimizes the embedding to best reflect the data's topology, represented using a neighborhood graph in Scanpy. The implementation of UMAP is based on umap-learn by McInnes (2018). For comparisons with tSNE, refer to Becht (2018).
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    min_dist: Optional[Any] = Field(
        0.5,
        description="The effective minimum distance between embedded points, influencing clustering or dispersion. Default value in `umap-learn` is 0.1. Original type annotation: float",
        title="Min Dist",
    )
    spread: Optional[Any] = Field(
        1.0,
        description="The scale of embedded points affecting clustering. Original type annotation: float",
        title="Spread",
    )
    n_components: Optional[Any] = Field(
        2,
        description="Number of dimensions in the embedding. Original type annotation: int",
        title="N Components",
    )
    maxiter: Any = Field(
        None,
        description="Number of optimization iterations. Original type annotation: int | None",
        title="Maxiter",
    )
    alpha: Optional[Any] = Field(
        1.0,
        description="Initial learning rate for embedding optimization. Original type annotation: float",
        title="Alpha",
    )
    gamma: Optional[Any] = Field(
        1.0,
        description="Weighting for negative samples in low-dimensional embedding. Original type annotation: float",
        title="Gamma",
    )
    negative_sample_rate: Optional[Any] = Field(
        5,
        description="Number of negative samples per positive sample for optimization. Original type annotation: int",
        title="Negative Sample Rate",
    )
    init_pos: Optional[Any] = Field(
        "spectral",
        description="Initialization method for low-dimensional embedding. Options include 'paga', 'spectral', 'random'. Original type annotation: _InitPos | np.ndarray | None",
        title="Init Pos",
    )
    random_state: Optional[Any] = Field(
        0,
        description="Seed for random number generation. Original type annotation: _LegacyRandom",
        title="Random State",
    )
    a: Any = Field(
        None,
        description="Specific parameter for embedding control, automatically set if None. Original type annotation: float | None",
        title="A",
    )
    b: Any = Field(
        None,
        description="Specific parameter for embedding control, automatically set if None. Original type annotation: float | None",
        title="B",
    )
    method: Optional[Any] = Field(
        "umap",
        description="Chosen implementation between 'umap' and 'rapids'. Original type annotation: Literal['umap', 'rapids']",
        title="Method",
    )
    key_added: Any = Field(
        None,
        description="Specifies storage location for embedding and parameters. Original type annotation: str | None",
        title="Key Added",
    )
    neighbors_key: Optional[Any] = Field(
        "neighbors",
        description="Key to access neighbors settings and connectivities for Umap. Original type annotation: str",
        title="Neighbors Key",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Option to return a copy instead of modifying the input data. Original type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.tl.umap")
    _products_original: list[str] = PrivateAttr(
        default=['data.obsm["X_umap"]', 'data.uns["umap"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlTsne(BaseAPI):
    """
    t-SNE was proposed for visualizating single-cell data and can be implemented using scikit-learn or Multicore-tSNE.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    n_pcs: Any = Field(
        None,
        description="Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`. Original type annotation: int | None",
        title="N Pcs",
    )
    use_rep: Any = Field(
        None,
        description="Use the indicated representation. `'X'` or any key for `.obsm` is valid. If `None`, the representation is chosen automatically. Original type annotation: str | None",
        title="Use Rep",
    )
    perplexity: Optional[Any] = Field(
        30,
        description="The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Consider selecting a value between 5 and 50. Original type annotation: float",
        title="Perplexity",
    )
    metric: Optional[Any] = Field(
        "euclidean",
        description="Distance metric to calculate neighbors on. Original type annotation: str",
        title="Metric",
    )
    early_exaggeration: Optional[Any] = Field(
        12,
        description="Controls how tight natural clusters in the original space are in the embedded space. Original type annotation: float",
        title="Early Exaggeration",
    )
    learning_rate: Optional[Any] = Field(
        1000,
        description="The learning rate can be a critical parameter. It should be between 100 and 1000. Original type annotation: float",
        title="Learning Rate",
    )
    random_state: Optional[Any] = Field(
        0,
        description="Change this to use different initial states for the optimization. Original type annotation: _LegacyRandom",
        title="Random State",
    )
    use_fast_tsne: Optional[Any] = Field(
        False,
        description="No description available. Original type annotation: bool",
        title="Use Fast Tsne",
    )
    n_jobs: Any = Field(
        None,
        description="Number of jobs for parallel computation. Original type annotation: int | None",
        title="N Jobs",
    )
    key_added: Any = Field(
        None,
        description="Specifies where the embedding and parameters are stored. Original type annotation: str | None",
        title="Key Added",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Return a copy instead of writing to `adata`. Original type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.tl.tsne")
    _products_original: list[str] = PrivateAttr(
        default=['data.obsm["X_tsne"]', 'data.uns["tsne"]']
    )
    _data_name: str = PrivateAttr(default="adata")

# Jiahang(severe): not being correctly summarized.
class ScTlDiffmap(BaseAPI):
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
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix.\nOriginal type annotation: AnnData",
        title="Adata",
    )
    n_comps: Optional[Any] = Field(
        15,
        description="The number of dimensions of the representation.\nOriginal type annotation: int",
        title="N Comps",
    )
    neighbors_key: Any = Field(
        None,
        description="If not specified, diffmap looks in .uns['neighbors'] for neighbors settings\nand .obsp['connectivities'] and .obsp['distances'] for connectivities and\ndistances, respectively (default storage places for pp.neighbors).\nIf specified, diffmap looks in .uns[neighbors_key] for neighbors settings and\n.obsp[.uns[neighbors_key]['connectivities_key']] and\n.obsp[.uns[neighbors_key]['distances_key']] for connectivities and distances,\nrespectively.\nOriginal type annotation: str | None",
        title="Neighbors Key",
    )
    random_state: Optional[Any] = Field(
        0,
        description="A numpy random seed\nOriginal type annotation: _LegacyRandom",
        title="Random State",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Return a copy instead of writing to adata.\nOriginal type annotation: bool",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.tl.diffmap")
    _products_original: list[str] = PrivateAttr(
        default=['data.obsm["X_diffmap"]', 'data.uns["diffmap_evals"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlEmbeddingDensity(BaseAPI):
    """
    Calculate the density of cells in an embedding (per condition). Gaussian kernel density estimation is used to calculate the density of cells in an embedded space. This can be performed per category over a categorical cell annotation. The cell density can be plotted using the `pl.embedding_density` function. Note that density values are scaled to be between 0 and 1, making them comparable only within the same category. The KDE estimate used (`scipy.stats.gaussian_kde`) may become unreliable if there are not enough cells in a category. This function was written by Sophie Tritschler and implemented into Scanpy by Malte Luecken.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    basis: Optional[Any] = Field(
        "umap",
        description="The embedding over which the density will be calculated. This embedded representation is found in `adata.obsm['X_[basis]']`. Original type annotation: str",
        title="Basis",
    )
    groupby: Any = Field(
        None,
        description="Key for categorical observation/cell annotation for which densities are calculated per category. Original type annotation: str | None",
        title="Groupby",
    )
    key_added: Any = Field(
        None,
        description="Name of the `.obs` covariate that will be added with the density estimates. Original type annotation: str | None",
        title="Key Added",
    )
    components: Any = Field(
        None,
        description="The embedding dimensions over which the density should be calculated. This is limited to two components. Original type annotation: str | Sequence[str] | None",
        title="Components",
    )
    _api_name: str = PrivateAttr(default="sc.tl.embedding_density")
    _products_original: list[str] = PrivateAttr(
        default=['data.obs["umap_density"]', 'data.uns["umap_density_params"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlRankGenesGroups(BaseAPI):
    """
    Rank genes for characterizing groups. Expects logarithmized data.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix with original type annotation AnnData.",
        title="Adata",
    )
    groupby: Any = Field(
        ...,
        description="Key of the observations grouping to consider with original type annotation str.",
        title="Groupby",
    )
    mask_var: Any = Field(
        None,
        description="Subset of genes to use in statistical tests with original type annotation NDArray[np.bool_] | str | None.",
        title="Mask Var",
    )
    use_raw: Any = Field(
        None,
        description="Boolean flag to indicate whether to use raw attribute of adata if present.",
        title="Use Raw",
    )
    groups: Optional[Any] = Field(
        "all",
        description="Subset of groups for comparison with original type annotation Literal['all'] | Iterable[str].",
        title="Groups",
    )
    reference: Optional[Any] = Field(
        "rest",
        description="Identifier to specify which group to compare with in statistical tests.",
        title="Reference",
    )
    n_genes: Any = Field(
        None,
        description="Number of genes to include in the returned tables with default to all genes.",
        title="N Genes",
    )
    rankby_abs: Optional[Any] = Field(
        False,
        description="Flag to rank genes by the absolute value of the score rather than the score itself.",
        title="Rankby Abs",
    )
    pts: Optional[Any] = Field(
        False,
        description="Flag to compute the fraction of cells expressing the genes.",
        title="Pts",
    )
    key_added: Any = Field(
        None,
        description="Key in adata.uns where information is saved to.",
        title="Key Added",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Boolean flag to indicate whether to copy adata or modify it inplace.",
        title="Copy",
    )
    method: Any = Field(
        None,
        description="Method used for statistical testing with options like 't-test', 'wilcoxon', 'logreg', and more.",
        title="Method",
    )
    corr_method: Optional[Any] = Field(
        "benjamini-hochberg",
        description="Method for p-value correction used in certain statistical tests.",
        title="Corr Method",
    )
    tie_correct: Optional[Any] = Field(
        False,
        description="Flag to enable tie correction for 'wilcoxon' scores.",
        title="Tie Correct",
    )
    layer: Any = Field(
        None,
        description="Key from adata.layers used in performing statistical tests.",
        title="Layer",
    )
    kwds: Any = Field(
        ...,
        description="Parameters passed to test methods affecting certain sklearn.linear_model.LogisticRegression parameters.",
        title="Kwds",
    )
    _api_name: str = PrivateAttr(default="sc.tl.rank_genes_groups")
    _products_original: list[str] = PrivateAttr(
        default=['data.uns["rank_genes_groups"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlFilterRankGenesGroups(BaseAPI):
    """
    Filter out genes based on two criteria: log fold change and fraction of genes expressing the gene within and outside the 'groupby' categories. See the function scanpy.tl.rank_genes_groups. Results are stored in adata.uns[key_added] (default: 'rank_genes_groups_filtered'). To preserve the original structure of adata.uns['rank_genes_groups'], filtered genes are set to NaN.
    """

    adata: Any = Field(
        ...,
        description="AnnData is the original type annotation for this parameter.",
        title="Adata",
    )
    key: Any = Field(
        None,
        description="A string that may or may not be provided as an input.",
        title="Key",
    )
    groupby: Any = Field(
        None,
        description="A string that may or may not be provided as an input.",
        title="Groupby",
    )
    use_raw: Any = Field(
        None,
        description="A boolean value or None which determines whether to use raw data or not.",
        title="Use Raw",
    )
    key_added: Optional[Any] = Field(
        "rank_genes_groups_filtered",
        description="A string indicating the key to be added.",
        title="Key Added",
    )
    min_in_group_fraction: Optional[Any] = Field(
        0.25,
        description="A floating-point number representing a fraction.",
        title="Min In Group Fraction",
    )
    min_fold_change: Optional[Any] = Field(
        1,
        description="A floating-point number representing the minimum fold change.",
        title="Min Fold Change",
    )
    max_out_group_fraction: Optional[Any] = Field(
        0.5,
        description="A floating-point number representing the maximum out group fraction.",
        title="Max Out Group Fraction",
    )
    compare_abs: Optional[Any] = Field(
        False,
        description="A boolean parameter that decides whether to compare absolute values of log fold change.",
        title="Compare Abs",
    )
    _api_name: str = PrivateAttr(default="sc.tl.filter_rank_genes_groups")
    _products_original: list[str] = PrivateAttr(
        default=['data.uns["rank_genes_groups"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlMarkerGeneOverlap(BaseAPI):
    """
    Calculate an overlap score between data-derived marker genes and provided markers. Marker gene overlap scores can be quoted as overlap counts, overlap coefficients, or Jaccard indices. The method returns a pandas dataframe which can be used to annotate clusters based on marker gene overlaps. This function was written by Malte Luecken.
    """

    adata: Any = Field(..., description="The annotated data matrix.", title="Adata")
    reference_markers: Any = Field(
        ...,
        description="A marker gene dictionary object where keys are cell identity names and values are sets or lists of strings corresponding to adata.var_name.",
        title="Reference Markers",
    )
    key: Optional[Any] = Field(
        "rank_genes_groups",
        description="The key in adata.uns where the rank_genes_groups output is stored, default is 'rank_genes_groups'.",
        title="Key",
    )
    method: Optional[Any] = Field(
        "overlap_count",
        description="Method to calculate marker gene overlap, with options like 'overlap_count', 'overlap_coef', and 'jaccard'.",
        title="Method",
    )
    normalize: Any = Field(
        None,
        description="Normalization option for marker gene overlap output, dependent on the method chosen.",
        title="Normalize",
    )
    top_n_markers: Any = Field(
        None,
        description="Number of top data-derived marker genes to use, default is 100, can be overridden by adj_pval_threshold.",
        title="Top N Markers",
    )
    adj_pval_threshold: Any = Field(
        None,
        description="Significance threshold on adjusted p-values to select marker genes, applicable when p-values calculated by sc.tl.rank_genes_groups().",
        title="Adj Pval Threshold",
    )
    key_added: Optional[Any] = Field(
        "marker_gene_overlap",
        description="Name of the field in .uns that will store the marker overlap scores.",
        title="Key Added",
    )
    inplace: Optional[Any] = Field(
        False,
        description="Option to return a marker gene dataframe or store it directly in adata.uns.",
        title="Inplace",
    )
    _api_name: str = PrivateAttr(default="sc.tl.marker_gene_overlap")
    _products_original: list[str] = PrivateAttr(
        default=['data.uns["marker_gene_overlap"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlScoreGenes(BaseAPI):
    """
    Score a set of genes by calculating the average expression after subtracting the average expression of a reference set of genes, sampled from a gene pool. The method is based on the approach in Seurat by Satija (2015) and has been implemented for Scanpy by Davide Cittaro.
    """

    adata: Any = Field(
        ...,
        description="The annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    gene_list: Any = Field(
        ...,
        description="The list of gene names used for score calculation. Original type annotation: Sequence[str] | pd.Index[str]",
        title="Gene List",
    )
    ctrl_as_ref: Optional[Any] = Field(
        True,
        description="Allow the algorithm to use the control genes as reference. Will be changed to `False` in scanpy 2.0. Original type annotation: bool",
        title="Ctrl As Ref",
    )
    ctrl_size: Optional[Any] = Field(
        50,
        description="Number of reference genes to be sampled from each bin. If `len(gene_list)` is not too low, you can set `ctrl_size=len(gene_list)`. Original type annotation: int",
        title="Ctrl Size",
    )
    gene_pool: Any = Field(
        None,
        description="Genes for sampling the reference set. Default is all genes. Original type annotation: Sequence[str] | pd.Index[str] | None",
        title="Gene Pool",
    )
    n_bins: Optional[Any] = Field(
        25,
        description="Number of expression level bins for sampling. Original type annotation: int",
        title="N Bins",
    )
    score_name: Optional[Any] = Field(
        "score",
        description="Name of the field to be added in `.obs`. Original type annotation: str",
        title="Score Name",
    )
    random_state: Optional[Any] = Field(
        0,
        description="The random seed for sampling. Original type annotation: _LegacyRandom",
        title="Random State",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Copy `adata` or modify it inplace. Original type annotation: bool",
        title="Copy",
    )
    use_raw: Any = Field(
        None,
        description="Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present. Default value changed from `False` to `None` in version 1.4.5. Original type annotation: bool | None",
        title="Use Raw",
    )
    layer: Any = Field(
        None,
        description="Key from `adata.layers` whose value will be used to perform tests on. Original type annotation: str | None",
        title="Layer",
    )
    _api_name: str = PrivateAttr(default="sc.tl.score_genes")
    _products_original: list[str] = PrivateAttr(default=['data.obs["score"]'])
    _data_name: str = PrivateAttr(default="adata")


class ScTlScoreGenesCellCycle(BaseAPI):
    """
    Score cell cycle genes :cite:p:`Satija2015`. Given two lists of genes associated to S phase and G2M phase, calculates scores and assigns a cell cycle phase (G1, S or G2M). See :func:`~scanpy.tl.score_genes` for more explanation.
    """

    adata: Any = Field(..., description="The annotated data matrix.", title="Adata")
    s_genes: Any = Field(
        ..., description="List of genes associated with S phase.", title="S Genes"
    )
    g2m_genes: Any = Field(
        ..., description="List of genes associated with G2M phase.", title="G2M Genes"
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Copy `adata` or modify it inplace.",
        title="Copy",
    )
    _api_name: str = PrivateAttr(default="sc.tl.score_genes_cell_cycle")
    _products_original: list[str] = PrivateAttr(
        default=['data.obs["S_score"]', 'data.obs["G2M_score"]', 'data.obs["phase"]']
    )
    _data_name: str = PrivateAttr(default="adata")


class ScTlDrawGraph(BaseAPI):
    """
    Force-directed graph drawing: An alternative to tSNE that often preserves the topology of the data better. This method requires running scanpy.pp.neighbors first. The default layout uses the package 'fa2-modified', which can be installed via 'pip install fa2-modified'. Force-directed graph drawing describes a class of algorithms for visualizing graphs and has been suggested for visualizing single-cell data. Various layouts implemented in igraph are available, with similar approaches used by other studies.
    """

    adata: Any = Field(
        ...,
        description="Annotated data matrix. Original type annotation: AnnData",
        title="Adata",
    )
    layout: Optional[Any] = Field(
        "fa",
        description="'fa' (`ForceAtlas2`) or any valid `igraph layout`. Specify layout options like 'fr', 'grid_fr', 'kk', 'lgl', 'drl', and 'rt'. Original type annotation: _Layout",
        title="Layout",
    )
    init_pos: Any = Field(
        None,
        description="`'paga'`/`True`, `None`/`False`, or any valid 2d-`.obsm` key. Use precomputed coordinates for initialization. If `False`/`None`, initialize randomly. Original type annotation: str | bool | None",
        title="Init Pos",
    )
    root: Any = Field(
        None,
        description="Root for tree layouts. Original type annotation: int | None",
        title="Root",
    )
    random_state: Optional[Any] = Field(
        0,
        description="For layouts with random initialization like 'fr', change this to use different initial states for the optimization. If `None`, no seed is set. Original type annotation: _LegacyRandom",
        title="Random State",
    )
    n_jobs: Any = Field(
        None,
        description="No description available. Original type annotation: int | None",
        title="N Jobs",
    )
    adjacency: Any = Field(
        None,
        description="Sparse adjacency matrix of the graph, defaults to neighbors connectivities. Original type annotation: spmatrix | None",
        title="Adjacency",
    )
    key_added_ext: Any = Field(
        None,
        description="By default, append `layout`. Original type annotation: str | None",
        title="Key Added Ext",
    )
    neighbors_key: Any = Field(
        None,
        description="Specify where to look for connectivities. Original type annotation: str | None",
        title="Neighbors Key",
    )
    obsp: Any = Field(
        None,
        description="Use .obsp[obsp] as adjacency. Cannot specify both `obsp` and `neighbors_key` at the same time. Original type annotation: str | None",
        title="Obsp",
    )
    copy_: Optional[Any] = Field(
        False,
        alias="copy",
        description="Return a copy instead of writing to adata. Original type annotation: bool",
        title="Copy",
    )
    kwds: Any = Field(..., description="No description available.", title="Kwds")
    _api_name: str = PrivateAttr(default="sc.tl.draw_graph")
    _products_original: list[str] = PrivateAttr(
        default=['data.uns["draw_graph"]', 'data.obsm["X_draw_graph_fa"]']
    )
    _data_name: str = PrivateAttr(default="adata")


TOOLS_DICT = {
    "sc.tl.paga": ScTlPaga,
    "sc.tl.leiden": ScTlLeiden,
    "sc.tl.louvain": ScTlLouvain,
    "sc.tl.umap": ScTlUmap,
    "sc.tl.tsne": ScTlTsne,
    "sc.tl.diffmap": ScTlDiffmap,
    "sc.tl.embedding_density": ScTlEmbeddingDensity,
    "sc.tl.rank_genes_groups": ScTlRankGenesGroups,
    "sc.tl.filter_rank_genes_groups": ScTlFilterRankGenesGroups,
    "sc.tl.marker_gene_overlap": ScTlMarkerGeneOverlap,
    "sc.tl.score_genes": ScTlScoreGenes,
    "sc.tl.score_genes_cell_cycle": ScTlScoreGenesCellCycle,
    "sc.tl.draw_graph": ScTlDrawGraph,
}
