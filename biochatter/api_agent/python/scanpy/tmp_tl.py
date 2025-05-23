from __future__ import annotations

from typing import Any, Optional

from pydantic import ConfigDict, Field, PrivateAttr

from biochatter.api_agent.base.agent_abc import BaseAPI
from biochatter.api_agent.python.scanpy.info_hub import dep_graph_dict


class ScanpyToolsUmap(BaseAPI):
    
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
    _dep_graph_dict = PrivateAttr(default=dep_graph_dict)



class ScanpyToolsTsne(BaseAPI):
    
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
    _dep_graph_dict = PrivateAttr(default=dep_graph_dict)