from importlib import import_module

DATA = [
    {
        'query': 'Plot umap embedding of cells.',
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.umap(adata=adata)
        sc.pl.scatter(adata=adata, basis='umap')
        """
    },
    {
        'query': 'Plot umap embedding of cells, which are colored by leiden clustering.',
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.umap(adata=adata)
        sc.tl.leiden(adata=adata)
        sc.pl.umap(adata=adata, color='leiden')
        """
    },
    {
        
        'query': "Plot heatmap of gene expressions of genes ['TMSB4X', 'MALAT1', 'B2M'], where cells are clustered by louvain algorithm.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.louvain(adata=adata)
        sc.pl.heatmap(adata=adata, var_names=['TMSB4X', 'MALAT1', 'B2M'], groupby='louvain')
        """
    },
    {
        'query': "Plot dotplot with dendrogram of gene expressions of genes ['TMSB4X', 'MALAT1', 'B2M'], where cells are clustered by louvain algorithm.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.louvain(adata=adata)
        sc.pl.dotplot(adata=adata, var_names=['TMSB4X', 'MALAT1', 'B2M'], groupby='louvain')
        """
    },
    {
        'query': "Plot violinplot of gene expressions of genes ['TMSB4X', 'MALAT1', 'B2M'], where cells are clustered by louvain algorithm.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.louvain(adata=adata)
        sc.pl.violin(adata=adata, var_names=['TMSB4X', 'MALAT1', 'B2M'], groupby='louvain')
        """
    },
    {
        'query': "visualize dendrogram of clusters defined by louvain algorithm on cells.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.louvain(adata=adata)
        sc.pl.dendrogram(adata=adata, groupby='louvain')
        """
    },
    {
        'query': "visualize diffusion map embedding of cells which are clustered by leiden algorithm.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.diffmap(adata=adata)
        sc.tl.leiden(adata=adata)
        sc.pl.diffmap(adata=adata, color='leiden')
        """
    },
    {
        'query': "visualize dispersions versus mean expressions of genes in scatter plot.",
        'codes': \
        """
        sc.pp.log1p(data=adata)
        sc.pp.highly_variable_genes(adata=adata)
        sc.pl.highly_variable_genes(adata=adata)
        """
    },
    {
        'query': "visualize PCA embedding of cells.",
        'codes': \
        """
        sc.pp.pca(data=adata)
        sc.pl.pca(adata=adata)
        """
    },
    {
        'query': "visualize PCA embedding of cells which are clustered by louvain algorithm.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.louvain(adata=adata)
        sc.pl.pca(adata=adata, color='louvain')
        """
    },
    {
        'query': "visualize tSNE embedding of cells which are clustered by leiden algorithm.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.tsne(adata=adata)
        sc.tl.leiden(adata=adata)
        sc.pl.tsne(adata=adata, color='leiden')
        """
    },
    {
        'query': "visualize umap embedding density of cells.",
        'codes': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.umap(adata=adata)
        sc.tl.embedding_density(adata=adata, basis='umap')
        sc.pl.embedding_density(adata=adata, basis='umap')
        """
    },
    {
        'query': None,
        'codes': \
        """
        sc.pp.pca(data=adata)
        sc.pl.pca_variance_ratio(adata=adata)
        """
    },
    {
        'query': None,
        'code': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.leiden(adata=adata)
        sc.tl.rank_genes_groups(adata=adata, 'louvain')
        sc.pl.rank_genes_groups_dotplot(adata=adata, groupby='leiden')
        """
    },
    {
        'query': None,
        'code': \
        """
        sc.pp.neighbors(adata=adata)
        sc.tl.leiden(adata=adata)
        sc.pl.correlation_matrix(adata=adata, groupby='leiden')        
        """
        
    }
]

META = {
    'state_name': 'scanpy',
    'sc': import_module('scanpy')
}