from biochatter.api_agent.gen_data_model import apis_to_data_models

def test_generated_data_model():
    file_name = "biochatter/api_agent/python/scanpy/tmp_tl.py"
    import scanpy as sc
    api_dict = {
        "meta": {
            "package": sc,
            "module": sc.tl,
        },
        "api_list": [
            {
                "api": sc.tl.umap,
                "products": [
                    "data.obsm[\"X_umap\"]",
                    "data.uns[\"umap\"]"
                ],
                "data_name": "adata"
            },
            {
                "api": sc.tl.tsne,
                "products": [
                    "data.obsm[\"X_tsne\"]",
                    "data.uns[\"tsne\"]"
                ],  
                "data_name": "adata"
            }
        ]
    }
    data_models, codes_list = apis_to_data_models(api_dict)

    codes = "\n\n".join(codes_list)

    with open(file_name, "w") as f:
        f.write(codes)

    import os
    from biochatter.api_agent.python.scanpy import tmp_tl

    sc.settings.datasetdir = os.environ.get("DATA")
    data = sc.datasets.pbmc3k()
    umap = tmp_tl.ScanpyToolsUmap(
        adata=data,
        min_dist=0.5,
        spread=1.0,
        n_components=2,
        maxiter=None,
        alpha=1.0,
    )
    tsne = tmp_tl.ScanpyToolsTsne(
        adata=data,
        n_pcs=None,
        use_rep=None,
        perplexity=30,
        metric='euclidean',
        early_exaggeration=12,
        learning_rate=1000,
    )

    pass

if __name__ == "__main__":
    test_generated_data_model()