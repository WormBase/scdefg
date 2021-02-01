# scvi-de-flask
A proper Flask app to do the things on https://github.com/munfred/scdefg


The app takes in a pretrained scvi model and the corresponding anndata with embeddings stored followomg the scanpy convention of starting with X_, e.g., `anndata.obsm['X_tsne'], anndata.obsm['X_umap']`. 

To avoid storing data on git while making the app easy to test and deploy, trained models and the adata are uploaded to the releases page, and may be retried with `wget -nc` (to avoid downloading multiple times). Or to unzip and then remove the downloaded file:

```
wget -q -O tmp.zip https://github.com/Munfred/scvi-de-flask/releases/xxxxxx && unzip tmp.zip && rm tmp.zip
```

https://github.com/Munfred/scvi-de-flask/releases
