# scvi-de-flask
A proper Flask app to do the things on https://github.com/munfred/scdefg


The app takes in a pretrained scvi model and the corresponding anndata with embeddings stored followomg the scanpy convention of starting with X_, e.g., `anndata.obsm['X_tsne'], anndata.obsm['X_umap']`. 

To avoid storing data on git while making the app easy to test and deploy, trained models and the adata are uploaded to the releases page, and may be retried with `wget -nc` (to avoid downloading multiple times). Or to unzip and then remove the downloaded file:

the two lines below will download the taylor model and rename the folder to be called just `model`

```
wget -q -O tmp.zip https://github.com/Munfred/scvi-de-flask/releases/download/taylor2020/taylor2020_100955cells_11569genes_20210129_scvi.zip && unzip tmp.zip && rm tmp.zip

mv taylor2020_100955cells_11569genes_20210129_scvi movel
```

https://github.com/Munfred/scvi-de-flask/releases


**Note on loading other adata files**

At the moment the code relies on the specific names of the adata.obs for cells and experiments, which here are `cell_type` and `experiment_code` so if the new model adata file you load has different column names, you'll need to change the code to match


## TODOS:
- Fix the DE results table, which currently breaks with more than ~2000 entries. Something to do with datatables...
- Make the results page tables show up side by side
- Add "loading" page while ressults are processing (5-10s) - flask currently redirects immediately for some reason, and it will break it nothing has been computed, or load the old results
- Figure if the p-values can be a little finer to avoid discretization on the volcano plot...
