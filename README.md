# scvi-de-flask
This app offer a single page Flask UI that allows you to quickly select cell groups and perform differential expression on single cell RNA sequencing data using [scvi-tools](https://scvi-tools.org)

## Select some cell types and click submit:
<img width="1094" alt="Screen Shot 2021-02-10 at 02 18 14" src="https://user-images.githubusercontent.com/12504176/107468161-44bc3580-6b46-11eb-9175-d10e9749f747.png">





## Results are returned in ~15s with an interactive volcano plot and sortable tables:
<img width="1214" alt="Screen Shot 2021-02-10 at 02 16 37" src="https://user-images.githubusercontent.com/12504176/107468037-07f03e80-6b46-11eb-8b27-9bccc3b5e9a6.png">


# How to deploy


The app takes in a pretrained scvi model and the corresponding anndata with embeddings stored followomg the scanpy convention of starting with X_, e.g., `anndata.obsm['X_tsne'], anndata.obsm['X_umap']`. 

To avoid storing data on git while making the app easy to test and deploy, trained models and the adata are uploaded to the releases page, and may be retried with `wget -nc` (to avoid downloading multiple times). Or to unzip and then remove the downloaded file:

the two lines below will download the taylor model and rename the folder to be called just `model`

This is data from the CeNGEN project 2020 data release, described in the preprint Molecular topography of an entire nervous system. Just select cell types and experiments to compare and some genes to highlight. It will produce an interactive volcano plot and table with DE results in ~15s.

```
wget -q -O tmp.zip https://github.com/Munfred/scvi-de-flask/releases/download/taylor2020/taylor2020_100955cells_11569genes_20210129_scvi.zip && unzip tmp.zip && rm tmp.zip

mv taylor2020_100955cells_11569genes_20210129_scvi model
```

https://github.com/Munfred/scvi-de-flask/releases


**Note on loading other adata files**

At the moment the code relies on the specific names of the adata.obs for cells and experiments, which here are `cell_type` and `experiment_code` so if the new model adata file you load has different column names, you'll need to change the code to match


### Simple deployment

If you don't care about production robustness and just want a few users to access the app concurrently,
an easy deployment can be done following this tutorial up to and including section 3: 
https://www.digitalocean.com/community/tutorials/how-to-serve-flask-applications-with-gunicorn-and-nginx-on-ubuntu-18-04

Then start the app using gunicorn and specify the number of workers and timeout. 5000 is the default flask port.

```
gunicorn --workers 3 --timeout 120  --bind 0.0.0.0:5000 wsgi:app 
```

If you are doing a production deployment, then follow the guide until the end.