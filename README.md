# This repo was moved to https://github.com/WormBase/scdefg
Please refer to that repo instead.

========

# scdefg: scvi-tools Differential Expression Flask GUI

This app offer a single page Flask UI that allows you to quickly select cell groups and perform differential expression
on single cell RNA sequencing data using [scvi-tools](https://scvi-tools.org)

## Select some cell types and click submit:

![](https://user-images.githubusercontent.com/12504176/107468161-44bc3580-6b46-11eb-9175-d10e9749f747.png)

## Results are returned in ~15s with an interactive volcano plot and sortable tables:

![](https://user-images.githubusercontent.com/12504176/107468037-07f03e80-6b46-11eb-8b27-9bccc3b5e9a6.png)

# How to deploy

The app takes in a pretrained scvi model and the corresponding anndata. At a minimum, the adata file should have the
fields `adata.var.gene_name` and `adata.var.gene_id` as strings containing the respective gene names and gene IDs. If
the field
`adata.var.gene_description` is also present, they will be shown during the mousover on the interactive plots. Gene
descriptions are available at the site of the Alliance for Genome Resources: https://www.alliancegenome.org/downloads

Additionally, at least one field should be present in `adata.obs`. The data can be stratified according to any 
number of columns that are present in the adata file.

Columns by which to stratify data may be provided in the field `adata.uns['selection_columns']`
in a list, e.g. `adata.uns['selection_columns']=['sample','cell_type']`. Alternatively you may provide the name
of one column at a time with the -s flag , eg `-s cell_type`.

You _should_ try provide the names of tha adata.obs columns that contain relevant conditions to stratify the data by.
Otherwise the app defaults to using all columns, and the selection tables will become extremely slow with 100k+ rows,
so you have been warned.

Below is an example deploy showing how to stratify the data using the columns `adata.obs['tissue']`, 
`adata.obs['cell_type']` and `adata.obs['experiment_code']`. 
```
## Example deploy with data from the CeNGEN project 2020 data release (www.cengen.org)
git clone https://github.com/Munfred/scdefg.git
cd scdefg
wget -q -O tmp.zip https://github.com/Munfred/scvi-de-flask/releases/download/taylor2020/taylor2020_100955cells_11569genes_20210129_scvi.zip && unzip tmp.zip && rm tmp.zip
mv taylor2020_100955cells_11569genes_20210129_scvi model
scdefg/app.py ./model -s tissue -s cell_type -s experiment_code
```
