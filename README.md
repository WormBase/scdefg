# scdefg: scvi-tools Differential Expression Flask GUI

This app offer a single page Flask UI that allows you to quickly select cell groups and perform differential expression
on single cell RNA sequencing data using [scvi-tools](https://scvi-tools.org)

You can see a deployment with 3 _C. elegans_ datasets at [scdefg.textpressolab.com](https://scdefg.textpressolab.com/). 
For more information on this see [single-cell.wormbase.org](https://single-cell.wormbase.org).



# How to launch

The app takes in a pretrained [scVI model](https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html) saved with
and the corresponding anndata (using the option `save_anndata = True`). You need to specify the path where the model is
saved at when launching the app. 

The `adata` anndata should have the
fields `adata.var.gene_name` and `adata.var.gene_id` as strings containing the respective gene names and gene IDs, so that results 
will be properly displayed. If the field `adata.var.gene_description` is also present, they will be shown during the mousover on the interactive plots.
Text files with gene descriptions can be downloaded [here](https://www.alliancegenome.org/downloads)

Additionally, at least one column should be present in the `adata.obs`, such as for example `cell_type`. 
The cell selection menu can be stratified according to any 
number of columns that are present in the adata file. When launching the app, you just need to provide the name
of each column with the -s flag , eg `-s cell_type`. So for example, to offer the user the option to group cells
by `cell_type`, `tissue` and `experiment_code` you would launch the app using the following arguments:
```
scdefg/app.py ./model -s tissue -s cell_type -s experiment_code
```

Alternatively you could provide the columns by which to stratify data may be provided in the field `adata.uns['selection_columns']`
in a list, e.g. `adata.uns['selection_columns']=['sample','cell_type']`. 

You _should_ try provide the names of tha adata.obs columns that contain relevant conditions to stratify the data by.
Otherwise the app defaults to using all columns, and the selection tables will become **extremely** slow with 100k+ rows.

Below is an example for a trained model that you can try to download and launch. It show how to stratify the data using the columns `adata.obs['tissue']`, 
`adata.obs['cell_type']` and `adata.obs['experiment_code']`. T

```
## Example deploy with data from the CeNGEN project 2020 data release (www.cengen.org)

## clone the repo
git clone https://github.com/Munfred/scdefg.git
cd scdefg
## download example trained model and rename the folder to just `model`
wget -q -O tmp.zip https://github.com/Munfred/scvi-de-flask/releases/download/taylor2020/taylor2020_100955cells_11569genes_20210129_scvi.zip && unzip tmp.zip && rm tmp.zip
mv taylor2020_100955cells_11569genes_20210129_scvi model

## launch the app and stratify the selection menu with 3 columns: `cell_type`, `tissue`, `experiment_code`
scdefg/app.py ./model -s tissue -s cell_type -s experiment_code
```


### Here is what the selection menu would look like in this case

![](https://user-images.githubusercontent.com/12504176/107468161-44bc3580-6b46-11eb-9175-d10e9749f747.png)

### Here is a view of the full page after results are returned 

![](https://user-images.githubusercontent.com/12504176/107468037-07f03e80-6b46-11eb-8b27-9bccc3b5e9a6.png)
