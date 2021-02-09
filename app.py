print('starting imports....')

import warnings
# ignore annoying pandas future warning
warnings.simplefilter(action='ignore', category=FutureWarning)
from flask import Flask, Response, jsonify, request, render_template, Blueprint, send_file, redirect, url_for
import logging
import pandas as pd
import json
from io import StringIO
import scvi
import numpy as np
import plotly.graph_objects as go
import time

print('scvi-tools version:', scvi.__version__)


logging.basicConfig(level=logging.INFO)
logger=logging.getLogger(__name__)
logger.info('Starting scvi-de-flask...')

app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['JSONIFY_PRETTYPRINT_REGULAR'] = False
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
tables = Blueprint('tables', __name__, url_prefix='/tables')
de_enriched_tables = Blueprint('de_enriched_tables', __name__, url_prefix='/de_enriched_tables')
de_depleted_tables = Blueprint('de_depleted_tables', __name__, url_prefix='/de_depleted_tables')

############################ SCVI AND PANDAS PREP ##############################
# this will load the pretrained model
# if stored in a differently named folder, change the model_name variable to match
# the model must have been saved with the anndata
model_name='model'
model = scvi.model.SCVI.load(model_name, use_cuda=False)
adata = model.adata.copy()
# our adata has annotations by experiment, by cell type and by tissue
# we use the experiment code and cell type to generate a dataframe
# that counts how many cells of each type are present in each experiment
# this dataframe is needed to create the table that allows the user to select
# which cells to compare on the app
# cell type does on the index, experiment on the columns, entries are number of
# cells of that type per experiment
# if interested in selecting using another parameter (eg 'patient') change the 'experiment_code' string to match
experiments= np.sort(adata.obs['experiment_code'].unique())
unique_cell_types=np.sort(adata.obs['cell_type'].unique())
census=pd.DataFrame(index= unique_cell_types)
# census['ct'] = unique_cell_types
# for experiment in experiments:
#     census[experiment] = census.index.map(adata.obs[adata.obs['experiment_code']==experiment]['cell_type'].value_counts())
census.index=census.index.rename('Cell Type')


# the index of the dataframe is ignored when rendering the tables,
# so we do reset_index to put the cell names in the first column
census = census.reset_index()
df_nice_names = census.copy()
df_nice_names.columns = df_nice_names.columns.str.replace('_',' ')
df_nice_names.columns = df_nice_names.columns.str.replace('-','&#8209;')

# convert df to dict for sending as json to datatables, replace underscore with space for aesthetics

dict_df = df_nice_names.replace('_',' ', regex=True).to_dict(orient='records')
# convert column names into dict for sending as json to datatables
columns = [{"data": item, "title": item} for item in df_nice_names.columns]

#### datatables to render the selection tables ####
@tables.route("/", methods=['GET', 'POST'])
def clientside_table_content():
    return jsonify({'data': dict_df, 'columns': columns})
app.register_blueprint(tables)


#### datatables to render the DE results table ####
@de_enriched_tables.route("/", methods=['GET','POST'])
def de_enriched_table_content():
    # convert df to dict for sending as json to datatables
    de_df=pd.read_csv('results.csv', index_col=0).reset_index()

    de_df=de_df[['index','gene_name','minuslog10pval','lfc_mean']].fillna('-')
    de_df.columns=['Gene ID', 'Gene Name', '-log10 p-value','mean log2 fold change' ]
    de_df = de_df.copy()
    # convert df to dict for sending as json to datatables
    de_dict_df = de_df.to_dict(orient='records')
    # convert column names into dict for sending as json to datatables
    columns = [{"data": item, "title": item} for item in de_df.columns]
    return jsonify({'data': de_dict_df, 'columns': columns})
app.register_blueprint(de_enriched_tables)

#### datatables to render the DE results DEPLEPTED table ####
@de_depleted_tables.route("/", methods=['GET','POST'])
def de_depleted_table_content():
    # convert df to dict for sending as json to datatables
    de_df=pd.read_csv('results.csv', index_col=0).reset_index()

    de_df=de_df[['index','gene_name','minuslog10pval','lfc_mean']].fillna('-')
    de_df.columns=['Gene ID', 'Gene Name', '-log10 p-value','mean log2 fold change' ]
    de_df = de_df.copy()
    # convert df to dict for sending as json to datatables
    de_dict_df = de_df.to_dict(orient='records')
    # convert column names into dict for sending as json to datatables
    columns = [{"data": item, "title": item} for item in de_df.columns]
    return jsonify({'data': de_dict_df, 'columns': columns})
app.register_blueprint(de_depleted_tables)




# this is the landing page
@app.route("/", methods=['GET', 'POST'])
def home():
    return render_template("home.html")

@app.route("/results", methods=['POST', 'GET'])
def results():
    return render_template("results.html")

@app.route('/submit', methods=['POST', 'GET'])
def receive_submission():
    logger.info('Got a submission!')
    timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")

    # answer is a dict of json strings containing selected row and column index numbers
    answer = request.form.to_dict(flat=False)
    # need to convert the json strings to dict, then to a data frame
    # data1 is the selection for the first group, data2 for the second
    data1 = json.loads(answer['data1'][0])
    data1_df = pd.DataFrame.from_dict(data1[0])
    data2 = json.loads(answer['data2'][0])
    data2_df = pd.DataFrame.from_dict(data2[0])

    print(data1)
    print(data1_df)

    # now map the index number to experiment name and cell type name
    group1 = pd.DataFrame()
    group1['cell_type1'] = data1_df['row'].map(census['Cell Type'])
    print(group1)

    group2 = pd.DataFrame()
    group2['cell_type2'] = data2_df['row'].map(census['Cell Type'])

    genes = StringIO(json.loads(answer['genes'][0]))
    genes_df = pd.read_csv(genes, names=['selected_genes'])

    selected_groups_df = pd.concat([group1, group2, genes_df], axis=1)
    selected_groups_df.to_csv('selected_groups.csv')


#### Creates the masks for the selected cell types

    # first create the mask as an array of all false
    # then for each group in the data add them to the mask
    group1_mask = adata.obs.index != adata.obs.index
    for idx, row in group1.iterrows():
        mask = adata.obs['cell_type']==row['cell_type1']
        group1_mask = group1_mask | mask

    group2_mask = adata.obs.index != adata.obs.index
    for idx, row in group2.iterrows():
        mask = adata.obs['cell_type']==row['cell_type2']
        group2_mask = group2_mask | mask

    # the masks then define the two groups of cells on which to perform DE
    fdr_target=0.005
    de = model.differential_expression( adata,
                                       idx1=group1_mask,
                                       idx2=group2_mask,
                                       fdr_target=fdr_target)

#### Wrangles the DE results dataframe a bit

    # first we create these variables to customize the hover text in plotly's heatmap
    # the text needs to be arranged in a matrix the same shape as the heatmap
    #try to add gene descriptions and gene names if the adata has those, otherwise add a blank

    try:
        de['gene_description']=de.index.map(adata.var['gene_description'])
        # for the gene descriptions text, which can be long, we add line breaks
        de['gene_description_html'] = de['gene_description'].str.wrap(80).str.replace('\n','<br>')
        de['gene_name']=de.index.map(adata.var['gene_name']).astype(str)
        de['gene_id'] = de.index.astype(str)
        # de['gene_id'] = de.index.map(adata.var['gene_id'])

    except:
        de['gene_description_html'] = 'gene description here'
        de['gene_name']='gene name here'
        de['gene_id']='gene ID here'

    de['gene_name']=de['gene_name'].fillna('-')
    # calculate the -log10(p-value) for the volcano
    de['minuslog10pval']=-np.log10(de['proba_not_de'] + 0.00001)

    # all genes are initially colored black
    de['color'] = 'black'
    # uncomment line below to color genes by FDR significance
    # de['color'] = de['is_de_fdr_'+str(fdr_target)].map({True:'steelblue',False:'gray'})

    # then we loops through the list of genes provided to color some red
    # gene ids should be a perfect match
    pd.set_option('mode.chained_assignment',None) #supress warning
    de['color'][de['gene_id'].isin(genes_df['selected_genes'].values)] = 'red'
    # gene names should be a partial match
    for partial_string in genes_df['selected_genes'].values:
        de['color'][de['gene_name'].str.contains(partial_string)] = 'red'


    group1_str = ''
    for cell1 in group1.cell_type1.values:
        if group1_str != '':  group1_str = group1_str + ', '
        group1_str = group1_str + str(cell1)
    print(group1_str)

    group2_str = ''
    for cell2 in group2.cell_type2.values:
        if group2_str != '':  group2_str = group2_str + ', '
        group2_str = group2_str + str(cell2)
    print(group2_str)
#### This makes the volcano plot using plotly
    fig = go.Figure(
                    data=go.Scatter(
                            x=de["lfc_mean"].round(3)
                            , y=de["minuslog10pval"].round(3)
                            , mode='markers'
                            , marker=dict(
                                        color=de['color'],
                                        opacity=0.5)
                            , hoverinfo='text'
                            , text=de['gene_description_html']
                            , customdata=de.gene_name + '<br>' + de.gene_id
                            , hovertemplate='%{customdata} <br>' +
                                            '-log10 p-value: %{y}<br>' +
                                            'Mean log fold change: %{x}' +
                                            '<extra>%{text}</extra>'
                            )
                            , layout= {
                                    "title": {"text":
                                            group1_str + ' cells versus  <br>' + group2_str + " <br> differential expression results, (dashes mark p = 0.01)"
                                            , 'x':0.5
                                            }
                                    , 'xaxis': {'title': {"text": "Mean log fold change"}}
                                    , 'yaxis': {'title': {"text": "-log10 p-value"}}
        #                             , "height": 600,
        #                             , "width":1000
                            }
                )
    fig.update_layout(hovermode='closest', template='none')
    fig.add_shape(type="line", x0=-6, y0=2, x1=6, y1=2, line=dict(color="lightsalmon", width=2, dash="dash"))

    # overwrites the last figure in order to serve it in the results page
    fig.write_html('./static/fig.html')
    de_csv = de[['gene_name','minuslog10pval','lfc_mean','lfc_std','proba_not_de',
                 # 'is_de_fdr_'+str(fdr_target)
                 ]].round(2)
    de_csv.to_csv('results.csv')

    # if you want to keep all results uncomment the lines below
    # de_csv.to_csv(timestamp+'results.csv')
    # fig.write_html( timestamp + '_fig.html')
    # selected_groups_df.to_csv(timestamp + 'selected_groups.csv')


    return redirect(url_for('results'))

@app.route("/get_results_csv")
def get_results_csv():
    with open("results.csv") as fp:
        csv = fp.read()
    return Response( csv, mimetype="text/csv",
        headers={"Content-disposition":"attachment; filename=DE_results.csv"})

@app.route("/get_groups_csv")
def get_groups_csv():
    with open("selected_groups.csv") as fp:
        csv = fp.read()
    return Response( csv, mimetype="text/csv",
        headers={"Content-disposition":"attachment; filename=selected_groups.csv"})

@app.route("/get_plot")
def get_plot():
    with open("./static/fig.html") as fp:
        csv = fp.read()
    return Response( csv, mimetype="text/html",
        headers={"Content-disposition":"attachment; filename=fig.html"})

if __name__ == "__main__":
    app.run()