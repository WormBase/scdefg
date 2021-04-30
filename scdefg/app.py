print('starting imports....')
# the import below is the user configs
from scdefg import config
import warnings
# ignore annoying pandas future warning
warnings.simplefilter(action='ignore', category=FutureWarning)
from flask import Flask, jsonify, request, render_template, Blueprint
import logging
import pandas as pd
import json
from io import StringIO
import scvi
import numpy as np
import plotly.graph_objects as go

print('scvi-tools version:', scvi.__version__)
logging.basicConfig(level=logging.INFO)
logger=logging.getLogger(__name__)
logger.info('Starting scdefg...')
app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['JSONIFY_PRETTYPRINT_REGULAR'] = False
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
tables = Blueprint('tables', __name__, url_prefix='/tables')
user_configs = Blueprint('user_configs', __name__, url_prefix='/user_configs')
de_tables = Blueprint('de_tables', __name__, url_prefix='/de_tables')
de_depleted_tables = Blueprint('de_depleted_tables', __name__, url_prefix='/de_depleted_tables')



introduction_html_text='default'
table_selection_columns=['experiment_code', 'tissue']
model_file='../model'


############################ SCVI AND PANDAS PREP ##############################
# this will load the pretrained model
# you specify where the model is stored in the config file
# the model must have been saved with the anndata .h5ad format

model = scvi.model.SCVI.load(model_file, use_gpu=False)
adata = model.adata.copy()
print(adata.obs.columns)

#crates the selection table using the columns specified in the config file
selection_table_df = adata.obs.groupby(table_selection_columns, observed = True).size().rename('#cells').reset_index()
selection_table_dict = selection_table_df.to_dict(orient='records')
# convert column names into dict for sending as json to datatables
columns = [{"data": item, "title": item} for item in selection_table_df.columns]

#### datatables to render the selection tables ####
@tables.route("/", methods=['GET', 'POST'])
def clientside_table_content():
    return jsonify({'data': selection_table_dict, 'columns': columns})
app.register_blueprint(tables)

#### datatables to render the selection tables ####
@user_configs.route("/", methods=['GET', 'POST'])
def send_user_configs():
    return jsonify({'intro':introduction_html_text,
                    'table_selection_columns':table_selection_columns})
app.register_blueprint(user_configs)

# this is the landing page
@app.route("/", methods=['GET', 'POST'])
def home():
    return render_template("home.html")

@app.route('/submit', methods=['POST', 'GET'])
def receive_submission():
    logger.info('Got a submission!')
    timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")

    # answer is a dict of json strings containing selected row and column index numbers
    answer = request.form.to_dict(flat=False)
    # need to convert the json strings to dict, then to a data frame
    # data1 is the selection for the first group, data2 for the second
    data1 = json.loads(answer['data1'][0])
    data1_df = pd.DataFrame.from_records(data1)
    data2 = json.loads(answer['data2'][0])
    data2_df = pd.DataFrame.from_records(data2)

    genes = StringIO(json.loads(answer['genes'][0]))
    genes_df = pd.read_csv(genes, names=['selected_genes'])

#### Creates the masks for the selected cell types

    # first create the mask as an array of all false
    # then for each group in the data add them to the mask
    group1_mask = adata.obs.index != adata.obs.index
    for idx, row in data1_df.iterrows():
        # create an array of booleans that starts with all of them true
        current_mask = adata.obs.index == adata.obs.index
        # need to do the _ replace for matching with the pretty names
        for col in table_selection_columns:
            # perform successive AND operations to create the mask
            # that keeps only the cells matching the user selection
            partial_mask = adata.obs[col] == row[col]
            current_mask = current_mask & partial_mask
        group1_mask = group1_mask | current_mask

    # first create the mask as an array of all false
    # then for each group in the data add them to the mask
    group2_mask = adata.obs.index != adata.obs.index
    for idx, row in data2_df.iterrows():
        # create an array of booleans that starts with all of them true
        current_mask = adata.obs.index == adata.obs.index
        for col in table_selection_columns:
            # perform successsive AND operations to create the mask
            # that keeps only the cells matching the user selection
            partial_mask = adata.obs[col] == row[col]
            current_mask = current_mask & partial_mask
        group2_mask = group2_mask | current_mask

    # the masks then define the two groups of cells on which to perform DE
    de = model.differential_expression( adata,
                                       idx1=group1_mask,
                                       idx2=group2_mask)

#### Wrangles the DE results dataframe a bit

    # first we create these variables to customize the hover text in plotly's heatmap
    # the text needs to be arranged in a matrix the same shape as the heatmap
    # try to add gene descriptions and gene names if the adata has those, otherwise add a blank
    # the adata.var fields should include a field named gene_id and gene_name, otherwise they will be filled with a blank default

    try:
        de['gene_description']=de.index.map(adata.var['gene_description'])
        # for the gene descriptions text, which can be long, we add line breaks
        de['gene_description_html'] = de['gene_description'].str.wrap(80).str.replace('\n','<br>')
        de['gene_name']=de.index.map(adata.var['gene_name']).astype(str)
        de['gene_id'] = de.index.astype(str)
        # de['gene_id'] = de.index.map(adata.var['gene_id'])

    except:
        de['gene_description_html'] = 'warning: adata.var["gene_description"] does not exist, filled with blank'
        de['gene_name']='warning: adata.var["gene_name"] does not exist, filled with blank'
        de['gene_id']='warning: adata.var["gene_id"] does not exist, filled with blank'

    de['gene_name']=de['gene_name'].fillna('-')
    # calculate the -log10(p-value) for the volcano
    de['minuslog10pval']=-np.log10(de['proba_not_de'] + 0.00001)

    de['log10mean_expression']=np.log10((de['scale1'] + de['scale2'])/2)


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


    group1_str = 'group1str'
    # for cell1 in group1.cell_type1.values:
    #     if group1_str != '':  group1_str = group1_str + ', '
    #     group1_str = group1_str + str(cell1)

    group2_str = 'group2str'
    # for cell2 in group2.cell_type2.values:
    #     if group2_str != '':  group2_str = group2_str + ', '
    #     group2_str = group2_str + str(cell2)

#### This makes the volcano plot using plotly
    defig = go.Figure(
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
                                            'Mean log2 fold change: %{x}' +
                                            '<extra>%{text}</extra>'
                            )
                            , layout= {
                                    "title": {"text":
                                            group1_str + ' cells versus  <br>' + group2_str + " <br> dashes mark p = 0.01"
                                            , 'x':0.5
                                            }
                                    , 'xaxis': {'title': {"text": "Mean log fold change"}}
                                    , 'yaxis': {'title': {"text": "-log10 p-value"}}
                                    , "height": 700,
        #                             , "width":1000
                            }
                )
    defig.update_layout(hovermode='closest', template='none')
    defig.add_shape(type="line", x0=-6, y0=2, x1=6, y1=2, line=dict(color="lightsalmon", width=2, dash="dash"))
    # overwrites the last figure in order to serve it in the results page
    defig = defig.to_html()

#### This makes the MA plot using plotly
    mafig = go.Figure(
        data=go.Scatter(
            x=de["log10mean_expression"].round(3)
            , y=de["lfc_mean"].round(3)
            , mode='markers'
            , marker=dict(
                color=-de['minuslog10pval'],
                colorscale='Viridis',
                opacity=1)
            , hoverinfo='text'
            , text=de['gene_description_html']
            , customdata=de.gene_name + '<br>' + de.gene_id + '<br>-log10 p-value: ' + de.minuslog10pval.round(2).astype(str)
            , hovertemplate='%{customdata} <br>' +
                            'log10 mean expression: %{x}<br>' +
                            'mean log2 fold change: %{y}' +
                            '<extra>%{text}</extra>'
        )
        , layout={
            "title": {"text":
                          group1_str + ' cells versus  <br>' + group2_str
                , 'x': 0.5
                      }
            , 'xaxis': {'title': {"text": "log10 scvi normalized expression"}}
            , 'yaxis': {'title': {"text": "Mean log2 fold change"}}
            , "height": 700,
            #                             , "width":1000
        }
    )
    mafig.update_layout(hovermode='closest', template='none')
    # overwrites the last figure in order to serve it in the results page
    mafig = mafig.to_html()

    ### creates the results tables in a dataframe that is converted to json
    de_df = de[['gene_name','minuslog10pval','lfc_mean','lfc_std','proba_not_de','log10mean_expression']]
    print(de_df)
    de_df.index.rename('gene_id', inplace=True)
    de_df=de_df.reset_index()
    de_df=de_df[['gene_id','gene_name','minuslog10pval','lfc_mean','log10mean_expression']].fillna('-')
    de_df.columns=['Gene ID', 'Gene Name', '-log10 p-value','mean log2 fold change', 'log 10 mean expression' ]
    de_df['mean log2 fold change']=de_df['mean log2 fold change'].astype(float).round(2)
    de_df['-log10 p-value'] = de_df['-log10 p-value'].astype(float).round(2)
    de_df['log 10 mean expression'] = de_df['log 10 mean expression'].astype(float).round(2)
    # convert df to dict for sending as json to datatables
    de_dict_df = de_df.to_dict(orient='records')
    # convert column names into dict for sending as json to datatables
    columns = [{"data": item, "title": item} for item in de_df.columns]

    title = group1_str + ' versus ' + group2_str

    return jsonify({'deplothtml':defig ,'maplothtml':mafig ,'dejsondata':{'data': de_dict_df, 'columns': columns}, 'title':title})

if __name__ == "__main__":
    app.run()