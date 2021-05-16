
import time
import warnings
from flask import Flask, jsonify, request, render_template, Blueprint
import logging
import pandas as pd
import json
from io import StringIO
import scvi
import numpy as np
import plotly.graph_objects as go
import click


# ignore annoying pandas future warning
warnings.simplefilter(action='ignore', category=FutureWarning)

print('🧮 🧬 📊 🧫 🧪 scdefg: scvi-tools Differential Expression Flask GUI 📱 🎲 📈 🦠 📉 🎰')
print('You are using scvi-tools version:', scvi.__version__)

# set up logs and flask blueprints
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
tables = Blueprint('tables', __name__, url_prefix='/tables')
user_configs = Blueprint('user_configs', __name__, url_prefix='/user_configs')
de_tables = Blueprint('de_tables', __name__, url_prefix='/de_tables')
de_depleted_tables = Blueprint('de_depleted_tables', __name__, url_prefix='/de_depleted_tables')


### declare command line arguments and help texts
@click.command()
@click.version_option("Development version 2021-04-30 - Eduardo's branch")
@click.argument('scvi_tools_model_path', type=click.Path(exists=True))
@click.option('-s', '--selection_columns', multiple=True,
              help="""Name of one of the adata.obs columns to be used in the selection tables.
                    Can be passed multiple times for using multiple columns.
                    A list of column names may also be provided on adata.uns['selection_columns'], e.g.
                    adata.uns['selection_columns']=['sample','cell_type'] """)
@click.option('-i', '--intro_text_html', type=click.Path(exists=True),
              help="""The path to a file containing an html snippet to be displayed at the introduction div of the app
                   Can also be provided at adata.uns['intro_text_html'] as a multiline string.
                   If not provided a default introduction is used. """)
# BUG: FOR SOME REASON PASSING 0.0.0.0 AS AN ARGUMENT DOESN'T WORK AND BREAKS
@click.option('-h', '--host', default='localhost',
              help="Host on which to launch the scdefg app. Default: localhost")
@click.option('-p', '--port', default='1337',
              help="Port on which to launch the scdefg app. Default: 1337")
def launch(scvi_tools_model_path, selection_columns, intro_text_html, host, port):
    """
    🧮 🧬 📊 🧫 🧪 scdefg: scvi-tools Differential Expression Flask GUI 📱 🎲 📈 🦠 📉 🎰

    Provide the path to the folder of a saved scvi-tools model and the scdefg app will be launched.
    The trained model must have been saved together with the adata.h5ad file.

    Columns by which to stratify data may be provided in the field adata.uns['selection_columns']
    in a list, e.g. `adata.uns['selection_columns']=['sample','cell_type']`. Alternatively you may provide the name
    of one column at a time with the -s flag , eg `-s cell_type`.

    You _should_ try provide the names of tha adata.obs columns that contain relevant conditions to stratify the data by.
    Otherwise the app defaults to using all columns, and the selection tables will become extremely slow with 100k+ rows,
    so you have been warned.
    """
    # the string above is used by the click library to provide the description of the --help text command line argument

    app = Flask(__name__, template_folder="templates", static_folder="static")
    app.config["JSONIFY_PRETTYPRINT_REGULAR"] = False
    app.config["TEMPLATES_AUTO_RELOAD"] = True
    app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0

    ############################ SCVI AND PANDAS PREP ##############################
    # this will load the pretrained model
    # you specify where the model is stored in the config file
    # the model must have been saved with the anndata .h5ad format

    model = scvi.model.SCVI.load(scvi_tools_model_path, use_gpu=False)
    adata = model.adata.copy()

    print(adata.obs.columns)
    print('❇️ ❇️ ❇️ ❇️ ❇️ ')
    ### if a file was provided, load intro text from it, if not try to load from adata and finally use the default
    print(intro_text_html)
    if intro_text_html is not None:
        with open(intro_text_html) as file:
            intro_text_html = file.read()
        print(intro_text_html)
    # if no intro text was provided, check the adata.obs.uns['intro_text_html']
    if intro_text_html is None:
        try:
            intro_text_html = adata.obs.uns["intro_text_html"]
            print('Loading introduction text from adata.obs.uns["intro_text_html"] ...')
        except Exception as e:
            print('No custom introduction text provided. Using just the default introduction text...')

    # if the user provided at least one selection column with the argument -s --selection_columns
    if selection_columns:  # if the selection_columns tuple is not empty
        # strips eventual trailing quotes because click's selection is bad
        # NOTE: if this is just trailing, use strip or rstrip with the character
        selection_columns = [item.replace('"', '') for item in selection_columns]
        selection_columns = [item.replace("'", '') for item in selection_columns]
        # check that user provided columns are valid
        for item in selection_columns:
            # raise an error if one of the used provided columns are not in adata.obs.columns
            if item not in adata.obs.columns:
                error_msg = item + """
                is not in adata.obs.columns. Valid options for this adata:

                """ + str(list(adata.obs.columns.values))
                raise Exception(error_msg)
    else:
        # try to load from adata.uns['selection_columns']
        try:
            print("Attempting to load columns from adata.uns['selection_columns']...")
            selection_columns = adata.uns['selection_columns']
            # check that user provided columns are valid
            for item in selection_columns:
                # raise an error if one of the values in adata.uns['selection_columns'] are not in adata.obs.columns
                if item not in adata.obs.columns:
                    error_msg = "Warning: adata.uns['selection_columns'] contains invalid values: " + item + """
                    is not in adata.obs.columns. Defaulting to using all columns. Valid options for this adata:
                    """ + str(list(adata.obs.columns.values))
                    print(error_msg)
                    raise Exception(error_msg)
        except Exception as e:
            print("Attempting to load all columns from adata.obs...")
            selection_columns = list(adata.obs.columns[:-4])

    # crates the selection table using the columns specified in the config file
    selection_table_df = adata.obs.groupby(selection_columns, observed=True).size().rename('#cells').reset_index()
    selection_table_dict = selection_table_df.to_dict(orient='records')
    print(selection_table_df)
    # convert column names into dict for sending as json to datatables
    columns = [{"data": item, "title": item} for item in selection_table_df.columns]

    @tables.route("/", methods=['GET', 'POST'])
    def clientside_table_content():
        return jsonify({"data": selection_table_dict, "columns": columns})

    #### datatables to render the selection tables ####
    @user_configs.route("/", methods=["GET", "POST"])
    def send_user_configs():
        return jsonify({"intro": intro_text_html,
                        "table_selection_columns": selection_columns})
    if selection_columns is None:
        print("SELECTION COLUMN IS NONE👛👛👛👛👛👛")

    app.register_blueprint(tables)
    app.register_blueprint(user_configs)

    # this is the landing page
    @app.route("/", methods=['GET', 'POST'])
    def home():
        return render_template("home.html")

    @app.route('/submit', methods=['POST', 'GET'])
    def receive_submission():
        logger.info('Got a submission!')
        # timestamp = time.strftime("%Y-%m-%d_%H:%M:%S")

        # answer is a dict of json strings containing selected row and column index numbers
        answer = request.form.to_dict(flat=False)
        # need to convert the json strings to dict, then to a data frame
        # data1 is the selection for the first group, data2 for the second
        data1 = json.loads(answer["data1"][0])
        data1_df = pd.DataFrame.from_records(data1)
        data2 = json.loads(answer["data2"][0])
        data2_df = pd.DataFrame.from_records(data2)

        genes = StringIO(json.loads(answer["genes"][0]))
        genes_df = pd.read_csv(genes, names=["selected_genes"])
        jobname = json.loads(answer["jobname"][0])
        print("    🈶 🈶   Received job: ", jobname)

        #### Creates the masks for the selected cell types
        # first create the mask as an array of all false
        # then for each group in the data add them to the mask
        group1_mask = adata.obs.index != adata.obs.index
        for idx, row in data1_df.iterrows():
            # create an array of booleans that starts with all of them true
            current_mask = adata.obs.index == adata.obs.index
            # need to do the _ replace for matching with the pretty names
            for col in selection_columns:
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
            for col in selection_columns:
                # perform successsive AND operations to create the mask
                # that keeps only the cells matching the user selection
                partial_mask = adata.obs[col] == row[col]
                current_mask = current_mask & partial_mask
            group2_mask = group2_mask | current_mask

        # the masks then define the two groups of cells on which to perform DE
        de = model.differential_expression(adata,
                                           idx1=group1_mask,
                                           idx2=group2_mask)

        #### Wrangles the DE results dataframe a bit
        # first we create these variables to customize the hover text in plotly's heatmap
        # the text needs to be arranged in a matrix the same shape as the heatmap
        # try to add gene descriptions and gene names if the adata has those, otherwise add a blank
        # the adata.var fields should include a field named gene_id and gene_name, otherwise they will be filled with a blank default
        try:
            de["gene_description"] = de.index.map(adata.var["gene_description"])
            # for the gene descriptions text, which can be long, we add line breaks
            de["gene_description_html"] = de['gene_description'].str.wrap(80).str.replace('\n', '<br>')
            de["gene_name"] = de.index.map(adata.var['gene_name']).astype(str)
            de["gene_id"] = de.index.astype(str)
            # de['gene_id'] = de.index.map(adata.var['gene_id'])

        except Exception as e:
            de["gene_description_html"] = 'warning: adata.var["gene_description"] does not exist, filled with blank'
            de["gene_name"] = 'warning: adata.var["gene_name"] does not exist, filled with blank'
            de["gene_id"] = 'warning: adata.var["gene_id"] does not exist, filled with blank'

        de["gene_name"] = de["gene_name"].fillna("-")

        # calculate the -log10(p-value) for the volcano
        de["minuslog10pval"] = -np.log10(de["proba_not_de"] + 0.00001)
        de["log10mean_expression"] = np.log10((de["scale1"] + de['scale2']) / 2)

        # all genes are initially colored black
        de["color"] = "black"
        # uncomment line below to color genes by FDR significance
        # de['color'] = de['is_de_fdr_'+str(fdr_target)].map({True:'steelblue',False:'gray'})

        # then we loops through the list of genes provided to color some red
        # gene ids should be a perfect match
        pd.set_option('mode.chained_assignment', None)  # supress warning
        de['color'][de['gene_id'].isin(genes_df['selected_genes'].values)] = "red"

        # gene names should be a partial match
        for partial_string in genes_df['selected_genes'].values:
            de["color"][de["gene_name"].str.contains(partial_string)] = "red"

        #### This makes the volcano plot using plotly
        defig = go.Figure(
                        data=go.Scatter(
                                x=de["lfc_mean"].round(3),
                                y=de["minuslog10pval"].round(3),
                                mode="markers",
                                marker=dict(
                                            color=de['color'],
                                            opacity=0.5),
                                hoverinfo='text',
                                text=de['gene_description_html'],
                                customdata=de.gene_name + '<br>' + de.gene_id,
                                hovertemplate='%{customdata} <br>' + '-log10 p-value: %{y}<br>' + 'Mean log2 fold change: %{x}' + '<extra>%{text}</extra>'),
                                layout={
                                        "title": {"text": str(jobname) + " <br> dashes mark p = 0.01", "x": 0.5},
                                        "xaxis": {"title": {"text": "Mean log fold change"}},
                                        "yaxis": {"title": {"text": "-log10 p-value"}},
                                        "height": 700,
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
                x=de["log10mean_expression"].round(3),
                y=de["lfc_mean"].round(3),
                mode='markers',
                marker=dict(
                    color=-de['minuslog10pval'],
                    colorscale='Viridis',
                    opacity=1),
                hoverinfo='text',
                text=de['gene_description_html'],
                customdata=de.gene_name + '<br>' + de.gene_id + '<br>-log10 p-value: ' + de.minuslog10pval.round(2).astype(str),
                hovertemplate='%{customdata} <br>' + 'log10 mean expression: %{x}<br>' + 'mean log2 fold change: %{y}' + '<extra>%{text}</extra>'
            ),
            layout={
                "title": {"text": str(jobname), "x": 0.5},
                'xaxis': {'title': {"text": "log10 scvi normalized expression"}},
                'yaxis': {'title': {"text": "Mean log2 fold change"}},
                "height": 700,
            }
        )
        mafig.update_layout(hovermode='closest', template='none')
        # overwrites the last figure in order to serve it in the results page
        mafig = mafig.to_html()

        ### creates the results tables in a dataframe that is converted to json
        de_df = de[['gene_name', 'minuslog10pval', 'lfc_mean', 'lfc_std', 'proba_not_de', 'log10mean_expression']]
        print(de_df)
        de_df.index.rename('gene_id', inplace=True)
        de_df = de_df.reset_index()
        de_df = de_df[['gene_id', 'gene_name', 'minuslog10pval', 'lfc_mean', 'log10mean_expression']].fillna('-')
        de_df.columns = ['Gene ID', 'Gene Name', '-log10 p-value', 'mean log2 fold change', 'log 10 mean expression']
        de_df['mean log2 fold change'] = de_df['mean log2 fold change'].astype(float).round(2)
        de_df['-log10 p-value'] = de_df['-log10 p-value'].astype(float).round(2)
        de_df['log 10 mean expression'] = de_df['log 10 mean expression'].astype(float).round(2)

        # convert df to dict for sending as json to datatables
        de_csv_df = de_df.to_csv()
        # convert column names into dict for sending as json to datatables
        columns = [{"data": item, "title": item} for item in de_df.columns]

        return jsonify({'deplothtml': defig, 'maplothtml': mafig, 'decsv': {'data': de_csv_df, 'columns': columns}, 'title': str(jobname)})

    ######### END OF FUNCTION DEFS #########

    print('🔵🔵🔵🔵🔵🔵   GOING TO RUN THE APP NOW    🔵🔵🔵🔵🔵')
    app.run(host=host, port=str(port))

    print('📙 📙 ENDED APP.RUN LOOP 📙📙')
    return app


if __name__ == '__main__':
    print('Starting launch function...')
    launch()
