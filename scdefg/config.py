### model_file should be the path to the saved scvi model
model_file='../model'

# these are the columns that will be used to slow filtering for viewing only desired groups prior to selection
filter_columns=['experiment_code','tissue']

#table_selection_columns should be the solumns in adata.obs for which you want to stratify groups
# for for example it could be table_selection_columns=['batch','cell_type']
# so that users may select cells of a given time for specific batches, not all at once
table_selection_columns=['experiment_code', 'tissue']

introduction_html_text="""
<h3> Differential expression on single cell RNA sequencing data</h3>

<ul>
    <li>This app allows you to quickly perform differential expression on single cell RNA
        sequencing data.
    </li>
    <li>Just select cell types and experiments to compare and some genes to
        highlight.
    </li>
    <li>Hold shift to select multiple rows.</li>
    <li>Results take about 15s to compute.</li>
    <li>You will get an interactive volcano plot, MA plot and tables with DE values.</li>

    <li>This app uses <a href="https://scvi-tools.org">scvi-tools</a>, code
        is available on
        <a href="https://github.com/Munfred/scdefg/">GitHub.</a>
</ul>
<p>
To highlight genes add one gene per line below. 
Gene names (<samp>bus-1</samp>) or ID are accepted (<samp>WBGene00018223</samp>)
</p>
"""