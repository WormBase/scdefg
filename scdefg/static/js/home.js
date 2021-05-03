/*jslint browser: true*/
/*global $*/
var table1
var table2

const getInputTableContent = (data) => (
    {
        data: data.data,
        dom: 'Plfrtip',
        lengthMenu: [[15, 30, 50, 100, -1], [15, 30, 50, 100, "All"]],
        scrollY: "70em",
        scrollX: true,
        scrollCollapse: true,
        fixedColumns: {
            leftColumns: 0,
            rightColumns: 0
        },
        columnDefs: [{
            orderable: false,
            className: 'select-checkbox',
            targets: -1
        }],
        columns: data.columns,
        select: {
            style: 'multi+shift',
        },
    });

const convertResultCSVToJson = csvString => {
    let lines = csvString.split('\n');
    let headers = lines.shift().split(',');
    headers[0] = "rowNum";
    let allLines = lines.map(line => {
        let linearr = line.split(',');
        let singeLine = Object.assign({}, ...headers.map((header, idx) => {
            let headerObj = {};
            headerObj[header] = linearr[idx];
            return headerObj;
        }));
        if (Object.values(singeLine).every(value => value !== undefined)) {
            return singeLine;
        }
    });
    return allLines.filter(line => line !== undefined);
}

const getOutputTableContent = (data, columns, title, asc = true) => ({
    data: data,
    "paging": true,
    "ordering": true,
    "lengthMenu": [[100, 250, 500, -1], [100, 250, 500, "All"]],
    "order": [[3, asc ? "asc" : "desc"]],
    "info": true,
    "searching": true,
    dom: 'Bfrtipl',
    buttons: [
        {
            extend: 'csvHtml5',
            title: title,
            text: 'Download csv',
        },
        {
            extend: 'excelHtml5',
            title: title,
            text: 'Download Excel',
        },
        {
            extend: 'copyHtml5',
            text: 'Copy all',
            title: title,
        },
        {
            extend: 'copyHtml5',
            text: 'Copy current page',
            title: title,
            exportOptions: {
                modifier: {
                    page: 'current'
                }
            }
        }],
    scrollY: "70em",
    columns: columns,
    // select: {
    //     style: 'multi+shift',
    //     items: 'cell'
    // }
});

$(document).ready(function () {
    $("#results-div").hide();
    $("#spinner-div").hide();
    $("#newjob-div").hide();
    $.get('/user_configs/', function (user_configs) {
        if (user_configs.intro != "default") {
            $("#intro-div").html(user_configs.intro)
        }
    });

    $.get('/tables/', function (data) {
        table1 = $('#FIRST_TABLE').DataTable(getInputTableContent(data));
        table2 = $('#SECOND_TABLE').DataTable(getInputTableContent(data));
    });

    $('buttonbutton').click(function () {
        console.log(table1.rows('.selected').data())
        alert(table1.rows('.selected').data() + ' row(s) selected');
    });
    $('#submitButton').click(function () {
        var data1 = table1.rows(['.selected']).data().toArray();
        var json1 = JSON.stringify(data1);
        var data2 = table2.rows(['.selected']).data().toArray();
        var json2 = JSON.stringify(data2);
        var nrows1 = table1.rows('.selected').data().length
        var nrows2 = table2.rows('.selected').data().length
        var genes = $('#genesList').val()
        var jobname = $('#jobname').val()
        console.log(genes)
        // console.log(jobname)
        if (table1.rows('.selected').data().length == 0) {
            alert(' You did not select any cells for group 1')
        }
        if (table2.rows('.selected').data().length == 0) {
            alert(' You did not select any cells for group 2')
        }

        if (nrows1 != 0 && nrows2 != 0) {

            var confirmation = true;

            if (confirmation == true) {
                json_genes = JSON.stringify(genes);
                json_jobname = JSON.stringify(jobname);
                $("#button-div").hide();
                $("#spinner-div").show();

                $.post("/submit", {
                    // "contentType": "application/json",
                    'data1': json1,
                    'data2': json2,
                    'genes': json_genes,
                    'jobname': json_jobname
                })
                    .done(function (data) {
                        $("#spinner-div").hide();
                        $("#newjob-div").show();
                        $("#results-div").show();
                        $("#de-plot-display-div").html(data.deplothtml)
                        $("#ma-plot-display-div").html(data.maplothtml)
                        table3 = $('#DE_ENRICHED_TABLE').DataTable(getOutputTableContent(
                            convertResultCSVToJson(data.decsv.data), data.decsv.columns, data.title, false));
                        table4 = $('#DE_DEPLETED_TABLE').DataTable(getOutputTableContent(
                            convertResultCSVToJson(data.decsv.data), data.decsv.columns, data.title, true));
                        document.getElementById("results-div").scrollIntoView();
                    })
                    .fail(function () {
                            alert("Something went wrong. Refresh the page and try again. If it keeps happening email eduardo@wormbase.org")
                        });
            }
        }
    });
});

