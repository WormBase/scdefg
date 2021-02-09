/*jslint browser: true*/
/*global $*/
var table1
var table2

$(document).ready(function () {
    $.get('/de_tables/', function (data) {
        table1 = $('#DE_ENRICHED_TABLE').DataTable({
            data: data.data,
            "paging": false,
            "ordering": true,
            "order": [[3, "desc"]],
            "info": true,
            "searching": true,
            dom: 'Bfrtip',
            buttons: ['copy', 'csv', 'excel'],
            scrollY: "70em",
            columns: data.columns,
            select: {
                style: 'multi+shift',
                items: 'cell'
            }
        });
    });


    $.get('/de_tables/', function (data) {
        table2 = $('#DE_DEPLETED_TABLE').DataTable({
            data: data.data,
            "paging": false,
            "ordering": true,
            "order": [[3, "asc"]],
            "info": true,
            "searching": true,
            dom: 'Bfrtip',
            buttons: ['copy', 'csv', 'excel'],
            scrollY: "70em",
            columns: data.columns,
            select: {
                style: 'multi+shift',
                items: 'cell'
            }
        });
    });


});

