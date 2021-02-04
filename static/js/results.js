/*jslint browser: true*/
/*global $*/
var table1
var table2

$(document).ready(function () {
    $.get('/de_results_tables/', function (data) {
        table1 = $('#DE_RESULTS_TABLE').DataTable({
            data: data.data,
            "paging": true,
            "ordering": true,
            "info": true,
            "searching": true,
            "pageLength": 25,
            dom: 'Bfrtip',
            buttons: ['copy', 'csv', 'excel', 'pdf', 'print'],
            scrollY: "70em",
            scrollX: true,
            scrollCollapse: true,
            fixedColumns: {
                leftColumns: 1,
                rightColumns: 0
            },
            columns: data.columns,
            select: {
                style: 'multi+shift',
                items: 'cell',
                selector: 'td:not(:first-child)'
            },
            "columnDefs": [
                {"width": "20em", "targets": 0}
            ],
        });
    });


    $.get('/tables/', function (data) {
        table2 = $('#SELECTED_GROUPS_TABLE').DataTable({
            data: data.data,
            "paging": true,
            "ordering": true,
            "info": true,
            "searching": false,
            "pageLength": 20,
            dom: 'frtip',
            scrollY: "50em",
            scrollX: true,
            scrollCollapse: true,
            fixedColumns: {
                leftColumns: 1,
                rightColumns: 0
            },
            columns: data.columns,
            select: {
                style: 'multi+shift',
                items: 'cell',
                selector: 'td:not(:first-child)'
            },
            "columnDefs": [
                {"width": "20em", "targets": 0}
            ],
        });
    });


});

