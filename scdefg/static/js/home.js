/*jslint browser: true*/
/*global $*/
var table1
var table2

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
        table1 = $('#FIRST_TABLE').DataTable({
            data: data.data,
            dom: 'Plfrtip',
            scrollY: "30em",
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
    });

    $.get('/tables/', function (data) {
        table2 = $('#SECOND_TABLE').DataTable({
            data: data.data,
            dom: 'Plfrtip',
            scrollY: "30em",
            scrollX: true,
            scrollCollapse: true,
            fixedColumns: {
                leftColumns: 1,
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
    });

    $('buttonbutton').click(function () {
        console.log(table1.rows('.selected').data())
        alert(table1.rows('.selected').data() + ' row(s) selected');
    });
    $('button').click(function () {
        var data1 = table1.rows(['.selected']).data().toArray();
        var json1 = JSON.stringify(data1);
        var data2 = table2.rows(['.selected']).data().toArray();
        var json2 = JSON.stringify(data2);
        var nrows1 = table1.rows('.selected').data().length
        var nrows2 = table2.rows('.selected').data().length
        var form_data = $('form').serializeArray()
        var genes = form_data[0].value
        console.log(genes)

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
                $("#button-div").hide();
                $("#spinner-div").show();

                $.post("/submit", {
                    // "contentType": "application/json",
                    'data1': json1,
                    'data2': json2,
                    'genes': json_genes,
                })
                    .done(function (data) {
                        $("#spinner-div").hide();
                        $("#newjob-div").show();
                        $("#results-div").show();
                        $("#de-plot-display-div").html(data.deplothtml)
                        $("#ma-plot-display-div").html(data.maplothtml)
                        table3 = $('#DE_ENRICHED_TABLE').DataTable({
                            data: data.dejsondata.data,
                            "paging": true,
                            "ordering": true,
                            "lengthMenu": [[100, 250, 500, -1], [100, 250, 500, "All"]],
                            "order": [[3, "desc"]],
                            "info": true,
                            "searching": true,
                            dom: 'Bfrtipl',
                            buttons: [
                                {
                                    extend: 'csvHtml5',
                                    title: data.title,
                                    text: 'Download csv',
                                },
                                {
                                    extend: 'excelHtml5',
                                    title: data.title,
                                    text: 'Download Excel',
                                },
                                {
                                    extend: 'copyHtml5',
                                    text: 'Copy all',
                                    title: data.title,
                                },
                                {
                                    extend: 'copyHtml5',
                                    text: 'Copy current page',
                                    title: data.title,
                                    exportOptions: {
                                        modifier: {
                                            page: 'current'
                                        }
                                    }
                                }],
                            scrollY: "70em",
                            columns: data.dejsondata.columns,
                            // select: {
                            //     style: 'multi+shift',
                            //     items: 'cell'
                            // }
                        });

                        table4 = $('#DE_DEPLETED_TABLE').DataTable({
                            data: data.dejsondata.data,
                            "paging": true,
                            "ordering": true,
                            "lengthMenu": [[100, 250, 500, -1], [100, 250, 500, "All"]],
                            "order": [[3, "asc"]],
                            "info": true,
                            "searching": true,
                            dom: 'Bfrtipl',
                            buttons: [
                                {
                                    extend: 'csvHtml5',
                                    title: data.title,
                                    text: 'Download csv',
                                },
                                {
                                    extend: 'excelHtml5',
                                    title: data.title,
                                    text: 'Download Excel',
                                },
                                {
                                    extend: 'copyHtml5',
                                    text: 'Copy all',
                                    title: data.title,
                                },
                                {
                                    extend: 'copyHtml5',
                                    text: 'Copy current page',
                                    title: data.title,
                                    exportOptions: {
                                        modifier: {
                                            page: 'current'
                                        }
                                    }
                                }],
                            scrollY: "70em",
                            columns: data.dejsondata.columns,
                            // select: {
                            //     style: 'multi+shift',
                            //     items: 'cell'
                            // }
                        });

                    })
                    .fail(function () {
                            alert("Something went wrong. Refresh the page and try again. If it keeps happening email eduardo@wormbase.org")
                        }
                    );
            }
        }
        ;
    });


});

