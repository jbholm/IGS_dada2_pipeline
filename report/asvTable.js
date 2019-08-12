
(function () {
    var x = [];
    var y = [];
    var z = [[${ zmin }, ${ zmax }]];

    var trace1 = {
        x: [0, 0],
        y: [0, 0],
        mode: 'markers',
        marker: {
            size: 0,
            color: [0, 0.931],
            colorscale: 'Viridis',
            colorbar: {
                title: 'Relative abundance',
                x: -0.1 // Also shift y in pos. direction, dependent on nSamples
            },
        }

    };
    var data = [trace1];


    var layout = {
        width: 150,
        height: 500,
        margin: {
            l: 0,
            b: 250,
            t: 80
        },
        xaxis: {
            color: "#FFF",
            range: [1, 2],
            tickfont: { color: "#FFF" },
            fixedrange: true

        },
        yaxis: {
            color: "#FFF",
            range: [1, 2],
            tickfont: { color: "#FFF" },
            fixedrange: true

        }
    };
    Plotly.newPlot('heatmapScale-${hmNbr}', data, layout, { displaylogo: false, responsive: true, displayModeBar: false });
})();


(function () {
    var select = function (arr, indexes) {
        var ans = new Array(indexes.length);
        for (var i = 0; i < indexes.length; i++) {
            ans[i] = arr[indexes[i]]
        }
        return arr;
    }

    var divname = 'heatmapInner-${hmNbr}';

    var zValues = ${ abundances };

    var xBotValues = [${ taxa }];
    var xTopValues = [${ asvs }];

    var yValues = [
        ${ samples }
    ];


    var data = ${ data };

    var xTickPos = [];

    for (var i = 0; i < xBotValues.length; i++) {
        xTickPos.push(i);
    }

    var layout = {
        //annotations: [],
        width: Math.ceil(${ nCol } * 14.485 + 277.57),
        height: 500,
        margin: {
            b: 250,
        },
        yaxis: {
            title: "<b>Sample</b>",
            autorange: 'reversed', // Because plotly uses Cartesian coordinates
            showgrid: false,
            type: 'category',
            automargin: true,
            fixedrange: true,
            linecolor: '#333',
            linewidth: 1,
        },
        ${ xAxes }
        };
Plotly.newPlot(divname, data, layout, { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['sendDataToCloud', 'toggleSpikelines', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'zoom2d', 'pan2d', 'select2d', 'lasso2d'] });

var d = document.getElementById(divname);
d.highlight = makeHighlight(zValues, xTopValues, "#" + divname);
d.highlightMax = highlightMax(zValues);
d.highlightElems = [];

d.on('plotly_click', function (e) {
    point = e.points[0];
    d.highlight({ sample: point.pointIndex[0], asvID: point.pointIndex[1], type: "index" });
});

})();
