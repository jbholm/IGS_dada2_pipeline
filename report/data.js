
var makeHighlight = function (zValues, xTopValues, divSelector) {
    z = zValues;
    x = xTopValues;

    // FUNCTION FACTORY. We create an enclosing environment to which the returned
    // function retains access, so each plot can have its own tailored self.highlight.
    var xtick = xTopValues.length > 0 ? "x2tick" : "xtick";
    var nrow = zValues.length;
    var ncol = zValues[0].length;
    var canvas = $(divSelector + " .nsewdrag");
    var svg = $(divSelector + " .main-svg")[0];
    var container = $(divSelector + " .plot-container.plotly");
    
    return function ({ sample, asvID, scroll = false, type = "names" } = {}) {
        removeHighlight(this.highlightElems);
        
        var height = parseInt(canvas.attr("height")) / nrow;
        var width = parseInt(canvas.attr("width")) / ncol;

        var x;
        var y;
        // if (type == 'names') {
        //     var xtick = $("#all-samples .heatmap.plotly-container ." + xtick + " > text[data-unformatted='" + asvID + "']");
        //     var ytick = $("#all-samples .heatmap.plotly-container .subplot.xy .yaxislayer-above .ytick > text[data-unformatted='" + sample + "']");
            
        //     var xtickdata = parseTf(xtick.attr("transform"));
        //     var ytickdata = parseTf(ytick.attr("transform"));
    
        //     x = xtickdata['translate'][0];
        //     y = ytickdata['translate'][1];
        // } else
            if (type == 'index') {
            x = parseInt(asvID) * width + 0.5 * width + parseInt(canvas.attr('x'));
            y = parseInt(sample) * height + 0.5 * height + parseInt(canvas.attr('y'));
        }
    
        var highlight = document.createElementNS("http://www.w3.org/2000/svg", "rect");
        highlight.setAttributeNS(null, "id", "sample-asv-highlighter-black");
        highlight.setAttributeNS(null, "x", x - Math.round(width / 2));
        highlight.setAttributeNS(null, "y", y - Math.round(height / 2));
        highlight.setAttributeNS(null, "height", height);
        highlight.setAttributeNS(null, "width", width);
        highlight.setAttributeNS(null, "fill", "transparent");
        highlight.setAttributeNS(null, "stroke", "rgb(0, 0, 0)");
        highlight.setAttributeNS(null, "stroke-width", "1");
    
        var highlightblack = highlight.cloneNode(true);
        highlightblack.setAttributeNS(null, "stroke", "rgb(255, 255, 255)");
        highlightblack.setAttributeNS(null, "stroke-width", "3");
        highlightblack.setAttributeNS(null, "id", "sample-asv-highlighter-white");
        svg.appendChild(highlightblack);
        svg.appendChild(highlight);
    
        var use = document.createElementNS("http://www.w3.org/2000/svg", "use");
        use.setAttributeNS('http://www.w3.org/1999/xlink', "xlink:href", "#sample-asv-highlighter-black");
        use.id = "sample-asv-highlighter-use";
        svg.appendChild(use);
        if (scroll) {
            container.scrollLeft(x - Math.round(width / 2) - $(window).width() / 2);
            container.scrollTop(y - Math.round(height / 2) - $(window).height() / 2);
        }

        this.highlightElems = [use, highlight, highlightblack];
        return false;
    };
};
var highlightMax = function (zValues) {
    
    var lookup = zValues.map(row => row.reduce((iMax, x, i, arr) => x > arr[iMax] ? i : iMax, 0) );
    return function(sampleI) {
        console.log(lookup);
        this.highlight({ sample: sampleI, asvID: lookup[sampleI], scroll: true, type: "index" });
    };
};
var asvHeatmapMaxHighlight = function (sampleI, divname) {
    var d = document.getElementById(divname);
    d.highlightMax(sampleI);
};
    

    ${asvtablesjs}


    $(".heatmap").click(function(e) {
        var x = e.offsetX;
        var y = e.offsetY;
        var canvas = $('#all-samples .heatmap.plotly-container rect.nsewdrag');
        if(x > canvas.attr('X') && x < parseInt(canvas.attr('X')) + parseInt(canvas.attr("width")) && (y > canvas.attr('Y')) && y < parseInt(canvas.attr('Y')) + parseInt(canvas.attr("height"))) {
        } else {
            var d = document.getElementById('heatmapInner-1');
            removeHighlight(d.highlightElems);
        }

    })

    $("div.plotly-container > div").css({
        'height': '100%', // Should take up the height of its parent (which fills the screen). 
        // same for width, but we defined this property in the <head> stylesheet
        'overflow': 'auto' // Causes the div immediately outside svg to have scrollbars automatically
    });
    var bodyFontSize = window.getComputedStyle(document.body, null).fontSize;
    bodyFontSize = parseInt(bodyFontSize, 10) * 1.35;

    var fillRed = function (n, v) {
        var styleArr = v.split(';');
        var newstyle = [];
        for (i = 0; i < styleArr.length; i++) {
            if (styleArr[i].split(':')[0].toLowerCase().trim() != 'fill') {
                newstyle.push(styleArr[i]);
            } else {
                newstyle.push('fill: #f05f5e');
            }
        }
        return newstyle.join(';');
    }
    var contaminants = [${ contaminants }];

    $("#all-samples .x2tick > text").filter(function () {
        return contaminants.indexOf($(this).data('unformatted')) > -1
    }).attr("style", fillRed);



function parseTf (a)
        {
            var b={};
            for (var i in a = a.match(/(\w+\((\-?\d+\.?\d*e?\-?\d*,?)+\))+/g))
            {
                var c = a[i].match(/[\w\.\-]+/g);
                b[c.shift()] = c;
            }
            return b;
        }
        function removeHighlight(removables) {
            removables.map(e => e.remove());
    }
    
    function renderCtrPlot(div, data, sems, cols, sample, taxaLabs, asvLabs) {
        var select = function (arr, indexes) {
        var ans = new Array(indexes.length);
        for (var i = 0; i < indexes.length; i++) {
            ans[i] = arr[indexes[i]]
        }
        return ans;
        }
        
        
        var xTickPos = Array(data.length);
        for (var i = 1; i <= xTickPos.length; i++) {
            xTickPos[i - 1] = i;
        }
          var ys = data;
        var type;
        var semVisible;
        var dotColor;
        if(sample == 0) {
          type = "box";
          semVisible = true;
        } else {
          type = "scatter"
          semVisible = false;
          ys = ys.map(arr => [arr[sample - 1]]);
          dotColor = cols[sample - 1];
        }
        console.log(ys);
        ysFlat = ys.reduce((prevVal, curVal) => prevVal.concat(curVal));
        var dsorted = ysFlat.slice(0).sort((a, b) => a - b);
        var min = 0;
        while(min == 0) {
            min = dsorted.shift();
        }
        var diffs = Array(ysFlat.length - 1);
        for (var i = 1; i < ysFlat.length; i++)  diffs[i - 1] = Math.abs((ysFlat[i] - ysFlat[i - 1]))
        var minDiff = Math.min(...diffs);
        ysFlat = ysFlat.map(i => i + Math.min(min, minDiff) / 10);
        
        console.log(xTickPos.map(val => Array(ys[0].length).fill(val)).reduce((prevVal, curVal) => prevVal.concat(curVal)));
        console.log(ysFlat);
        
        Plotly.newPlot(div, {
        data: [{
        
        "mode": "markers",
          "type": type,
          "marker": {
            size: 8,
            color: dotColor,
        },
        "x": xTickPos.map(val => Array(ys[0].length).fill(val)).reduce((prevVal, curVal) => prevVal.concat(curVal)),
        "y": ysFlat,
        "error_y": {
        type: 'data',
        array: sems,
        visible: semVisible,
        },
        "orientation": "v",
          "showscale": false,
          boxpoints: 'Outliers',
          jitter: 0.3,
          "hoverinfo": "y+z",
          "line": {
            "color": "#D55E00",
            "width": 1,
        },
        }, {
        
        "type": "box",
        "x": asvLabs,
        "y": Array(xTickPos.length).fill(Math.min(min, minDiff) / 10),
        "xaxis": "x2",
        "orientation": "v",
        "hoverinfo": "none",
        }],
        layout: {
        "xaxis": {
          "ticktext": taxaLabs,
          "tickvals": xTickPos,
          tickfont: {color: "#333"},
        "type": "category",
        "range": [-0.5, xTickPos.length + 0.5],
          "domain": [0, 1],
          fixedrange: true,
          title: "<b>Taxa</b>",
          "overlaying": "x2",
          linecolor: '#333',
          linewidth: 1,
        },
        "yaxis": {
        "type": "log",
        title: "<b>Abundance</b>",
        "domain": [0, 1],
        fixedrange: true,
        linecolor: '#333',
        linewidth: 1,
        },
        "xaxis2": {
        "side": "top",
        title: "<b>ASVs</b>",
        tickfont: {color: "#333"},
        "type": "category",
        "range": [-0.5, xTickPos.length + 0.5],
        "anchor": "y",
        fixedrange: true,
        },
        width: Math.ceil( 20.013 * xTickPos.length + 259.91),
        // 1324 bars: 26758 wide. 7 bars: 400 wide (bars are on the wide side).
        height: 500,
        margin: {
          b: 250,
        },
        showlegend: false,
        
        }, config: { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['sendDataToCloud', 'toggleSpikelines', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'zoom2d', 'pan2d', 'select2d', 'lasso2d'] }
        });
        }
        
        ${ pcrNegCtrlJs }
        ${pcrPosCtrlJs}
        ${ extNegCtrlJs }

        
    