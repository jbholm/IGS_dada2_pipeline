function renderCtrPlot(data, sems, cols, sample, taxaLabs, asvLabs) {
  var select = function (arr, indexes) {
    var ans = new Array(indexes.length);
    for (var i = 0; i < indexes.length; i++) {
        ans[i] = arr[indexes[i]]
    }
    return ans;
}

  
    var xTickPos = Array(length(data));
    for (var i = 1; i <= length(xTickPos); i++) {
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
    console.log(min);
    var diffs = Array(ysFlat.length - 1);
    for (var i = 1; i < ysFlat.length; i++)  diffs[i - 1] = Math.abs((ysFlat[i] - ysFlat[i - 1]))
    console.log(diffs);
    var minDiff = Math.min(...diffs);
    console.log(minDiff);
  ysFlat = ysFlat.map(i => i + Math.min(min, minDiff) / 10);
  
  console.log(xTickPos.map(val => Array(ys[0].length).fill(val)).reduce((prevVal, curVal) => prevVal.concat(curVal)));
  console.log(ysFlat);
  
  Plotly.newPlot('negCtrlBar', {
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
