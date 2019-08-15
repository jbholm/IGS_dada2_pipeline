function renderPcrNegCtrlPlot (sample) {
  var xTickPos = [${ xs }];
  var ys = ${ ys };
  var sems = [${ sems }];
  var cols = [${ colors }];
  var taxaLabs = [${ taxa }];
  var asvLabs = [${ asvs }];

  renderCtrPlot('pcrNegCtrlPlot', ys, sems, cols, sample, taxaLabs, asvLabs);

};
renderPcrNegCtrlPlot(0);

$('#pcrNegCtrlSelect').on('change', function(e) {
  renderPcrNegCtrlPlot($('#pcrNegCtrlSelect').prop('selectedIndex'));
});
