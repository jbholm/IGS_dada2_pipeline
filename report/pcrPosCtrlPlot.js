function renderPcrPosCtrlPlot (sample) {
  var xTickPos = [${ xs }];
  var ys = ${ ys };
  var sems = [${ sems }];
  var cols = [${ colors }];
  var taxaLabs = [${ taxa }];
  var asvLabs = [${ asvs }];

  renderCtrPlot('pcrPosCtrlPlot', ys, sems, cols, sample, taxaLabs, asvLabs);

};
renderPcrPosCtrlPlot(0);

$('#pcrPosCtrlSelect').on('change', function(e) {
  renderPcrPosCtrlPlot($('#pcrPosCtrlSelect').prop('selectedIndex'));
});
