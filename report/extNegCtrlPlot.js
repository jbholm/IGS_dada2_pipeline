function renderExtNegCtrlPlot (sample) {
  var xTickPos = [${ xs }];
  var ys = ${ ys };
  var sems = [${ sems }];
  var cols = [${ colors }];
  var taxaLabs = [${ taxa }];
  var asvLabs = [${ asvs }];

  renderCtrPlot('extNegCtrlPlot', ys, sems, cols, sample, taxaLabs, asvLabs);

};
renderExtNegCtrlPlot(0);

$('#extNegCtrlSelect').on('change', function(e) {
  renderExtNegCtrlPlot($('#extNegCtrlSelect').prop('selectedIndex'));
});
