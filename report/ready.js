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
