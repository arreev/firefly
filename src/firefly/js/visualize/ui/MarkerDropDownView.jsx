/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React, {PropTypes} from 'react';
import {SingleColumnMenu} from '../../ui/DropDownMenu.jsx';
import MarkerTool from '../../drawingLayers/MarkerTool.js';
import {ToolbarButton,
        DropDownVerticalSeparator} from '../../ui/ToolbarButton.jsx';
import {dispatchCreateDrawLayer,
    dispatchAttachMarkerLayerToPlot} from '../DrawLayerCntlr.js';
import {visRoot} from '../ImagePlotCntlr.js';
import {has, get} from 'lodash';

var idCnt=0;

const markerItem = {
    marker: {label: 'Marker', drawLayerType: MarkerTool.TYPE_ID}
};

function displayItemText(itemName) {
    if (has(markerItem, itemName)) {
        return `Add ${markerItem[itemName].label}`;
    }
}

function addNewDrawLayer(pv, itemName) {
    if (!has(markerItem, itemName)) return;

    var plotId = pv ? pv.plotId : get(visRoot(), 'activePlotId');  // pv should be true, otherwise marker is disabled.
    var drawLayerId = `${markerItem[itemName].label}-${idCnt++}`;

    if (plotId) {
        dispatchCreateDrawLayer(markerItem[itemName].drawLayerType, {Title: drawLayerId, drawLayerId});
        dispatchAttachMarkerLayerToPlot(drawLayerId, plotId, true);
    }
 }

export function MarkerDropDownView({plotView:pv}) {
    var enabled = !!pv;

    return (
        <SingleColumnMenu>
            <ToolbarButton text={displayItemText('marker')}
                           enabled={enabled} horizontal={false}
                           onClick={()=> addNewDrawLayer(pv, 'marker')}/>
            <DropDownVerticalSeparator/>
        </SingleColumnMenu>
    );
}

MarkerDropDownView.propTypes= {
    plotView : PropTypes.object
};
