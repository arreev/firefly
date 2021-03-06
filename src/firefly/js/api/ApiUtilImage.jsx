/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React from 'react';
import {take,race,call} from 'redux-saga/effects';
import {has} from 'lodash';
import {MouseState} from '../visualize/VisMouseSync.js';
import ImagePlotCntlr, {visRoot, ExpandType} from '../visualize/ImagePlotCntlr.js';
import {primePlot} from '../visualize/PlotViewUtil.js';
import {dispatchAddSaga} from '../core/MasterSaga.js';
import {DefaultApiReadout} from '../visualize/ui/DefaultApiReadout.jsx';
import {reduxFlux} from '../core/ReduxFlux.js';
import {PopupMouseReadoutFull} from '../visualize/ui/PopupMouseReadoutFull.jsx';
import DialogRootContainer from '../ui/DialogRootContainer.jsx';
import {PopupPanel, LayoutType} from '../ui/PopupPanel.jsx';
import {dispatchShowDialog,dispatchHideDialog, isDialogVisible} from '../core/ComponentCntlr.js';
import {readoutRoot,isAutoReadIsLocked, isLockByClick,STANDARD_READOUT} from '../visualize/MouseReadoutCntlr.js';
import {mouseUpdatePromise} from '../visualize/VisMouseSync.js';
import {RangeValues} from '../visualize/RangeValues.js';



const API_READOUT= 'apiReadout';

// NOTE 
// NOTE 
//----------------------------------------------------------------
// Anything that is exported here becomes part of the lowlevel API
// It should have a jsdoc
//----------------------------------------------------------------
// NOTE 
// NOTE

/**
 * @public
 *
 */
export {RangeValues} from '../visualize/RangeValues.js';
export {WPConst, WebPlotRequest, findInvalidWPRKeys, confirmPlotRequest} from '../visualize/WebPlotRequest.js';
export {RequestType} from '../visualize/RequestType';
export {ExpandType, dispatchApiToolsView} from '../visualize/ImagePlotCntlr.js';

export {CysConverter} from '../visualize/CsysConverter.js';
export {CCUtil} from '../visualize/CsysConverter.js';
export {startCoverageWatcher} from '../visualize/saga/CoverageWatcher.js';
export {startImageMetadataWatcher} from '../visualize/saga/ImageMetaDataWatcher.js';
export {getSelectedRegion} from '../drawingLayers/RegionPlot.js';

export {extensionAdd, extensionRemove} from '../core/ExternalAccessUtils.js';
export {makeWorldPt, makeScreenPt, makeImagePt, parsePt} from '../visualize/Point.js';

export {IMAGE, NewPlotMode} from '../visualize/MultiViewCntlr';


/**
 * Get plot object with the given plot id, when plotId is not included, active plot is returned.
 * @param {string} [plotId] the plotId, optional
 * @returns {WebPlot}
 * @function getPrimePlot
 * @public
 * @memberof firefly.util.image
 *
 */
export function getPrimePlot(plotId) {
    return primePlot(visRoot(), plotId);
}



var isInit= false;
/**
 * Set a defaults object on for a draw layer type.
 * The following draw layers are supported: 'ACTIVE_TARGET_TYPE', 'CATALOG_TYPE'
 * @param {string} drawLayerTypeId
 * @param {DrawingDef} defaults
 * @public
 * @function setDrawLayerDefaults
 * @memberof firefly.util.image
 *
 * @example
 * firefly.util.image.setDrawLayerDefaults('ACTIVE_TARGET_TYPE', {symbol:'x', color:'pink', size:15});
 * or
 * firefly util.image.setDrawLayerDefaults('CATALOG_TYPE', {symbol:'cross', color:'red'});
 *
 * @see DrawingDef
 * @see DrawSymbol
 */
export function setDrawLayerDefaults(drawLayerTypeId, defaults) {
    reduxFlux.setDrawLayerDefaults(drawLayerTypeId, defaults);
}

/**
 * @summary  initialize the auto readout. Must be call once at the begging to get the popup readout running.
 * @param {object} ReadoutComponent - either a PopupMouseReadoutMinimal or PopupMouseReadoutFull
 * @param {object} props - a list of the properties
 * @public
 * @function initAutoReadout
 * @memberof firefly.util.image
 */
export function initAutoReadout(ReadoutComponent= DefaultApiReadout,
      props={MouseReadoutComponent:PopupMouseReadoutFull, showThumb:true,showMag:true } ){
    if (isInit) return;
    dispatchAddSaga(autoReadoutVisibility, {ReadoutComponent,props});
    isInit= true;
}


/**
 *
 * @param stretchType the type of stretch may be 'Percent', 'Absolute', 'Sigma'
 * @param lowerValue lower value of stretch, based on stretchType
 * @param upperValue upper value of stretch, based on stretchType
 * @param algorithm the stretch algorithm to use, may be 'Linear', 'Log', 'LogLog', 'Equal', 'Squared', 'Sqrt'
 * @public
 * @function serializeSimpleRangeValues
 * @memberof firefly.util.image
 */
export function serializeSimpleRangeValues(stretchType,lowerValue,upperValue,algorithm) {
    const rv= RangeValues.makeSimple(stretchType,lowerValue,upperValue,algorithm);
    return rv.toJSON();
}


//========== Private ================================
//========== Private ================================
//========== Private ================================
//========== Private ================================

const delay = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function *autoReadoutVisibility({ReadoutComponent,props}) {
    var inDialog= false;
    var showing;
    var mouseState;
    // var action;
    var doYield= true;
    var winner;
    while (true) {
        // if (doYield) action= yield take([MOUSE_STATE_CHANGE, ImagePlotCntlr.CHANGE_EXPANDED_MODE]);

        if (doYield) {
            winner = yield race({
                expandedChange: take([ImagePlotCntlr.CHANGE_EXPANDED_MODE]),
                mouse: call(mouseUpdatePromise)
            });
        }


        doYield= true;
        if (winner.mouse) {
            mouseState = winner.mouse.mouseState;
        }
        if (visRoot().expandedMode!==ExpandType.COLLAPSE) {
            hideReadout();
            continue;
        }
        if (!mouseState) continue;
        showing= isDialogVisible(API_READOUT);
        if (mouseState!==MouseState.EXIT && !showing) {
            showReadout(ReadoutComponent,props, (inD) => {
                inDialog= inD;
            });
        }
        else if (mouseState===MouseState.EXIT && showing) {
             winner = yield race({
                           expandedChange: take([ImagePlotCntlr.CHANGE_EXPANDED_MODE]),
                           mouse: call(mouseUpdatePromise),
                           timer: call(delay, 3000)
                         });
            if ((has(winner,'timer')) && !inDialog && !isLockByClick(readoutRoot()) && !isAutoReadIsLocked(readoutRoot())) {
                hideReadout();
            }
            else {
                doYield= false;
            }
        }
    }
}



function showReadout(DefaultApiReadout, props={},  mouseInDialog) {

    const  readout=readoutRoot();
    const title = readout[STANDARD_READOUT] && readout[STANDARD_READOUT].readoutItems  && readout[STANDARD_READOUT].readoutItems.title?readout[STANDARD_READOUT].readoutItems.title.value:'';
    const popup= (
        <PopupPanel title={title} layoutPosition={LayoutType.TOP_RIGHT} mouseInDialog={mouseInDialog} >
            <DefaultApiReadout {...props} />
        </PopupPanel>
    );
    DialogRootContainer.defineDialog(API_READOUT, popup);
    dispatchShowDialog(API_READOUT);
}

function hideReadout() {
    if (isDialogVisible(API_READOUT)) dispatchHideDialog(API_READOUT);
}

