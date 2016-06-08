/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React from 'react';
import {take,race,call} from 'redux-saga/effects';
import {MouseState} from '../visualize/VisMouseSync.js';
import ImagePlotCntlr, {visRoot, ExpandType} from '../visualize/ImagePlotCntlr.js';
import {dispatchAddSaga} from '../core/MasterSaga.js';
import  {DefaultApiReadout} from '../visualize/ui/DefaultApiReadout.jsx';
//import  {PopupMouseReadoutMinimal} from '../visualize/ui/PopupMouseReadoutMinimal.jsx';
import  {PopupMouseReadoutFull} from '../visualize/ui/PopupMouseReadoutFull.jsx';
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

export {RangeValues} from '../visualize/RangeValues.js';
export {WPConst, WebPlotRequest, findInvalidWPRKeys, confirmPlotRequest} from '../visualize/WebPlotRequest.js';
export {RequestType} from '../visualize/RequestType';
export {ExpandType, dispatchApiToolsView} from '../visualize/ImagePlotCntlr.js';


/**
 * initialize the auto readout. Must be call once at the begging to get the popup readout running.
 * @param ReadoutComponent
 * @param props
 */
export function initAutoReadout(ReadoutComponent= DefaultApiReadout,
         //   props={MouseReadoutComponent:PopupMouseReadoutMinimal, showThumb:false,showMag:false}){
      props={MouseReadoutComponent:PopupMouseReadoutFull, showThumb:false,showMag:false } ){


    dispatchAddSaga(autoReadoutVisibility, {ReadoutComponent,props});
}


/**
 *
 * @param stretchType the type of stretch may be 'Percent', 'Absolute', 'Sigma'
 * @param lowerValue lower value of stretch, based on stretchType
 * @param upperValue upper value of stretch, based on stretchType
 * @param algorithm the stretch algorithm to use, may be 'Linear', 'Log', 'LogLog', 'Equal', 'Squared', 'Sqrt'
 */
export function serializeSimpleRangeValues(stretchType,lowerValue,upperValue,algorithm) {
    const rv= RangeValues.makeSimple(stretchType,lowerValue,upperValue,algorithm);
    return rv.serialize();
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
            if ((!winner.expandedChange || !winner.mouse) && !inDialog && !isLockByClick(readoutRoot()) && !isAutoReadIsLocked(readoutRoot())) {
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
