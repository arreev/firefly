/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React from 'react';
import PropTypes from 'prop-types';
import {ThumbnailView} from './ThumbnailView.jsx';
import {MagnifiedView} from './MagnifiedView.jsx';
import {getActivePlotView, getPlotViewById} from '../PlotViewUtil.js';
import {MouseReadout} from './MouseReadout.jsx';
import InputFieldLabel from '../../ui/InputFieldLabel.jsx';


/**
 *
 * @param visRoot visualization store root
 * @param {object} currMouseState  the current state of the mouse
 * @return {XML}
 */
export function VisHeaderView({readout}) {

    var rS= {
        width: 700,
        minWidth:660,
        height: 32,
        minHeight:32,
        display: 'inline-block',
        position: 'relative',
        verticalAlign: 'top',
        cursor:'pointer',
        whiteSpace : 'nowrap',
        overflow : 'hidden'
    };

    return (
        <div style={{display:'inline-block', float:'right', whiteSpace:'nowrap'}}>
            <div style={rS}>
                <div style={{position:'absolute', color:'white'}}>
                    <MouseReadout readout={readout}/>
                </div>
            </div>
        </div>
    );
}

VisHeaderView.propTypes= {
    readout:  PropTypes.object.isRequired,
};


/**
 * @param {Object} p        React props object.
 * @param {Object} p.visRoot
 * @param {Object} p.currMouseState
 * @returns {HTML}
 * @constructor
 */
export function VisPreview({visRoot,currMouseState}) {
    var pv= getActivePlotView(visRoot);
    var mousePv= getPlotViewById(visRoot,currMouseState.plotId);
    return (
        <div>
            <ThumbnailView plotView={pv}/>
            <MagnifiedView plotView={mousePv} size={70} mouseState={currMouseState} />
        </div>
    );
}

VisPreview.propTypes= {
    visRoot : PropTypes.object.isRequired,
    currMouseState :PropTypes.object,
};
