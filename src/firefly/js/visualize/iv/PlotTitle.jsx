/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */


import React from 'react';
import PropTypes from 'prop-types';
import Enum from 'enum';
import numeral from 'numeral';
import {getZoomDesc} from '../ZoomUtil.js';
import {primePlot} from '../PlotViewUtil.js';

import './PlotTitle.css';
import LOADING from 'html/images/gxt/loading.gif';

export const TitleType= new Enum(['INLINE', 'HEAD', 'EXPANDED']);

export function PlotTitle({plotView:pv, titleType, brief, working}) {
    let styleName= '';
    const plot= primePlot(pv);
    switch (titleType) {
        case TitleType.INLINE:
            styleName= 'plot-title-inline-title-container';
        break;
        case TitleType.HEAD:
            styleName= 'plot-title-header-title-container';
            break;
        case TitleType.EXPANDED:
            styleName= 'plot-title-expanded-title-container';
            break;

    }
    let zlStr= getZoomDesc(pv);
    let rotString= null;
    if (pv.rotation) {
        if (pv.plotViewCtx.rotateNorthLock) {
            rotString= 'North';
        } else {
            const angleStr= numeral(360-pv.rotation).format('#');
            rotString= angleStr + String.fromCharCode(176);
        }
        zlStr+=',';
    }

    return (
        <div className={styleName} title={plot.title}>
            <div className='plot-title-title'>{plot.title}</div>
            {!brief ? <div className='plot-title-zoom'><div dangerouslySetInnerHTML={{__html:zlStr}}/> </div> : ''}
            {!brief && rotString ? <div className='plot-title-rotation'>{rotString}</div> : ''}
            {working ?<img style={{width:14,height:14,padding:'0 3px 0 5px'}} src={LOADING}/> : ''}
        </div>
    );
}

PlotTitle.propTypes= {
    plotView : PropTypes.object,
    titleType: PropTypes.object.isRequired,
    annotationOps : PropTypes.object,
    brief : PropTypes.bool.isRequired
};

