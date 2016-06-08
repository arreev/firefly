/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import {isString} from 'lodash';
import {dispatchOnAppReady} from '../core/AppDataCntlr.js';

// Used for dispatch and action type constants
import * as TableStatsCntlr from '../visualize/TableStatsCntlr.js';
import * as HistogramCntlr from '../visualize/HistogramCntlr.js';
import * as XYPlotCntlr from '../visualize/XYPlotCntlr.js';
import * as TablesCntlr from '../tables/TablesCntlr.js';
import * as ReadoutCntlr from '../visualize/MouseReadoutCntlr.js';
import * as ImPlotCntlr from '../visualize/ImagePlotCntlr.js';
import * as MultiViewCntlr from '../visualize/MultiViewCntlr.js';
import * as AppDataCntlr from '../core/AppDataCntlr.js';
import * as DrawLayerCntlr from '../visualize/DrawLayerCntlr.js';
import {ApiExpandedView} from './ApiExpandedView.jsx';

// Parts of the lowlevel api
import * as ApiUtil from './ApiUtil.js';
import  * as ApiUtilChart from './ApiUtilChart.jsx';
import  * as ApiUtilImage from './ApiUtilImage.jsx';
import  * as ApiUtilTable from './ApiUtilTable.jsx';

// UI component
import {MultiImageViewer} from '../visualize/ui/MultiImageViewer.jsx';
import {ImageViewer} from '../visualize/iv/ImageViewer.jsx';
import {ImageMetaDataToolbar} from '../visualize/ui/ImageMetaDataToolbar.jsx';
import {MultiViewStandardToolbar} from '../visualize/ui/MultiViewStandardToolbar.jsx';
import {ExpandedModeDisplay} from '../visualize/iv/ExpandedModeDisplay.jsx';
import {ApiExpandedDisplay} from '../visualize/ui/ApiExpandedDisplay.jsx';
import {TablesContainer} from '../tables/ui/TablesContainer.jsx';
import {TablePanel} from '../tables/ui/TablePanel.jsx';
import {ChartsContainer} from '../visualize/ChartsContainer.jsx';
import {ChartsTableViewPanel} from '../visualize/ChartsTableViewPanel.jsx';

// builds the highlevel api
import {buildHighLevelApi} from './ApiHighlevelBuild.js';
import {buildViewerApi} from './ApiViewer.js';



// CSS
import './ApiStyle.css';


/**
 * Start in api mode. Will create the api and call window.onFireflyLoaded(firefly)
 */
export function initApi() {
    const lowlevelApi= buildLowlevelAPI();
    const viewInterface= buildViewerApi();
    const highLevelApi= buildHighLevelApi(lowlevelApi);
    if (!window.firefly) window.firefly= {};
    Object.assign(window.firefly, lowlevelApi, highLevelApi, viewInterface);
    const firefly= window.firefly;
    dispatchOnAppReady(() => {
        window.onFireflyLoaded && window.onFireflyLoaded(firefly);
    });
    initExpandedView();
    window.firefly.ignoreHistory = true;
}



/**
Structure of API
   {
              //--- High level API , all high level api are in the root
     all high level functions....
     
            //--- Low level API, lowlevel api are under action, ui, util
     action : { all dispatch functions...
               type: {all action type constants}
              }
     ui: { high level react components }
     util { renderDom, unrenderDom, isDebug, debug       // built by ApiUtil.js
            image : {image utility routines}                // imported from ApiUtilImage.js
            xyplot : {xyplot utility routines}               // imported from ApiUtilXYPlot.js //todo
            table : {table utility routines}                // imported from ApiUtilTable.js //todo
            data : {data utility routines???? }            // //todo do we need this?????
     }
   }
*/



/**
 * Return the api object.
 * 
 * @return {{action:{},ui:{},util:{}}}
 *
 */
export function buildLowlevelAPI() {


    const type= Object.assign({},
        findActionType(TableStatsCntlr,TableStatsCntlr.TBLSTATS_DATA_KEY),
        findActionType(HistogramCntlr, HistogramCntlr.HISTOGRAM_DATA_KEY),
        findActionType(XYPlotCntlr, XYPlotCntlr.XYPLOT_DATA_KEY ),
        findActionType(TablesCntlr, TablesCntlr.DATA_PREFIX),
        findActionType(TablesCntlr, TablesCntlr.RESULTS_PREFIX),
        findActionType(TablesCntlr, TablesCntlr.UI_PREFIX),
        findActionType(ReadoutCntlr, ReadoutCntlr.READOUT_PREFIX),
        findActionType(MultiViewCntlr, MultiViewCntlr.IMAGE_MULTI_VIEW_PREFIX),
        findActionType(ImPlotCntlr.default, ImPlotCntlr.PLOTS_PREFIX),
        findActionType(AppDataCntlr, AppDataCntlr.APP_DATA_PATH),
        findActionType(DrawLayerCntlr.default, DrawLayerCntlr.DRAWLAYER_PREFIX)
    );


    const action= Object.assign({},
        {type},
        findDispatch(TableStatsCntlr),
        findDispatch(HistogramCntlr),
        findDispatch(XYPlotCntlr),
        findDispatch(TablesCntlr),
        findDispatch(ReadoutCntlr),
        findDispatch(MultiViewCntlr),
        findDispatch(ImPlotCntlr),
        findDispatch(AppDataCntlr),
        findDispatch(DrawLayerCntlr)
    );

    const ui= {
        ImageViewer,
        MultiImageViewer,
        MultiViewStandardToolbar,
        ApiExpandedDisplay,
        ExpandedModeDisplay,
        ImageMetaDataToolbar,
        TablesContainer,
        TablePanel,
        ChartsContainer,
        ChartsTableViewPanel
    };
    
    const util= Object.assign({}, ApiUtil, {image:ApiUtilImage}, {chart:ApiUtilChart}, {table:ApiUtilTable}, {data:{}} );

    return { action, ui, util };
}


/**
 * pull all the dispatch functions out of the object
 * @param obj
 * @return {*}
 */
function findDispatch(obj) {
   return Object.keys(obj).reduce( (res,key) => {
        if (key.startsWith('dispatch')) res[key]= obj[key];
        return res;
    },{} );
}


/**
 * pull all the action type constants out of the object
 * @param obj
 * @param prefix
 * @return {*}
 */
function findActionType(obj,prefix) {
    return Object.keys(obj).reduce( (res,key) => {
        if (isString(obj[key]) && obj[key].startsWith(prefix) && obj[key].length>prefix.length) {
            res[key]= obj[key];
        }
        return res;
    },{} );
}


function initExpandedView(div){

    const EXPANDED_DIV= 'expandedArea';
    var expandedDivEl;
    if (div) {
        expandedDivEl= isString(div) ? document.getElementById(div) : div;
    } else {
        expandedDivEl= document.createElement('div');
        document.body.appendChild(expandedDivEl);
        expandedDivEl.id= EXPANDED_DIV;
    }
    
    ApiUtil.renderDOM(expandedDivEl, ApiExpandedView);
}