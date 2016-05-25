/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
import {flux} from '../Firefly.js';

import {updateSet, updateMerge} from '../util/WebUtil.js';
import {get, has, omit, omitBy, isUndefined, isString} from 'lodash';

import {doFetchTable, getTblById} from '../tables/TableUtil.js';
import * as TablesCntlr from '../tables/TablesCntlr.js';
import {serializeDecimateInfo} from '../tables/Decimate.js';
import {logError} from '../util/WebUtil.js';

export const XYPLOT_DATA_KEY = 'xyplot';
export const LOAD_PLOT_DATA = `${XYPLOT_DATA_KEY}/LOAD_COL_DATA`;
export const UPDATE_PLOT_DATA = `${XYPLOT_DATA_KEY}/UPDATE_COL_DATA`;
export const SET_SELECTION = `${XYPLOT_DATA_KEY}/SET_SELECTION`;
const SET_ZOOM = `${XYPLOT_DATA_KEY}/SET_ZOOM`;
const RESET_ZOOM = `${XYPLOT_DATA_KEY}/RESET_ZOOM`;

/*
 Possible structure of store:
  /xyplot
    chartId: Object - the name of this node matches chart id
    {
         // tblXYPlotData
         tblId: string // table id
         isPlotDataReady: boolean
         decimatedUnzoomed: boolean // tells that unzoomed data are decimated
         xyPlotData: {
                    rows: [[x: string, y: string, rowIdx: string]*] ,
                    decimateKey: string,
                    xMin: string,
                    xMax: string,
                    yMin: string,
                    yMax: string,
                    weightMin: string,
                    weightMax: string,
                    idStr: string
         }
         xyPlotParams: {
           title: string
           xyRatio: number
           selection: {xMin, xMax, yMin, yMax} // currently selected rectangle
           zoom: {xMin, xMax, yMin, yMax} // currently zoomed rectangle
           stretch: string (fit|fill)
           x: {
                columnOrExpr
                label
                unit
                options: [grid,log,flip]
                nbins
                min
                max
              }
           y: {
                columnOrExpr
                label
                unit
                options: [grid,log,flip]
                nbins
                min
                max
           }
     }
 */

/*
 * Load xy plot data
 * @param {string} chartId - if no chart id is specified table id is used as chart id
 * @param {Object} xyPlotParams - XY plot options (column names, etc.)
 * @param {ServerRequest} searchRequest - table search request
 * @param {function} dispatcher only for special dispatching uses such as remote
 */
export function dispatchLoadPlotData(chartId, xyPlotParams, searchRequest, dispatcher= flux.process) {
    dispatcher({type: LOAD_PLOT_DATA, payload: {chartId: (chartId||searchRequest.tbl_id), xyPlotParams, searchRequest}});
}

/*
 * Set selection to give user choice of actions on selection (zoom, filter, or select points)
 * @param {String} chartId - chart id
 * @param {Object} selection - {xMin, xMax, yMin, yMax}
 */
export function dispatchSetSelection(chartId, selection) {
    flux.process({type: SET_SELECTION, payload: {chartId, selection}});
}

/*
 * Zoom XY plot to a given selection or reset zoom if no selection is given
 * @param {String} chartId - chart id
 * @param {String} tblId - table id
 * @param {Object} selection - {xMin, xMax, yMin, yMax}
 */
export function dispatchZoom(chartId, tblId, selection) {
    const {xyPlotData, xyPlotParams, decimatedUnzoomed} = get(flux.getState(), [XYPLOT_DATA_KEY,chartId], {});
    if (xyPlotData && xyPlotParams) {
        if (selection) {
            // zoom to selection
            if (xyPlotData.decimateKey) {
                const tableModel = getTblById(tblId);
                if (tableModel) {
                    const paramsWithZoom = Object.assign({}, xyPlotParams, {zoom: xyPlotParams.selection});
                    dispatchLoadPlotData(chartId, paramsWithZoom, tableModel.request);
                }
            } else {
                dispatchSetZoom(chartId, selection);
            }
        } else {
            // reset zoom
            if (decimatedUnzoomed || isUndefined(decimatedUnzoomed)) {
                const tableModel = getTblById(tblId);
                if (tableModel) {
                    const paramsWithoutZoom = Object.assign({}, omit(xyPlotParams, 'zoom'));
                    dispatchLoadPlotData(chartId, paramsWithoutZoom, tableModel.request);
                }
            } else {
                dispatchResetZoom(chartId);
            }

        }
    }
}

function dispatchSetZoom(chartId, selection) {
    flux.process({type: SET_ZOOM, payload: {chartId, selection}});
}

function dispatchResetZoom(chartId) {
    flux.process({type: RESET_ZOOM, payload: {chartId}});
}


/*
 * @param rawAction (its payload should contain searchRequest to get source table and histogram parameters)
 * @returns function which loads plot data (x, y, rowIdx, etc.)
 */
export function loadPlotData (rawAction) {
    return (dispatch) => {
        const {searchRequest, xyPlotParams, chartId} = rawAction.payload;
        dispatch({ type : LOAD_PLOT_DATA, payload : rawAction.payload });
        if (searchRequest && xyPlotParams) {
            fetchPlotData(dispatch, searchRequest, xyPlotParams, chartId);
        }

    };
}

/**
 * The data is an object with
 * chartId - string, chart id,
 * isPlotDataReady - boolean, flags that xy plot data are available
 * xyPlotData - an array of data rows
 * xyPlotParams - plot parameters
 * decimatedUnzoomed - tells if unzoomed data are decimated
 * @param data {Object} the data to merge with the xyplot branch
 * @returns {{type: string, payload: object}}
 */
function updatePlotData(data) {
    return { type : UPDATE_PLOT_DATA, payload: data };
}

export function reducer(state={}, action={}) {
    switch (action.type) {
        case (TablesCntlr.TABLE_REMOVE)  :
        {
            const tbl_id = action.payload.tbl_id;
            const chartsToDelete = [];
            Object.keys(state).forEach((cid) => {
                if (state[cid].tblId === tbl_id) {
                    chartsToDelete.push(cid);
                }
            });
            return (chartsToDelete.length > 0) ?
                Object.assign({}, omit(state, chartsToDelete)) : state;
        }
        case (LOAD_PLOT_DATA)  :
        {
            const {chartId, xyPlotParams, searchRequest} = action.payload;
            return updateSet(state, chartId,
                { tblId: searchRequest.tbl_id, isPlotDataReady: false, xyPlotParams, decimatedUnzoomed: get(state, [chartId,'decimatedUnzoomed'])});
        }
        case (UPDATE_PLOT_DATA)  :
        {
            const {isPlotDataReady, decimatedUnzoomed, xyPlotData, chartId, xyPlotParams} = action.payload;
            if (state[chartId].xyPlotParams === xyPlotParams) {
                const decimatedUnzoomedNext = isUndefined(decimatedUnzoomed) ? state[chartId].decimatedUnzoomed : decimatedUnzoomed;
                const newParams = xyPlotParams.selection ?  updateSet(xyPlotParams, 'selection', undefined) : xyPlotParams;
                return updateMerge(state, chartId,
                    {isPlotDataReady, decimatedUnzoomed: decimatedUnzoomedNext, xyPlotData, xyPlotParams: newParams});
            }
            return state;
        }
        case (SET_SELECTION) :
        {
            const {chartId, selection} = action.payload;
            return updateSet(state, [chartId,'xyPlotParams','selection'], selection);
        }
        case (SET_ZOOM) :
        {
            const {chartId, selection} = action.payload;
            const newState = updateSet(state, [chartId,'xyPlotParams','zoom'], selection);
            Reflect.deleteProperty(newState[chartId].xyPlotParams, 'selection');
            return newState;
        }
        case (RESET_ZOOM) :
        {
            const chartId = action.payload.chartId;
            const newParams = Object.assign({}, omit(state[chartId].xyPlotParams, ['selection', 'zoom']));
            return updateSet(state, [chartId,'xyPlotParams'], newParams);
        }
        case (TablesCntlr.TABLE_SELECT) :
        {
            const tbl_id = action.payload.tbl_id; //also has selectInfo
            let newState = state;
            Object.keys(state).forEach((cid) => {
                if (state[cid].tblId === tbl_id || has(state[cid], ['xyPlotParams','selection'])) {
                    newState = updateSet(newState, [cid,'xyPlotParams','selection'], undefined);
                }
            });
            return newState;
        }
        default:
            return state;
    }
}



/**
 * fetches xy plot data
 * set isColStatsReady to true once done.
 * @param dispatch
 * @param activeTableServerRequest table search request to obtain source table
 * @param xyPlotParams object, which contains xy plot parameters
 * @param {string} chartId  - chart id
 */
function fetchPlotData(dispatch, activeTableServerRequest, xyPlotParams, chartId) {

    if (!xyPlotParams) { return; }

    let limits = [];
    if (xyPlotParams.zoom) {
        const {xMin, xMax, yMin, yMax}  = xyPlotParams.zoom;
        limits = [xMin, xMax, yMin, yMax];
    }

    const req = Object.assign({}, omit(activeTableServerRequest, ['tbl_id', 'META_INFO']), {
        'startIdx' : 0,
        'pageSize' : 1000000,
        //'inclCols' : `${xyPlotParams.x.columnOrExpr},${xyPlotParams.y.columnOrExpr}`, // ignored if 'decimate' is present
        'decimate' : serializeDecimateInfo(xyPlotParams.x.columnOrExpr, xyPlotParams.y.columnOrExpr, 10000, 1.0, ...limits)
    });

    req.tbl_id = `xy-${chartId}`;

    doFetchTable(req).then(
        (tableModel) => {
            if (tableModel.tableData && tableModel.tableData.data) {
                const {tableMeta} = tableModel;
                const xyPlotData = omitBy({
                    rows: tableModel.tableData.data,
                    decimateKey: tableMeta['decimate_key'],
                    xMin: tableMeta['decimate.X-MIN'],
                    xMax: tableMeta['decimate.X-MAX'],
                    yMin: tableMeta['decimate.Y-MIN'],
                    yMax: tableMeta['decimate.Y-MAX'],
                    weightMin: tableMeta['decimate.WEIGHT-MIN'],
                    weightMax: tableMeta['decimate.WEIGHT-MAX'],
                    idStr: tableMeta['tbl_id']
                }, isUndefined);

                // convert strings with numbers into numbers
                Object.keys(xyPlotData).forEach( (prop) => {
                    const val = xyPlotData[prop];
                    if (isString(val) && isFinite(val)) { xyPlotData[prop] = Number(val); }
                });

                dispatch(updatePlotData(
                    {
                        isPlotDataReady : true,
                        // when zoomed, we don't know if the unzoomed data are decimated or not
                        decimatedUnzoomed: Boolean(tableMeta['decimate_key']) || (xyPlotParams.zoom ? undefined : false),
                        xyPlotParams,
                        xyPlotData,
                        chartId
                    }));
            }
        }
    ).catch(
        (reason) => {
            logError(`Failed to fetch XY plot data: ${reason}`);
        }
    );

}


