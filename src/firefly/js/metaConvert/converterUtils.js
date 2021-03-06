/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import {isEmpty, get} from 'lodash';
import {ServerRequest} from '../data/ServerRequest.js';
import {WebPlotRequest} from '../visualize/WebPlotRequest.js';
import {ZoomType} from '../visualize/ZoomType.js';
import {getTblById,getTblInfo, getCellValue} from '../tables/TableUtil.js';
import {getCenterColumns} from '../tables/TableInfoUtil.js';
import {converterFactory} from './ConverterFactory.js';
import {MetaConst} from '../data/MetaConst.js';

const dataSourceUpper= 'DATASOURCE';

const getSetInSrByRow= (table,sr,rowNum) => (col) => {
    sr.setSafeParam(col.name, getCellValue(table,rowNum,col.name));
};


/**
 *
 * @param table table data
 * @param {Array.<string>} colToUse columns from table
 * @param {Array.<string>} headerParams meta data parameters
 * @param {RangeValues} rv rangeValues
 * @param {number} colorTableId color table id
 * @return {function} see below, function takes plotId, reqKey,title, rowNum, extranParams and returns a WebPlotRequest
 *
 */
export function makeServerRequestBuilder(table, colToUse, headerParams, rv=null, colorTableId=0) {
    /**
     * @param plotId - the plot id for the request
     * @param reqKey - search processor request key
     * @param title - title of plot
     * @param rowNum - get the row number of data in the table
     * @param extraParams can be an object with single key or an array of objects with single key
     * @return {WebPlotRequest}
     */
    return (plotId, reqKey, title, rowNum, extraParams) => {
        const sr= new ServerRequest(reqKey);
        if (typeof extraParams === 'object') {
            if (!Array.isArray(extraParams)) extraParams= [extraParams];
            extraParams.forEach( (p) => sr.setParam(p));
        }
        const {columns}= table.tableData;
        const {tableMeta:meta}= table;
        const setInSr= getSetInSrByRow(table,sr,rowNum);

        if (!Array.isArray(colToUse) && typeof colToUse === 'string') colToUse= [colToUse];
        if (!Array.isArray(headerParams) && typeof headerParams=== 'string') headerParams= [headerParams];


        if (isEmpty(colToUse) || colToUse[0].toUpperCase()==='ALL') {
            columns.forEach(setInSr);
        }
        else {
            columns.filter((c) => colToUse.includes(c.name)).forEach( setInSr);
        }

        if (!isEmpty(headerParams)) {
            if (headerParams[0].toUpperCase()==='ALL') {
                Object.keys(meta).forEach( (metaKey) => sr.setSafeParam(metaKey, meta[metaKey]) );
                
            }
            else {
                Object.keys(meta).filter( (m) => headerParams.includes(m) )
                    .forEach( (metaKey) => sr.setSafeParam(metaKey, meta[metaKey]) );
            }
        }
        const wpReq= WebPlotRequest.makeProcessorRequest(sr,title);
        // wpReq.setZoomType(ZoomType.FULL_SCREEN);
        wpReq.setInitialColorTable(colorTableId);
        wpReq.setTitle(title);
        wpReq.setPlotId(plotId);
        wpReq.setZoomType(ZoomType.TO_WIDTH);
        if (rv) wpReq.setInitialRangeValues(rv);
        return wpReq;
    };
}


export const computePlotId= (plotIdRoot ,plotIdx) => `${plotIdRoot}-row-${plotIdx}`;

/**
 * helper function to make a list of table row that should be used to plot in grid mode
 * @param table
 * @param maxRows
 * @param plotIdRoot
 * @return {Array}
 */
export function findGridTableRows(table,maxRows, plotIdRoot) {

    const {startIdx, endIdx, highlightedRow}= getTblInfo(table, maxRows);

    let j= 0;
    const retval= [];

    for(let i=startIdx; (i<endIdx );i++) {
        retval[j++] = {plotId: computePlotId(plotIdRoot, i), row: i, highlight: i === highlightedRow};
    }
    return retval;
}


/**
 * Guess if this table contains image meta data
 * @param tbl_id
 * @return {boolean} true if there is image meta data
 */
export function isMetaDataTable(tbl_id) {
    const table= getTblById(tbl_id);
    if (isEmpty(table)) return false;
    const {tableMeta, totalRows} = table;
    if (!tableMeta || !totalRows) return false;

    const hasDsCol= Boolean(Object.keys(tableMeta).find( (key) => key.toUpperCase()===dataSourceUpper));

    return Boolean(tableMeta[MetaConst.DATASET_CONVERTER] || hasDsCol);
}

/**
 * Guess if this table contains catalog data
 * @param tbl_id
 * @return {boolean} true if this is catalog data
 */
export function isCatalogTable(tbl_id) {
    const table= getTblById(tbl_id);
    return get(table, 'totalRows') && !isEmpty(getCenterColumns(table));
}

