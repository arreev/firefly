/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import {get, set, unset, has, isEmpty, uniqueId, cloneDeep, omit, omitBy, isNil, isPlainObject, isArray} from 'lodash';
import * as TblCntlr from './TablesCntlr.js';
import {SortInfo, SORT_ASC, UNSORTED} from './SortInfo.js';
import {FilterInfo} from './FilterInfo.js';
import {SelectInfo} from './SelectInfo.js';
import {flux} from '../Firefly.js';
import {encodeServerUrl} from '../util/WebUtil.js';
import {fetchTable} from '../rpc/SearchServicesJson.js';
import {DEF_BASE_URL} from '../core/JsonUtils.js';
import {ServerParams} from '../data/ServerParams.js';

export const MAX_ROW = Math.pow(2,31) - 1;
/* TABLE_REQUEST should match QueryUtil on the server-side */

const LSSTQueryPID = 'LSSTCataLogSearch';
/**
 *  @public
 */
/*----------------------------< creator functions ----------------------------*/


/**
 * Creates a table request object for the given id.
 * @param {string} id       required.  SearchProcessor ID.
 * @param {string} [title]  title to display with this table.
 * @param {object} [params] the parameters to include with this request.
 * @param {TableRequest} [options] more options.  see TableRequest for details.
 * @returns {TableRequest}
 * @pubic
 * @func  makeTblRequest
 * @memberof firefly.util.table
 */
export function makeTblRequest(id, title, params={}, options={}) {
    var req = {startIdx: 0, pageSize: 100};
    title = title || id;
    const tbl_id = options.tbl_id || uniqueTblId();
    var META_INFO = Object.assign(options.META_INFO || {}, {title, tbl_id});
    options = omit(options, 'tbl_id');
    return omitBy(Object.assign(req, options, params, {META_INFO, tbl_id, id}), isNil);
}

/**
 * Creates a table request for tabular data from a file.  Source of file may be
 * from a url or an absolute path on the server.
 * @param {string} [title]      title to display with this table.
 * @param {string} source       required; location of the ipac table. url or file path.
 * @param {string} [alt_source] use this if source does not exists.
 * @param {TableRequest} [options]  more options.  see TableRequest for details.
 * @returns {TableRequest}
 * @pubic
 * @func makeFileRequest
 * @memberof firefly.util.table
 */
export function makeFileRequest(title, source, alt_source, options={}) {
    const id = 'IpacTableFromSource';
    var req = {startIdx: 0, pageSize: 100};
    title = title || source;
    const tbl_id = options.tbl_id || uniqueTblId();
    options = omit(options, 'tbl_id');
    var META_INFO = Object.assign(options.META_INFO || {}, {title, tbl_id});
    return omitBy(Object.assign(req, options, {source, alt_source, META_INFO, tbl_id, id}), isNil);
}


/**
 * Parameters for cone search
 * @typedef {object} ConeParams
 * @global
 * @prop {string} SearchMethod  'Cone'.
 * @prop {string} position  name or coordinates of the search
 * @prop {string} radius    radius of the search in arcsec
 *
 */

/**
 * Parameters for eliptical search
 * @typedef {object} ElipParams
 * @global
 * @prop {string} SearchMethod  'Eliptical'.
 * @prop {string} position  name or coordinates of the search
 * @prop {string} radius    radius of the search in arcsec
 * @prop {string} radunits  the units for the radius or side, must be arcsec,arcmin,degree, default arcsec
 * @prop {string} ratio     ratio for elliptical request
 * @prop {string} posang    pa for elliptical request
 *
 */

/**
 * Parameters for box search
 * @typedef {object} BoxParams
 * @global
 * @prop {string} SearchMethod 'Eliptical'.
 * @prop {string} position  name or coordinates of the search
 * @prop {string} size      the length of a side for a box search
 *
 */

/**
 * creates the request to query IRSA catalogs.
 * @param {string} title    title to be displayed with this table result
 * @param {string} project
 * @param {string} catalog  the catalog name to search
 * @param {ConeParams|BoxParams|ElipParams} params   one of 'Cone','Eliptical','Box','Polygon','Table','AllSky'.
 * @param {TableRequest} [options]
 * @returns {TableRequest}
 * @access public
 * @func makeIrsaCatalogRequest
 *  @memberof firefly.util.table
 */
export function makeIrsaCatalogRequest(title, project, catalog, params={}, options={}) {
    var req = {startIdx: 0, pageSize: 100};
    title = title || catalog;
    options.use = options.use || 'catalog_overlay';
    const tbl_id = options.tbl_id || uniqueTblId();
    const id = 'GatorQuery';
    const UserTargetWorldPt = params.UserTargetWorldPt || params.position;  // may need to convert to worldpt.
    const catalogProject = project;
    var META_INFO = Object.assign(options.META_INFO || {}, {title, tbl_id});

    options = omit(options, 'tbl_id');
    params = omit(params, 'position');

    return omitBy(Object.assign(req, options, params, {id, tbl_id, META_INFO, UserTargetWorldPt, catalogProject, catalog}), isNil);
}

/**
 * creates the request to query LSST catalogs.  // TODO: more detail to be updated based on the LSST catalog DD content
 * @param {string} title    title to be displayed with this table result
 * @param {string} project
 * @param {string} catalog  the catalog name to search
 * @param {ConeParams|BoxParams|ElipParams} params   one of 'Cone','Eliptical','Box','Polygon','Table','AllSky'.
 * @param {TableRequest} [options]
 * @returns {TableRequest}
 * @access public
 * @func makeLsstCatalogRequest
 *  @memberof firefly.util.table
 */
export function makeLsstCatalogRequest(title, project, catalog, params={}, options={}) {
    var req = {startIdx: 0, pageSize: 100};

    title = title || catalog;
    options.use = options.use || 'lsst_catalog_overlay';
    const tbl_id = options.tbl_id || uniqueTblId();
    const id = LSSTQueryPID;
    const UserTargetWorldPt = params.UserTargetWorldPt || params.position;  // may need to convert to worldpt.
    const table_name = 'RunDeepForcedSource_limit100';    //catalog+'_queryresult';
    const meta_table = catalog;
    var META_INFO = Object.assign(options.META_INFO || {}, {title, tbl_id});


    options = omit(options, 'tbl_id');
    params = omit(params, 'position');

    return omitBy(Object.assign(req, options, params,
                                {id, tbl_id, META_INFO, UserTargetWorldPt, table_name, meta_table, project}), isNil);
}

/**
 * creates the request to query VO catalog
 * @param {string} title    title to be displayed with this table result
 * @param {ConeParams|BoxParams|ElipParams} params   one of 'Cone','Eliptical','Box','Polygon','Table','AllSky'.
 * @param {TableRequest} [options]
 * @returns {TableRequest}
 * @public
 * @func makeVOCatalogRequest
 *  @memberof firefly.util.table
 */
export function makeVOCatalogRequest(title, params={}, options={}) {
    var req = {startIdx: 0, pageSize: 100};
    options.use = options.use || 'catalog_overlay';
    const tbl_id = options.tbl_id || uniqueTblId();
    const id = 'ConeSearchByURL';
    const UserTargetWorldPt = params.UserTargetWorldPt || params.position;  // may need to convert to worldpt.
    var META_INFO = Object.assign(options.META_INFO || {}, {title, tbl_id});

    options = omit(options, 'tbl_id');
    params = omit(params, 'position');

    return omitBy(Object.assign(req, options, params, {id, tbl_id, META_INFO, UserTargetWorldPt}), isNil);
}

/**
 * create a deep clone of the given request.  tbl_id is removed from the cloned request.
 * @param {TableRequest} request  the original request to clone
 * @param {Object} params   additional parameters to add to the cloned request
 * @returns {TableRequest}
 * @public
 * @func cloneRequest
 * @memberof firefly.util.table
 */
export function cloneRequest(request, params = {}) {
    const req = cloneDeep(omit(request, 'tbl_id'));
    unset(req, 'META_INFO.tbl_id');
    return Object.assign(req, params);
}

/*---------------------------- creator functions >----------------------------*/


/**
 * tableRequest will be sent to the server as a json string.
 * @param {TableRequest} tableRequest is a table request params object
 * @param {number} [hlRowIdx] set the highlightedRow.  default to startIdx.
 * @returns {Promise.<TableModel>}
 * @public
 * @func doFetchTable
 * @memberof firefly.util.table
 * @deprecated  use firefly.util.table.fetchTable instead
 */
export function doFetchTable(tableRequest, hlRowIdx) {
    return fetchTable(tableRequest, hlRowIdx);
}

export function doValidate(type, action) {
    if (type !== action.type) {
        error(action, `Incorrect type:${action.type} was sent to a ${type} actionCreator.`);
    }
    if (!action.payload) {
        error(action, 'Invalid action.  Payload is missing.');
    }
    var {request} = action.payload;
    if (type === TblCntlr.TABLE_FETCH ) {
        if (request.id) {
            error(action, 'Required "id" field is missing.');
        }
        if (request.tbl_id) {
            error(action, 'Required "tbl_id" field is missing.');
        }
    } else if(type === TblCntlr.TABLE_HIGHLIGHT) {
        const idx = action.payload.highlightedRow;
        if (!idx || idx<0) {
            error(action, 'highlightedRow must be a positive number.');
        }
    }
    return action;
}

/**
 * updates the given action with a new error given by cause.
 * action.err is stored as an array of errors.  Errors may be a String or an Error type.
 * @param {Object} action  the action to update
 * @param {string} cause  the error to be added.
 * @public
 * @func error
 * @memberof firefly.util.table
 */
export function error(action, cause) {
    (action.err = action.err || []).push(cause);
}

/**
 * returns true is there is data within the given range.  this is needed because
 * of paging table not loading the full table.
 * @param {number} startIdx
 * @param {number} endIdx
 * @param {TableModel} tableModel
 * @returns {boolean}
 * @public
 * @func isTblDataAvail
 * @memberof firefly.util.table
 */
export function isTblDataAvail(startIdx, endIdx, tableModel) {
    if (!tableModel) return false;
    endIdx =  endIdx >0 ? Math.min( endIdx, tableModel.totalRows) : startIdx;
    if (startIdx >=0 && endIdx > startIdx) {
        const data = get(tableModel, 'tableData.data', []);
        const dataCount = Object.keys(data.slice(startIdx, endIdx)).length;
        return dataCount === (endIdx-startIdx);
    } else return false;
}

/**
 * returns the table model with the given tbl_id
 * @param tbl_id
 * @returns {TableModel}
 * @func getTblById
 * @memberof firefly.util.table
 */
export function getTblById(tbl_id) {
    return get(flux.getState(),[TblCntlr.TABLE_SPACE_PATH, 'data', tbl_id]);
}

/**
 * returns the table group information
 * @param {string} tbl_group    the group name to look for
 * @returns {TableGroup}
 * @public
 * @memberof firefly.util.table
 * @func getTableGroup
 */
export function getTableGroup(tbl_group='main') {
    return get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH, 'results', tbl_group]);
}

/**
 * returns the table group name given a tbl_id.  it will return undefined if
 * the given tbl_id is not in a group.
 * @param {string} tbl_id    table id
 * @returns {TableGroup}
 * @public
 * @memberof firefly.util.table
 * @func findGroupByTblId
 */
export function findGroupByTblId(tbl_id) {
    const resultsRoot = get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH, 'results'], {});
    const groupName = Object.keys(resultsRoot).find( (tbl_grp_id) => {
        return has(resultsRoot, [tbl_grp_id, 'tables', tbl_id]);
    });
    return groupName;
}

/**
 * returns an array of tbl_id for the given tbl_group_id
 * @param {string} tbl_group_id    table group name.  defaults to 'main' if not given
 * @returns {String[]} array of tbl_id
 * @public
 * @func getTblIdsByGroup
 * @memberof firefly.util.table
 */
export function getTblIdsByGroup(tbl_group_id = 'main') {
    const tableGroup = get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH, 'results', tbl_group_id]);
    return Object.keys(get(tableGroup, 'tables', {}));
}

/**
 * returns the table information for the given id and group.
 * @param {string} tbl_id       table id.
 * @param {string} tbl_group    table group name.  defaults to 'main' if not given
 * @returns {TableModel}
 * @public
 * @func getTableInGroup
 * @memberof firefly/util/table
 */
export function getTableInGroup(tbl_id, tbl_group='main') {
    return get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH, 'results', tbl_group, 'tables',  tbl_id]);
}

/**
 * get the table working state by tbl_ui_id
 * @param {string} tbl_ui_id     table UI id.
 * @returns {Object}
 * @public
 *  @memberof firefly.util.table
 *  @func  getTableUiById
 */
export function getTableUiById(tbl_ui_id) {
    return get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH, 'ui', tbl_ui_id]);
}

/**
 * returns the first table working state for the given tbl_id
 * @param {string} tbl_id
 * @returns {Object}
 * @memberof firefly.util.table
 * @func getTableUiByTblId
 */
export function getTableUiByTblId(tbl_id) {
    const uiRoot = get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH, 'ui'], {});
    const tbl_ui_id = Object.keys(uiRoot).find( (ui_id) => {
        return get(uiRoot, [ui_id, 'tbl_id']) === tbl_id;
    });
    return tbl_ui_id || uiRoot[tbl_ui_id];
}

/**
 * returns the working state of the currently expanded table.
 * @returns {Object}
 * @public
 * @memberof firefly.util.table
 * @func getTblExpandedInfo
 */
export function getTblExpandedInfo() {
    return get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH, 'ui', 'expanded'], {});
}

/**
 * returns true if the table referenced by the given tbl_id is fully loaded.
 * @param {string} tbl_id
 * @returns {boolean}
 * @public
 * @memberof firefly.util.table
 * @func  isFullyLoaded
 */
export function isFullyLoaded(tbl_id) {
    return isTableLoaded(getTblById(tbl_id));
}

/**
 * Returns the column index with the given name; otherwise, -1.
 * @param {TableModel} tableModel
 * @param {string} colName
 * @returns {number}
 * @memberof firefly.util.table
 * @func getColumnIdx
 */
export function getColumnIdx(tableModel, colName) {
    const cols = get(tableModel, 'tableData.columns', []);
    return cols.findIndex((col) => {
        return col.name === colName;
    });
}

/**
 * returns column information for the given name.
 * @param {TableModel} tableModel
 * @param {string} colName
 * @returns {TableColumn}
 * @memberof firefly.util.table
 * @func getColumn
 */
export function getColumn(tableModel, colName) {
    const colIdx = getColumnIdx(tableModel, colName);
    if (colIdx >= 0) {
        return get(tableModel, `tableData.columns.${colIdx}`);
    }
}

export function getFilterCount(tableModel) {
    const filterInfo = get(tableModel, 'request.filters');
    const filterCount = filterInfo ? filterInfo.split(';').length : 0;
    return filterCount;
}

export function clearFilters(tableModel) {
    const request = get(tableModel, 'request');
    if (request && request.filters) {
        const newRequest = Object.assign({}, request, {filters: ''});
        TblCntlr.dispatchTableFilter(newRequest, 0);
    }
}

export function clearSelection(tableModel) {
    if (tableModel) {
        const selectInfoCls = SelectInfo.newInstance({rowCount: tableModel.totalRows});
        TblCntlr.dispatchTableSelect(tableModel.tbl_id, selectInfoCls.data);
    }
}

/**
 * return the tbl_id of the active table for the given group.
 * @param {string} tbl_group group name; defaults to 'main' if not given.
 * @returns {string}
 * @public
 * @memberof firefly.util.table
 * @func getActiveTableId
 */
export function getActiveTableId(tbl_group='main') {
    return get(flux.getState(), [TblCntlr.TABLE_SPACE_PATH,'results',tbl_group,'active']);
}

/**
 *
 * @param {TableModel} tableModel
 * @param {number} rowIdx
 * @param {string} colName
 * @return {string}
 * @memberof firefly.util.table
 * @func getCellValue
 */
export function getCellValue(tableModel, rowIdx, colName) {
    if (get(tableModel, 'tableData.data.length', 0) > 0) {
        const colIdx = getColumnIdx(tableModel, colName);
        // might be undefined if row is not loaded
        return get(tableModel.tableData.data, [rowIdx, colIdx]);
    }
}


/**
 * return true if the given table is fully loaded.
 * @param {TableModel} tableModel
 * @returns {boolean}
 * @public
 *  @memberof ffirefly.util.table
 * @func isTableLoaded
 */
export function isTableLoaded(tableModel) {
    const status = tableModel && !tableModel.isFetching && get(tableModel, 'tableMeta.Loading-Status', 'COMPLETED');
    return status === 'COMPLETED';
}

/**
 * This function transform the json data from the server to fit the need of the UI.
 * For instance, the column's name is repeated after transform.  This is good for the UI.
 * But, it's more efficient to not include it during data transfer from the server.
 * @param {TableModel} tableModel
 * @returns {*}
 * @public
 * @memberof firefly.util.table
 * @func transform
 */
export function transform(tableModel) {

    if (tableModel.tableData && tableModel.tableData.data) {
        const cols = tableModel.tableData.columns;
        // change row data from [ [val] ] to [ {cname:val} ]
        tableModel.tableData.data = tableModel.tableData.data.map( (row) => {
            return cols.reduce( (nrow, col, cidx) => {
                nrow[col.name] = row[cidx];
                return nrow;
            }, {});
        });
    }
}

/**
 * This function merges the source object into the target object
 * by traversing and comparing every like path.  If a value was
 * merged at any data node in the data graph, the node and all of its
 * parent nodes will be shallow cloned and returned.  Otherwise, the target's value
 * will be returned.
 * @param {Object} target
 * @param {Object} source
 * @returns {Object}
 * @public
 * @memberof firefly.util.table
 * @func smartMerge
 */
export function smartMerge(target, source) {
    if (!target) return source;

    if (isPlainObject(source)) {
        const objChanges = {};
        Object.keys(source).forEach((k) => {
            const nval = smartMerge(target[k], source[k]);
            if (nval !== target[k]) {
                objChanges[k] = nval;
            }
        });
        return (isEmpty(objChanges)) ? target : Object.assign({}, target, objChanges);
    } else if (isArray(source)){
        const aryChanges = [];
        source.forEach((v, idx) => {
            const nval = smartMerge(target[idx], source[idx]);
            if (nval !== target[idx]) {
                aryChanges[idx] = nval;
            }
        });
        if (isEmpty(aryChanges)) return target;
        else {
            const nAry = target.slice();
            aryChanges.forEach((v, idx) => nAry[idx] = v);
            return nAry;
        }
    } else {
        return (target === source) ? target : source;
    }
}

/**
 * sort the given tableModel based on the given request
 * @param {TableModel} origTableModel original table model.  this is returned when direction is UNSORTED.
 * @param {string} sortInfoStr
 * @returns {TableModel}
 * @public
 * @func sortTable
 * @memberof firefly.util.table
 */
export function sortTable(origTableModel, sortInfoStr) {
    const tableModel = cloneDeep(origTableModel);
    set(tableModel, 'request.sortInfo', sortInfoStr);
    const {data, columns} = tableModel.tableData;
    sortTableData(data, columns, sortInfoStr);
    return tableModel;
}

/**
 * sort table data in place.
 * @param {TableData} tableData
 * @param {TableColumn[]} columns
 * @param {string} sortInfoStr
 * @returns {TableData}
 * @public
 * @memberof firefly.util.table
 * @func sortTableData
 */
export function sortTableData(tableData, columns, sortInfoStr) {
    const sortInfoCls = SortInfo.parse(sortInfoStr);
    const colName = get(sortInfoCls, 'sortColumns.0');
    const dir = get(sortInfoCls, 'direction', UNSORTED);
    if (dir === UNSORTED || get(tableData, 'length', 0) === 0) return tableData;

    const multiplier = dir === SORT_ASC ? 1 : -1;
    const colIdx = columns.findIndex( (col) => {return col.name === colName;} );
    const col = columns[colIdx];
    if (!col) return tableData;

    var comparator;
    if (!col.type || ['char', 'c'].includes(col.type) ) {
        comparator = (r1, r2) => {
            const [s1, s2] = [r1[colIdx], r2[colIdx]];
            return multiplier * (s1 > s2 ? 1 : -1);
        };
    } else {
        comparator = (r1, r2) => {
            const [v1, v2] = [r1[colIdx], r2[colIdx]];
            return multiplier * (Number(v1) - Number(v2));
        };
    }
    tableData.sort(comparator);
    return tableData;
}

/**
 * return a new table after applying the given filters.
 * @param {TableModel} tableModel
 * @param {string} filterInfoStr filters are separated by comma(',').
 * @returns {TableModel}
 * @memberof firefly.util.table
 * @func filterTable
 */
export function filterTable(tableModel, filterInfoStr) {
    const filtered = cloneDeep(tableModel);
    set(filtered, 'request.filters', filterInfoStr);
    const comparators = filterInfoStr.split(';').map((s) => s.trim()).map((s) => FilterInfo.createComparator(s, tableModel));
    filtered.tableData.data = tableModel.tableData.data.filter((row) => {
        return comparators.reduce( (rval, match) => rval && match(row), true);
    } );
    filtered.totalRows = tableModel.tableData.data.length;
    return filtered;
}

/**
 * collects all available table information given the tbl_id
 * @param {string} tbl_id
 * @param {number} aPageSize  use this pageSize instead of the one in the request.
 * @returns {{tableModel, tbl_id, title, totalRows, request, startIdx, endIdx, hlRowIdx, currentPage, pageSize, totalPages, highlightedRow, selectInfo, error}}
 * @public
 * @memberof firefly.util.table
 * @func getTblInfoById
 */
export function getTblInfoById(tbl_id, aPageSize) {
    const tableModel = getTblById(tbl_id);
    return getTblInfo(tableModel, aPageSize);
}

/**
 * collects all available table information given the tableModel.
 * @param {TableModel} tableModel
 * @param {number} aPageSize  use this pageSize instead of the one in the request.
 * @returns {{tableModel, tbl_id, title, totalRows, request, startIdx, endIdx, hlRowIdx, currentPage, pageSize, totalPages, highlightedRow, selectInfo, error}}
 * @public
 * @memberof firefly.util.table
 * @func getTblInfo
 */
export function getTblInfo(tableModel, aPageSize) {
    if (!tableModel) return {};
    var {tbl_id, request, highlightedRow=0, totalRows=0, tableMeta={}, selectInfo, error} = tableModel;
    const {title} = tableMeta;
    const pageSize = aPageSize || get(request, 'pageSize', MAX_ROW);  // there should be a pageSize.. default to 1 in case of error.  pageSize cannot be 0 because it'll overflow.
    if (highlightedRow < 0 ) {
        highlightedRow = 0;
    } else  if (highlightedRow >= totalRows-1) {
        highlightedRow = totalRows-1;
    }
    const currentPage = highlightedRow >= 0 ? Math.floor(highlightedRow / pageSize)+1 : 1;
    const hlRowIdx = highlightedRow >= 0 ? highlightedRow % pageSize : 0;
    const startIdx = (currentPage-1) * pageSize;
    const endIdx = Math.min(startIdx+pageSize, totalRows) || get(tableModel,'tableData.data.length', startIdx) ;
    var totalPages = Math.ceil((totalRows || 0)/pageSize);
    return { tableModel, tbl_id, title, totalRows, request, startIdx, endIdx, hlRowIdx, currentPage, pageSize,totalPages, highlightedRow, selectInfo, error};
}

/**
 * returns the url to download a snapshot of the current table data.
 * @param {string} tbl_ui_id  UI id of the table
 * @returns {string}
 * @public
 * @memberof firefly.util.table
 * @func getTableSourceUrl
 */
export function getTableSourceUrl(tbl_ui_id) {
    const def = {
        startIdx: 0,
        pageSize : MAX_ROW
    };
    const {columns, request} = getTableUiById(tbl_ui_id) || {};
    const tableRequest = Object.assign(def, cloneDeep(request));
    const visiCols = columns.filter( (col) => {
                return isNil(col) || col.visibility === 'show';
            }).map( (col) => {
                return col.name;
            } );
    if (visiCols.length !== columns.length) {
        tableRequest['inclCols'] = visiCols.toString();
    }
    Reflect.deleteProperty(tableRequest, 'tbl_id');
    const params = {
        [ServerParams.COMMAND]: ServerParams.TABLE_SAVE,
        [ServerParams.REQUEST]: JSON.stringify(tableRequest),
        file_name: tableRequest.file_name
    };
    return encodeServerUrl(DEF_BASE_URL, params);
}

/**
 * returns an object map of the column name and its width.
 * The width is the number of characters needed to display
 * the header and the data in a table given columns and dataAry.
 * @param {TableColumn[]} columns  array of column object
 * @param {TableData} dataAry  array of array.
 * @returns {Object.<string,number>} a map of cname -> width
 * @public
 * @memberof firefly.util.table
 * @func calcColumnWidths
 */
export function calcColumnWidths(columns, dataAry) {
    return columns.reduce( (pv, cv, idx) => {
        const cname = cv.label || cv.name;
        var width = Math.max(cname.length, get(cv, 'units.length', 0),  get(cv, 'type.length', 0));
        width = dataAry.reduce( (maxWidth, row) => {
            return Math.max(maxWidth, get(row, [idx, 'length'], 0));
        }, width);  // max width of data
        pv[idx] = width;
        return pv;
    }, {ROWID: 8});
}

/**
 * create a unique table id (tbl_id)
 * @returns {string}
 * @public
 * @memberof firefly.util.table
 * @func uniqueTblId
 */
export function uniqueTblId() {
    const id = uniqueId('tbl_id-');
    if (getTblById(id)) {
        return uniqueTblId();
    } else {
        return id;
    }
}

/**
 * create a unique table UI id (tbl_ui_id)
 * @returns {string}
 * @public
 * @memberof firefly.util.table
 * @func uniqueTblUiId
 */
export function uniqueTblUiId() {
    return uniqueId('tbl_ui_id-');
}
/**
 *  This function provides a patch until we can reliably determine that the ra/dec columns use radians or degrees.
 * @param tableOrMeta the table object or the tableMeta object
 * @public
 * @memberof firefly.util.table
 * @func isTableUsingRadians
 *
 */
export function isTableUsingRadians(tableOrMeta) {
    if (!tableOrMeta) return false;
    const tableMeta= tableOrMeta.tableMeta || tableOrMeta;
    return has(tableMeta, 'HIERARCH.AFW_TABLE_VERSION');
}

export function createErrorTbl(tbl_id, error) {
    return {tbl_id, error};
}