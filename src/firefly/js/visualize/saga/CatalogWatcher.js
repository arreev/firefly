/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import {take} from 'redux-saga/effects';
import {isEmpty, isNil, get} from 'lodash';
import {TABLE_LOADED, TABLE_SELECT,TABLE_HIGHLIGHT,TABLE_REMOVE,TABLE_UPDATE} from '../../tables/TablesCntlr.js';
import {getDlRoot, SUBGROUP, dispatchAttachLayerToPlot, dispatchChangeVisibility, dispatchCreateDrawLayer,
        dispatchDestroyDrawLayer, dispatchModifyCustomField} from '../DrawLayerCntlr.js';
import ImagePlotCntlr, {visRoot} from '../ImagePlotCntlr.js';
import {getTblById, doFetchTable, getTableGroup, isTableUsingRadians} from '../../tables/TableUtil.js';
import {cloneRequest, makeTableFunctionRequest, MAX_ROW} from '../../tables/TableRequestUtil.js';
import {serializeDecimateInfo} from '../../tables/Decimate.js';
import {getDrawLayerById, getPlotViewById} from '../PlotViewUtil.js';
import {dlRoot} from '../DrawLayerCntlr.js';
import {MetaConst} from '../../data/MetaConst.js';
import Catalog from '../../drawingLayers/Catalog.js';
import {CoordinateSys} from '../CoordSys.js';
import {logError} from '../../util/WebUtil.js';


/**
 * this saga does the following:
 * <ul>
 *     <li>Waits until first fits image is plotted
 *     <li>loads all the table that are catalogs
 *     <li>Then loops:
 *     <ul>
 *         <li>waits for a table new table, update, highlight or select change and then updates the drawing layer
 *         <li>waits for a new plot and adds any catalog
 *     </ul>
 * </ul>
 */
export function* watchCatalogs() {


    yield take([ImagePlotCntlr.PLOT_IMAGE, ImagePlotCntlr.PLOT_HIPS]);

    const tableGroup= get(getTableGroup(), 'tables', {});  // get the main table group.
    if (!isEmpty(tableGroup)) {
        Object.keys(tableGroup).forEach( (tbl_id) => handleCatalogUpdate(tbl_id) );
    }


    while (true) {
        const action= yield take([TABLE_LOADED, TABLE_SELECT,TABLE_HIGHLIGHT, TABLE_UPDATE,
                                  TABLE_REMOVE, ImagePlotCntlr.PLOT_IMAGE, ImagePlotCntlr.PLOT_HIPS]);
        const {tbl_id}= action.payload;
        switch (action.type) {
            case TABLE_LOADED:
                handleCatalogUpdate(tbl_id);
                break;
            
            case TABLE_SELECT:
                dispatchModifyCustomField(tbl_id, {selectInfo:action.payload.selectInfo});
                break;
            
            case TABLE_HIGHLIGHT:
            case TABLE_UPDATE:
                dispatchModifyCustomField(tbl_id, {highlightedRow:action.payload.highlightedRow});
                break;
                
            case TABLE_REMOVE:
                dispatchDestroyDrawLayer(tbl_id);
                break;

            case ImagePlotCntlr.PLOT_HIPS:
            case ImagePlotCntlr.PLOT_IMAGE:
                attachToAllCatalogs(action.payload.pvNewPlotInfoAry);
                break;
        }
    }
}


const isCName = (name) => (c) => c.name===name;

//todo - this fucntion should start using TableInfoUtil.getCenterColumns
function handleCatalogUpdate(tbl_id) {
    const sourceTable= getTblById(tbl_id);

    
    const {tableMeta,totalRows,tableData, request, highlightedRow,selectInfo}= sourceTable;
    



    if (!totalRows ||
        !tableMeta[MetaConst.CATALOG_OVERLAY_TYPE] ||
        (!tableMeta[MetaConst.CATALOG_COORD_COLS] && !tableMeta[MetaConst.CENTER_COLUMN])) {
        return; 
    }

    const {ignoreTables}=  getDlRoot();
    if (ignoreTables.some( (obj) => obj.tableId===tbl_id && obj.drawLayerTypeId===Catalog.TYPE_ID)) {
        return;
    }

    const cenData= tableMeta[MetaConst.CATALOG_COORD_COLS] || tableMeta[MetaConst.CENTER_COLUMN];

    var s;
    if (cenData) s= cenData.split(';');
    if (!s || s.length!== 3) return;

    const columns= {
        lonCol: s[0],
        latCol: s[1],
        csys : CoordinateSys.parse(s[2])
    };

    if (!tableData.columns.find( isCName(columns.lonCol)) && !tableData.columns.find(isCName(columns.latCol))) {
        return;
    }

    const params= {
        startIdx : 0,
        pageSize : MAX_ROW,
        inclCols : `"${columns.lonCol}","${columns.latCol}","ROW_IDX"`        // column names should be in quotes
    };

    let req = cloneRequest(sourceTable.request, params);
    var dataTooBigForSelection= false;
    if (totalRows>50000) {
        const sreq = cloneRequest(sourceTable.request, {inclCols: `"${columns.lonCol}","${columns.latCol}"`});
        req = makeTableFunctionRequest(sreq, 'DecimateTable', 'heatmap',
            {decimate: serializeDecimateInfo(columns.lonCol, columns.latCol, 10000), pageSize: MAX_ROW});
        dataTooBigForSelection= true;
    }

    req.tbl_id = `cat-${tbl_id}`;

    doFetchTable(req).then(
        (tableModel) => {
            if (tableModel.tableData && tableModel.tableData.data) {
                updateDrawingLayer(tbl_id, tableModel.title,
                    tableModel.tableData, tableModel.tableMeta,
                    request, highlightedRow, selectInfo, columns, dataTooBigForSelection);
            }
        }
    ).catch(
        (reason) => {
            logError(`Failed to catalog plot data: ${reason}`, reason);
        }
    );
}

function updateDrawingLayer(tbl_id, title, tableData, tableMeta, tableRequest,
                            highlightedRow, selectInfo, columns, dataTooBigForSelection) {

    const plotIdAry= visRoot().plotViewAry.map( (pv) => pv.plotId);

    const dl= getDrawLayerById(dlRoot(),tbl_id);
    if (dl) { // update drawing layer
        dispatchModifyCustomField(tbl_id, {title, tableData, tableMeta, tableRequest,
                                           highlightedRow, selectInfo, columns,
                                           dataTooBigForSelection});
    }
    else { // new drawing layer
        const angleInRadian= isTableUsingRadians(tableMeta);
        dispatchCreateDrawLayer(Catalog.TYPE_ID,
            {catalogId:tbl_id, title, tableData, tableMeta, tableRequest, highlightedRow,
                                selectInfo, columns, dataTooBigForSelection, catalog:true,
                                angleInRadian});
        dispatchAttachLayerToPlot(tbl_id, plotIdAry);
        const dl= getDrawLayerById(dlRoot(),tbl_id);
        if (dl.supportSubgroups  &&  dl.tableMeta[SUBGROUP]) {
            plotIdAry.map( (plotId) =>  getPlotViewById(visRoot(), plotId))
                .filter( (pv) => dl.tableMeta[SUBGROUP]!==get(pv, 'drawingSubGroupId'))
                .forEach( (pv) => pv && dispatchChangeVisibility({id:dl.drawLayerId, visible:false,
                                                                  plotId:pv.plotId, useGroup:false}));
        }
    }
}


function attachToAllCatalogs(pvNewPlotInfoAry=[]) {
    dlRoot().drawLayerAry.forEach( (dl) => {
        if (dl.drawLayerTypeId===Catalog.TYPE_ID && dl.catalog) {
            pvNewPlotInfoAry.forEach( (info) => {
                dispatchAttachLayerToPlot(dl.drawLayerId, info.plotId);
                const pv= getPlotViewById(visRoot(), info.plotId);
                const pvSubGroup= get(pv, 'drawingSubGroupId');
                const tableSubGroup= dl.tableMeta[SUBGROUP];
                if (!isNil(pvSubGroup) && !isNil(tableSubGroup)  && pvSubGroup!==tableSubGroup) {
                    pv && dispatchChangeVisibility({id:dl.drawLayerId, visible:false, plotId:pv.plotId, useGroup:false});
                }
            });
        }
    });
}


