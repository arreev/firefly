/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */


import Point, {makeScreenPt,makeImagePt,pointEquals} from '../Point.js';
import * as AppDataCntlr from '../../core/AppDataCntlr.js';
import BrowserInfo, {Browser} from '../../util/BrowserInfo.js';
import DrawUtil from './DrawUtil.js';
import Color from '../../util/Color.js';
import CsysConverter, {CCUtil} from '../CsysConverter.js';
import {POINT_DATA_OBJ} from './PointDataObj.js';
import DrawOp from './DrawOp.js';
import {get, isFunction} from 'lodash';


const ENABLE_COLORMAP= false;
var drawerCnt=0;



class Drawer {


    constructor() {
        this.drawerID;
        this.drawingDef= null;
        this.data;
        this.highlightData;
        this.drawConnect= null;
        this.selectedIndexes=[];

        this.plot= null;
        this.primaryCanvas= null;
        this.selectCanvas= null;
        this.highlightCanvas= null;
        this.drawingCanceler= null;
        this.plotTaskId= null;
        //this.dataUpdater= null;
        this.isPointData= false;
        this.decimate= false;

        this.decimatedData= null;
        this.decimateDim= null;
        this.lastDecimationPt= null;
        this.lastDecimationColor= null;
        this.drawTextAry= [];
        this.textUpdateCallback= null;
        //this.highPriorityLayer= false;
        this.drawerId= drawerCnt++;
        this.deferredDrawingCompletedCB= null;
    }




    //setDataTypeHint(dataTypeHint) { this.dataTypeHint= dataTypeHint; }

    //setHighPriorityLayer(highPriorityLayer) { this.highPriorityLayer = highPriorityLayer; }


    /*
     * TODO - figure this out in the new system
     * when the image resize, like a zoom, then fire events to redraw
     * Certain types of data will need to recompute the data when the image size changes so this
     * methods disables the default automactic handling
     * By default, this property is true
     * @param h handle image changes
     */
    //setHandleImageChanges(h) { this.handleImagesChanges = h; }



    dispose() {
        this.primaryCanvas= null;
        this.selectCanvas= null;
        this.highlightCanvas= null;
        this.data= null;
        this.decimatedData= null;
    }


    setPointConnector(connector) { this.drawConnect= connector; }

    setEnableDecimationDrawing(d) { this.decimate= d; }

    setPlotChangeDataUpdater(dataUpdater) { this.dataUpdater= dataUpdater; }


    cancelRedraw() {
        if (this.drawingCanceler) {
            this.drawingCanceler();
            this.drawingCanceler= null;
        }
    }


    /**
     *
     * @param {Array} data the list of DataObj
     * @param {number[]} selectedIndexes
     * @param {WebPlot} plot
     * @param {number} width
     * @param {number} height
     * @param {DrawingDef} drawingDef
     * @param {boolean} forceUpdate
     */
    setData(data,selectedIndexes,plot,width,height,drawingDef,forceUpdate= false) {
        if (data && !Array.isArray(data)) data= [data];
        var cWidth, cHeight, dWidth, oldDWidth, dHeight;
        var oldDHeight, zfact, oldZfact, oldTestPtStr, testPtStr, pt;

        width= Math.floor(width);
        height= Math.floor(height);
        if (this.primaryCanvas) {
            cWidth= this.primaryCanvas.width;
            cHeight= this.primaryCanvas.height;
        }

        if (plot) {
            dWidth= plot.dataWidth;
            dHeight= plot.dataHeight;
            zfact= Math.round(plot.zoomFactor*100000)/100000;
            pt= CCUtil.getWorldCoords(plot,makeImagePt(1,1));
            testPtStr= pt ? pt.toString() : '';
        }

        if (this.plot) {
            oldDWidth= this.plot.dataWidth;
            oldDHeight= this.plot.dataHeight;
            oldZfact= Math.round(this.plot.zoomFactor*100000)/100000;
            pt= CCUtil.getWorldCoords(this.plot,makeImagePt(1,1));
            oldTestPtStr= pt ? pt.toString() : '';
        }

        var viewUpdated= true;

        if (drawingDef===this.drawingDef &&
            cWidth===width && cHeight===height &&
            dWidth===oldDWidth && dHeight===oldDHeight  &&
            zfact===oldZfact  && testPtStr===oldTestPtStr ) {
            viewUpdated= false;
        }
        
        
        const primaryUpdated= (data && data!==this.data) || viewUpdated;

        const selectedUpdated= (selectedIndexes!==this.selectedIndexes) || viewUpdated;


        if (data!==this.data || this.plot!==plot) {
            this.decimatedData= null;
        }

        if (!primaryUpdated && !selectedUpdated  && !forceUpdate) return;

        this.plot= plot;
        this.data = data;
        this.selectedIndexes = selectedIndexes;
        this.drawingDef = drawingDef;
        

        if (primaryUpdated || forceUpdate) {
            this.dataUpdated(width,height);
        }

        if (selectedUpdated || forceUpdate) {
            this.updateDataSelectLayer();
        }




        //======== DEBUG =============================
        //var changes= [this.primaryCanvas ? '' : 'no canvas'];
        // var changes= [];
        // if (data!==this.data) changes.push('data');
        // //if (plot!==this.plot) changes.push('plot');
        // if (cWidth!==width ) changes.push(`width: ${width}, ${cWidth}`);
        // if (cHeight!==height ) changes.push(`height: ${height}, ${cHeight}`);
        // if (dWidth!==oldDWidth ) changes.push(`data width: ${oldDWidth}, ${dWidth}`);
        // if (dHeight!==oldDHeight ) changes.push(`data height: ${oldDHeight}, ${dHeight}`);
        // if (zfact!==oldZfact ) changes.push(`zoom factor ${oldZfact}, ${zfact}`);
        // if (testPtStr!==oldTestPtStr ) changes.push('test pt');
        // if (oldvpY!==vpY ) changes.push(`vpY: ${oldvpY}, ${vpY}`);
        // if (drawingDef!==this.drawingDef ) changes.push('drawingDef');
        // var changeStr= join(', ',...changes);
        // if (false && primaryUpdated) console.log(`Drawer ${this.drawerId}: redraw- changes: ${changeStr}`);
        //=====================================
        
    }



    setPrimCanvas(c, width, height) {
        if (c && c!==this.primaryCanvas) {
            this.primaryCanvas= c;
            //console.log(`Drawer ${this.drawerId}: redraw primary- canvas update`);
            this.dataUpdated(width,height);
        }
    }

    setHighlightCanvas(c, width, height) {
        if (c && c!==this.highlightCanvas) {
            this.highlightCanvas= c;
            //console.log(`Drawer ${this.drawerId}: redraw highlight- canvas update`);
            updateCanvasSize(width,height,c);
            this.updateDataHighlightLayer(this.highlightData, width, height);
        }
    }

    setSelectCanvas(c, width, height) {
        if (c && c!==this.selectCanvas) {
            this.selectCanvas= c;
            //console.log(`Drawer ${this.drawerId}: redraw select- canvas update`);
            updateCanvasSize(width,height,c);
            this.updateDataSelectLayer();
        }
    }

    dataUpdated(width,height) {
        if (!this.primaryCanvas) return;
        this.cancelRedraw();
        updateCanvasSize(width,height,this.primaryCanvas,this.selectCanvas,this.highlightCanvas);
        if (this.data && this.data.length>0) {
            this.redraw();
        }
        else {
            this.clear();
        }
    }

    updateDataSelectLayer() {
        var {plot,selectCanvas,selectedIndexes,data}= this;
        var sCtx= selectCanvas? selectCanvas.getContext('2d') : null;

        if (sCtx) {
            var cc= CsysConverter.make(plot);
            this.redrawSelected(selectCanvas, sCtx, cc, data, selectedIndexes,
                selectCanvas.width, selectCanvas.height);
        }
    }

    /**
     *
     * @param highlightData
     * @param width
     * @param height
     */
    updateDataHighlightLayer(highlightData,width,height) {
        var {highlightCanvas, drawingDef}= this;
        this.highlightData=highlightData;
        var hCtx= highlightCanvas? highlightCanvas.getContext('2d') : null;
        if (hCtx)  {
            var cc= CsysConverter.make(this.plot);
            this.redrawHighlight(hCtx, cc, highlightData,drawingDef,width,height);
        }
    }



    clearSelectLayer() { DrawUtil.clearCanvas(this.selectCanvas); }

    clearHighlightLayer() { DrawUtil.clearCanvas(this.highlightCanvas); }

    clear() {
        this.cancelRedraw();
        var {primaryCanvas,selectCanvas,highlightCanvas}= this;

        DrawUtil.clearCanvas(primaryCanvas);
        DrawUtil.clearCanvas(selectCanvas);
        DrawUtil.clearCanvas(highlightCanvas);
        this.removeTask();
    }


    redraw() {
        var {primaryCanvas,selectCanvas,highlightCanvas}= this;
        if (!primaryCanvas) return;

        var pCtx= primaryCanvas ? primaryCanvas.getContext('2d') : null;
        var sCtx= selectCanvas? selectCanvas.getContext('2d') : null;
        var hCtx= highlightCanvas? highlightCanvas.getContext('2d') : null;

        var w= primaryCanvas.width;
        var h= primaryCanvas.height;

        var cc= CsysConverter.make(this.plot);
        this.redrawPrimary(primaryCanvas, cc, pCtx, this.data, this.drawingDef, w,h);
        this.redrawHighlight(hCtx, cc, this.highlightData, this.drawingDef, w,h);
        this.redrawSelected(selectCanvas, sCtx, cc, this.data, this.selectedIndexes, w, h);
    }


    redrawSelected(selectCanvas, ctx, cc, data, selectedIndexes,w,h) {
        if (!ctx) return;
        var selectedData;
        DrawUtil.clear(ctx,w,h);
        if (canDraw(ctx,data)) {
            if (!selectedIndexes) return;
            if (Array.isArray(selectedIndexes)) {
                selectedData= this.decimateData(this.decimate,
                             selectedIndexes.map( (dataIdx)=> data[dataIdx]), cc,false,null);
            }
            else if (typeof selectedIndexes ==='function') {
                selectedData= this.decimateData(this.decimate,
                    data.filter( (d,idx)=> selectedIndexes(idx)), cc,false,null);
            }
            else {
                return;
            }
            const selDrawDef= Object.assign({}, this.drawingDef, {color:this.drawingDef.selectedColor});
            const params= makeDrawingParams(selectCanvas, null, selDrawDef, cc,selectedData,null, Number.MAX_SAFE_INTEGER);
            this.doDrawing(params);
        }
    }

    /**
     * 
     * @param ctx
     * @param cc
     * @param highlightData
     * @param drawingDef
     * @param w
     * @param h
     */
    redrawHighlight(ctx, cc, highlightData,drawingDef,w,h) {
        if (!ctx) return;
        DrawUtil.clear(ctx,w,h);
        if (!highlightData || !highlightData.length) return;

        if (canDraw(ctx,highlightData)) {
            var sPtM= makeScreenPt(0,0);
            highlightData.forEach( (pt) => drawObj(ctx, [], drawingDef, cc, pt, sPtM, false) );
        }
    }



    redrawPrimary(canvas, cc, ctx, data, drawingDef,w,h) {
        if (!ctx) return;
        this.clear();
        var params;
        if (canDraw(ctx,data)) {
            this.decimatedData= this.decimateData(this.decimate, data, cc,true,this.decimatedData);
            var drawData= this.decimatedData;
            this.drawTextAry= [];
            if (drawData.length>500) {
                params= makeDrawingParams(canvas, this.drawTextAry, drawingDef,cc,drawData,
                                         this.drawConnect, getMaxChunk(drawData,this.isPointData),
                                         this.deferredDrawingCompletedCB);
                this.cancelRedraw();
                this.drawingCanceler= makeDrawingDeferred(this,params);
                this.removeTask();
                if (drawData.length>15000) this.addTask();
            }
            else {
                params= makeDrawingParams(canvas, this.drawTextAry, drawingDef,
                                          cc,drawData,this.drawConnect, Number.MAX_SAFE_INTEGER);
                this.doDrawing(params);
                if (this.textUpdateCallback) this.textUpdateCallback(this.drawTextAry);
            }
        }
        else {
            this.removeTask();
        }
    }


    //decimateData(decimate, inData, cc, useColormap, oldDecimatedData) {
    //    if (decimate && inData.length>150) {
    //        return this.doDecimateData(inData,oldDecimatedData,cc,useColormap);
    //    }
    //    else {
    //        return inData;
    //    }
    //}

    /**
     *
     * @param decimate
     * @param inData
     * @param cc
     * @param useColormap
     * @param oldDecimatedData
     * @return {*}
     */
    decimateData(decimate, inData, cc, useColormap, oldDecimatedData) {
        if (!decimate || inData.length<=150) return inData;

        var retData= inData;
        var dim = cc.viewDim;
        var spt= cc.getScreenCoords(makeScreenPt(0,0));
        var defCol= this.drawingDef.color;
        if (!oldDecimatedData ||
            dim.width!==this.decimateDim.width ||
            dim.height!==this.decimateDim.height ||
            defCol!==this.lastDecimationColor ||
            !pointEquals(spt,this.lastDecimationPt))  {
            retData= doDecimation(inData, cc, useColormap);
            this.lastDecimationColor= defCol;
            this.lastDecimationPt=spt;
            this.decimateDim= dim;
        }
        else if (oldDecimatedData) {
            retData= oldDecimatedData;
        }
        return retData;
    }



    doDrawing(params) {
        if (params.begin) {
            params.begin= false;
            params.canvas.style.visibility= 'hidden';
            params.done= false;
        }
        else {
            params.deferCnt++;
        }


        if (!params.done) {
            if (params.drawConnect) params.drawConnect.beginDrawing();
            var nextChunk= getNextChuck(params);
            if (nextChunk.optimize) {
                drawChunkOptimized(nextChunk.drawList, params, params.ctx);
                params.opCnt++;
            }
            else {
                drawChunkNormal(nextChunk.drawList, params, params.ctx);
            }
            if (params.next.done) { //loop finished
                params.done= true;
                this.removeTask();
            }
        }

        if (params.done ) {
            params.canvas.style.visibility= 'visible';
            if (params.drawConnect) {
                params.drawConnect.endDrawing();
            }
            if (params.deferCnt && isFunction(params.deferredDrawingCompletedCB)) params.deferredDrawingCompletedCB();
        }

    }


    removeTask() {
        var {plot,plotTaskId}= this;
        if (plot && plotTaskId) {
            setTimeout( () => AppDataCntlr.dispatchRemoveTaskCount(plot.plotId,plotTaskId) ,0);
            this.plotTaskId= null;
        }
    }

    addTask() {
        var {plot}= this;
        if (plot) {
            this.plotTaskId= AppDataCntlr.makeTaskId();
            setTimeout( () => AppDataCntlr.dispatchAddTaskCount(plot.plotId,this.plotTaskId) ,0);
        }
    }




    static makeDrawer() {
        return new Drawer();
    }

}



//=======================================================================
//------------------ private static functions
//------------------ private static functions
//------------------ private static functions
//------------------ private static functions
//------------------ private static functions
//=======================================================================


function nextPt(i,fuzzLevel, max) {
    i= Math.trunc(i);
    var remainder= i%fuzzLevel;
    var retval= (remainder===0) ? i : i+(fuzzLevel-remainder);
    if (retval===max) retval= max-1;
    return retval;
}

/**
 *
 * @param canvas
 * @param drawTextAry
 * @param drawingDef
 * @param csysConv
 * @param data
 * @param drawConnect
 * @param maxChunk
 * @param {Function} deferredDrawingCompletedCB - called when drawing has completed
 * @return {Object}
 */
function makeDrawingParams(canvas, drawTextAry, drawingDef, csysConv, data, drawConnect, maxChunk, deferredDrawingCompletedCB) {
    var params= {
        canvas,    //const
        ctx : canvas.getContext('2d'),    //const
        drawTextAry,
        drawingDef,    //const
        csysConv,    //const
        data,    //const
        maxChunk,    //const
        drawConnect,    //const
        iterator : data[Symbol.iterator](),
        startTime: Date.now(),    //const
        opCnt: 0, //only for debug
        done : false,
        begin : true,
        deferCnt: 0,
        vpPtM : makeScreenPt(0,0), //const
        deferredDrawingCompletedCB //const
    };
    params.next= params.iterator.next();
    return params;
}

/**
 *
 * @param {{x:number,y:number,type:string}} pt
 * @param {{x:number,y:number,type:string}}mSpPt
 * @param {CysConverter} cc
 * @return {*}
 */
function getScreenCoords(pt, mSpPt, cc) {
    var retval;
    if (pt.type===Point.W_PT) {
        var success= cc.getScreenCoordsOptimize(pt,mSpPt);
        retval= success ? mSpPt : null;
    }
    else {
        retval= cc.getScreenCoords(pt);
    }
    return retval;
}



function makeDrawingDeferred(drawer,params) {
    var id= window.setInterval( () => {
        if (params.done) {
            window.clearInterval(id);
            if (drawer.textUpdateCallback) drawer.textUpdateCallback(drawer.drawTextAry);
        }
        drawer.doDrawing(params);
    },0);
    return () => window.clearInterval(id);
}

/**
 *
 * @param ctx canvas context object
 * @param drawTextAry
 * @param def DrawingDef
 * @param csysConv web csysConv
 * @param obj DrawObj
 * @param vpPtM mutable viewport point
 * @param {boolean} onlyAddToPath
 */
function drawObj(ctx, drawTextAry, def, csysConv, obj, vpPtM, onlyAddToPath) {
    DrawOp.draw(obj, ctx, drawTextAry, csysConv, def, vpPtM,onlyAddToPath);
}

/**
 * An optimization of drawing.  Check is the Object is a PointDataObj (most common and simple) and then checks
 * it is draw on a WebPlot, and if it is in the drawing area.
 * Otherwise it will always return true
 * @param {CysConverter} csysConv the WebPlot to draw on
 * @param obj the DrawingObj to check
 * @return {boolean} true is it should be drawn
 */
function shouldDrawObj(csysConv, obj) {
    if (!obj) return false;
    let retval= true;
    if (csysConv && obj.pt && obj.type===POINT_DATA_OBJ) {
        retval= csysConv.pointOnDisplay(obj.pt);
    }
    return retval;
}

function canDraw(ctx,data) {
    return (ctx && data && data.length);
}

/**
 *
 * @param ctx canvas object
 * @param def DrawingDef
 * @param csysConv web csysConv
 * @param dc drawConnector
 * @param obj DrawObj
 * @param lastObj DrawObj
 */
function drawConnector(ctx, def, csysConv, dc, obj, lastObj) {
    if (!obj && !lastObj) return;
    if (csysConv) {
        var wp1= csysConv.getWorldCoords(DrawOp.getCenterPt(lastObj));
        var wp2= csysConv.getWorldCoords(DrawOp.getCenterPt(obj));
        if (!csysConv.coordsWrap(wp1,wp2)) {
            dc.draw(ctx,csysConv,def, wp1,wp2);
        }
    }
    else {
        if (DrawOp.getCenterPt(lastObj).type===Point.SPT && DrawOp.getCenterPt(obj).type===Point.SPT) {
            dc.draw(ctx,def, DrawOp.getCenterPt(lastObj), DrawOp.getCenterPt(obj));
        }
    }
}

function drawChunkOptimized(drawList, params, ctx) {
    if (!drawList.length) return;
    DrawUtil.beginPath(ctx,params.drawingDef.color,params.drawingDef.lineWidth);
    for(var obj of drawList) {
        drawObj(ctx, params.drawTextAry, params.drawingDef, params.csysConv, obj,params.vpPtM, true);
    }
    DrawUtil.stroke(ctx);
}

function drawChunkNormal(drawList, params, ctx) {
    var lastObj= null;
    var {drawingDef,drawConnect,csysConv,vpPtM}= params;
    for(var obj of drawList) {
        if (drawConnect) { // in this case doDraw was already called
            drawObj(ctx, params.drawTextAry, drawingDef, csysConv, obj,vpPtM, false);
        }
        else  {
            if (shouldDrawObj(csysConv,obj)) { // doDraw must be call when there is a connector
                drawObj(ctx, params.drawTextAry, drawingDef, csysConv, obj,vpPtM, false);
                if (drawConnect) {
                    drawConnector(ctx,drawingDef,csysConv,drawConnect,obj,lastObj);
                }
            }
            lastObj= obj;
        }
    }
}


function getNextChuck(params) {
    var drawList= [];
    var optimize= params.drawConnect?false:true;
    var objLineWidth;
    var objColor;
    var {drawingDef}= params;
    var i;


    var obj= params.next.value;
    var color= drawingDef.color;
    var lineWidth=  get(obj,'lineWidth',false) || drawingDef.lineWidth || 1;

    for(i= 0; (!params.next.done && i<params.maxChunk ); ) {
        obj= params.next.value;
        params.next= params.iterator.next();
        if (!params.drawConnect) {
            if (shouldDrawObj(params.csysConv, obj)) {
                drawList.push(obj);
                if (optimize) {
                    objLineWidth= obj.lineWidth || lineWidth;
                    objColor= obj.color || color;
                    optimize= (DrawOp.usePathOptimization(obj,drawingDef) &&
                               lineWidth===objLineWidth &&
                               color===objColor);
                }
                i++;
            }
        }
        else if (obj) {
            drawList.push(obj);
            i++;
        }
    }
    return {drawList,optimize};
}




function getMaxChunk(drawData,isPointData) {
    var maxChunk= 1;
    if (!drawData.length) return maxChunk;
    if (isPointData) {
        maxChunk= BrowserInfo.isBrowser(Browser.SAFARI) || BrowserInfo.isBrowser(Browser.CHROME) ? 2000 : 500;
    }
    else {
        maxChunk= BrowserInfo.isBrowser(Browser.SAFARI) || BrowserInfo.isBrowser(Browser.CHROME) ? 1000 : 200;
    }
    return maxChunk;
}


function updateCanvasSize(w,h,...cAry) {
    cAry.forEach( (c) => {
        if (!c) return;
        c.width= w;
        c.height= h;
    });
}

function makeColorMap(mapSize,color) {
    return Color.makeSimpleColorMap(color,mapSize,true);
}


function setupColorMap(data, maxEntry) {
    var colorMap= makeColorMap(maxEntry);
    if (colorMap)  {
        var cnt;
        var obj= null;
        var idx;
        if (maxEntry>colorMap.length) {
            var maxCnt = maxEntry+1; // to include draw obj with cnt==maxEntry into the last color band
            for(obj of data) {
                cnt= obj.representCnt || 1;
                idx = cnt*colorMap.length/maxCnt;
                obj.color=colorMap[idx];
            }
        }  else {
            for(obj of data) {
                cnt= obj.representCnt || 1;
                //if (cnt>colorMap.length) cnt=colorMap.length;
                obj.color=colorMap[cnt-1];
            }
        }
    }
}

function doDecimation(inData, cc, useColormap) {
    var i,j;
    var dim = cc.viewDim;

    var supportCmap= useColormap && ENABLE_COLORMAP;

    //var drawArea= dim.width*dim.height;
    //var percentCov= inData.length/drawArea;

    var fuzzLevel= 5;
    //var start = Date.now();

    var {width,height}= dim;

    var decimateObs= new Array(width);
    for(i=0; (i<decimateObs.length);i++) decimateObs[i]= new Array(height);

    var seedPt= makeScreenPt(0,0);
    var sPt;
    var pt;
    var maxEntry= -1;
    var entryCnt;

//        GwtUtil.getClientLogger().log(Level.INFO,"doDecimation: " + (enterCnt++) + ",data.size= "+ _data.size() +
//                ",drawID="+drawerID+
//                ",data="+Integer.toHexString(_data.hashCode()));

    var first200= [];
    var decimatedAddCnt= 0;
    var totalInViewPortCnt= 0;

    for(var obj of inData) {
        if (obj) {
            pt= DrawOp.getCenterPt(obj);
            if (pt.type===Point.W_PT) {
                sPt= cc.pointInPlotRoughGuess(pt) ? getScreenCoords(pt,seedPt,cc) : null;
            }
            else {
                sPt= getScreenCoords(pt,seedPt,cc);
            }

        }
        else {
            sPt= null;
        }
        if (sPt) {
            i= nextPt(sPt.x,fuzzLevel,width);
            j= nextPt(sPt.y, fuzzLevel,height);
            if (i>=0 && j>=0 && i<width && j<height) {
                if (!decimateObs[i][j]) {
                    decimateObs[i][j]= supportCmap ? Object.assign({},obj) : obj;
                    if (supportCmap) {
                        decimateObs[i][j].representCnt= obj.representCnt;
                        entryCnt= decimateObs[i][j].representCnt;
                        if (entryCnt>maxEntry) maxEntry= entryCnt;
                    }
                    decimatedAddCnt++;
                }
                else {
                    if (supportCmap) {
                        decimateObs[i][j].representCnt+=(obj.representCnt||1);
                        entryCnt= decimateObs[i][j].representCnt;
                        if (entryCnt>maxEntry) maxEntry= entryCnt;
                    }
                }
                if (totalInViewPortCnt<200) first200.push(obj);
                totalInViewPortCnt++;
            }
        }
    }

    var retData;
    if (totalInViewPortCnt<200) {
        retData= first200;
    }
    else {
        retData= [];
        for(i= 0; (i<decimateObs.length); i++) {
            for(j= 0; (j<decimateObs[i].length); j++) {
                if (decimateObs[i][j]) retData.push(decimateObs[i][j]);
            }
        }
    }



    if (supportCmap) setupColorMap(retData,maxEntry);

    return retData;
}

export default Drawer;
