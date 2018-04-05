/**
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 * Lijun
 * 4/14/16
 *
 * 05/02/17
 * DM-6501
 *   Add the Regrid and LinearInterpolator in the utility
 *   Use the Regrid to adding more gride lines and the number of the points in the line
 *   Rotate the labels based on the angles
 *   Cleaned up the codes
 *
 *  02/13/18
 *  IRSA-1391
 *  Add the Grid lines for HiPs map
 *  The HiPs map is very different from the regular image map.  The HiPs map has a huge range in longitude.
 *  In order to get finer grid calculation, I adjust the algorithm, instead of using 4 for the intervals when finding lines,
 *  I used interval = range/0.5.
 *
 *  For HiPs map, the range is calculated differently.  The range calculation for regular image is using the data width
 *  data height (which are the naxis1/naxis2).  But for HiPs, there are no data width and height.  Intead of using
 *  data width and height by walking four corners, I used the field view angle and ceter point.  If the angle is >=180,
 *  the image is sphere.  If the fov < 180, the ranges are calculated using the cell values.
 *
 *  For the grid lines calculation, all valid points are in the range, there is no need to check with the screen width.
 *
 *  For Label location, for the latitude, the labels lay on the longitude=0.
 *
 *  3/2/18
 *  It takes a lot longer to draw the grid lines when zoom out.  The reason is that the points search taken longer and
 *  and the number of points gotten larger and larger.  To solve this issue, I changed the algorithm as described belowL
 *  1. Use larger delta=10 (deg) to calculate line intervals and points.
 *  2. Use the view port as a frame to test if the grid lines need to be reduced.  When the image is zoomed out, the ranges
 *  are getter smaller and smaller, the grid lines will fall beyond the boundary.  Those lines no longer need to be calculated.
 *  Thus, there will be some time saved.
 *  However, the view port does not always return the correct border because it depends on the image's projection etc.
 *     a. Four corner is calculated.  If four corners found, the image is in zoom out status, the four corners are on the
 *        image.  In this case, the grid lines are filtered.  Only the lines are inside the image will be calculated.
 *
 *     b. If the image has discontinuity i.e (0, 360 next each other in longitude), I test in the range 0- min view port
 *        and max viewport - 360.  The way to test the discontinuity is to see 0 is in the levels array.
 *
 *     c. If the image is smaller than the view port, there is no grid line filtering needed.
 *
 *
 *
 *   3. Since the HIPS map has a larger range,  the levels are always regrided to finer grid lines.
 *
 *   4. Filter the grid lines
 *     a.  For HIPS map, I checked if the longitude=0 is in the image, that is if the pole is inside the viewable area.
 *     If it is, the algorithm is different from the no pole case.
 *
 *     b. In order to get more than one grid lines,  the original levels are regrided to finer levels if only one
 *     line is found, using the new finer levels to look for the grid lines.  The maximum tries is 5 times.
 *
 *
 */

import { makeWorldPt, makeImagePt,makeDevicePt,makeImageWorkSpacePt} from '../visualize/Point.js';
import ShapeDataObj from '../visualize/draw/ShapeDataObj.js';
import CoordinateSys from '../visualize/CoordSys.js';
import CoordUtil from '../visualize/CoordUtil.js';
import numeral from 'numeral';
import { getDrawLayerParameters} from './WebGrid.js';
import { getPointMaxSide} from '../visualize/HiPSUtil.js';
import {Regrid} from '../util/Interp/Regrid.js';
import {getCenterOfProjection} from '../visualize/PlotViewUtil';

const precision3Digit = '0.000';
const RANGE_THRESHOLD = 1.02;
const minUserDistance= 0.25;   // user defined max dist. (deg)
const maxUserDistance= 3.00;   // user defined min dist. (deg)
var userDefinedDistance = false;
const angleStepForHipsMap=4.0;

/**
 *
 * @param plot
 * @param cc
 * @param useLabels
 * @param numOfGridLines
 */
export function makeGridDrawData (plot,  cc, useLabels, numOfGridLines=11){


    const {width,height, screenWidth, csys, labelFormat} = getDrawLayerParameters(plot);

    const wpt = cc.getWorldCoords(makeImageWorkSpacePt(1, 1), csys);
    const aitoff = (!wpt);

    const centerWpt = getCenterOfProjection(plot);

    //TODO the center point does not work for image, need to find out why
    //TODO the angle does not look good using the center project as a label location.
    //The rotation angle can be used in that case...
    //TODO there was small issue about lines in centain cases, there should be something wrong with the logic
    const centerScreenPt = plot.plotType==='hips'?cc. getImageWorkSpaceCoords(centerWpt , plot.zoom):null;

    if (width > 0 && height >0) {
        const bounds = new Rectangle(0, 0, width, height);

        //for hips, the range is always the full sky ranges
        var range = plot.plotType==='hips'?[[0, 360], [-90, 90]]: getRange(csys, width, height, cc);
        const {xLines, yLines, labels} = plot.plotType==='hips'?computeHipGridLines(cc, csys, range,screenWidth, numOfGridLines, labelFormat,plot)
        :computeImageGridLines(cc, csys, range,screenWidth, numOfGridLines, labelFormat,plot);

        return  drawLines(bounds, labels, xLines, yLines, aitoff, screenWidth, useLabels, cc,plot.plotType,  centerScreenPt);
    }
}

function getViewPortInfo(plot, csys,  cc,ranges) {

    const {width, height} = plot.viewDim;
    const corners = [
    cc.getWorldCoords(makeDevicePt(0, 0), csys),
    cc.getWorldCoords(makeDevicePt(0, height), csys),
    cc.getWorldCoords(makeDevicePt(width, height), csys),
    cc.getWorldCoords(makeDevicePt(width, 0), csys)

    ];

    if( corners.indexOf(null)>-1)  return {corners:undefined, ranges};

    var vals = [];
    var viewBorder = [[1.e20, -1.e20], [1.e20, -1.e20]];
    for (let i = 0; i < corners.length; i++) {
        if (corners[i]) {
            vals[0] = corners[i].getLon();
            vals[1] = corners[i].getLat();
            //assign the new lower and upper longitude if found
            if (vals[0] < viewBorder[0][0]) viewBorder[0][0] = vals[0];
            if (vals[0] > viewBorder[0][1]) viewBorder[0][1] = vals[0];

            //assign the new lower and upper latitude if found
            if (vals[1] < viewBorder[1][0]) viewBorder[1][0] = vals[1];
            if (vals[1] > viewBorder[1][1]) viewBorder[1][1] = vals[1];
        }
    }
    return {corners, viewBorder};
}


function getLonInRightRange(lon){
    if (lon<0) {
        return 360+lon;
    }
    else if(lon>360) {
        return lon-360;
    }
    else return lon;
}


function getLatInRightRange(lat){
    if (lat<-90) {
        return 180+lat;
    }
    else  if(lat>90){
        return 180-90;
    }
    else return lat;
}

function getPoleInfo(cc, csys){


  var poles=0;


  var  sharedLon= 0.0;
  var  sharedLat= 90.0;
  if (cc.pointOnDisplay(makeWorldPt(sharedLon, sharedLat, csys))){
    poles += 1;

  }

  sharedLon= 0.0;
  sharedLat= -90.0;
  if (cc.pointOnDisplay(makeWorldPt(sharedLon, sharedLat, csys))){
    poles += 2;

  }

  return poles;
}

/**
 * Define a rectangle object
 * @param x - the x coordinate
 * @param y - the y coordinate
 * @param width - the width of the rectangle
 * @param height - the height of the rectangle
 * @constructor
 */
function Rectangle(x, y, width, height){
    this.x= x;
    this.y= y;
    this.width= width;
    this.height= height;
}
/**
 * walk around four corners to find the proper ranges
 * @param intervals
 * @param x0 - the starting x value of the corner
 * @param y0 - the starting y value of the corner
 * @param dx - the increment
 * @param dy - the increment
 * @param range - the x and y rangs
 * @param csys - the coordinate system the grid is drawing with
 * @param wrap - boolean value, true or false
 * @param cc - the CoordinateSys object
 */
function edgeRun (intervals,x0,  y0,dx, dy, range, csys, wrap, cc) {


    var  x = x0;
    var  y = y0;

    var i = 0;
    var vals = [];

    while (i <= intervals) {

        var wpt =cc.getWorldCoords(makeImageWorkSpacePt(x,y),csys);
        //look for lower and upper longitude and latitude
        if (wpt) {
            vals[0] = wpt.getLon();
            vals[1] = wpt.getLat();

            if (wrap && vals[0] > 180) vals[0] = vals[0]-360;
            if (wrap!=null //true and false
                || wrap==null && cc.pointInPlot(makeWorldPt(vals[0], vals[1], csys))//temporary solution to solve the problem reported.
                || csys===CoordinateSys.EQ_B2000
                || csys===CoordinateSys.EQ_J2000 ) {
                //assign the new lower and upper longitude if found
                if (vals[0] < range[0][0]) range[0][0] = vals[0];
                if (vals[0] > range[0][1]) range[0][1] = vals[0];

                //assign the new lower and upper latitude if found
                if (vals[1] < range[1][0]) range[1][0] = vals[1];
                if (vals[1] > range[1][1]) range[1][1] = vals[1];
            }
        }
        x += dx;
        y += dy;

        ++i;
    }

}
/**
 * get the ege value
 * @param intervals - the number of interval to look for the values
 * @param width - the width of the image
 * @param height - the height of the image
 * @param csys - the coordinate system the grid is drawing with
 * @param wrap - boolean value, true or false
 * @param cc - the CoordinateSys object
 * @returns the value at the edge of the image
 */
function  edgeVals(intervals, width, height, csys, wrap, cc) {

    var  range = [[1.e20, -1.e20],[1.e20, -1.e20]];
    var xmin=0;
    var ymin=0;
    var xmax= width;
    var ymax= height;
    var xdelta, ydelta, x, y;


    // four corners.
    // point a[xmin, ymax], the top left point, from here toward to right top point b[xmax, ymax], ie the line
    //   a-b
    y = ymax;
    x = xmin;
    xdelta = (width / intervals) - 1; //define an interval of the point in line a-b
    ydelta = 0; //no change in the y direction, so ydelta is 0, thus the points should be alone line a-b
    edgeRun(intervals, x, y, xdelta, ydelta, range, csys,wrap, cc);

    // Bottom: right to left
    y = ymin;
    x = xmax;
    xdelta = -xdelta;
    edgeRun(intervals, x, y, xdelta, ydelta, range, csys, wrap, cc);

    // Left.  Bottom to top.
    xdelta = 0;
    ydelta = (height / intervals) - 1;
    y = ymin;
    x = xmin;
    edgeRun(intervals, x, y, xdelta, ydelta, range, csys, wrap, cc);

    // Right. Top to bottom.
    ydelta = -ydelta;
    y = ymax;
    x = xmax;
    edgeRun(intervals, x, y, xdelta, ydelta, range, csys, wrap,cc);

    // grid in the middle
    xdelta = (width / intervals) - 1;
    ydelta = (height / intervals) - 1;
    x = xmin;
    y = ymin;
    edgeRun(intervals, x, y, xdelta, ydelta, range, csys, wrap,cc);


    return range;
}
/**
 * Test to see if the edge is a real one
 * @param xrange
 * @param trange
 * @returns {boolean}
 */
function  testEdge( xrange, trange)
{
    /* This routine checks if the experimental minima and maxima
     * are significantly changed from the old test minima and
     * maxima.  xrange and trange are assumed to be multidimensional
     * extrema of the form double[ndim][2] with the minimum
     * in the first element and the maximum in the second.
     *
     * Note that xrange is modified to have the most extreme
     * value of the test or old set of data.
     */

    /* Find the differences between the old data */
    var delta =trange.map( (t)=>{
        return Math.abs(t[1]-t[0]);
    });

    for (let i=0; i<trange.length; i += 1) {
        var  ndelta = Math.abs(xrange[i][1]-xrange[i][0]);
        /* If both sets have nil ranges ignore this dimension */
        if (ndelta <= 0. && delta[i] <= 0.){
            continue;
        }

        /* If the old set was more extreme, use that value. */
        if (xrange[i][0] > trange[i][0]) {
            xrange[i][0] = trange[i][0];
        }

        if (xrange[i][1] < trange[i][1])  {
            xrange[i][1] = trange[i][1];
        }

        /* If the previous range was 0 then we've got a
         * substantial change [but see above if both have nil range]
         */
        if (!delta[i]) {
            return false;
        }

        /* If the range has increased by more than 2% than */
        if ( Math.abs(xrange[i][1]-xrange[i][0])/delta[i] >
            RANGE_THRESHOLD) {
            return false;
        }
    }

    return true;
}
/**
 * Get the line ranges
 * @param  {object} csys - the coordinate system the grid is drawing with
 * @param  {double} width - the width of the image
 * @param  {double} height - the height of the image
 * @param  {object} cc - the CoordinateSys object
 * @returns the range array
 */
function getRange( csys, width, height, cc) {


    //using three state boolean wrap=true/false, the null, null for temporary solution to solve the problem reported.
    var range=[[0.,0.],[0.,0.]];
    var poles=0;  // 0 = no poles, 1=north pole, 2=south pole, 3=both

    /* Get the range of values for the grid. */
    /* Check for the poles.  We allow the poles to
     * be a pixel outside of the image and still consider
     * them to be included.
     */

    var wrap=false;	/* Does the image wrap from 360-0. */
    var  sharedLon= 0.0;
    var  sharedLat= 90.0;



    if (cc.pointInPlot(makeWorldPt(sharedLon, sharedLat, csys))){
        range[0][0] = -179.999;
        range[0][1] =  179.999;
        range[1][1] = 90;
        poles += 1;
        wrap = true;
    }

    sharedLon= 0.0;
    sharedLat= -90.0;
    if (cc.pointInPlot(makeWorldPt(sharedLon, sharedLat, csys))){
        range[0][0] = -179.999;
        range[0][1] =  179.999;
        range[1][0] = -90;
        poles += 2;
        wrap = true;
    }

    /* If both poles are included we can just return */
    if (poles == 3){
        return range;
    }

    /* Now we have to go around the edges to find the remaining
     * minima and maxima.
     */
    var  trange = edgeVals(1, width, height, csys, wrap,cc);
    if (!wrap) {
        /* If we don't have a pole inside the map, check
         * to see if the image wraps around.  We do this
         * by checking to see if the point at [0, averageDec]
         * is in the map.
         */

        sharedLon= 0.0;
        sharedLat= (trange[1][0] + trange[1][1])/2;

        //this block is modified to fix the issue reported that only one line is drawn in some case
        if (cc.pointInPlot(makeWorldPt(sharedLon, sharedLat, csys)))
        {
            wrap = true;

            // Redo min/max
            trange = edgeVals(1, width, height, csys,wrap,cc);
        }
        //this block was a temporary solution to solve the problem once found.  Comment it out for now in js version
        else if (csys===CoordinateSys.GALACTIC){
            trange=edgeVals(1, width, height, csys,null,cc);
            sharedLon = 0.0;
            sharedLat = (trange[1][0] + trange[1][1]) / 2;
            if (cc.pointInPlot(makeWorldPt(sharedLon, sharedLat,csys))) {
                wrap=true;
                // Redo min/max
                trange =edgeVals(1, width, height, csys,wrap,cc);
            }
        }

    }

    var  xrange = trange;
    var xmin= 0;
    var  adder=2;
    for (let  intervals = xmin+adder;
         intervals < width; intervals+= adder) {

        xrange = edgeVals(intervals, width, height, csys,wrap,cc);
        if (testEdge(xrange, trange)) {
            break;
        }
        trange = xrange;
        adder*= 2;
    }

    if (!poles && wrap){
        xrange[0][0] += 360;
    }
    else if (poles ===1) {
        range[1][0] = xrange[1][0];
        return range;
    }
    else if (poles === 2)  {
        range[1][1] = xrange[1][1];
        return range;
    }

    return xrange;

}

function lookup(val, factor){

    const conditions=[val < 1,val > 90,val > 60 ,val > 30,val > 23,val > 18,val > 6, val > 3];
    const index = conditions.indexOf(true);
    var values =[val, 30, 20, 10, 6, 5, 2, 1] ;
    var retval = (index && index>=0)? values[index] :0.5;
    if (factor >=4.0) {
        retval = retval/2.0;
    }
    return retval;



}

function calculateDelta(min, max,factor){
    var delta = (max-min)<0? max-min+360:max-min;
    var  qDelta;
    if (delta > 1) { // more than 1 degree
        qDelta = lookup(delta,factor);
    }
    else if (60*delta > 1) {// more than one arc minute
        qDelta = lookup(60*delta,factor)/60;
    }
    else if (3600*delta > 1) {// more than one arc second
        qDelta= lookup(3600*delta,factor)/3600;
    }
    else{
        qDelta= Math.log(3600*delta)/Math.log(10);
        qDelta= Math.pow(10,Math.floor(qDelta));
    }
    if (userDefinedDistance && !(minUserDistance < qDelta && qDelta < maxUserDistance)){
        var minTry= Math.abs(minUserDistance-qDelta);
        var maxTry= Math.abs(maxUserDistance-qDelta);
        qDelta= (minTry<maxTry) ? minUserDistance :
            maxUserDistance;
    }

    return qDelta;

}


/**
 * Find the labels according to the coordinates.
 * @param {Array.<Array.<number>>} levels
 * @param {CoordinateSys} csys
 * @param {String} labelFormat  pass 'hms' for sexigesimal
 * @returns {Array.<String>}
 */
function getLabels(levels,csys, labelFormat) {

    const labels = [];
    const isHms= labelFormat === 'hms' && (csys===CoordinateSys.EQ_J2000 || csys===CoordinateSys.EQ_B1950);

    for (let i=0; i < 2; i++){
        const toHms= i===0 ? CoordUtil.convertLonToString : CoordUtil.convertLatToString;
        for (let j=0; j < levels[i].length; j++) {
            const value= levels[i][j];
            labels.push(isHms ? toHms(value, csys) : numeral(value).format(precision3Digit));
        }
    }
    return labels;
}

function findLine(cc,csys, direction, value, range, screenWidth, type='image', centerWp=null, doCorrecction=false, hasPM=false){

    var intervals;
    var x, dx, y, dy;

    const dLength=direction===0?range[1][1]-range[1][0]:range[0][1]-range[0][0];

    var nInterval = 4;
    if (type==='hips'){
       const n = parseInt(dLength/angleStepForHipsMap);
       nInterval = n>4?n:4;
    }

    if (!direction )  {// longitude lines
        x  = value;
        dx = 0;
        y  = range[1][0];
        dy = (range[1][1]-range[1][0])/nInterval;
    }
    else { // latitude lines

            y = value;
            dy = 0;
            x = range[0][0];
            dx = (range[0][1] - range[0][0]);
            if (doCorrecction && hasPM){
              dx = 2* getLonRangeForPoints(cc, csys, centerWp, range);
            }
            dx = dx < 0 ? dx + 360 : dx;
            dx /= nInterval;

    }
    var opoints = findPoints(cc, csys,nInterval, x, y, dx, dy, null,doCorrecction,hasPM);
    if (type==='hips') return fixPoints(opoints);

    //NO need to do this, but left here since it was here originally
    var  straight = isStraight(opoints);
    var npoints = opoints;
    intervals = 2* nInterval;
    var nstraight;
    var count=1;
    while (intervals < screenWidth  && count<10) { //limit longer loop
        dx /= 2;
        dy /= 2;
        npoints = findPoints(cc, csys, intervals, x, y, dx, dy, opoints, hasPM);
        nstraight = isStraight(npoints);
        if (straight && nstraight) {
            break;
        }
        straight = nstraight;
        opoints = npoints;
        intervals *= 2;
        count++;
    }

    return fixPoints(npoints);


}


function isStraight(points){

    /* This function returns a boolean value depending
     * upon whether the points do not bend too rapidly.
     */
    const len = points[0].length;
    if (len < 3) return true;

    var dx0,  dy0,  len0;
    var  crossp;

    var dx1 = points[0][1]-points[0][0];
    var dy1 = points[1][1]-points[1][0];

    var len1 = (dx1*dx1) + (dy1*dy1);

    for (let i=1; i < len-1; i += 1)   {

        dx0 = dx1;
        dy0 = dy1;
        len0 = len1;
        dx1 = points[0][i+1]-points[0][i];
        dy1 = points[1][i+1]-points[1][i];

        if (dx1>=1.e20 || dy1>=1.e20) continue;
        len1 = (dx1*dx1) + (dy1*dy1);
        if (!len0  || !len1 ){
            continue;
        }
        crossp = (dx0*dx1 + dy0*dy1);
        var  cos_sq = (crossp*crossp)/(len0*len1);
        if (!cos_sq) return false;
        if (cos_sq >= 1) continue;

        var  tan_sq = (1-cos_sq)/cos_sq;
        if (tan_sq*(len0+len1) > 1) {
            return false;
        }
    }
    return true;
}

/**
 * For image map, the interval is hard coded as 4 in the original version (java).  I think that
 * since each image only covers a small stripe of the sky, the range of longitude usually less than 1 degree.
 * Howver for HiPs map, the range is whole sky (0-360).  The hard coded interval 4 is not good enough to find the
 * good points.
 *
 * For HiPs map, I use the interval = longitude-range/0.5, and interval = latitude-range/0.5.  Thus, more points
 * are checked and found for each line.
 *
 * @param cc
 * @param csys
 * @param intervals
 * @param x0
 * @param y0
 * @param dx
 * @param dy
 * @param opoints
 * @param doCorrecction
 * @param hasPM
 * @returns {[null,null]}
 */
function findPoints(cc,csys, intervals, x0, y0,dx, dy,  opoints, doCorrecction, hasPM){



    var  xpoints = [[],[]];
    var lon=[], lat=[];
    var i0, di;
    if (opoints)  {
        i0 = 1;
        di = 2;
        for (let i=0; i <= intervals; i += 2) {
            xpoints[0][i] = opoints[0][Math.trunc(i/2)];
            xpoints[1][i] = opoints[1][Math.trunc(i/2)];
        }
    }
    else {
        i0 = 0;
        di = 1;
    }

    var sharedLon, wpt, ip, xy,sharedLat,tx, ty;
    for (let i=i0; i <= intervals; i += di) {
        tx= doCorrecction && hasPM ? x0-i*dx: x0+i*dx;
        tx = tx > 360?tx-360:tx;
        tx=tx<0?tx+360:tx;
        ty=y0+i*dy;
        ty=ty>90?ty-180:ty;
        ty=ty<-90?ty+180:ty;
        sharedLon= tx;
        sharedLat= ty;
        wpt= makeWorldPt(sharedLon, sharedLat, csys);
        ip = cc.getImageWorkSpaceCoords(wpt);
        if (ip) {

            xy = makeImagePt(ip.x, ip.y);
        }
        else {

            xy=makeImagePt(1.e20,1.e20);
        }
        lon[i]= sharedLon;
        lat[i]=sharedLat;
        xpoints[0][i] = xy.x;
        xpoints[1][i] = xy.y;

    }
    return xpoints;
}

function fixPoints(points){

    // Convert points to fixed values.
    var len = points[0].length;
    for (let i=0; i < len; i += 1){
        if (points[0][i] < 1.e10) continue;
        points[0][i] = -10000;
        points[1][i] = -10000;

    }

    //points.map( (point)=> point.filter( (item) => item < 1.e20) );

    return points;
}

/**
 *
 * @param drawData
 * @param bounds
 * @param label
 * @param x
 * @param y
 * @param aitoff
 * @param screenWidth
 * @param useLabels
 * @param cc
 * @param isRaLine
 * @param plotType
 * @param centerScreenPt
 */
function drawLabeledPolyLine (drawData, bounds,  label,  x, y, aitoff,screenWidth, useLabels,cc, isRaLine,  plotType='image',  centerScreenPt=null){


    //add the  draw line data to the drawData
    var ipt0, ipt1;
    var slopAngle;
    var labelPoint;
    if(!x) return;

    for (let i=0; i<x.length-1; i+=1) {
        //check the x[i] and y[i] are inside the image screen
        if (x[i] > -1000 && x[i+1] > -1000 &&
            ((x[i] >= bounds.x) &&
            ((x[i] - bounds.x) < bounds.width) &&
            (y[i] >= bounds.y) &&
            ((y[i]-bounds.y) < bounds.height) ||
            // bounds check on x[i+1], y[i+1]
            (x[i+1] >= bounds.x) &&
            ((x[i+1] - bounds.x) < bounds.width) &&
            (y[i+1] >= bounds.y) &&
            ((y[i+1]-bounds.y) < bounds.height))) {
            ipt0= makeImageWorkSpacePt(x[i],y[i]);
            ipt1= makeImageWorkSpacePt(x[i+1], y[i+1]);
            //For image, the ra/dec interval is 8, so the points needed to be checked if they are located within the interval
            //For hips, the range for ra is 360, so no check is needed.
            if ( plotType==='hips' ||
                 plotType==='image' && (!aitoff  ||  ((Math.abs(ipt1.x-ipt0.x)<screenWidth /8 ) && (aitoff))) ) {

                 drawData.push(ShapeDataObj.makeLine(ipt0, ipt1));

                 //find the middle point of the line, index from 0, so minus 1
                 if ( (plotType==='hips' && (!isRaLine && i===0 || isRaLine && i===Math.round(x.length / 2) - 1) )  ||
                    (  plotType==='image' &&  i===Math.round(x.length / 2) - 1 )  ) {
                        var wpt1 = cc.getScreenCoords(ipt0);
                        var wpt2 = cc.getScreenCoords(ipt1);
                        const slope = (wpt2.y - wpt1.y) / (wpt2.x - wpt1.x);
                        slopAngle = Math.atan(slope) * 180 / Math.PI;
                        //since atan is multi-value function, the slopAngle is unique, for raLine, I set it is in the range of 0-180
                        if (isRaLine && slopAngle < 0) {
                            slopAngle += 180;
                        }

                        //for dec line, I set it to -90 to 90
                        if (!isRaLine && slopAngle > 90) {
                            slopAngle = 180 - slopAngle;
                        }
                        if (!isRaLine && slopAngle < -90) {
                            slopAngle = 180 + slopAngle;
                        }
                        //put the label in the center of the image line
                        if (centerScreenPt){
                            if (isRaLine) {
                                labelPoint=makeImageWorkSpacePt(x[i],centerScreenPt.y);
                            }
                            else {
                                labelPoint=makeImageWorkSpacePt(centerScreenPt.x, y[i]);
                            }
                        }
                        else {
                            labelPoint =wpt1;
                        }
                        // both screen coordinates or ImageWorkSpacePt are OK.

                }
            }
        } //if
    } // for


    // draw the label.
    if (useLabels  ){
        drawData.push(ShapeDataObj.makeText(labelPoint, label, slopAngle+'deg'));
    }
}

function drawLines(bounds, labels, xLines,yLines, aitoff,screenWidth, useLabels,cc, plotType,  centerScreenPt) {
    // Draw the lines previously computed.
    //get the locations where to put the labels
    var drawData=[];

    var  lineCount = xLines.length;


    for (let i=0; i<lineCount; i++) {
            drawLabeledPolyLine(drawData, bounds, labels[i] ,
            xLines[i], yLines[i], aitoff,screenWidth, useLabels,cc, i<lineCount/2, plotType,  centerScreenPt);
    }
    return drawData;

}

function getAdjustRange(levels, viewBorder, hasPM){
    var range=[[],[]];
    var  min = Math.min(Math.min(...levels[0]), viewBorder[0][0]);
    var  max = Math.max(Math.max(...levels[0]), viewBorder[0][1]);
    if (hasPM){
        min = Math.max(Math.min(...levels[0]), viewBorder[0][0]);
        max = Math.min(Math.max(...levels[0]), viewBorder[0][1]);
    }
    range[0]=[min, max];
    min = Math.min(Math.min(...levels[1]), viewBorder[1][0]);
    max = Math.max(Math.max(...levels[1]), viewBorder[1][1]);
    range[1]=[min, max];
    return range;
}

function computeImageGridLines(cc, csys, range,  screenWidth, numOfGridLines, labelFormat, plot) {


    const factor = plot.zoomFactor<1?1:plot.zoomFactor;

    //get levels for the whole longitude and latitude range
    var levels = getLevels(range, factor, numOfGridLines, plot.plotType);


    /* This is where we do all the work. */
    /* range and levels have a first dimension indicating x or y
     * and a second dimension for the different values (2 for range)
     * and a possibly variable number for levels.
     */
    const labels = getLabels(levels, csys, labelFormat);

    var xLines = [];
    var yLines = [];
    var offset = 0;
    var points=[];
    for (let i=0; i<2; i++) {
        for (let j=0; j<levels[i].length; j++) {
            points= findLine(cc, csys, i, levels[i][j],range, screenWidth, plot);
            xLines[offset] = points[0];
            yLines[offset] = points[1];
            offset += 1;

        }
    }
    return {xLines, yLines, labels};
}


function computeHipGridLines(cc, csys, range,  screenWidth, nGridLines, labelFormat, plot) {


    const factor = plot.zoomFactor<1?1:plot.zoomFactor;

    const {fov, centerWp}= getPointMaxSide(plot, plot.viewDim);

    var numOfGridLines=nGridLines;

    if (fov>150){
        numOfGridLines = 25;
    }
    else {
        numOfGridLines=35;
    }

    //get levels for the whole longitude and latitude range
    var levels = getLevels(range, factor, numOfGridLines, plot.plotType);

    /*get the view border, NOTE, the border does not mean the maximum and minimum of the range. If there the Prime meridian
    is in the image, the minimum of border is the value in the range from 0-minBorder, and the maximum border is the
    range [maximum border - 360]
    */
    const  {corners, viewBorder} = getViewPortInfo(plot, csys,  cc, range);

    const poles= getPoleInfo(cc, csys);

    //isPlotNorth does not work??
    /*
     * If the image is NorthnUp, the algorithm is different, however, the isPlotNorth does not work.
     *   const isNorthUp = isPlotNorth(plot);
     */


    const {zeroLonWpt, hasPM}= isPMinInView(cc, csys);

    /*
     If no corner is null, there is no pole in side, the levels can be recalculated based on the viewBorder and Prime Meridian
     information. If there is a pole, the longitude has to be from 0-360.
     If the latitude range is changed, the levels no need to be re-calculated, because it is taking care when
     the points are calculated.
     */
    var doCorrecction=false;
    if (corners && corners.indexOf(null)===-1 && poles===0  && plot.plotType==='hips'){
        levels = getLevelsHips(viewBorder, centerWp,zeroLonWpt, hasPM, numOfGridLines);

        doCorrecction=true;
    }
    /*

        const minLat = Math.min(Math.min(...levels[1]), viewBorder[1][0]);
        const maxLat = Math.max(Math.max(...levels[1]), viewBorder[1][1]);
    */

    // const r = doCorrecction ?[viewBorder[0], [minLat, maxLat]] :range;

    const vRange = doCorrecction ?getAdjustRange(levels, viewBorder, hasPM):range;

    /* This is where we do all the work. */
    /* range and levels have a first dimension indicating x or y
     * and a second dimension for the different values (2 for range)
     * and a possibly variable number for levels.
     */
    const labels = getLabels(levels, csys, labelFormat);

    var xLines = [];
    var yLines = [];
    var offset = 0;
    var points=[];
    for (let i=0; i<2; i++) {
        for (let j=0; j<levels[i].length; j++) {
            points= findLine(cc, csys, i, levels[i][j], vRange, screenWidth, plot.plotType,centerWp,doCorrecction, hasPM);
            xLines[offset] = points[0];
            yLines[offset] = points[1];
            offset += 1;

        }
    }
    return {xLines, yLines, labels};
}

function isPMinInView(cc, csys){

    var zeroLonWpt=null, lat;
    var hasPM=false;

    const det = 0.1;
    const n=180/0.1;
    for (let i=0; i<=n; i++){
        lat=-90+i*det;
        zeroLonWpt=makeWorldPt(0, lat, csys);
        if (cc.pointOnDisplay(zeroLonWpt)){
            hasPM=true;
            break;
        }
    }
    return {zeroLonWpt, hasPM};
}
function getLonRangeForPoints(cc, csys, centerWpt, viewBorder){

    const {zeroLonWpt, hasPM}= isPMinInView(cc, csys);

    var lonRange;
    if (hasPM) {//centerWpt is the PM
        if (zeroLonWpt.x.toFixed(1)===centerWpt.x.toFixed(1)){
            lonRange=viewBorder[0][0]-centerWpt.x;
        }
        else if (centerWpt.x >viewBorder[0][0]) { //PM is on the left hand of the center
            lonRange = viewBorder[0][1]-centerWpt.x;
        }
        else if (centerWpt.x < viewBorder[0][0]){ //PM is on the right of center
            lonRange = viewBorder[0][0]-centerWpt.x;

        }
    }
    else {
        lonRange =( viewBorder[0][1]-viewBorder[0][0])/2;
    }
    return lonRange;
}


function getLevelsHips(viewBorder, centerWpt,zeroLonWpt, hasPM, maxLines){

    var lonLevels=[];
    var latLevels=[];

    lonLevels.push(centerWpt.x);

    var deltaLon;
    const nHalfLines =4;
    var pt1, pt2;
    if (hasPM) {//centerWpt is the PM
       if (zeroLonWpt.x.toFixed(1)===centerWpt.x.toFixed(1)){
             deltaLon = (viewBorder[0][0]-centerWpt.x)/nHalfLines;
             for (let i=1; i<=nHalfLines; i++){
                 lonLevels.push(centerWpt.x+deltaLon*i); //on the left side of the center
                 lonLevels.push(360-deltaLon*i); //on the right side of the center

             }
       }
       else if (centerWpt.x >viewBorder[0][0]) { //PM is on the left hand of the center
           deltaLon =( viewBorder[0][1]-centerWpt.x)/nHalfLines;
           for (let i=1; i<=nHalfLines; i++){
               lonLevels.push(centerWpt.x-deltaLon*i); //on the right side of the center
               lonLevels.push(getLonInRightRange(centerWpt.x + deltaLon*i)); //on the left side of the center
           }


       }
       else if (centerWpt.x < viewBorder[0][0]){ //PM is on the right of center
           deltaLon =( viewBorder[0][0]-centerWpt.x)/nHalfLines;
           for (let i=1; i<=nHalfLines; i++){
               lonLevels.push(centerWpt.x+deltaLon*i); //on the left side of the center
               lonLevels.push(getLonInRightRange(centerWpt.x - deltaLon*i)); //on the left side of the center
           }
       }
    }
    else {
        deltaLon =( viewBorder[0][1]-viewBorder[0][0])/(2*nHalfLines);
        for (let i=1; i<=nHalfLines; i++){
            pt1=centerWpt.x + i * deltaLon;
            pt2=centerWpt.x - i * deltaLon;
            lonLevels.push(getLonInRightRange(pt1));
            lonLevels.push(getLonInRightRange(pt2));
        }
    }

    //lat lines
    const deltaLat =( viewBorder[1][1]-viewBorder[1][0])/(2*nHalfLines);
    latLevels.push(centerWpt.y);
    for (let i=1; i<=nHalfLines; i++){
            pt1=centerWpt.y + i * deltaLat;
            pt2=centerWpt.y - i * deltaLat;
            latLevels.push(getLatInRightRange(pt1));
            latLevels.push(getLatInRightRange(pt2));
    }
    return [[...new Set(lonLevels.sort( (a,b)=>{return a-b; }))],[...new Set(latLevels.sort( (a,b)=>{return a-b; }))] ];
}

/**
 *
 * @param ranges
 * @param factor
 * @param maxLines
 * @param plotType
 * @returns {Array}
 */
function getLevels(ranges,factor, maxLines, plotType){

    var levels=[];
    var  min, max, delta;
    for (let i=0; i<ranges.length; i++){
        /* Expect max and min for each dimension */
        if (ranges[i].length!==2){
            levels[i]=[];
        }
        else {
            min = ranges[i][0];
            max =ranges[i][1];
            if (min===max){
                levels[i]=[];
            }
            else if ( Math.abs(min - (-90.0))  < 0.1 && Math.abs(max - 90.0) <0.1){ //include both poles
                levels[i]= [-75.,-60., -45., -30., -15., 0., 15., 30.,  45., 60.,  75.];
            }
            else {
                /* LZ DM-10491: introduced this simple algorithm to calculate the intervals.  The previous one
                commented below caused line missing. For example, 45,0 wise, 45, 90 wise etc.

                The algorithm previous used (commented ) missed one line. I don't understand the purpose of
                 the algorithm.  The levels can be simply defined as the loop below
                 */
                levels[i] = [];

                delta =calculateDelta (min, max,factor);


                var count = Math.ceil ( (max -min)/delta);
                if (count<=2){
                    delta=delta/2.0;
                    count=2*count;
                }
                for (let j=0; j<count; j++){
                    levels[i][j] = j*delta + min;
                  if (!i && levels[i][j] > 360){
                    levels[i][j] -= 360;
                  }
                  else if (!i && levels[i][j] < 0){
                    levels[i][j] += 360;
                  }

                }


                /* We've now got the increment between levels.
                 * Now find all the levels themselves.
                 */

                //LZ comment out the original algorithm to calculate the intervals
                /* min=(max<min)?min-360:min;
                 val = min<0? min-min%delta : min + delta-min%delta;
                 count=0;
                 while (val + count*delta <= max){
                     count++;
                 }
                 if (count<=2){
                     delta=delta/2.0;
                     count=2*count;
                 }
                 levels[i] = [];
                 for (let j=0; j<count; j++){
                     levels[i][j] = j*delta + val;
                     if (!i && levels[i][j] > 360){
                         levels[i][j] -= 360;
                     }
                     else if (!i && levels[i][j] < 0){
                         levels[i][j] += 360;
                     }

                 }*/
            }
        }

    }

    return levels.map( (row)=>{
        if (row.length<maxLines){
            return Regrid(row,  maxLines, true);
        }
        else {
            return row;
        }

     });

}

