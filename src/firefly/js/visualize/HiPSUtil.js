/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */


/**
 * Many utility functions for working with HiPS.  Several of the more computationally intensive ones use cached results
 *
 * Here is some documentation on HiPS.
 * https://aladin.u-strasbg.fr/hips/hipsdoc.pdf
 * http:www.ivoa.net/documents/HiPS/20170519/REC-HIPS-1.0-20170519.pdf
 */


import {get, isUndefined} from 'lodash';
import {clone} from '../util/WebUtil.js';
import {CysConverter} from './CsysConverter.js';
import {makeWorldPt, makeDevicePt} from './Point.js';
import {SpatialVector, HealpixIndex} from '../externalSource/aladinProj/HealpixIndex.js';
import {convert, computeDistance} from './VisUtil.js';
import {replaceHeader, getScreenPixScaleArcSec} from './WebPlot.js';
import {primePlot} from './PlotViewUtil.js';
import CoordinateSys from './CoordSys';
import {encodeServerUrl} from '../util/WebUtil.js';
import {getRootURL} from '../util/BrowserUtil.js';
import {toDegrees} from './VisUtil.js';


export const MAX_SUPPORTED_HIPS_LEVEL= 14;

/**
 *
 * @param {WebPlot} plot
 * @param {WorldPt} wp new center of projection
 */
export function changeProjectionCenter(plot, wp) {
    if (!plot) return undefined;
    wp= convert(wp, plot.projection.coordSys);
    const header= clone(plot.projection.header, {crval1:wp.x, crval2:wp.y});
    return replaceHeader(plot,header);
}


/**
 * Determine how deep we show a grid on the current Hips map.  This will change depending on how deeply it is zoomed.
 * @param plot
 * @return {number}
 */
export function getMaxDisplayableHiPSGridLevel(plot) {
    let {norder}= getHiPSNorderlevel(plot);
    norder = norder>3 ? norder+3 : norder+2;
    if (norder>14) norder= 14;
    return norder;
}



const  angSizeCacheMap = new WeakMap(); // use week map since we don't want to keep old plot objects

/**
 * Return the angular size of the pixel at the nOrder level for this plot.
 * Results are cached so this function is very efficient.
 * @param plot
 * @return {*}
 */
export function getPlotTilePixelAngSize(plot) {
    if (!plot) return 0;
    let size= angSizeCacheMap.get(plot);
    if (isUndefined(size)) {
        size= getTilePixelAngSize(getHiPSNorderlevel(plot,true).norder);
        angSizeCacheMap.set(plot, size);
    }
    return size;
}



const tilePixelAngSizeCacheMap= {};

/**
 * Return the angular size of the pixel of a nOrder level. this function assumes 512x512 tiles sizes
 * Results are cached so this function is very efficient.
 * @param {Number} nOrder
 * @return {Number} the angular size of a pixel in a HiPS tile
 */
export function getTilePixelAngSize(nOrder) {
    nOrder= Math.trunc(nOrder);
    if (tilePixelAngSizeCacheMap[nOrder]) return tilePixelAngSizeCacheMap[nOrder];
    const rad= Math.sqrt(4*Math.PI / (12*Math.pow(512*Math.pow(2,nOrder) , 2)));
    tilePixelAngSizeCacheMap[nOrder]= toDegrees(rad);
    return tilePixelAngSizeCacheMap[nOrder];
}

/**
 *
 * @param {WebPlot} plot
 * @param {boolean} [limitToImageDepth] When true, do not return a number that is greater than what this HiPS map
 * can display.  Use hipsProperties.hips_order to determine.
 * @return {{useAllSky:boolean, norder:number}} norder is the result, useAllSky true when the norder is 2 or 3 but
 * is zoom out so much that the full norder 3 tiles are not necessary, if all sky map is available the is should be used
 * for drawing.
 */
export function getHiPSNorderlevel(plot, limitToImageDepth= false) {
    if (!plot) return {norder:-1, useAllSky:false};

    const screenPixArcsecSize= getScreenPixScaleArcSec(plot);
    if (screenPixArcsecSize> 130) return  {useAllSky:true, norder:2};
    if (screenPixArcsecSize> 100) return  {useAllSky:true, norder:3};

    let norder= getNOrderForPixArcSecSize(screenPixArcsecSize);

    if (limitToImageDepth) {
        const maxOrder= Number(get(plot, 'hipsProperties.hips_order', '3'));
        norder= Math.min(norder, maxOrder);
    }
    return {norder, useAllSky:false};

}

const nOrderForPixAsSizeCacheMap= {};

/**
 * Return the best norder for a given screen pixel angular size in arc seconds
 * Results are cached so this function is very efficient.
 * @param {Number} sizeInArcSec - pixel size in arc seconds
 * @return {Number} the best norder for the pixel
 */
function getNOrderForPixArcSecSize(sizeInArcSec) {
    const sizeInArcSecKey= Math.trunc(sizeInArcSec*10000)+'';
    let norder= nOrderForPixAsSizeCacheMap[sizeInArcSecKey];
    if (isUndefined(norder)) {
        let nside = HealpixIndex.calculateNSide(sizeInArcSec*512); // 512 size tiles hardcoded, should fix
        if (nside===8192 && sizeInArcSec<.04) nside=16384 ;

        norder = Math.log(nside)/Math.log(2); // convert to a base 2 log - 	logb(x) = logc(x) / logc(b)
        norder= Math.max(3, norder);
        nOrderForPixAsSizeCacheMap[sizeInArcSecKey]= norder;
    }
    return norder;
}


export function makeHiPSTileUrl(plot, nOrder, tileNumber) {
    if (!plot) return null;
    const dir= Math.floor(tileNumber/10000)*10000;
    const exts= get(plot, 'hipsProperties.hips_tile_format', 'jpg');
    const cubeExt= plot.cubeDepth>1 && plot.cubeIdx>0 ? '_'+plot.cubeIdx : '';
    return makeHipsUrl(`${plot.hipsUrlRoot}/Norder${nOrder}/Dir${dir}/Npix${tileNumber}${cubeExt}.${getHiPSTileExt(exts)}`, plot.proxyHips);
}

/**
 * @param urlRoot
 * @param exts
 * @param cubeIdx
 * @param proxy
 * @return {*}
 */
export function makeHiPSAllSkyUrl(urlRoot,exts,cubeIdx= 0, proxy= false) {
    if (!urlRoot || !exts) return null;
    const cubeExt= cubeIdx? '_'+cubeIdx : '';
    return makeHipsUrl(`${urlRoot}/Norder3/Allsky${cubeExt}.${getHiPSTileExt(exts)}`, proxy);
}

export function makeHiPSAllSkyUrlFromPlot(plot) {
    if (!plot) return null;
    const exts= get(plot, 'hipsProperties.hips_tile_format', 'jpg');
    const cubeIdx= plot.cubeDepth>1 && plot.cubeIdx>0 ? plot.cubeIdx : 0;
    return makeHiPSAllSkyUrl(plot.hipsUrlRoot, exts, cubeIdx, plot.proxyHips);

}


/**
 * Build the HiPS url
 * @param {String} url
 * @param {boolean} proxy - make a URL the uses proxying
 * @return {String} the modified url
 */
export function makeHipsUrl(url, proxy) {
    if (proxy) {
        const params= {
            hipsUrl : url
        };
        return encodeServerUrl(getRootURL() + 'servlet/Download', params);
    }
    else {
        return url;
    }

}


/**
 * Choose an extension to use from the available extensions.
 * @param {Array.<String>} exts
 * @return {string}
 */
export function getHiPSTileExt(exts) {
    if (!exts) return null;
    if (exts.includes('png')) return'png';
    else if (exts.includes('jpeg') || exts.includes('jpg')) return 'jpg';
    else return 'jpg';

}

/**
 *
 * @param {WebPlot} plot
 * @param {Dimension} viewDim
 * @return {{centerWp:WorldPt, fov: number, centerDevPt: DevicePt}}
 */
export function getPointMaxSide(plot, viewDim) {

    const cc= CysConverter.make(plot);
    const {width,height}= viewDim;
    const centerDevPt= makeDevicePt(width/2, height/2);
    const centerWp= cc.getWorldCoords( centerDevPt, plot.imageCoordSys);

    const pt0= cc.getWorldCoords( makeDevicePt(0,0), plot.imageCoordSys);
    const ptWidth= cc.getWorldCoords( makeDevicePt(width,0), plot.imageCoordSys);
    const ptHeight= cc.getWorldCoords( makeDevicePt(0,height), plot.imageCoordSys);


    if (!ptWidth || !ptHeight) {
        return {centerWp,centerDevPt, fov:180};
    }

    const widthInDeg= computeDistance(pt0, ptWidth);
    const heightInDeg= computeDistance(pt0, ptHeight);
    return {centerWp, centerDevPt, fov:Math.max(widthInDeg,heightInDeg)};
}


/**
 *
 * @param {PlotView} pv
 * @param {number} size in degrees
 * @return {number} a zoom level
 */
export function getHiPSZoomLevelToFit(pv,size) {
    const {width,height}=pv.viewDim;
    const plot= primePlot(pv);
    if (!plot || !width || !height) return 1;

    // make version of plot centered at 0,0
    const tmpPlot= changeProjectionCenter(plot, makeWorldPt(0,0, plot.imageCoordSys));
    const cc= CysConverter.make(tmpPlot);
    const pt1= cc.getImageCoords( makeWorldPt(0,0, plot.imageCoordSys));
    const pt2= cc.getImageCoords( makeWorldPt(size,0, plot.imageCoordSys));
    return Math.min(width, height)/Math.abs(pt2.x-pt1.x);
}

/**
 *
 * @param {PlotView} pv
 * @return {number} fov in degrees
 */
export function getHiPSFoV(pv) {
    const cc= CysConverter.make(primePlot(pv));
    const {width,height}=pv.viewDim;
    if (!cc || !width || !height) return;

    const pt1= cc.getWorldCoords( makeDevicePt(1,height/2));
    const pt2= cc.getWorldCoords( makeDevicePt(width-1,height/2));
    return (pt1 && pt2) ? computeDistance(pt1,pt2) : 180;
}


function makeCorners(hpIdx, pixList, coordSys) {
    const spVec = new SpatialVector();
    return pixList.map( (ipix) => {
        const corners = hpIdx.corners_nest(ipix, 1);
        const wpCorners= corners.map( (c) => {
            spVec.setXYZ(c.x, c.y, c.z);
            return makeWorldPt(spVec.ra(), spVec.dec(), coordSys);
        });
        return { ipix, wpCorners };
    });
}




/**
 * This function make an object (with functions) to cache allsky type computations.
 * @return {*}
 */
function makeSimpleHpxCornerCache() {
    const tmpHealpixIdx3 = new HealpixIndex(8);
    const tmpHealpixIdx2 = new HealpixIndex(4);
    const npix = HealpixIndex.nside2Npix(8);
    const cachedCorners8= [];
    for (let ipix=0; ipix<npix; ipix++) {
        cachedCorners8.push(tmpHealpixIdx3.corners_nest(ipix, 1));
    }


    const level3pixelCnt = HealpixIndex.nside2Npix(8);
    const level2pixelCnt = HealpixIndex.nside2Npix(4);
    const cachedLevel3FullPixelList= [];
    const cachedLevel2FullPixelList= [];
    for (let ipix=0; ipix<level3pixelCnt; ipix++) cachedLevel3FullPixelList[ipix]= ipix;
    const j2000Leve3Corners= makeCorners(tmpHealpixIdx3, cachedLevel3FullPixelList, CoordinateSys.EQ_J2000);
    const galLevel3Corners= makeCorners(tmpHealpixIdx3, cachedLevel3FullPixelList, CoordinateSys.GALACTIC);


    for (let ipix=0; ipix<level2pixelCnt; ipix++) cachedLevel2FullPixelList[ipix]= ipix;
    const j2000Leve2Corners= makeCorners(tmpHealpixIdx2, cachedLevel2FullPixelList, CoordinateSys.EQ_J2000);
    const galLevel2Corners= makeCorners(tmpHealpixIdx2, cachedLevel2FullPixelList, CoordinateSys.GALACTIC);


    return {
        cornersNest(ipix,nside, healpixIdx) {
            return nside === 8 ? cachedCorners8[ipix] : healpixIdx.corners_nest(ipix, 1);
        },
        getFullLevel3CornerList(coordSys)  {
            switch (coordSys) {
                case CoordinateSys.EQ_J2000: return j2000Leve3Corners;
                case CoordinateSys.GALACTIC: return galLevel3Corners;
                default: return null;
            }
        },
        getFullLevel2CornerList(coordSys)  {
            switch (coordSys) {
                case CoordinateSys.EQ_J2000: return j2000Leve2Corners;
                case CoordinateSys.GALACTIC: return galLevel2Corners;
                default: return null;
            }
        }

    };

}

/**
 * @Function
 * lazily defined.
 * @param {number} ipix
 * @param {number} nside
 * @param {HealpixIndex} healpixIdx
 */
let healpixCache;



/**
 *
 * @param norder
 * @param {WorldPt} centerWp - center of visible area, coorindate system of this point should be same as the projection
 * @param {number} fov - Math.max(width, height) of the field of view in degrees (i think)
 * @param {CoordinateSys} dataCoordSys
 * @return {Array.<{ipix:string, wpCorners:Array.<WorldPt>}>} an array of objects the contain the healpix
 *            pixel number and a worldPt array of corners
 */
export function getVisibleHiPSCells (norder, centerWp, fov, dataCoordSys) {
    if (!healpixCache) healpixCache= makeSimpleHpxCornerCache();
    const nside = Math.pow(2, norder);
    const dataCenterWp= convert(centerWp, dataCoordSys);
    let hpxIdx;

         // ------------------------------
         // first, find the Healpix pixels for the field of view, if the fov is large just return all of them
         // ------------------------------
    let pixList;
    if (fov>80 && norder===3) { // this case if so common, don't recompute, use cache
        return filterAllSky(dataCenterWp, healpixCache.getFullLevel3CornerList(dataCoordSys));
    }
    else if (fov>80 && norder===2) { // this case if so common, don't recompute, use cache
        return filterAllSky(dataCenterWp, healpixCache.getFullLevel2CornerList(dataCoordSys));
    }
    else if (fov>80) { // with norder 1 or 2
        hpxIdx = new HealpixIndex(nside);
        pixList= [];
        const npix = HealpixIndex.nside2Npix(nside);
        for (let ipix=0; ipix<npix; ipix++) pixList[ipix]= ipix;
    }
    else {
        hpxIdx = new HealpixIndex(nside);
        const spatialVector = new SpatialVector();
        spatialVector.set(dataCenterWp.x, dataCenterWp.y);
        let radius = fov/2;
                          // we need to extend the radius (suggestion from Aladin)
        if (fov>60) radius *= 1.6;
        else if (fov>12) radius *=1.45;
        else radius *= 1.1;

        pixList = hpxIdx.queryDisc(spatialVector, radius*Math.PI/180.0, true, true);
    }

         // ------------------------------
         // second, find the 4 corners for every Healpix pixel
         // ------------------------------
    const spVec = new SpatialVector();
    let cells= pixList.map( (ipix) => {
        const corners = healpixCache.cornersNest(ipix, nside, hpxIdx);
        const wpCorners= corners.map( (c) => {
            spVec.setXYZ(c.x, c.y, c.z);
            return makeWorldPt(spVec.ra(), spVec.dec(), dataCoordSys);
        });
        return { ipix, wpCorners };
    });

    if (fov>80) {
        cells= filterAllSky(dataCenterWp, cells);
    }

    return cells;
}



function filterAllSky(centerWp, cells) {
    return cells.filter( (cell) =>{
        const {wpCorners}= cell;
        return (computeDistance(centerWp, wpCorners[0]) <90);
    });
}

export const API_HIPS_CONSTANTS= {
    TWO_MASS: 'http://alasky.u-strasbg.fr/2MASS/Color',
    DSS_COLORED: 'http://alasky.u-strasbg.fr/DSS/DSSColor',
    DSS2_RED: 'http://alasky.u-strasbg.fr/DSS/DSS2Merged',
    AllWISE_COLOR: 'http://alasky.u-strasbg.fr/AllWISE/RGB-W4-W2-W1/',
    IRAC_COLOR: 'http://alasky.u-strasbg.fr/SpitzerI1I2I4color',
};

/**
 * return the value of the constant or if not found return the given value.
 * @param c
 * @return {*}
 */
export function resolveHiPSConstant(c) {
    return API_HIPS_CONSTANTS[c] || c;

}

