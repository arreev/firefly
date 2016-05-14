/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

/**
 * for Julian precession.
 * fk5prc to have simpler arg list when no proper motions are present.
 * calls fk5prp.  See fk5prp for comments on method.
 *
 * @param {number} eqx1
 * @param {number} rain
 * @param {number} decin
 * @param {number} eqx2
 * @return {number} raou
 * @return {number} decou
 * @return {{raou:number,decou:number}}
 */
export function fk5prc(eqx1, rain, decin, eqx2) {

    const {raou,decou}= fk5prp(eqx1,rain,decin,0.0,0.0,0.0,0.0, eqx2);
    return {raou,decou};
}



var rtod, dtor, delt,f,leqx1=-1.0,leqx2=-1.0;
var p= [ new Array(3), new Array(3), new Array(3) ];


/**
 *  j.bennett (1992) - for Julian precession:
 *              using IAU (1976) System of Astronomical Constants
 *
 *  allows precession of proper motion values (also will make use of
 *    parallax and radial velocity if not zero.
 *
 * All arguments are double precision:
 *  @param {number} eqx1 is equinox of input position (in year - e.g. 2000.0d0)
 *  @param {number} rain, decin is input position in decimal degrees at eqx1.
 *  @param {number} decin
 *  @param {number} pmain is proper motion in ra in seconds of time per Julian century.
 *  @param {number} pmdin is proper motion in dec in seconds of arc per Julian century.
 *  @param {number} pin is parallax in arc seconds (0.0d0 if none)
 *  @param {number} vin is radial velocity in km/sec (0.0d0 if none)
 *  @param {number} eqx2 is eqinox of output position (in year - e.g. 1992.0d0)
 *                 (or 1949.99979d0 if quick and dirty approx of B1950 desired)
 *  @return {number} raou, decou is position in decimal degrees at eqx2.
 *  @return {number} decou
 *  @return {number} pmaou is proper motion in ra in seconds of time per Julian century for eqx2.
 *  @return {number} pmdou is proper motion in dec in seconds of arc per Julian century for eqx2.
 *  @return {{raou:number,decou:number,pmaou:number,pmdou:number}}
 *
 *******************************************************************************
 * method from Fricke,W.,Schwan,H.,Lederle,T.,1988,Fifth Fundamental Catalogue
 *     (FK5) Part I. (pages 10-11)
 *******************************************************************************
 */
export function fk5prp( eqx1,  rain, decin, pmain, pmdin, pin, vin,
                        eqx2) {


    var raou, decou, pmaou, pmdou;
    var i;
    var zetar,zr,thetar,zeta,z,theta;
    var tau, t;
    var r0= [],rdot0= [];
    var r= [], rdot= [];
    var pivelf;
    var czet,cthet,cz,szet,sthet,sz;
    var cosa,sina,cosd,sind,rar,decr;
    var cosao,sinao,cosdo,sindo;
    var pmas,pmds,rdiv,raour,decour,pmaous,pmdous;
    var duda= [],dudd= [];

    if(eqx1===eqx2) {
        raou  = rain;
        decou = decin;
        pmaou = pmain;
        pmdou = pmdin;
        return {raou,decou,pmaou,pmdou};
    }


    if(eqx1 !==leqx1 || eqx2 !==leqx2) {
        dtor = Math.atan(1.0) / 45.0;
        rtod = 1.0 /dtor;
        /* tau, t, delt in Julian centuries: */
        tau =  (eqx1 - 2000.0) * 0.01;
        t    = (eqx2 - eqx1) * 0.01;
        delt = t;
        f    = 4.8481368e-06;
        leqx1 = eqx1;
        leqx2 = eqx2;

        /*
         *  zeta,theta,z from FK5 Catalogue
         *  zeta, theta, z in seconds of arc:
         */
        zeta = (2306.2181 + 1.39656*tau - 0.000139*tau*tau)*t +
            (0.30188  - 0.000344*tau)*t*t + 0.017998*t*t*t;
        z = (2306.2181 + 1.39656*tau - 0.000139*tau*tau)*t +
            (1.09468 + 0.000066*tau)*t*t + 0.018203*t*t*t;
        theta = (2004.3109 - 0.85330*tau - 0.000217*tau*tau)*t -
            (0.42665 + 0.000217*tau)*t*t - 0.041833*t*t*t;
        zetar = (zeta/3600.0)*dtor;
        zr = (z/3600.0)*dtor;
        thetar = (theta/3600.0)*dtor;
        czet = Math.cos(zetar);
        szet = Math.sin(zetar);
        cz   = Math.cos(zr);
        sz   = Math.sin(zr);
        cthet= Math.cos(thetar);
        sthet= Math.sin(thetar);

        /*  p matrix from Green,Robin M.,1985,Spherical Astronomy, Cambridge University
         *   Press, p. 221-222.  See also: Lieske,J.,1979,Astron.Astrophys.,Vol.73,
         *   p. 282-284.
         *
         */
        p[0][0] = czet*cthet*cz - szet*sz;
        p[1][0] = czet*cthet*sz + szet*cz;
        p[2][0] = czet*sthet;
        p[0][1] = -szet*cthet*cz - czet*sz;
        p[1][1] = -szet*cthet*sz + czet*cz;
        p[2][1] = -szet*sthet;
        p[0][2] = -sthet*cz;
        p[1][2] = -sthet*sz;
        p[2][2] = cthet;
    }
    rar  = dtor*rain;
    decr = dtor*decin;
    cosa = Math.cos(rar);
    sina = Math.sin(rar);
    cosd = Math.cos(decr);
    sind = Math.sin(decr);
    r0[0] = cosd*cosa;
    r0[1] = cosd*sina;
    r0[2] = sind;

    pmas = pmain * 15.0;
    pmds = pmdin;
    if(vin=== 0 || pin ===0)
    {
        rdot0[0] = f*(pmas*(-cosd)*sina + pmds*(-sind)*cosa);
        rdot0[1] = f*(pmas*cosd*cosa + pmds*(-sind)*sina);
        rdot0[2] = f*(pmds*cosd);
    }
    else
    {
        pivelf = 21.094953 * pin *vin;
        rdot0[0] = f*(pmas*(-cosd)*sina + pmds*(-sind)*cosa+pivelf*r0[0]);
        rdot0[1] = f*(pmas*cosd*cosa + pmds*(-sind)*sina + pivelf*r0[1]);
        rdot0[2] = f*(pmds*cosd + pivelf*r0[2]);
    }

    for(i=0; i<3; i++) {
        rdot[i] = p[i][0]*rdot0[0] + p[i][1]*rdot0[1] + p[i][2]*rdot0[2];
        r[i] = p[i][0]*(r0[0]+rdot0[0]*delt) +
            p[i][1]*(r0[1]+rdot0[1]*delt) +
            p[i][2]*(r0[2]+rdot0[2]*delt);
    }

    raour  =Math.atan2(r[1],r[0]);
    decour =Math.atan2(r[2],Math.sqrt(r[0]*r[0] + r[1]*r[1]));
    rdiv  = Math.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    cosdo = Math.cos(decour);
    sindo = Math.sin(decour);
    cosao = Math.cos(raour);
    sinao = Math.sin(raour);

    duda[0] = -cosdo*sinao;
    duda[1] =  cosdo*cosao;
    duda[2] =  0.0;
    dudd[0] = -sindo*cosao;
    dudd[1] = -sindo*sinao;
    dudd[2] =  cosdo;
    rdot[0] = rdot[0] / rdiv;
    rdot[1] = rdot[1] / rdiv;
    rdot[2] = rdot[2] / rdiv;

    pmaous = (rdot[0]*duda[0] + rdot[1]*duda[1] + rdot[2]*duda[2])/
        (f*cosdo*cosdo);
    pmdous = (rdot[0]*dudd[0] + rdot[1]*dudd[1] + rdot[2]*dudd[2])/f;

    raou = raour * rtod;
    if(raou > 360.0) {
        raou = raou - 360.0;
    }
    else if (raou < 0.0) {
        raou = raou + 360.0;
    }

    decou = decour * rtod;
/*
* since decour is output of atan2(y,x) with x always positive, if dec
*  .gt. 90. it is due to rounding when being converted to degrees.
*/
    if(decou > 90.0) {
        decou = 90.0;
    }
    else if (decou < -90.0) {
        decou = -90.0;
    }

    pmaou = pmaous / 15.0;
    pmdou = pmdous;

    return {raou,decou,pmaou,pmdou};
}
