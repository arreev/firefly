/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

/**
 * Date: Mar 13, 2009
 * @author loi
 */
/* eslint prefer-template:0 */

import Enum from 'enum';
import { getRootURL } from './BrowserUtil.js';

const ParamType= new Enum(['POUND', 'QUESTION_MARK']);

const saveAsIpacUrl = getRootURL() + 'servlet/SaveAsIpacTable';


/**
 * Returns a string where all characters that are not valid for a complete URL have been escaped.
 * Also, it will do URL rewriting for session tracking if necessary.
 * Fires SESSION_MISMATCH if the seesion ID on the client is different from the one on the server.
 *
 * @param url    this could be a full or partial url.  Delimiter characters will be preserved.
 * @param paramType  if the the parameters are for the server use QUESTION_MARK if the client use POUND
 * @param {array|Object} params parameters to be appended to the url.  These parameters may contain
 *               delimiter characters.  Unlike url, delimiter characters will be encoded as well.
 * @return {string} encoded url
 */
export const encodeUrl= function(url, paramType, params) {
    var paramChar= paramType===ParamType.QUESTION_MARK ? '?': '#';
    var parts = url.split('\\'+paramChar, 2);
    var baseUrl = parts[0];
    var queryStr = encodeURI(parts.length===2 ? parts[1] : '');

    var paramAry= params || [];

    if (params.length===1 && Object.keys(params[0]).length>0) {
        paramAry= Object.keys(params[0]).reduce( (ary,key) => {
            ary.push({name:key, value : params[0][key]});
            return ary;
        },[]);
    }

    queryStr= paramAry.reduce((str,param,idx) => {
        if (param && param.name) {
            var key = encodeURI(param.name.trim());
            var val = param.value ? encodeURIComponent(param.value.trim()) : '';
            str += val.length ? key + '=' + val + (idx < paramAry.length ? '&' : '') : key;
            return str;
        }
    },'');

    return encodeURI(baseUrl) + (queryStr.length ? paramChar + queryStr : '');
};


/**
 * Returns a string where all characters that are not valid for a complete URL have been escaped.
 * Also, it will do URL rewriting for session tracking if necessary.
 * Fires SESSION_MISMATCH if the seesion ID on the client is different from the one on the server.
 *
 * @param url    this could be a full or partial url.  Delimiter characters will be preserved.
 * @param params parameters to be appended to the url.  These parameters may contain
 *               delimiter characters.  Unlike url, delimiter characters will be encoded as well.
 * @return encoded url
 */
export const encodeServerUrl= function(url, params) {
    return encodeUrl(url, ParamType.QUESTION_MARK,params);
};




/**
 *
 * @param {ServerRequest} request
 * @return {string} encoded
 */
export const getTableSourceUrl= function(request) {
    request.setStartIndex(0);
    request.setPageSize(Number.MAX_SAFE_INTEGER);
    var source = { name : 'request', value : request.toString()};  //todo : i don't think I got this line right
    var filename = request.getParam('file_name');
    if (!filename) filename = request.getRequestId();
    var fn = { name: 'file_name', value : filename};
    return encodeServerUrl(saveAsIpacUrl, source, fn);
};



/**
 * A wrapper for the underlying window.fetch function.
 * see https://github.com/github/fetch for usage.
 * see https://developer.mozilla.org/en-US/docs/Web/API/Fetch_API for current status on the API
 * see https://fetch.spec.whatwg.org/ for official standard

 * This function applies default behaviors before fetching.
 * options.params is a custom property used to carry a set of parameters.  It does not need to
 *                be encoded.  Base on the method used, it will be handled internally.
 *
 * @param url
 * @param options
 * @return a promise of the response when successful, or reject with an Error.
 */
export function fetchUrl(url, options) {

    if (!url) return;
    url = url.trim();
    options = options || {};

    if ( !(url.startsWith('http') || url.startsWith('/')) ) {
        url = getRootURL() + url;
    }

    // define defaults request options
    const req = { method: 'GET',
            mode: 'cors',
            credentials: 'include',
            cache: 'default'
        };
    options = Object.assign(req, options);

    const headers = {};
    if (options.method.toUpperCase() === 'POST') {
        // add default content-type header when method is 'post'
        headers['Content-type'] = 'application/x-www-form-urlencoded; charset=UTF-8';
    }
    options.headers = Object.assign(headers, options.headers);

    if (options.params) {
        if (options.method.toUpperCase() === 'GET') {
            url = encodeUrl(url, ParamType.QUESTION_MARK, [options.params]);
        } else {
            if (!options.body) {
                // if 'post' but, body is not provided, add the parameters into the body.
                options.body = new FormData();
                Object.keys(options.params).forEach( (key) => {
                    options.body.append(key, options.params[key]);
                });
            }
        }
    }

    // do the actually fetch, then return a promise.
    return fetch(url, options)
        .then( (response) => {
            if (response.ok) {
                return Promise.resolve(response);
            } else {
                return Promise.reject(new Error(response.statusText));
            }
        }).catch( (error) => {
            return Promise.reject(new Error(`Request failed: ${url}`, error));
        });
}

