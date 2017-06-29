import React from 'react';
import PropTypes from 'prop-types';
import {get} from 'lodash';

import {InputFieldView} from './InputFieldView.jsx';
import {fieldGroupConnector} from './FieldGroupConnector.jsx';
import {fetchUrl} from '../util/WebUtil.js';
import {getRootURL} from '../util/BrowserUtil.js';
import {ServerParams} from '../data/ServerParams.js';

import LOADING from 'html/images/gxt/loading.gif';
const UL_URL = `${getRootURL()}sticky/CmdSrv?${ServerParams.COMMAND}=${ServerParams.UPLOAD}`;


function FileUploadView({fileType, isLoading, label, valid, wrapperStyle, message, onChange, value}) {
    var style = {color: 'transparent', border: 'none', background: 'none'};
    var fileName = value ?  value.split(/(\\|\/)/g).pop() : 'No file chosen';

    return (
        <div>
            <InputFieldView
                valid = {valid}
                visible = {true}
                message = {message}
                onChange = {onChange}
                type = 'file'
                label = {label}
                value = {value}
                tooltip = {value}
                labelWidth = {0}
                inline={true}
                style = {style}
                wrapperStyle = {wrapperStyle}
            />
            {fileName && <div style={{display:'inline-block', marginLeft: -150}}>{fileName}</div>}
            {isLoading && <img style={{position: 'inline-block', marginLeft: 10, width:14,height:14}} src={LOADING}/> }
        </div>
    );
}


FileUploadView.propTypes = {
    fileType   : PropTypes.string.isRequired,
    isLoading: PropTypes.bool.isRequired,
    message: PropTypes.string.isRequired,
    onChange: PropTypes.func.isRequired,
    value : PropTypes.string.isRequired
};

FileUploadView.defaultProps = {
    fileType: 'TABLE',
    isLoading: false
};

export const FileUpload = fieldGroupConnector(FileUploadView, getProps);




/*---------------------------- private ----------------------------*/

function getProps(params, fireValueChange) {
    return Object.assign({}, params,
        {
            value: params.displayValue,
            onChange: (ev) => handleChange(ev, fireValueChange, params.fileType, params.fileAnalysis)
        });
}

function handleChange(ev, fireValueChange, type, fileAnalysis) {
    var file = get(ev, 'target.files.0');
    var displayValue = get(ev, 'target.value');

    fireValueChange({
        displayValue,
        value: !fileAnalysis ? makeDoUpload(file, type) : makeDoUpload(file, type, fileAnalysis)()
    });
}

function makeDoUpload(file, type, fileAnalysis) {
    return () => {
        return doUpload(file, {type, fileAnalysis}).then(({status, message, cacheKey, analysisResult}) => {
            const valid = status === '200';
            return {isLoading: false, valid, message, value: cacheKey, analysisResult};
        }).catch((err) => {
            return {isLoading: false, valid: false, message: `Unable to upload file: ${get(file, 'name')}`};
        });
    };
}

/**
 * post the data in
 * @param {File} file
 * @param {Object} params additional parameters if any
 * @returns {Promise.<{Object}>}  The returned object is : {status:string, message:string, cacheKey:string}
 */
export function doUpload(file, params={}) {
    if (!file) return Promise.reject('Required file parameter not given');
    params = Object.assign(params, {file});   // file should be the last param due to AnyFileUpload limitation
    const options = {method: 'multipart', params};

    return fetchUrl(UL_URL, options).then( (response) => {
        return response.text().then( (text) => {
            // text is in format ${status}::${message}::${cacheKey}
            const [status, message, cacheKey, analysisResult] =  text.split('::');
            return {status, message, cacheKey, analysisResult};
        });
    });
}
