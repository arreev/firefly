import {flux} from '../Firefly.js';
import {get, isArray, isEmpty, set} from 'lodash';
import {fetchUrl} from '../util/WebUtil.js';
import {hipsSURVEYS} from './HiPSUtil.js';
import Enum from 'enum';

export const HIPSSURVEY_PREFIX = 'HIPSCntlr';
export const HIPSSURVEY_PATH = 'hipssurveys';
export const HIPSSURVEY_CREATE_PATH = `${HIPSSURVEY_PREFIX}.createPath`;
export const HIPSSURVEY_IN_LOADING = `${HIPSSURVEY_PREFIX}.inLoading`;
export const HiPSSurveyTableColumm = new Enum(['url', 'title', 'type', 'order', 'sky_fraction']);
export const HiPSSurveysURL = 'http://alasky.unistra.fr/MocServer/query?hips_service_url=*&get=record';
export const HiPSPopular = 'popular';
export const _HiPSPopular = '_popular';

export default {actionCreators, reducers };

function actionCreators() {
    return {
        [HIPSSURVEY_CREATE_PATH]: updateHiPSSurveys
    };
}

function reducers() {
    return {
        [HIPSSURVEY_PATH]: reducer
    };
}

export const HiPSDataType= new Enum([ 'image', 'cube', 'catalog'], { ignoreCase: true });

/**
 * update HiPS surveys id based on if it is for popular surveys
 * @param id
 * @param isPopular
 * @returns {*}
 */
export function updateHiPSId(id, isPopular=false) {

    const sId = id.endsWith(_HiPSPopular) ? id.substring(0, (id.length - _HiPSPopular.length)) : id;
    return isPopular ? sId+_HiPSPopular : sId;
}

export function dispatchHiPSSurveySearch(dataType, id, callback, dispatcher=flux.process) {
    dispatcher({type: HIPSSURVEY_CREATE_PATH, payload: {dataType, id, callback}});
}

/**
 * get HiPS list, used as a callback for HiPS selection
 * @param dataType
 * @param id
 * @param isPopular
 */
export function onHiPSSurveys(dataType, id, isPopular) {
    if (!getHiPSSurveys(updateHiPSId(id, isPopular))) {
        dispatchHiPSSurveySearch(dataType, updateHiPSId(id, false));
    }
}

function startLoadHiPSSurvery(dispatch, id) {
    dispatch({type: HIPSSURVEY_IN_LOADING, payload: {isLoading: true, id}});
}

function endLoadHiPSSurvey(dispatch, id, message) {
    dispatch({type: HIPSSURVEY_IN_LOADING, payload: {isLoading: false, id, message}});
}

/**
 * update HiPS surveys set including popular surveys set
 * @param action
 * @returns {Function}
 */

function updateHiPSSurveys(action) {
    return (dispatch) => {
        const {id} = get(action, 'payload');

        startLoadHiPSSurvery(dispatch, id);

        let dataType = get(action, ['payload', 'dataType'], Object.keys(HiPSDataType));


        if (!isArray(dataType)) {
            dataType = [dataType];
        }

        const dataTypeStr = dataType.map((oneType) => {
            if (typeof oneType === 'string') {
                return oneType;
            } else {
                return oneType.key;
            }
        });

        const dataProduct = HiPSDataType.enums.reduce((prev, aType) => {
            if (!dataTypeStr.includes(aType.key)) {
                prev += 'dataproduct_type=!' + aType.key;
            }
            return prev;
        }, '');

        const hipsListQuery = (dp) => {
            return dp ? `${HiPSSurveysURL}&${dp}` : HiPSSurveysURL;
        };


        fetchUrl(hipsListQuery(dataProduct), {}, true, false)
            .then((result)=>result.text())
            .then((s)=>parseHiPSList(s))
            .then((HiPSSet) => {
                set(action, ['payload', 'data'], HiPSSet);
                set(action, ['payload', 'popularData'], getPopularSurveys(HiPSSet));
                dispatch(action);
            })
            .catch( (error) => {
                console.log(error.message);
                endLoadHiPSSurvey(dispatch, id, error.message||'HiPS loading error');
            } );

    };
}

function getPopularSurveys(HiPSList) {
    const defaultUrls =  hipsSURVEYS.map((oneSurvey) => oneSurvey.url);

    return HiPSList.filter((oneSurvey) => defaultUrls.includes(oneSurvey.url));
}

function parseHiPSList(str) {
    let currentRec = {};
    const hipsColumnMap = {obs_title:        HiPSSurveyTableColumm.title.key,
                           hips_service_url: HiPSSurveyTableColumm.url.key,
                           dataproduct_type: HiPSSurveyTableColumm.type.key,
                           hips_order:       HiPSSurveyTableColumm.order.key,
                           moc_sky_fraction: HiPSSurveyTableColumm.sky_fraction.key};
    const testRow = {
        ID: 'test/row',
        [HiPSSurveyTableColumm.title.key]: 'test_record\\(%for test only%)',
        [HiPSSurveyTableColumm.url.key]: 'http:\\\\www.caltech.edu\\ipac',
        [HiPSSurveyTableColumm.type.key]: 'test.1',
        [HiPSSurveyTableColumm.order.key]: '3',
        [HiPSSurveyTableColumm.sky_fraction.key]: '0.5'
    };

    return str.split('\n')
                    .map( (s) => s.trim())
                    .filter( (s) => !s.startsWith('#') && s)
                    .map( (s) => s.split('='))
                    .reduce( (preHiPS, sAry) => {
                        if (sAry[0].trim() === 'ID') {
                            currentRec = Object.assign({}, {ID: sAry[1].trim()});
                            preHiPS.push(currentRec);
                        } else if (!isEmpty(currentRec)) {
                            const prop = sAry[0].trim();

                            if (Object.keys(hipsColumnMap).includes(prop)) {
                                Object.assign(currentRec, {[hipsColumnMap[prop]]: sAry[1].trim()});
                            }
                        }
                        return preHiPS;
                    }, [])
                    .reduce( (preHiPS, row) => {
                        if (row[HiPSSurveyTableColumm.type.key]) {    // remove the item with undefined type
                            preHiPS.push(row);
                            if (!row[HiPSSurveyTableColumm.title.key]) {
                                row[HiPSSurveyTableColumm.title.key] = '';   // reset title to empty string if undefined
                            }
                        }
                        return preHiPS;
                    }, []);

}

/**
 * get HiPS survey list from the store
 * @param id
 * @returns {*}
 */
export function getHiPSSurveys(id) {
    return get(flux.getState(), [HIPSSURVEY_PATH, id, 'data']);
}



/**
 * check if HiPS survey list is under loading
 * @param id
 * @returns {*}
 */
export function isLoadingHiPSSurverys(id) {
    return get(flux.getState(), [HIPSSURVEY_PATH, id, 'isLoading'], false);
}

/**
 * get loading message if there is
 * @param id
 * @returns {*}
 */
export function getHiPSLoadingMessage(id) {
    return get(flux.getState(), [HIPSSURVEY_PATH, id, 'message']);
}

function reducer(state={}, action={}) {
    const {type} = action;
    if (!type || !action.payload) return state;

    const {data, popularData, isLoading, id, message} = action.payload;
    let   retState = state;

    switch(type) {
        case HIPSSURVEY_CREATE_PATH:
            if (data) {
                retState = addSurveys(data, popularData, id,  state);
            }
            set(retState, [id, 'isLoading'], false);
            break;
        case HIPSSURVEY_IN_LOADING:
            retState = setLoading(isLoading, id, state);
            if (message) {
                set(retState, [id, 'message'], message);
            }
            break;
    }
    return retState;
}

function addSurveys(data, popularData, id, state ) {
    const newState = Object.assign({}, state);

    set(newState, [id, 'data'], data);
    set(newState, [id+'_popular', 'data'], popularData);
    set(newState, [id, 'valid'], true);
    return newState;
}

function setLoading(isLoading, id, state) {
    const newState = Object.assign({}, state);

    set(newState, [id, 'isLoading'], isLoading);
    set(newState, [id, 'message'], '');
    return newState;
}


