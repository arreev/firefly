/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React, {PureComponent} from 'react';
import PropTypes from 'prop-types';
import {get,omit} from 'lodash';
import {dispatchMountComponent,dispatchValueChange} from '../fieldGroup/FieldGroupCntlr.js';
import FieldGroupUtils, {getFieldGroupState} from '../fieldGroup/FieldGroupUtils.js';
import {flux} from '../Firefly.js';

const defaultConfirmInitValue= (v) => v;

export class FieldGroupEnable extends PureComponent {

    constructor(props, context) {
        super(props, context);
        this.fireValueChange = this.fireValueChange.bind(this);
        const {fieldKey} = props;
        const groupKey= getGroupKey(props,context);

        let fieldState;
        if (props.forceReinit) {
            fieldState = props.initialState || FieldGroupUtils.getGroupFields(groupKey)[fieldKey];
        }
        else {
            fieldState = FieldGroupUtils.getGroupFields(groupKey)[fieldKey] || props.initialState || {};
        }
        this.state = {fieldState};
    }

    fireValueChange(payload) {
        const {fieldKey}= this.props;
        const modPayload= Object.assign({},payload, {fieldKey,groupKey:getGroupKey(this.props,this.context)});
        dispatchValueChange(modPayload);
    }

    componentWillReceiveProps(nextProps, context) {
        const {fieldKey}= nextProps;
        const groupKey= getGroupKey(this.props,this.context);
        const nGroupKey= getGroupKey(nextProps, context);
        if (this.props.fieldKey!==fieldKey || nGroupKey !== groupKey) {
            const fieldState= get(FieldGroupUtils.getGroupFields(nGroupKey),fieldKey);
            this.storeUnmount(this.props.fieldKey,groupKey);
            this.reinit(nextProps,fieldState, context);
        }
        else {
            this.updateFieldState(nextProps, context);
        }
    }

    updateFieldState(props, context) {
        const {fieldKey}= props;
        const groupState = getFieldGroupState(getGroupKey(props,context));
        if (get(groupState,'mounted') && get(groupState,['fields',fieldKey])) {
            if (this.iAmMounted && this.state.fieldState !== groupState.fields[fieldKey]) {
                this.setState({fieldState: groupState.fields[fieldKey]});
            }
        }
    }


    storeUnmount(fieldKey,groupKey) {
        if (this.storeListenerRemove) this.storeListenerRemove();
        this.storeListenerRemove= null;
        dispatchMountComponent(
            groupKey, fieldKey, false, get(this.state, 'fieldState.value')
        );
    }


    reinit(props,fieldState, context) {
        let value= get(fieldState, 'value');
        const {confirmInitialValue= defaultConfirmInitValue} = props;
        if (props.forceReinit) {
            value= props.initialState.value;
        }
        if (this.storeListenerRemove) this.storeListenerRemove();
        this.storeListenerRemove = flux.addListener(()=> this.updateFieldState(props, context));
        dispatchMountComponent(
            getGroupKey(props, context),
            props.fieldKey,
            true,
            confirmInitialValue(value,props,fieldState),
            props.initialState
        );
    }

    componentDidMount() {
        this.reinit(this.props,this.state.fieldState, this.context);
        this.iAmMounted= true;
    }

    componentWillUnmount() {
        this.storeUnmount(this.props.fieldKey, getGroupKey(this.props,this.context));
        this.iAmMounted= false;
    }

    render() {
        const {children}= this.props;
        const propFromStore= getPropsFromStore(this.state.fieldState,this.props);
        return children(propFromStore, this.fireValueChange);
    }
}

FieldGroupEnable.contextTypes = { groupKey: PropTypes.string };

FieldGroupEnable.propTypes = {
    fieldKey: PropTypes.string.isRequired,
    groupKey: PropTypes.string, // usually comes from context but this is a fallback
    initialState: PropTypes.object,
    forceReinit: PropTypes.bool,
    confirmInitialValue: PropTypes.func
};


const STORE_OMIT_LIST= ['fieldKey', 'groupKey', 'initialState', 'fireReducer',
                             'confirmInitialValue', 'forceReinit', 'children', 'mounted'];
const defValidatorFunc= () => ({valid:true,message:''});

function getPropsFromStore(fieldState,props) {
    const {message= '', valid= true, visible= true, value= '', displayValue= '',
           tooltip= '', validator= defValidatorFunc,
           ...theRestOfFieldState}= fieldState;

    const tmpProps= Object.assign({ message, valid, visible, value, displayValue, tooltip, validator},
                                    theRestOfFieldState, props);

    return omit(tmpProps, STORE_OMIT_LIST);
}


const getGroupKey= (props,context)=> props.groupKey || get(context,'groupKey');
