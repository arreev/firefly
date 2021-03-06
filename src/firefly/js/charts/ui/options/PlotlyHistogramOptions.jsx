import React from 'react';
import {get, isUndefined} from 'lodash';

import {FieldGroup} from '../../../ui/FieldGroup.jsx';
import {ValidationField} from '../../../ui/ValidationField.jsx';
import {ListBoxInputField} from '../../../ui/ListBoxInputField.jsx';
import {ColumnOrExpression} from '../ColumnOrExpression.jsx';
import {BasicOptionFields, basicFieldReducer} from './BasicOptions.jsx';
import {getChartData} from '../../ChartsCntlr.js';

import {SimpleComponent} from '../../../ui/SimpleComponent.jsx';
import {getColValStats} from '../../TableStatsCntlr.js';

const fieldProps = {labelWidth: 62, size: 15};

/**
 * This are the options for Plotly histogram chart
 * Plotly histogram does not display well with firefly histogram, which is using Plotly bar chart
 * Known issues: shen shown with firefly histogram, Plotly histogram shows as lines or does not show at all
 */
export class HistogramOptions extends SimpleComponent {

    getNextState() {
        const {chartId} = this.props;
        const {activeTrace:cActiveTrace=0} = getChartData(chartId);
        // activeTrace is passed via property, when used from NewTracePanel
        const activeTrace = isUndefined(this.props.activeTrace) ? cActiveTrace : this.props.activeTrace;
        return {activeTrace};
    }

    render() {
        const {chartId, tbl_id:tblIdProp} = this.props;
        const {tablesources, activeTrace:cActiveTrace=0} = getChartData(chartId);
        const activeTrace = isUndefined(this.props.activeTrace) ? cActiveTrace : this.props.activeTrace;
        const groupKey = this.props.groupKey || `${chartId}-ffhist-${activeTrace}`;
        const tablesource = get(tablesources, [cActiveTrace], tblIdProp && {tbl_id: tblIdProp});
        const tbl_id = get(tablesource, 'tbl_id');
        const colValStats = getColValStats(tbl_id);
        const xProps = {fldPath:`_tables.data.${activeTrace}.x`, label: 'X:', name: 'X', nullAllowed: false, colValStats, groupKey, labelWidth: 62};
        return (
            <FieldGroup className='FieldGroup__vertical' keepState={false} groupKey={groupKey} reducerFunc={fieldReducer({chartId, activeTrace})}>
                    <ColumnOrExpression {...xProps}/>
                    <ListBoxInputField fieldKey={`data.${activeTrace}.histfunc`} options={[{value:'count'}, {value:'sum'}, {value:'avg'}, {value:'min'}, {value:'max'}]}/>
                    <ValidationField fieldKey={`data.${activeTrace}.nbinsx`}/>
                    <ValidationField fieldKey={`data.${activeTrace}.xbins.size`}/>
                    <ValidationField fieldKey={`data.${activeTrace}.xbins.start`}/>
                    <ValidationField fieldKey={`data.${activeTrace}.xbins.end`}/>
                    <br/>
                    <BasicOptionFields {...{activeTrace, groupKey}}/>
                </FieldGroup>
        );
    }
}

export function fieldReducer({chartId, activeTrace}) {
    const basicReducer = basicFieldReducer({chartId, activeTrace});

    const getFields = () => {
        const {data, tablesources = {}} = getChartData(chartId);
        const tablesourceMappings = get(tablesources[activeTrace], 'mappings');

        const fields = {
            [`_tables.data.${activeTrace}.x`]: {
                fieldKey: `_tables.data.${activeTrace}.x`,
                value: get(tablesourceMappings, 'x', ''),
                tooltip: 'X axis',
                label: 'X:',
                ...fieldProps
            },
            [`data.${activeTrace}.histfunc`]: {
                fieldKey: `data.${activeTrace}.histfunc`,
                value: get(data, `${activeTrace}.histfunc`),
                tooltip: 'Binning function used for this histogram trace',
                label: 'Function:',
                ...fieldProps
            },
            [`data.${activeTrace}.nbinsx`]: {
                fieldKey: `data.${activeTrace}.nbinsx`,
                value: get(data, `${activeTrace}.nbinsx`),
                tooltip: 'Maximum number of desired bins',
                label: 'Num of bins:',
                ...fieldProps
            },
            [`data.${activeTrace}.xbins.size`]: {
                fieldKey: `data.${activeTrace}.xbins.size`,
                value: get(data, `${activeTrace}.xbins.size`),
                tooltip: 'Step in between value each x axis bin',
                label: 'Bin width:',
                ...fieldProps
            },
            [`data.${activeTrace}.xbins.start`]: {
                fieldKey: `data.${activeTrace}.xbins.start`,
                value: get(data, `${activeTrace}.xbins.start`),
                tooltip: 'Sets the starting value for the x axis bins',
                label: 'Min:',
                ...fieldProps
            },
            [`data.${activeTrace}.xbins.end`]: {
                fieldKey: `data.${activeTrace}.xbins.end`,
                value: get(data, `${activeTrace}.xbins.end`),
                tooltip: 'Sets the end value for the x axis bins.',
                label: 'Max:',
                ...fieldProps
            },
            ...basicReducer(null)
        };
        return fields;
    };
    return (inFields, action) => {
        if (!inFields) {
            return getFields();
        } else {
            return inFields;
        }
    };
}
