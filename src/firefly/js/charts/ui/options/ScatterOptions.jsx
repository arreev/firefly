import React from 'react';
import {get, isUndefined} from 'lodash';

import {getChartData} from '../../ChartsCntlr.js';
import {FieldGroup} from '../../../ui/FieldGroup.jsx';
import {VALUE_CHANGE} from '../../../fieldGroup/FieldGroupCntlr.js';

//import {ValidationField} from '../../../ui/ValidationField.jsx';
import {ListBoxInputField} from '../../../ui/ListBoxInputField.jsx';
import {BasicOptionFields, OptionTopBar, basicFieldReducer, submitChanges} from './BasicOptions.jsx';
import {updateSet} from '../../../util/WebUtil.js';
import {SimpleComponent} from '../../../ui/SimpleComponent.jsx';
import {getColValStats} from '../../TableStatsCntlr.js';
import {ColumnOrExpression} from '../ColumnOrExpression.jsx';
import {Errors, errorTypeFieldKey, errorFieldKey, errorMinusFieldKey} from './Errors.jsx';

const fieldProps = {labelWidth: 62, size: 15};

export class ScatterOptions extends SimpleComponent {

    getNextState() {
        const {chartId} = this.props;
        const {activeTrace:cActiveTrace} = getChartData(chartId);
        // activeTrace is passed via property, when used from NewTracePanel
        const activeTrace = isUndefined(this.props.activeTrace) ? cActiveTrace : this.props.activeTrace;
        return {activeTrace};
    }

    render() {
        const {chartId} = this.props;
        //const {activeTrace=0} = this.state;
        const {tablesources, data, layout, activeTrace:cActiveTrace=0} = getChartData(chartId);
        const activeTrace = isUndefined(this.props.activeTrace) ? cActiveTrace : this.props.activeTrace;
        const groupKey = this.props.groupKey || `${chartId}-scatter-${activeTrace}`;
        const tablesource = get(tablesources, [cActiveTrace]);
        const tbl_id = get(tablesource, 'tbl_id');

        return (
            <div style={{padding:'0 5px 7px'}}>
                {isUndefined(this.props.activeTrace) && <OptionTopBar {...{groupKey, activeTrace, chartId, tbl_id, submitChangesFunc: submitChangesScatter}}/>}
                <FieldGroup className='FieldGroup__vertical' keepState={false} groupKey={groupKey} reducerFunc={fieldReducer({data, layout, activeTrace, tablesources})}>
                    <ListBoxInputField fieldKey={`data.${activeTrace}.mode`} options={[{value:'markers'}, {value:'lines'}, {value:'lines+markers'}]}/>
                    <ListBoxInputField fieldKey={`data.${activeTrace}.marker.symbol`}
                                       options={[{value:'circle'}, {value:'circle-open'}, {value:'square'}, {value:'square-open'}, {value:'diamond'}, {value:'diamond-open'},
                                                 {value:'cross'}, {value:'x'}, {value:'triangle-up'}, {value:'hexagon'}, {value:'star'}]}/>
                    {tablesource && <TableSourcesOptions {...{tablesource, activeTrace, groupKey}}/>}
                    <br/>
                    <BasicOptionFields {...{layout, data, activeTrace}}/>
                </FieldGroup>
            </div>
        );
    }
}

export function fieldReducer({data, layout, activeTrace, tablesources={}}) {
    const tablesourceMappings = get(tablesources[activeTrace], 'mappings');
    const basicReducer = basicFieldReducer({data, layout, activeTrace, tablesources});
    const fields = {
        [`data.${activeTrace}.mode`]: {
            fieldKey: `data.${activeTrace}.mode`,
            value: get(data, `${activeTrace}.mode`),
            tooltip: 'Select plot style',
            label: 'Plot Style:',
            ...fieldProps
        },
        [`data.${activeTrace}.marker.symbol`]: {
            fieldKey: `data.${activeTrace}.marker.symbol`,
            value: get(data, `${activeTrace}.marker.symbol`),
            tooltip: 'Select marker symbol',
            label: 'Symbol:',
            ...fieldProps
        },
        [errorTypeFieldKey(activeTrace, 'x')]: {
            fieldKey: errorTypeFieldKey(activeTrace, 'x'),
            value: get(data, `${activeTrace}.fireflyData.options.error_x.errorsType`, 'none')
        },
        [errorTypeFieldKey(activeTrace, 'y')]: {
            fieldKey: errorTypeFieldKey(activeTrace, 'y'),
            value: get(data, `${activeTrace}.fireflyData.options.error_y.errorsType`, 'none')
        },
        ...basicReducer(null)
    };
    const tblRelFields = {
        [`_tables.data.${activeTrace}.x`]: {
            fieldKey: `_tables.data.${activeTrace}.x`,
            value: get(tablesourceMappings, 'x', ''),
            //tooltip: 'X axis',
            label: 'X:',
            ...fieldProps
        },
        [`_tables.data.${activeTrace}.y`]: {
            fieldKey: `_tables.data.${activeTrace}.y`,
            value: get(tablesourceMappings, 'y', ''),
            //tooltip: 'Y axis',
            label: 'Y:',
            ...fieldProps
        },
        [errorFieldKey(activeTrace, 'x')]: {
            fieldKey: errorFieldKey(activeTrace, 'x'),
            value: get(tablesourceMappings, ['error_x.array'], ''),
            //tooltip: 'X error',
            label: 'X error\u2191:',
            ...fieldProps
        },
        [errorMinusFieldKey(activeTrace, 'x')]: {
            fieldKey: errorMinusFieldKey(activeTrace, 'x'),
            value: get(tablesourceMappings, ['error_x.arrayminus'], ''),
            //tooltip: 'X error',
            label: 'X Err\u2193:',
            ...fieldProps
        },
        [errorFieldKey(activeTrace, 'y')]: {
            fieldKey: errorFieldKey(activeTrace, 'y'),
            value: get(tablesourceMappings, ['error_y.array'], ''),
            //tooltip: '',
            label: 'Y error\u2191:',
            ...fieldProps
        },
        [errorMinusFieldKey(activeTrace, 'y')]: {
            fieldKey: errorMinusFieldKey(activeTrace, 'y'),
            value: get(tablesourceMappings, ['error_y.arrayminus'], ''),
            //tooltip: 'Y error',
            label: 'Y error\u2193:',
            ...fieldProps
        },
        [`_tables.data.${activeTrace}.marker.color`]: {
            fieldKey: `_tables.data.${activeTrace}.marker.color`,
            value: get(tablesourceMappings, 'marker.color', ''),
            //tooltip: 'Use a column for color map',
            label: 'Color Map:',
            ...fieldProps
        },
        [`_tables.data.${activeTrace}.marker.size`]: {
            fieldKey: `_tables.data.${activeTrace}.marker.size`,
            value: get(tablesourceMappings, 'marker.size', ''),
            //tooltip: 'Use a column for size map',
            label: 'Size Map:',
            ...fieldProps
        }
    };
    return (inFields, action) => {
        if (!inFields) {
            return tablesourceMappings? Object.assign({}, fields, tblRelFields) : fields;
        }

        const {payload:{fieldKey='', value=''}, type} = action;

        if (fieldKey.endsWith('marker.color') && type === VALUE_CHANGE && value.length === 1) {
            if (fieldKey.startsWith('_tables')) {
                const colorKey = Object.keys(inFields).find((k) => k.match(/data.+.marker.color/)) || '';
                if (colorKey) inFields = updateSet(inFields, [colorKey, 'value'], '');     // blanks out color when a color map is entered
            } else {
                const colorMapKey = Object.keys(inFields).find((k) => k.match(/_tables.+.marker.color/)) || '';
                if (colorMapKey) inFields = updateSet(inFields, [colorMapKey, 'value'], '');   // blanks out color map when a color is entered
            }
        }
        return inFields;

    };
}

export function TableSourcesOptions({tablesource={}, activeTrace, groupKey}) {
    // _tables.  is prefixed the fieldKey.  it will be replaced with 'tables::tbl_id,val' on submitChanges.
    const tbl_id = get(tablesource, 'tbl_id');
    const colValStats = getColValStats(tbl_id);
    const labelWidth = 30;
    const xProps = {fldPath:`_tables.data.${activeTrace}.x`, label: 'X:', name: 'X', nullAllowed: false, colValStats, groupKey, labelWidth};
    const yProps = {fldPath:`_tables.data.${activeTrace}.y`, label: 'Y:', name: 'Y', nullAllowed: false, colValStats, groupKey, labelWidth};

    const commonProps = {colValStats, groupKey, labelWidth: 62, nullAllowed: true};
    const flds = [
        {fldPath:`_tables.data.${activeTrace}.marker.color`,label: 'Color Map:', name: 'Color Map'},
        {fldPath:`_tables.data.${activeTrace}.marker.size`,label: 'Size Map:', name: 'Size Map'}
    ].map((e) => {return Object.assign(e, commonProps);});


    return (
        <div className='FieldGroup__vertical'>
            <br/>
            <ColumnOrExpression {...xProps}/>
            <Errors axis='x' {...{groupKey, colValStats, activeTrace, labelWidth}}/>
            <br/>
            <ColumnOrExpression {...yProps}/>
            <Errors axis='y' {...{groupKey, colValStats, activeTrace, labelWidth}}/>
            <br/>
            {
                flds.map((props,idx) => {
                    return <ColumnOrExpression key={idx} {...props}/>;
                })
            }
        </div>
    );
}

export function submitChangesScatter({chartId, activeTrace, fields, tbl_id}) {
    const colorMap = get(fields, `_tables.data.${activeTrace}.marker.color`);
    const sizeMap = get(fields, `_tables.data.${activeTrace}.marker.size`);

    const dataType = (!tbl_id || colorMap || sizeMap) ? 'scatter' : 'fireflyScatter';
    const changes = {[`fireflyData.${activeTrace}.dataType`] : dataType};
    if (dataType === 'fireflyScatter') {
        // add a mapping for rowIdx
        changes[`_tables.data.${activeTrace}.firefly.rowIdx`] = 'rowIdx'; // rowIdx is mapping table rows to data points
    }
    Object.assign(changes, fields);
    submitChanges({chartId, fields: changes, tbl_id});
}