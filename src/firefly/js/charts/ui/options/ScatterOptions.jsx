import React from 'react';
import PropTypes from 'prop-types';
import {isArray, get, isUndefined, omit} from 'lodash';

import {Expression} from '../../../util/expr/Expression.js';
import {getChartData, getTraceSymbol, hasUpperLimits} from '../../ChartsCntlr.js';
import {FieldGroup} from '../../../ui/FieldGroup.jsx';
import {VALUE_CHANGE} from '../../../fieldGroup/FieldGroupCntlr.js';

import {ListBoxInputField} from '../../../ui/ListBoxInputField.jsx';
import {RadioGroupInputField} from '../../../ui/RadioGroupInputField.jsx';
import {BasicOptionFields, basicFieldReducer, helpStyle, submitChanges} from './BasicOptions.jsx';
import {updateSet} from '../../../util/WebUtil.js';
import {SimpleComponent} from '../../../ui/SimpleComponent.jsx';
import {getColValStats} from '../../TableStatsCntlr.js';
import {ColumnOrExpression} from '../ColumnOrExpression.jsx';
import {Errors, errorTypeFieldKey, errorFieldKey, errorMinusFieldKey} from './Errors.jsx';
import {getAppOptions} from '../../../core/AppDataCntlr.js';

const fieldProps = {labelWidth: 62, size: 15};

/**
 * Should we display Upper Limit field under Y?
 * @returns {*}
 */
export function upperLimitUI() {
    return get(getAppOptions(), 'charts.upperLimitUI');
}

export class ScatterOptions extends SimpleComponent {

    getNextState() {
        const {chartId} = this.props;
        const {activeTrace:cActiveTrace=0} = getChartData(chartId);
        // activeTrace is passed via property, when used from NewTracePanel
        const activeTrace = isUndefined(this.props.activeTrace) ? cActiveTrace : this.props.activeTrace;
        return {activeTrace};
    }

    render() {
        const {chartId, groupKey:groupKeyProp, activeTrace:activeTraceProp, tbl_id:tblIdProp,showMultiTrace} = this.props;
        const {tablesources, activeTrace:cActiveTrace=0} = getChartData(chartId);
        const activeTrace = isUndefined(activeTraceProp) ? cActiveTrace : activeTraceProp;
        const groupKey = groupKeyProp || `${chartId}-scatter-${activeTrace}`;
        const tablesource = get(tablesources, [cActiveTrace], tblIdProp && {tbl_id: tblIdProp});

        const modeKey = `data.${activeTrace}.mode`;
        const modeOptions =[{label: 'points', value:'markers'},
            {label: 'connected points', value:'lines+markers'},
            {label: 'lines', value:'lines'}];

        return (
            <FieldGroup className='FieldGroup__vertical' keepState={false} groupKey={groupKey} reducerFunc={fieldReducer({chartId, activeTrace})}>

                {!showMultiTrace && <RadioGroupInputField alignment='horizontal' fieldKey={modeKey} options={modeOptions}/>}
                {showMultiTrace && <ListBoxInputField fieldKey={modeKey} options={modeOptions}/>}
                {showMultiTrace && <ListBoxInputField fieldKey={`data.${activeTrace}.marker.symbol`}
                                                      options={[{value:'circle'}, {value:'circle-open'}, {value:'square'}, {value:'square-open'}, {value:'diamond'}, {value:'diamond-open'},
                                                          {value:'cross'}, {value:'x'}, {value:'triangle-up'}, {value:'hexagon'}, {value:'star'}]}/>}
                {/* TODO: scattergl does not support 'open' symbols as of v1..28.2.  we'll add them back at a later time when they do.
                 options={[{value:'circle'}, {value:'square'}, {value:'diamond'},
                 {value:'cross'}, {value:'x'}, {value:'triangle-up'}, {value:'hexagon'}, {value:'star'}]}/>
                 */}

                {tablesource && <TableSourcesOptions {...{chartId, tablesource, activeTrace, groupKey,showMultiTrace}}/>}
                <BasicOptionFields {...{activeTrace, groupKey,showMultiTrace}}/>
            </FieldGroup>
        );
    }
}

export function fieldReducer({chartId, activeTrace}) {

    const basicReducer = basicFieldReducer({chartId, activeTrace});

    const getFields = () => {
        const {data, fireflyData, tablesources={}} = getChartData(chartId);
        const tablesourceMappings = get(tablesources[activeTrace], 'mappings');

        // when a symbol is substituted with an array,
        // the selected symbol is saved in fireflyData
        const symbol = getTraceSymbol(data, fireflyData, activeTrace);

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
                value: symbol,
                tooltip: 'Select marker symbol',
                label: 'Symbol:',
                ...fieldProps
            },
            [`data.${activeTrace}.marker.colorscale`]: {
                fieldKey: `data.${activeTrace}.marker.colorscale`,
                value: get(data, `${activeTrace}.marker.colorscale`),
                tooltip: 'Select colorscale for color map',
                label: 'Color Scale:',
                ...fieldProps
            },
            [errorTypeFieldKey(activeTrace, 'x')]: {
                fieldKey: errorTypeFieldKey(activeTrace, 'x'),
                value: get(data, errorTypeFieldKey(activeTrace, 'x').replace(/^data./, ''), 'none')
            },
            [errorTypeFieldKey(activeTrace, 'y')]: {
                fieldKey: errorTypeFieldKey(activeTrace, 'y'),
                value: get(data, errorTypeFieldKey(activeTrace, 'y').replace(/^data./, ''), 'none')
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
            [`_tables.fireflyData.${activeTrace}.yMax`]: {
                fieldKey: `_tables.fireflyData.${activeTrace}.yMax`,
                value: get(tablesourceMappings, `fireflyData.${activeTrace}.yMax`, ''),
                label: 'Limit:',
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
                label: 'Error\u2193:',
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
        return tablesourceMappings? Object.assign({}, fields, tblRelFields) : fields;
    };

    return (inFields, action) => {
        if (!inFields) {
            return getFields();
        }

        inFields = basicReducer(inFields, action);

        const {payload:{fieldKey='', value=''}, type} = action;

        if (type === VALUE_CHANGE) {
            if (fieldKey.endsWith('marker.color') && value.length === 1) {
                if (fieldKey.startsWith('_tables')) {
                    const colorKey = Object.keys(inFields).find((k) => k.match(/data.+.marker.color$/)) || '';
                    if (colorKey) inFields = updateSet(inFields, [colorKey, 'value'], '');     // blanks out color when a color map is entered
                } else {
                    const colorMapKey = Object.keys(inFields).find((k) => k.match(/_tables.+.marker.color$/)) || '';
                    if (colorMapKey) inFields = updateSet(inFields, [colorMapKey, 'value'], '');   // blanks out color map when a color is entered
                }
            }
            // when field changes, clear error fields
            ['x','y'].forEach((a) => {
                if (fieldKey === `_tables.data.${activeTrace}.${a}`) {
                    inFields = updateSet(inFields, [errorTypeFieldKey(activeTrace, `${a}`), 'value'], 'none');
                    inFields = updateSet(inFields, [errorFieldKey(activeTrace, `${a}`), 'value'], undefined);
                    inFields = updateSet(inFields, [errorMinusFieldKey(activeTrace, `${a}`), 'value'], undefined);
                }
            });
        }
        return inFields;

    };
}

export function TableSourcesOptions({chartId, tablesource={}, activeTrace, groupKey, showMultiTrace}) {
    // _tables.  is prefixed the fieldKey.  it will be replaced with 'tables::val' on submitChanges.
    const tbl_id = get(tablesource, 'tbl_id');
    const colValStats = getColValStats(tbl_id);
    if (!colValStats) { return null; }
    const labelWidth = 30;
    const xProps = {fldPath:`_tables.data.${activeTrace}.x`, label: 'X:', name: 'X', nullAllowed: false, colValStats, groupKey, labelWidth};
    const yProps = {fldPath:`_tables.data.${activeTrace}.y`, label: 'Y:', name: 'Y', nullAllowed: false, colValStats, groupKey, labelWidth};
    const yMaxProps = {fldPath:`_tables.fireflyData.${activeTrace}.yMax`, label: 'Limit:', name: 'Upper Limit', nullAllowed: true, colValStats, groupKey, labelWidth};
    
    const commonProps = {colValStats, groupKey, labelWidth: 62, nullAllowed: true};
    const sizemapTooltip = 'marker size. Please use expression to convert column value to valid pixels';
    const flds = [
        {key: 'colorMap', fldPath:`_tables.data.${activeTrace}.marker.color`,label: 'Color Map:', name: 'Color Map'},
        {key: 'sizeMap', fldPath:`_tables.data.${activeTrace}.marker.size`,label: 'Size Map:', name: 'Size Map', tooltip: sizemapTooltip}
    ].map((e) => {return Object.assign(e, commonProps);});
    const colorMapProps = flds[0];
    const sizeMapProps = flds[1];

    return (
        <div className='FieldGroup__vertical'>
            <br/>
            <div style={helpStyle}>
                For X and Y, enter a column or an expression<br/>
                ex. log(col); 100*col1/col2; col1-col2
            </div>
            <ColumnOrExpression {...xProps}/>
            <Errors axis='x' {...{groupKey, colValStats, activeTrace, labelWidth}}/>
            <br/>
            <ColumnOrExpression {...yProps}/>
            {(upperLimitUI() || hasUpperLimits(chartId, activeTrace)) && <ColumnOrExpression {...yMaxProps}/>}
            <Errors axis='y' {...{groupKey, colValStats, activeTrace, labelWidth}}/>
            {showMultiTrace &&  <div style={{paddingTop: 10}}>
                <ColumnOrExpression {...sizeMapProps}/>
                <ColumnOrExpression {...colorMapProps}/>
                <ListBoxInputField fieldKey={`data.${activeTrace}.marker.colorscale`}
                                   options={[{value: 'Default'}, {value: 'Bluered'}, {value: 'Blues'}, {value: 'Earth'}, {value: 'Electric'}, {value: 'Greens'},
                                       {value: 'Greys'}, {value: 'Hot'}, {value: 'Jet'}, {value: 'Picnic'}, {value: 'Portland'}, {value: 'Rainbow'},
                                       {value: 'RdBu'}, {value: 'Reds'}, {value: 'Viridis'}, {value: 'YlGnBu'}, {value: 'YlOrRd'}]}/>
            </div>
            }
        </div>
    );
}

TableSourcesOptions.propTypes = {
    tablesource: PropTypes.object,
    activeTrace: PropTypes.number,
    groupKey: PropTypes.string,
    showMultiTrace: PropTypes.bool
};

export function submitChangesScatter({chartId, activeTrace, fields, tbl_id}) {

    const changes = {[`data.${activeTrace}.type`] : getTraceType(chartId, activeTrace)};

    // check if size field is a constant
    const sizeMap = fields[`_tables.data.${activeTrace}.marker.size`];
    if (sizeMap) {
        const colValStats = getColValStats(tbl_id);
        const colNames = colValStats.map((colVal) => {return colVal.name;});
        const expr = new Expression(sizeMap, colNames);
        if (expr.isValid() && (expr.getParsedVariables().length === 0)) {
            const symSize = expr.getValue();
            changes[`data.${activeTrace}.marker.size`] = symSize;
            fields = omit(fields, `_tables.data.${activeTrace}.marker.size`);
        }
    }

    Object.assign(changes, fields);
    submitChanges({chartId, fields: changes, tbl_id});
}

/**
 * Returns gl or non-gl scatter type based on already used trace type
 * (GL and non-GL traces do not work well together)
 * @param chartId
 * @param activeTrace
 */
function getTraceType(chartId, activeTrace) {
    const chartData = getChartData(chartId);
    let type = get(chartData, `data.${activeTrace}.type`);
    if (isUndefined(type)) {
        if (activeTrace > 0) {
            // check previous trace type
            const isGL = get(chartData, `data.${activeTrace-1}.type`, 'scatter').endsWith('gl');
            type = isGL ? 'scattergl' : 'scatter';
        } else {
            type = 'scatter';
        }
    }
    return type;
}