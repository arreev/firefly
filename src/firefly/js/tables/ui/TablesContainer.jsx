/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React, {PureComponent} from 'react';
import PropTypes from 'prop-types';
import {isEmpty, get} from 'lodash';

import {flux} from '../../Firefly.js';
import * as TblUtil from '../TableUtil.js';
import {TablePanel} from './TablePanel.jsx';
import {Tabs, Tab} from '../../ui/panel/TabPanel.jsx';
import {dispatchTableRemove, dispatchActiveTableChanged} from '../TablesCntlr.js';


import {LO_VIEW, LO_MODE, dispatchSetLayoutMode, getExpandedMode} from '../../core/LayoutCntlr.js';
import {CloseButton} from '../../ui/CloseButton.jsx';


export class TablesContainer extends PureComponent {
    constructor(props) {
        super(props);
        this.state = this.nextState(props);
    }

    componentDidMount() {
        this.removeListener= flux.addListener(() => this.storeUpdate());
    }

    componentWillUnmount() {
        this.removeListener && this.removeListener();
        this.isUnmounted = true;
    }

    nextState(props) {
        var {mode, tbl_group, closeable, tableOptions} = props;
        const expandedMode = props.expandedMode || getExpandedMode() === LO_VIEW.tables;
        if (expandedMode && mode !== 'standard') {
            tbl_group = TblUtil.getTblExpandedInfo().tbl_group;
        }
        const {tables, layout, active} = TblUtil.getTableGroup(tbl_group) || {};

        return {closeable, tbl_group, expandedMode, tables, tableOptions, layout, active};
    }

    storeUpdate() {
        if (!this.isUnmounted) {
            this.setState(this.nextState(this.props));
        }
    }

    render() {
        const {closeable, tbl_group, expandedMode, tables, tableOptions, layout, active} = this.state;
        if (expandedMode) {
            return <ExpandedView {...{active, tables, tableOptions, layout, expandedMode, closeable, tbl_group}} />;
        } else {
            return isEmpty(tables) ? <div></div> : <StandardView {...{active, tables, tableOptions, expandedMode, tbl_group}} />;
        }
    }
}

TablesContainer.propTypes = {
    expandedMode: PropTypes.bool,
    closeable: PropTypes.bool,
    tbl_group: PropTypes.string,
    mode: PropTypes.oneOf(['expanded', 'standard', 'both'])
};
TablesContainer.defaultProps = {
    expandedMode: false,
    closeable: true,
    mode: 'standard'
};



function ExpandedView(props) {
    const {tables, closeable} = props;
    return (
        <div style={{ display: 'flex', height: '100%', flexGrow: 1, flexDirection: 'column', overflow: 'hidden'}}>
            <div style={{marginBottom: 3}}>
                {closeable && <CloseButton style={{display: 'inline-block', paddingLeft: 10}} onClick={() => dispatchSetLayoutMode(LO_MODE.expanded, LO_VIEW.none)}/>}
            </div>
            {!isEmpty(tables) &&
                <div style={{position: 'relative', flexGrow: 1}}>
                    <div style={{position: 'absolute', top:0, right:0, bottom:0, left:0}}>
                         <StandardView expandedMode={true} {...props} />
                    </div>
                </div>
            }
        </div>
    );

}


function StandardView(props) {
    const {tables, tableOptions, expandedMode, active, tbl_group} = props;

    var activeIdx = Object.keys(tables).findIndex( (tbl_ui_id) => get(tables,[tbl_ui_id,'tbl_id']) === active);
    activeIdx = activeIdx === -1 ? 0 : activeIdx;
    const onTabSelect = (idx) => {
        const tbl_id = get(tables, [Object.keys(tables)[idx], 'tbl_id']);
        tbl_id && dispatchActiveTableChanged(tbl_id, tbl_group);
    };
    const keys = Object.keys(tables);
    if (keys.length === 1) {
        return <SingleTable table={get(tables, [keys[0]])} expandedMode={expandedMode} tableOptions={tableOptions}/>;
    } else {
        return (
            <Tabs defaultSelected={activeIdx} onTabSelect={onTabSelect} resizable={true}>
                {tablesAsTab(tables, tableOptions, expandedMode)}
            </Tabs>
        );
    }
}

function SingleTable({table, tableOptions, expandedMode}) {
    var {tbl_id, title, removable, tbl_ui_id, options={}} = table;

    options = Object.assign({}, options, tableOptions);

    return  (
        <TablePanel key={tbl_id} border={true} {...{title, removable, tbl_id, tbl_ui_id, ...options, expandedMode}} />
    );
}

function tablesAsTab(tables, tableOptions, expandedMode) {

    return tables &&
        Object.keys(tables).map( (key) => {
            var {tbl_id, title, removable, tbl_ui_id, options={}} = tables[key];
            options = Object.assign({}, options, tableOptions);
            const onTabRemove = () => {
                dispatchTableRemove(tbl_id);
            };
            return  (
                <Tab key={tbl_ui_id} name={title} removable={removable} onTabRemove={onTabRemove}>
                    <TablePanel key={tbl_id} border={false} showTitle={false} {...{tbl_id, tbl_ui_id, ...options, expandedMode}} />
                </Tab>
            );
        } );
}

