/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React, {Component, PureComponent} from 'react';
import PropTypes from 'prop-types';
import shallowequal from 'shallowequal';
import {get, pick} from 'lodash';

import {getDropDownInfo} from '../core/LayoutCntlr.js';
import {flux, getVersion} from '../Firefly.js';
import {SearchPanel} from '../ui/SearchPanel.jsx';
import {ImageSearchDropDown} from '../visualize/ui/ImageSearchPanelV2.jsx';
import {TestSearchPanel} from '../ui/TestSearchPanel.jsx';
import {TestQueriesPanel} from '../ui/TestQueriesPanel.jsx';
import {ChartSelectDropdown} from '../ui/ChartSelectDropdown.jsx';
import {CatalogSelectViewPanel} from '../visualize/ui/CatalogSelectViewPanel.jsx';
import {LSSTCatalogSelectViewPanel} from '../visualize/ui/LSSTCatalogSelectViewPanel.jsx';
import {FileUploadDropdown} from '../ui/FileUploadDropdown.jsx';
import {WorkspaceDropdown} from '../ui/WorkspaceDropdown.jsx';
import {getAlerts} from '../core/AppDataCntlr.js';

import './DropDownContainer.css';

export const dropDownMap = {
    Search: <SearchPanel />,
    TestSearch: <TestSearchPanel />,
    TestSearches: <TestQueriesPanel />,
    ImageSearchPanelV2: <ImageSearchDropDown/>,
    ImageSelectDropDownCmd: <ImageSearchDropDown/>,
    ImageSelectDropDownSlateCmd: <ImageSearchDropDown gridSupport={true}/>,
    ChartSelectDropDownCmd: <ChartSelectDropdown />,
    IrsaCatalogDropDown: <CatalogSelectViewPanel/>,
    LsstCatalogDropDown: <LSSTCatalogSelectViewPanel/>,
    FileUploadDropDownCmd: <FileUploadDropdown />,
    WorkspaceDropDownCmd: <WorkspaceDropdown />
};




/**
 * The container for items appearing in the drop down panel.
 * This container mimic a card layout in which it will accept multiple cards.
 * However, only one selected card will be displayed at a time.
 * Items in this container must have a 'name' property.  It will be used to
 * compare to the selected card.
 */
export class DropDownContainer extends Component {
    constructor(props) {
        super(props);

        React.Children.forEach(this.props.children, (el) => {
            const key = get(el, 'props.name');
            if (key) dropDownMap[key] = el;
        });

        if (props.dropdownPanels) {
            props.dropdownPanels.forEach( (el) => {
                const key = get(el, 'props.name');
                if (key) dropDownMap[key] = el;
            } );
        }

        this.state = {
                visible: props.visible,
                selected: props.selected
            };
    }

    componentDidMount() {
        this.removeListener= flux.addListener(() => this.storeUpdate());
        this.iAmMounted = true;
    }

    componentWillUnmount() {
        this.iAmMounted = false;
        this.removeListener && this.removeListener();
    }
    
    shouldComponentUpdate(nProps, nState) {
        const check = ['visible','selected'];
        return !shallowequal(pick(nState, check), pick(this.state, check));
   }

    storeUpdate() {
        if (this.iAmMounted) {
            const {visible, view} = getDropDownInfo();
            if (visible !== this.state.visible || view !== this.state.selected) {
                this.setState({visible, selected: view});
            }
        }
    }

    render() {
        const {footer, alerts} = this.props;
        const { visible, selected }= this.state;
        const view = dropDownMap[selected];

        if (!visible) return <div/>;
        return (
            <div>
                <div className='DD-ToolBar'>
                    {alerts || <Alerts />}
                    <div style={{flexGrow: 1}}>
                        <div className='DD-ToolBar__content'>
                            {view}
                        </div>
                    </div>
                    <div id='footer' className='DD-ToolBar__footer'>
                        {footer}
                        <div className='DD-ToolBar__version'>
                            {getVersion()}
                        </div>
                    </div>
                </div>
            </div>
        );
    }
}

DropDownContainer.propTypes = {
    visible: PropTypes.bool,
    selected: PropTypes.string,
    dropdownPanels: PropTypes.arrayOf(PropTypes.element),
    footer: PropTypes.node,
    alerts: PropTypes.node
};
DropDownContainer.defaultProps = {
    visible: false
};

export class Alerts extends PureComponent {

    constructor(props) {
        super(props);
        this.state = Object.assign({}, props);
    }

    componentDidMount() {
        this.removeListener= flux.addListener(() => this.storeUpdate());
        this.iAmMounted = true;
    }

    componentWillUnmount() {
        this.iAmMounted = false;
        this.removeListener && this.removeListener();
    }

    storeUpdate() {
        if (this.iAmMounted) {
            this.setState(getAlerts());
        }
    }

    render() {
        const {msg, style} = this.state;
        if (msg) {
            /* eslint-disable react/no-danger */
            return (
                <div className='alerts__msg' style={style}>
                    <div dangerouslySetInnerHTML={{__html: msg}} />
                </div>
            );
        } else return <div/>;
    }
}

Alerts.propTypes = {
    msg: PropTypes.string,
    style: PropTypes.object
};
