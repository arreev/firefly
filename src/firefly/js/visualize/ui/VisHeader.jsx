/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React, {PureComponent} from 'react';
import PropTypes from 'prop-types';
import {visRoot} from '../ImagePlotCntlr.js';
import {flux} from '../../Firefly.js';
import {VisHeaderView, VisPreview} from './VisHeaderView.jsx';
import {addMouseListener, lastMouseCtx} from '../VisMouseSync.js';
import {readoutRoot} from '../../visualize/MouseReadoutCntlr.js';


export class VisHeader_old extends PureComponent {
    constructor(props) {
        super(props);
        this.state= {visRoot:visRoot(), currMouseState:lastMouseCtx()};
    }

    componentWillUnmount() {
        if (this.removeListener) this.removeListener();
        if (this.removeMouseListener) this.removeMouseListener();
    }


    componentDidMount() {
        this.removeListener= flux.addListener(() => this.storeUpdate());
        this.removeMouseListener= addMouseListener(() => this.storeUpdate());
    }

    storeUpdate() {
        if (visRoot()!==this.state.visRoot || lastMouseCtx() !==this.state.currMouseState) {
            this.setState({visRoot:visRoot(), currMouseState:lastMouseCtx()});
        }
    }

    render() {
        var {visRoot,currMouseState}= this.state;
        return <VisHeaderView visRoot={visRoot} currMouseState={currMouseState}/>;
    }
}



export class VisHeader extends PureComponent {
    constructor(props) {
        super(props);
        this.state= {visRoot:visRoot(), currMouseState:lastMouseCtx(), readout:readoutRoot()};
    }

    componentWillUnmount() {
        if (this.removeListener) this.removeListener();
        if (this.removeMouseListener) this.removeMouseListener();
    }


    componentDidMount() {
        this.removeListener= flux.addListener(() => this.storeUpdate());
        this.removeMouseListener= addMouseListener(() => this.storeUpdate());
    }

    storeUpdate() {
        const readout= readoutRoot();
        if (visRoot()!==this.state.visRoot || lastMouseCtx() !==this.state.currMouseState || readout!==this.state.readout) {
            this.setState({visRoot:visRoot(), currMouseState:lastMouseCtx(), readout:readout});
        }


    }

    render() {
        const {showHeader=true, showPreview=true} = this.props; 
        var {visRoot,currMouseState,readout}= this.state;
        return (
            <div>
                {showHeader && <VisHeaderView {...{showPreview, visRoot, currMouseState, readout}}/>}
                {showPreview && <VisPreview {...{showPreview, visRoot, currMouseState, readout}}/>}
            </div>
        );
    }
}

VisHeader.propTypes= {
    showHeader : PropTypes.bool,
    showPreview :PropTypes.bool,
};

