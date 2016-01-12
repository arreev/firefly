/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

import React from 'react';
import sCompare from 'react-addons-shallow-compare';
import {visRoot} from '../ImagePlotCntlr.js';
import {flux} from '../../Firefly.js';
import {VisHeaderView} from './VisHeaderView.jsx';
import {currMouseState} from '../VisMouseCntlr.js';



export class VisHeader extends React.Component {
    constructor(props) {
        super(props);
        this.state= {visRoot:visRoot()};
    }

    shouldComponentUpdate(np,ns) { return sCompare(this,np,ns); }

    componentWillUnmount() {
        if (this.removeListener) this.removeListener();
    }


    componentDidMount() {
        this.removeListener= flux.addListener(() => this.storeUpdate());
    }

    storeUpdate() {
        if (visRoot()!=this.state.visRoot || currMouseState() != this.state.currMouseState) {
            this.setState({visRoot:visRoot(), currMouseState:currMouseState()});
        }
    }

    render() {
        var {visRoot}= this.state;
        return (
            <VisHeaderView visRoot={visRoot} currMouseState={currMouseState()}/>
        );
    }
}
