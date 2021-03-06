/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

/**
 *  2/20/2018 LZ
 *  IRSA-663
 *
 *  Display meaningful downloaded file names
 *  Remove redundant codes
 *  Fixed the bug in DownloadDialog.jsx
 *
 */
import React, {PureComponent} from 'react';
import PropTypes from 'prop-types';
import {pick, get, isEmpty, cloneDeep} from 'lodash';
import SplitPane from 'react-split-pane';
import {flux} from '../../Firefly.js';
import {LO_VIEW, getLayouInfo, dispatchUpdateLayoutInfo} from '../../core/LayoutCntlr.js';
import {TablesContainer} from '../../tables/ui/TablesContainer.jsx';
import {ChartsContainer} from '../../charts/ui/ChartsContainer.jsx';
import {VisToolbar} from '../../visualize/ui/VisToolbar.jsx';
import {LcImageViewerContainer} from './LcImageViewerContainer.jsx';
import {SplitContent} from '../../ui/panel/DockLayoutPanel.jsx';
import {LC, getViewerGroupKey, updateLayoutDisplay} from './LcManager.js';
import FieldGroupUtils from '../../fieldGroup/FieldGroupUtils.js';
import {ValidationField} from '../../ui/ValidationField.jsx';
import {LcImageToolbar} from './LcImageToolbar.jsx';
import {DownloadOptionPanel, DownloadButton} from '../../ui/DownloadDialog.jsx';
import {showInfoPopup} from '../../ui/PopupUtil.jsx';
import CompleteButton from '../../ui/CompleteButton.jsx';
import {HelpIcon} from '../../ui/HelpIcon.jsx';
import {getTblById, doFetchTable, isTblDataAvail} from '../../tables/TableUtil.js';
import {MAX_ROW} from '../../tables/TableRequestUtil.js';
import {dispatchMultiValueChange, dispatchRestoreDefaults}  from '../../fieldGroup/FieldGroupCntlr.js';
import {logError} from '../../util/WebUtil.js';
import {getConverter, getMissionName,  DL_DATA_TAG} from './LcConverterFactory.js';
import {convertAngle} from '../../visualize/VisUtil.js';

const resultItems = ['title', 'mode', 'showTables', 'showImages', 'showXyPlots', 'searchDesc', 'images',
    LC.MISSION_DATA, LC.GENERAL_DATA, 'periodState'];


export class LcResult extends PureComponent {

    constructor(props) {
        super(props);

        this.state = Object.assign({}, pick(getLayouInfo(), resultItems));
    }

    componentDidMount() {
        this.iAmMounted = true;
        this.removeListener = flux.addListener(() => this.storeUpdate());
    }

    componentWillUnmount() {
        this.iAmMounted = false;
        this.removeListener && this.removeListener();
    }

    storeUpdate() {
        if (this.iAmMounted) {
            const nextState = pick(getLayouInfo(), resultItems);
            this.setState(nextState);
        }
    }

    render() {
        const {title, mode, showTables, showImages, showXyPlots, searchDesc, images,
            missionEntries, generalEntries, periodState} = this.state;
        var {expanded, standard} = mode || {};
        const content = {};
        var visToolbar;
        if (showImages) {
            visToolbar = <VisToolbar key='res-vis-tb'/>;
            content.imagePlot = (<LcImageViewerContainer key='res-images'
                                                         viewerId={LC.IMG_VIEWER_ID}
                                                         closeable={true}
                                                         forceRowSize={1}
                                                         imageExpandedMode={expanded===LO_VIEW.images}
                                                         Toolbar={LcImageToolbar}
                {...images}  />);
        }
        if (showXyPlots) {
            content.xyPlot = (<ChartsContainer key='res-charts'
                                               tbl_group='main'
                                               closeable={true}
                                               expandedMode={expanded===LO_VIEW.xyPlots}
            />);
        }
        if (showTables) {
            content.tables = (<TablesContainer key='res-tables'
                                               mode='both'
                                               closeable={true}
                                               expandedMode={expanded===LO_VIEW.tables}
                                               tableOptions={{help_id:'main1TSV.table'}}/>);
        }


        content.settingBox = (<SettingBox generalEntries={generalEntries}
                                          missionEntries={missionEntries}
                                          periodState={periodState}  />);


        expanded = LO_VIEW.get(expanded) || LO_VIEW.none;
        const expandedProps = {expanded, ...content};
        const standardProps = {visToolbar, title, searchDesc, standard, ...content};

        return (
            expanded === LO_VIEW.none
                ? <StandardView key='res-std-view' {...standardProps} />
                : <ExpandedView key='res-exp-view' {...expandedProps} />
        );
    }
}

const ExpandedView = ({expanded, imagePlot, xyPlot, tables}) => {
    const view = expanded === LO_VIEW.tables ? tables
        : expanded === LO_VIEW.xyPlots ? xyPlot
        : imagePlot;
    return (
        <div style={{ flex: 'auto', display: 'flex', flexFlow: 'column', overflow: 'hidden'}}>{view}</div>
    );
};

const buttonW = 650;

const StandardView = ({visToolbar, title, searchDesc, imagePlot, xyPlot, tables, settingBox}) => {

    const converterId = get(settingBox, ['props', 'missionEntries', LC.META_MISSION]);
    const convertData = getConverter(converterId);
    const cutoutSize = get(convertData, 'noImageCutout') ? undefined : get(settingBox, 'props.generalEntries.cutoutSize', '5');
    const mission = getMissionName(converterId) || 'Mission';
    const showImages = isEmpty(imagePlot);

    // convert the default Cutout size in arcmin to deg for WebPlotRequest, expected to be string in download panel
    const cutoutSizeInDeg = (convertAngle('arcmin','deg', cutoutSize)).toString();
    const currentTime =  (new Date()).toLocaleString('en-US', { hour12: false });
    const style = {width: 223};
    const defaultOptPanel = (m, c) => {
        return (
            <DownloadButton>
                <DownloadOptionPanel
                    groupKey = {mission}
                    dataTag = {DL_DATA_TAG}
                    cutoutSize={c}
                    title={'Image Download Options'}
                    style = {{width: 400}}
                    dlParams={{
                        MaxBundleSize: 200 * 1024 * 1024,    // set it to 200mb to make it easier to test multi-parts download.  each wise image is ~64mb
                        FilePrefix: `${m}_Files`,
                        BaseFileName: `${m}_Files`,
                        DataSource: `${m} images`,
                        FileGroupProcessor: 'LightCurveFileGroupsProcessor'
                    }}>
                    <ValidationField
                        style={style}
                        initialState={{
                            value: `${m}_Files: ${currentTime}`,
                            label: `${m}:`
                        }}
                        fieldKey='Title'
                        labelWidth={110}/>
                </DownloadOptionPanel>
            </DownloadButton>
        );
    };

    const downloaderOptPanel = convertData.downloadOptions || defaultOptPanel;

    let tsView = (err) => {

        if (!err) {
            return (
                <SplitPane split='horizontal' maxSize={-20} minSize={20} defaultSize={'60%'}>
                    <SplitPane split='vertical' maxSize={-20} minSize={20} defaultSize={buttonW}>
                        <SplitContent>
                            <div style={{display: 'flex', flexDirection: 'column', height: '100%'}}>
                                <div className='settingBox'>{settingBox}</div>
                                <div style={{flexGrow: 1, position: 'relative'}}>
                                    <div className='abs_full'>{tables}</div>
                                </div>
                            </div>
                        </SplitContent>
                        <SplitContent>{xyPlot}</SplitContent>
                    </SplitPane>
                    <SplitContent>{imagePlot}</SplitContent>
                </SplitPane>
            );
        } else {
            return (
                <SplitPane split='vertical' maxSize={-20} minSize={20} defaultSize={buttonW}>
                    <SplitContent>
                        <div style={{display: 'flex', flexDirection: 'column', height: '100%'}}>
                            <div className='settingBox'>{settingBox}</div>
                            <div style={{flexGrow: 1, position: 'relative'}}>
                                <div className='abs_full'>{tables}</div>
                            </div>
                        </div>
                    </SplitContent>
                    <SplitContent>{xyPlot}</SplitContent>
                </SplitPane>
            );
        }

    };
    return (
        <div style={{display: 'flex', flexDirection: 'column', flexGrow: 1, position: 'relative'}}>
            { visToolbar &&
            <div style={{display: 'inline-flex', justifyContent: 'space-between', alignItems: 'center'}}>
                <div>{visToolbar}</div>
                <div>
                    {downloaderOptPanel(mission, cutoutSizeInDeg)}
                </div>
            </div>
            }
            {searchDesc}
            {title && <h2 style={{textAlign: 'center'}}>{title}</h2>}
            <div style={{flexGrow: 1, position: 'relative'}}>
                <div style={{position: 'absolute', top: 0, right: 0, bottom: 0, left: 0}}>
                    {tsView(showImages)}
                </div>
            </div>
        </div>
    );
};


class SettingBox extends PureComponent {
    constructor(props) {
        super(props);
    }

    render() {
        var {generalEntries, missionEntries, periodState} = this.props;

        if (isEmpty(generalEntries) || isEmpty(missionEntries)) return false;

        const converterId = get(missionEntries, LC.META_MISSION);
        const converterData = converterId && getConverter(converterId);
        if (!converterId || !converterData) {
            return null;
        }
        const {MissionOptions} = converterData;

        const groupKey = getViewerGroupKey(missionEntries);
        return (
            <div>
                <div style={{position: 'relative', display: 'inline-flex', justifyContent: 'space-between', width: '100%'}}>
                  <div style={{alignSelf: 'flex-end'}}>
                      <MissionOptions {...{missionEntries, generalEntries}}/>
                  </div>

                  <div style={{display: 'flex', flexDirection: 'row-reverse'}}>
                      <HelpIcon helpId={'main1TSV.settings'}/>
                  </div>

              </div>
              <div >
                <CompleteButton
                     style={{ width:'1200px'}}
                     groupKey={groupKey}
                     onSuccess={setViewerSuccess(periodState)}
                     onFail={setViewerFail()}
                     text={'Period Finder...'}
                 />
              </div>
          </div>
        );
    }
}

SettingBox.propTypes = {
    generalEntries: PropTypes.object,
    missionEntries: PropTypes.object,
    periodState: PropTypes.string
};

/**
 * @summary callback to go to period finding page
 * @param {string} periodState
 * @returns {Function}
 */
function setViewerSuccess(periodState) {
    return (request) => {
        updateFullRawTable(()=>updateLayoutDisplay(periodState));
    };
}

function setViewerFail() {
    return (request) => {
        return showInfoPopup('Parameter setting error');
    };
}

function updateFullRawTable(callback) {
    var layoutInfo = getLayouInfo();
    const tableItems = ['tableData', 'tableMeta'];

    // fullRawTable for the derivation of other table, like phase folded table
    var setTableData = (tbl) => {
        const fullRawTable = pick(tbl, tableItems);

        // find tzero, tzeroMax, period min, period max from table data
        var {columns, data} = fullRawTable.tableData;
        var tIdx = columns.findIndex((col) => (col.name === get(layoutInfo, [LC.MISSION_DATA, LC.META_TIME_CNAME])));
        var arr = data.reduce((prev, e)=> {
            prev.push(parseFloat(e[tIdx]));
            return prev;
        }, []);

        var [tzero, tzeroMax] = arr.length > 0 ? [Math.min(...arr), Math.max(...arr)] : [0.0, 0.0];
        var max = 365;
        var min = Math.pow(10, -3);   // 0.001

        // var period = get( FieldGroupUtils.getGroupFields(LC.FG_PERIOD_FINDER), ['period', 'value'], '');

        const period = '';
        var fields = FieldGroupUtils.getGroupFields(LC.FG_PERIOD_FINDER);
        var initState;

        if (fields) {      // fields already exists and new table is loaded
            initState = [
                {fieldKey: 'time', value: get(layoutInfo, [LC.MISSION_DATA, LC.META_TIME_CNAME])},
                {fieldKey: 'flux', value: get(layoutInfo, [LC.MISSION_DATA, LC.META_FLUX_CNAME])},
                {fieldKey: 'periodMin', value: `${min}`},
                {fieldKey: 'periodMax', value: `${max}`},
                {fieldKey: 'period', value: `${period}`},
                {fieldKey: 'tzero', value: `${tzero}`},
                {fieldKey: 'tzeroMax', value: `${tzeroMax}`}];

            dispatchMultiValueChange(LC.FG_PERIOD_FINDER, initState);
        }
        fields = FieldGroupUtils.getGroupFields(LC.FG_PERIODOGRAM_FINDER);
        if (fields) {
            dispatchRestoreDefaults(LC.FG_PERIODOGRAM_FINDER);
        }

        //IRSA-464: the smartMerge does not replace array component.  It updates the array. Therefore, setting it
        //null first to force the new fullRawTable is updated to the layout
        dispatchUpdateLayoutInfo(Object.assign({}, layoutInfo, {
            fullRawTable:null,
            periodRange: {min, max, tzero, tzeroMax, period}
        }));
        //add the fullRawtable into the layout
        dispatchUpdateLayoutInfo(Object.assign({}, layoutInfo, {
            fullRawTable,
            periodRange: {min, max, tzero, tzeroMax, period}
        }));


        callback && callback();
    };

    var rawTable = getTblById(LC.RAW_TABLE);
    if (layoutInfo.fullRawTable && layoutInfo.fullRawTable.totalRows===rawTable) {
        callback && callback();
    } else {

        if (isTblDataAvail(0, rawTable.totalRows, rawTable)) {
            setTableData(rawTable);
        } else {
            var req = Object.assign(cloneDeep(rawTable.request), {pageSize: MAX_ROW});

            doFetchTable(req).then(
                (tableModel) => {
                    setTableData(tableModel);
                }
            ).catch(
                (reason) => {
                    logError(`Failed to get full raw table: ${reason}`, reason);
                }
            );
       }
    }
}
