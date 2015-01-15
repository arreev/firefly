/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
package edu.caltech.ipac.firefly.core.background;

import edu.caltech.ipac.firefly.data.packagedata.PackagedBundle;
import edu.caltech.ipac.firefly.data.packagedata.PackagedReport;
import edu.caltech.ipac.firefly.data.packagedata.SearchBundle;
import edu.caltech.ipac.util.HandSerialize;
import edu.caltech.ipac.util.StringUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
/**
 * User: roby
 * Date: Dec 15, 2009
 * Time: 10:27:30 AM
 */


/**
 * @author Trey Roby
 */
public class BackgroundReport implements BackgroundPart, Serializable, Iterable<BackgroundPart>, HandSerialize {

    private final static String SPLIT_TOKEN= "--BGReport--";
    private final static String MSG_SPLIT_TOKEN= "--MsgBGReport--";
    private final static String ATT_SPLIT_TOKEN= "--AttBGReport--";
    private final static String PART_SPLIT_TOKEN= "--AttBGReport--";
    private final static String PRO_SPLIT_TOKEN= "--ProBRReport--";

    public enum JobAttributes {Zipped, CanSendEmail, DownloadScript, EmailSent, LongQueue, Unknown, ClientActivated}


    public final static List<BackgroundPart> BLANK_LIST = Collections.emptyList();
    public final static List<BackgroundPart> CANCEL_LIST = dummyList(BackgroundState.CANCELED);
    public final static List<BackgroundPart> WAITING_LIST = dummyList(BackgroundState.WAITING);
    public final static List<BackgroundPart> FAIL_LIST = dummyList(BackgroundState.FAIL);


    public final static String NO_ID = "WARNING:_UNKNOWN_PACKAGE_ID";

    private List<BackgroundPart> _partList = new ArrayList<BackgroundPart>(1);//done
    private BackgroundState _state;//done
    private String _backgroundID;//done
    private boolean _done = false;//done
    private String _dataSource = "";//done
    private ArrayList<String> _messages;//done
    private ArrayList<PackageProgress> _progressList= new ArrayList<PackageProgress>(1);//done
    private Set<JobAttributes> _attributes = new HashSet<JobAttributes>();//done


//======================================================================
//----------------------- Constructors ---------------------------------
//======================================================================

    public BackgroundReport() {
    }

    public BackgroundReport(String packageID, BackgroundState state) {
        this(packageID, Arrays.asList((BackgroundPart) new DefaultBackgroundPart(state)), state);
    }

    public BackgroundReport(String packageID, List<BackgroundPart> parts, BackgroundState state) {
        _backgroundID = packageID;
        _state = state;
        _messages = null;
        if (parts != null) {
            for(BackgroundPart p : parts) addBackgroundPart(p);
        }
        _done = (state == BackgroundState.USER_ABORTED ||
                 state == BackgroundState.CANCELED ||
                 state == BackgroundState.FAIL ||
                 state == BackgroundState.SUCCESS ||
                 state == BackgroundState.UNKNOWN_PACKAGE_ID);
    }


    public static BackgroundReport createUnknownReport() {
        return new BackgroundReport(NO_ID, BLANK_LIST,
                                    BackgroundState.UNKNOWN_PACKAGE_ID);
    }

    public static BackgroundReport createCanceledReport(String packageID) {
        return new BackgroundReport(packageID, CANCEL_LIST,
                                    BackgroundState.CANCELED);
    }

    public static BackgroundReport createWaitingReport(String packageID) {
        return new BackgroundReport(packageID, WAITING_LIST,
                                    BackgroundState.WAITING);
    }

    public static BackgroundReport createFailReport(String packageID, String reason) {
        BackgroundReport report = new BackgroundReport(packageID, FAIL_LIST,
                                                       BackgroundState.FAIL);
        report.addMessage(reason);
        return report;
    }

//======================================================================
//----------------------- Public Methods -------------------------------
//======================================================================


    public Iterator<BackgroundPart> iterator() {
        List<BackgroundPart> l = new ArrayList<BackgroundPart>(_partList.size());
        l.addAll(_partList);
        return l.iterator();
    }

    public void addBackgroundPart(BackgroundPart part) {
        _partList.add(part);
    }


    public BackgroundState getState() {
        return _state;
    }

    public String getID() {
        return _backgroundID;
    }

    public boolean isDone() {
        return _done;
    }

    public int getPartCount() {
        return _partList.size();
    }


    public boolean isFail() {
        return (_state == BackgroundState.FAIL ||
                _state == BackgroundState.USER_ABORTED ||
                _state == BackgroundState.UNKNOWN_PACKAGE_ID ||
                _state == BackgroundState.CANCELED);
    }

    public boolean isSuccess() {
        return (_state == BackgroundState.SUCCESS);
    }


    public BackgroundState getPartState(int i) {
        if (_partList.size() == 0)
            throw new ArrayIndexOutOfBoundsException("Index " + i + " out of bounds. Report has no parts.");
        return _partList.get(i).getState();
    }

    public void setDataSource(String dataSource) {
        _dataSource = dataSource;
    }

    public String getDataSource() {
        return _dataSource;
    }

    public PackageProgress getPartProgress(int i) {
        PackageProgress retval= null;
        if (_progressList.size()>i) {
           retval= _progressList.get(i) ;
        }
        if (retval==null) retval= new PackageProgress();
        return retval;
    }

    public void setProgress(PackageProgress progress) {
        setPartProgress(progress,0);
    }

    public void setPartProgress(PackageProgress progress, int i) {
        if(_progressList.size()<i) { // ensure size
            for(int j=_progressList.size(); (j<=i);j++) {
                _progressList.add(null);
            }
        }
        _progressList.set(i,progress);
    }

    public BackgroundPart get(int idx) {
        return _partList.get(idx);
    }

    public void addMessage(String message) {
        if (_messages == null) {
            _messages = new ArrayList<String>();
        }
        _messages.add(message);
    }

    public int getNumMessages() {
        if (_messages == null)
            return 0;
        else
            return _messages.size();
    }

    public String getMessage(int idx) {
        if (_messages == null) throw new ArrayIndexOutOfBoundsException("No messages.");
        return _messages.get(idx);
    }

    public BackgroundReport cloneWithState(BackgroundState state) {
        BackgroundReport retval = new BackgroundReport(_backgroundID, state);
        retval.copyPartsFrom(this);
        retval.copyMessagesFrom(this);
        retval.copyAttributesFrom(this);
        retval._dataSource = _dataSource;
        return retval;
    }


    public String toString() {
        return "packageID= " + getID() + ", state= " + getState() + ", done= " + isDone();
    }

    public boolean hasFileKey() {
        return false;
    }

    public String getFileKey() {
        return null;
    }

    public boolean hasAttribute(JobAttributes a) {
        return _attributes.contains(a);
    }

    public Set<JobAttributes> getAllAttributes() {
        return _attributes;
    }

    public void addAttribute(JobAttributes a) {
        if (!_attributes.contains(a)) _attributes.add(a);
    }

    //======================================================================
//------------------ Private / Protected Methods -----------------------
//======================================================================

    protected void copyMessagesFrom(BackgroundReport rep) {
        if (rep._messages != null) {
            _messages = new ArrayList<String>(rep._messages.size());
            _messages.addAll(rep._messages);
        }
    }

    protected void copyPartsFrom(BackgroundReport rep) {
        _partList = new ArrayList<BackgroundPart>(rep._partList.size());
        for(BackgroundPart p : rep._partList) addBackgroundPart(p);
    }

    protected void copyAttributesFrom(BackgroundReport rep) {
        _attributes.clear();
        _attributes.addAll(rep._attributes);
    }

    private static List<BackgroundPart> dummyList(BackgroundState state) {
        BackgroundPart part = new DummyBackgroundPart(state);
        return Arrays.asList(part);
    }

    public String serialize() {

        StringBuilder sb= new StringBuilder(350);

        if (this instanceof BackgroundSearchReport) {
            sb.append("BackgroundSearchReport:");
        }
        else if (this instanceof PackagedReport) {
            sb.append("PackagedReport:");
        }
        else {
            sb.append("BackgroundReport:");
        }

        sb.append(_state.toString()).append(SPLIT_TOKEN);
        sb.append(_backgroundID).append(SPLIT_TOKEN);
        sb.append(_dataSource).append(SPLIT_TOKEN);
        if (_messages!=null) {
            sb.append(StringUtils.combineAry(MSG_SPLIT_TOKEN, _messages.toArray(new String[_messages.size()])));
        }
        else {
            sb.append("[]");
        }
        sb.append(SPLIT_TOKEN);
        String[] attStrAry= new String[_attributes.size()];
        int i= 0;
        for(JobAttributes att : _attributes) {
            attStrAry[i++]= att.toString();
        }
        sb.append(StringUtils.combineAry(ATT_SPLIT_TOKEN,attStrAry));
        sb.append(SPLIT_TOKEN);

        String[] partStrAry= new String[_partList.size()];
        String[] progressStrAry= new String[_progressList.size()];
        i= 0;
        for(BackgroundPart part : _partList) {
            partStrAry[i++]= BackgroundPartSerializer.serialize(part);
        }
        sb.append(StringUtils.combineAry(PART_SPLIT_TOKEN,partStrAry));
        sb.append(SPLIT_TOKEN);
        for(PackageProgress p : _progressList) {
            progressStrAry[i++]= p.serialize();
        }
        sb.append(StringUtils.combineAry(PRO_SPLIT_TOKEN,progressStrAry));
        sb.append(SPLIT_TOKEN);


        if (this instanceof BackgroundSearchReport) {
            sb.append(getFileKey()).append(SPLIT_TOKEN);
        }
        else if (this instanceof PackagedReport) {
            sb.append( ((PackagedReport)this).getTotalSizeInByte());
            sb.append(SPLIT_TOKEN);
        }

        return sb.toString();
    }


    public static BackgroundReport parse(String s) {

        BackgroundReport  retval= new BackgroundReport();
        if (s.startsWith("BackgroundSearchReport:")) {
            s= s.substring("BackgroundSearchReport:".length());
            retval= new BackgroundSearchReport(null,null,null,null);
        }
        else if (s.startsWith("CompositeJob:")) {
//            s= s.substring("CompositeJob:".length());
            throw new RuntimeException("not supported.");
        }
        else if (s.startsWith("PackagedReport:")) {
            s= s.substring("PackagedReport:".length());
            retval= new PackagedReport();
        }
        else if (s.startsWith("BackgroundReport:")) {
            s= s.substring("BackgroundReport:".length());
            retval= new BackgroundReport();
        }
        String sAry[]= StringUtils.parseHelper(s, 11, SPLIT_TOKEN);
        int i= 0;
        BackgroundState state= Enum.valueOf(BackgroundState.class, sAry[i++]);
        String backgroundID= StringUtils.checkNull(sAry[i++]);
        String dataSource= StringUtils.checkNull(sAry[i++]);
        List<String> messages= StringUtils.parseStringList(sAry[i++], MSG_SPLIT_TOKEN);
        List<String> attStringList= StringUtils.parseStringList(sAry[i++], ATT_SPLIT_TOKEN);
        List<String> partStringList= StringUtils.parseStringList(sAry[i++],PART_SPLIT_TOKEN);
        List<String> progStringList= StringUtils.parseStringList(sAry[i++],PRO_SPLIT_TOKEN);

        List<BackgroundPart> partList= new ArrayList<BackgroundPart>(partStringList.size());

        for(String pString : partStringList) {
            partList.add(BackgroundPartSerializer.parse(pString));
        }

        if (retval instanceof PackagedReport) {
            long sizeInBytes= StringUtils.getLong(sAry[i++]);
            PackagedBundle pbAry[]= partList.toArray(new PackagedBundle[partList.size()]);
            retval= new PackagedReport(backgroundID,pbAry,sizeInBytes,state);
        }
        else if (retval instanceof BackgroundSearchReport) {
            String fileKey= StringUtils.checkNull(sAry[i++]);
            SearchBundle sb= (SearchBundle)partList.get(0);
            retval= new BackgroundSearchReport(backgroundID,state,sb.getClientRequest(),sb.getServerRequest());
            ((BackgroundSearchReport)retval).setFilePath(fileKey);
        }
        else {
            retval= new BackgroundReport(backgroundID,partList,state);
            retval.setDataSource(dataSource);
            String ps;
            for(int j=0; (j<progStringList.size()); j++) {
                ps= StringUtils.checkNull(progStringList.get(j));
                retval.setPartProgress(PackageProgress.parse(ps),j);
            }
        }


        for(String m : messages) retval.addMessage(m);
        for(String aStr : attStringList)  retval.addAttribute(Enum.valueOf(JobAttributes.class, aStr));

        return retval;

    }

}
