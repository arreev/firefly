/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */


package edu.caltech.ipac.firefly.server.visualize.hips;

import edu.caltech.ipac.firefly.data.ServerParams;
import edu.caltech.ipac.firefly.server.util.Logger;
import edu.caltech.ipac.util.DataGroup;
import edu.caltech.ipac.util.DataObject;
import edu.caltech.ipac.util.DataType;
import edu.caltech.ipac.firefly.server.query.SearchProcessorImpl;
import edu.caltech.ipac.firefly.server.query.ParamDoc;
import edu.caltech.ipac.firefly.data.TableServerRequest;
import edu.caltech.ipac.firefly.server.query.DataAccessException;
import edu.caltech.ipac.firefly.server.query.EmbeddedDbProcessor;
import edu.caltech.ipac.firefly.server.db.EmbeddedDbUtil;
import edu.caltech.ipac.firefly.data.FileInfo;
import edu.caltech.ipac.firefly.server.db.DbAdapter;
import edu.caltech.ipac.firefly.server.visualize.hips.HiPSMasterListEntry.PARAMS;

import java.io.File;
import java.util.*;
import java.io.IOException;

import static edu.caltech.ipac.util.IpacTableUtil.*;

/**
 * @author Cindy Wang
 */
@SearchProcessorImpl(id = "HiPSSearch", params =
        {@ParamDoc(name = ServerParams.HIPS_DATATYPES, desc = "types of HiPS data to search"),
         @ParamDoc(name = ServerParams.HIPS_SOURCES, desc = "HiPS sources"),
         @ParamDoc(name = ServerParams.SORT_ORDER, desc = "HiPS order, source based")
        })
public class HiPSMasterList extends EmbeddedDbProcessor {
    public static String[] HiPSDataType = new String[]{ServerParams.IMAGE,
                                                       ServerParams.CUBE,
                                                       ServerParams.CATALOG};

    private static Map<String, HiPSMasterListSourceType> sources= new HashMap<>();
    public static String[] defaultSourceOrder = new String[]{ServerParams.IRSA, ServerParams.LSST, ServerParams.CDS,
                                                            ServerParams.EXTERNAL};

    static {
        sources.put(ServerParams.IRSA, new IrsaHiPSListSource());
        sources.put(ServerParams.CDS, new CDSHiPSListSource());
        sources.put(ServerParams.EXTERNAL, new ExternalHiPSListSource());
        sources.put(ServerParams.LSST, new LsstHiPSListSource());
    }

    private static final Logger.LoggerImpl _log= Logger.getLogger();
    private static final String errMsg = "HiPS Map search: no HiPS maps found";

    public FileInfo ingestDataIntoDb(TableServerRequest request, File dbFile) throws DataAccessException {
        String hipsSources = request.getParam(ServerParams.HIPS_SOURCES);
        String hipsDataTypes = request.getParam(ServerParams.HIPS_DATATYPES);
        String hipsMergePriority = request.getParam(ServerParams.HIPS_MERGE_PRIORITY);
        String workingSources[] = (hipsSources != null) ? hipsSources.split(",") : null;
        String workingTypes[] = (hipsDataTypes != null) ? hipsDataTypes.split(",") : null;
        String prioritySources[] = (hipsMergePriority != null) ? hipsMergePriority.split(",") : null;
        List<HiPSMasterListEntry> allSourceData = new ArrayList<>();

        DbAdapter dbAdapter = DbAdapter.getAdapter(request);


        if (workingSources == null || workingSources.length == 0 ||
                (workingSources.length == 1 && workingSources[0].equalsIgnoreCase(ServerParams.ALL))) {
            workingSources = defaultSourceOrder.clone();
        }
        if (workingTypes == null || workingTypes.length == 0 ||
                (workingTypes.length == 1 && workingTypes[0].equalsIgnoreCase(ServerParams.ALL))) {
            workingTypes = new String[]{ServerParams.IMAGE, ServerParams.CUBE};
        }

        try {

            for (String source : workingSources) {
                HiPSMasterListSourceType hipsls = sources.get(source);

                if (hipsls != null) {
                    List<HiPSMasterListEntry> hipsL = hipsls.getHiPSListData(workingTypes, source);
                    if (hipsL != null) {
                        allSourceData.addAll(hipsL);
                    }
                }
            }

            if (allSourceData.size() == 0) {
                throw new IOException(errMsg);
            }

            if (workingSources.length > 1 && prioritySources.length != 0) {
                allSourceData = mergeData(Arrays.asList(prioritySources), allSourceData);
            }

            DataGroup dg = createTableDataFromListEntry(allSourceData);
            dg.shrinkToFitData();

            setupMeta(dg, (workingSources.length > 1));

            return EmbeddedDbUtil.ingestDataGroup(dbFile, dg, dbAdapter, "data");
        } catch (Exception e) {
            _log.warn(e.getMessage());
            throw new DataAccessException(errMsg);
        }
    }

    private List<HiPSMasterListEntry> mergeData(List<String> prioritySources, List<HiPSMasterListEntry> allSourceData) {
        int totalS = allSourceData.size();
        HiPSMasterListEntry[] dataAry = new HiPSMasterListEntry[totalS];
        dataAry = allSourceData.toArray(dataAry);

        for (int i = 0; i < totalS; i++) {
            String ivo1 = dataAry[i].getMapInfo().get(PARAMS.IVOID.getKey());
            if (ivo1 == null) continue;

            String src1 = dataAry[i].getMapInfo().get(PARAMS.SOURCE.getKey());
            int n;

            // skip the entries from 'src1'
            for (n = i+1; n < totalS; n++) {
                if (!dataAry[n].getMapInfo().get(PARAMS.SOURCE.getKey()).equals(src1)) {
                    break;
                }
            }
            if (n >= totalS) break;

            // move to entries not from 'src1'
            for (int j = n; j < totalS; j++) {
                String ivo2 = dataAry[j].getMapInfo().get(PARAMS.IVOID.getKey());

                if (ivo2 == null || !ivo1.equals(ivo2)) continue;
                String src2 = dataAry[j].getMapInfo().get(PARAMS.SOURCE.getKey());
                String hsrc = getHigherPriority(src1, src2, prioritySources);

                if (hsrc == null) continue;
                if (hsrc.equals(src2)) {                           // skip to the one next to index i
                    dataAry[i].set(PARAMS.IVOID.getKey(), null);
                    break;
                } else {
                    dataAry[j].set(PARAMS.IVOID.getKey(), null);
                    for (int k = j+1; k < totalS; j = k, k++) {    // skip all entries from src2
                        if (!dataAry[k].getMapInfo().get(PARAMS.SOURCE.getKey()).equals(src2)) {
                            break;
                        }
                    }
                }
            }
        }

        List<HiPSMasterListEntry> newDataList = new ArrayList<>();
        for (HiPSMasterListEntry oneEntry : allSourceData) {
            if (oneEntry.getMapInfo().get(PARAMS.IVOID.getKey()) != null) {
                newDataList.add(oneEntry);
            }
        }

        return newDataList;
    }

    private String getHigherPriority(String src1, String src2, List<String> prioritySources) {
        if (src1.equals(src2)) return null;

        int idx1 = prioritySources.indexOf(src1);
        int idx2 = prioritySources.indexOf(src2);

        if (idx2 > idx1) {
            return idx1 == -1 ? src2 : src1;
        } else if (idx2 < idx1) {
            return idx2 == -1 ? src1 : src2;
        } else {
            return null;
        }
    }


    private void setupMeta(DataGroup dg, boolean bMulti) {
        int    sWidth = 30;

        for (DataType colDT : dg.getDataDefinitions()) {
            String colName = colDT.getKeyName();

            if (colDT.getDataType() != String.class) continue;

            int crtWidth = colDT.getFormatInfo().getWidth();
            if (crtWidth > sWidth && !colName.equals(PARAMS.PROPERTIES.getKey())) {
                dg.addAttribute(makeAttribKey(WIDTH_TAG, colName), Integer.toString(sWidth));
            }

            if ((!bMulti && colName.equals(PARAMS.SOURCE.getKey())) || colName.equals(PARAMS.URL.getKey())) {
                dg.addAttribute(makeAttribKey(VISI_TAG, colName), "hidden");
            }
        }
    }

    private static DataGroup createTableDataFromListEntry(List<HiPSMasterListEntry> hipsMaps) {
        List<DataType> cols = HiPSMasterListEntry.getHiPSEntryColumns();
        DataGroup dg = new DataGroup("HiPS Maps", cols);

        for (HiPSMasterListEntry entry : hipsMaps) {
            DataObject row = new DataObject(dg);
            Map<String, String> mapInfo = entry.getMapInfo();

            for (DataType col : cols) {
                String val = mapInfo.get(col.getKeyName());

                row.setDataElement(col, col.convertStringToData(val));
            }
            dg.add(row);
        }
        return dg;
    }
}


