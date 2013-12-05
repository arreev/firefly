package edu.caltech.ipac.hydra.server.download;


import edu.caltech.ipac.astro.IpacTableException;
import edu.caltech.ipac.firefly.data.DownloadRequest;
import edu.caltech.ipac.firefly.data.ReqConst;
import edu.caltech.ipac.firefly.data.ServerRequest;
import edu.caltech.ipac.firefly.server.ServerContext;
import edu.caltech.ipac.firefly.server.packagedata.FileGroup;
import edu.caltech.ipac.firefly.server.packagedata.FileInfo;
import edu.caltech.ipac.firefly.server.query.DataAccessException;
import edu.caltech.ipac.firefly.server.query.FileGroupsProcessor;
import edu.caltech.ipac.firefly.server.query.SearchManager;
import edu.caltech.ipac.firefly.server.query.SearchProcessorImpl;
import edu.caltech.ipac.firefly.server.util.Logger;
import edu.caltech.ipac.firefly.server.util.ipactable.DataGroupPart;
import edu.caltech.ipac.firefly.server.util.ipactable.IpacTableParser;
import edu.caltech.ipac.firefly.visualize.VisUtil;
import edu.caltech.ipac.hydra.data.PlanckTOIRequest;
import edu.caltech.ipac.hydra.server.query.PlanckTOIFileRetrieve;
import edu.caltech.ipac.util.StringUtils;
import edu.caltech.ipac.visualize.plot.WorldPt;

import java.io.File;
import java.io.IOException;
import java.util.*;


@SearchProcessorImpl(id = "planckTOIDownload")
public class PlanckTOIGroupsProcessor extends FileGroupsProcessor {

    private static final Logger.LoggerImpl logger = Logger.getLogger();

    
    public List<FileGroup> loadData(ServerRequest request) throws IOException, DataAccessException {
        assert (request instanceof DownloadRequest);
        try {
            return computeFileGroup((DownloadRequest) request);
        } catch (Exception e) {
            e.printStackTrace();
            throw new DataAccessException(e.getMessage());
        }
    }


    private List<FileGroup> computeFileGroup(DownloadRequest request) throws IOException, IpacTableException, DataAccessException {
        // create unique list of filesystem-based and url-based files


        Set<String> zipFiles = new HashSet<String>();

        Collection<Integer> selectedRows = request.getSelectedRows();
        DataGroupPart dgp = new SearchManager().getDataGroup(request.getSearchRequest());

        ArrayList<FileInfo> fiArr = new ArrayList<FileInfo>();
        long fgSize = 0;

        // values = folder or flat
        String zipType = request.getParam("zipType");
        boolean doFolders =  zipType != null && zipType.equalsIgnoreCase("folder");

        String type = request.getSafeParam("type");
        String Size = request.getSafeParam("sradius");
        String optBand = request.getSafeParam("planckfreq");
        String detector = request.getSafeParam("detc100");

        WorldPt pt;
        String pos = null;
        String userTargetWorldPt = request.getParam(ReqConst.USER_TARGET_WORLD_PT);
        if (userTargetWorldPt != null) {
            pt = WorldPt.parse(userTargetWorldPt);
            if (pt != null) {
                pt = VisUtil.convertToJ2000(pt);
                pos = pt.getLon() + "," + pt.getLat();
            }
        } else {
            throw new DataAccessException("No Name or Position found");
        }

        String baseFilename = "planck_toi_search_"+optBand+"-"+detector;
        IpacTableParser.MappedData dgData = IpacTableParser.getData(new File(dgp.getTableDef().getSource()),
                selectedRows, "t_begin", "t_end");

        String baseUrl = PlanckTOIFileRetrieve.getBaseURL(request);

        for (int rowIdx : selectedRows) {
            FileInfo fi = null;
            String fName = String.format("%d",rowIdx);
            String t_begin = String.format("%d",(Long)dgData.get(rowIdx, "t_begin"));
            String t_end = String.format("%d",(Long)dgData.get(rowIdx, "t_end"));
            int estSize = 30000;

            File fileName = new File(baseFilename);

            logger.briefInfo("filename=" +fileName);

            String extName = fName;

            String url = PlanckTOIFileRetrieve.createTOIURLString(baseUrl, pos, type, Size, optBand, detector,t_begin, t_end);
            // strip out filename when using file resolver
            logger.briefInfo("toi search url=" +url);

            // strip out filename when using file resolver
            if (doFolders) {
                int idx = extName.lastIndexOf("/");
                idx = idx < 0 ? 0 : idx;
                extName = extName.substring(0, idx) + "/";
            } else {
                extName = null;
            }

            String TOIFile = baseFilename+fName;

            fi = new FileInfo(url, TOIFile, estSize);

            if (fi != null) {
                fiArr.add(fi);
                fgSize += fi.getSizeInBytes();
            }
        }



        FileGroup fg = new FileGroup(fiArr, null, fgSize, "PlanckTOI Download Files");
        ArrayList<FileGroup> fgArr = new ArrayList<FileGroup>();
        fgArr.add(fg);
        return fgArr;


    }

}


/*
 * THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
 * INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. GOVERNMENT CONTRACT WITH
 * THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION (NASA). THE SOFTWARE
 * IS TECHNOLOGY AND SOFTWARE PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS
 * AND IS PROVIDED AS-IS TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND,
 * INCLUDING ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
 * A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC 2312- 2313)
 * OR FOR ANY PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS,
 * HOWEVER USED.
 *
 * IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA BE LIABLE
 * FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT LIMITED TO, INCIDENTAL
 * OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO
 * PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE
 * ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
 *
 * RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE SOFTWARE
 * AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH AND NASA FOR
 * ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE USE
 * OF THE SOFTWARE.
 */
