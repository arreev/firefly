package edu.caltech.ipac.firefly.server.query.lsst;


import edu.caltech.ipac.firefly.data.*;
import edu.caltech.ipac.firefly.server.query.*;
import edu.caltech.ipac.firefly.server.util.Logger;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;


/**
 * Created by zhang on 10/12/16.
 * This search processor is searching the MetaData (or Data Definition from DAX database, then save to a
 * IpacTable file.
 */
@SearchProcessorImpl(id = "LSSTImageSearch")
/**
 * Created by zhang on 11/3/16.
 */
public class LSSTImageSearch extends URLFileInfoProcessor {
    private static final Logger.LoggerImpl logger = Logger.getLogger();
    private static String DAX_URL="http://lsst-qserv-dax01.ncsa.illinois.edu:5000/image/v0/";

    /**
     * Implement the abstract method, "getURL"
     * @param sr - ServerRequest
     * @return a URL object
     * @throws MalformedURLException
     */
    @Override
    public URL getURL(ServerRequest sr) throws MalformedURLException {

        try {

            if (sr.getParam("tract") != null) {
                return getURLForDeepCoadd(sr);
            } else {
                return  getURLForCCDs(sr);
            }
        }
        catch (Exception e){
            throw new MalformedURLException(e.getMessage());
        }
    }

    /**
     * This method is using deepCoaddId to search the image.  It worked fine.  Comment it for now.  If it does not needed,
     * it will be deleted.
     * @param request - ServerRequest
     * @return URL Object
     * @throws IOException
     * @throws DataAccessException
     */
    private URL getURLForDeepCoadd(ServerRequest request) throws IOException, DataAccessException {
        String tract = request.getParam("tract");
        String patch = request.getParam("patch");
        String filterName = request.getParam("filterName");

        logger.info("create URL");
         return  new URL(createURLForDeepCoadd(tract, patch, filterName));
    }


    public static String createURLForDeepCoadd(String tract, String patch, String filterName) throws MalformedURLException {
        return  getBaseURL(true)+"tract="+tract+"&patch="+patch+"&filter="+filterName;
    }
    /**
     * This method uses a set of fields to search for image
     * @param request
     * @return
     * @throws IOException
     * @throws DataAccessException
     */
    private URL  getURLForCCDs(ServerRequest request)throws IOException, DataAccessException {
        logger.info("getting the parameters out from the request");
        String run = request.getParam("run");
        String camcol = request.getParam("camcol");
        String field = request.getParam("field");
        String filterName = request.getParam("filterName");
        logger.info("create URL");
        return new URL( createURLForScienceCCD(run, camcol,field, filterName) );

    }

    public static String  createURLForScienceCCD(String run, String camcol, String field, String filterName) throws MalformedURLException {
        return getBaseURL(false)+"run="+run+"&camcol="+camcol+"&field="+field+"&filter="+filterName;
    }

    private static  String getBaseURL(boolean isDeepCoadd){
        if (isDeepCoadd) {
            return DAX_URL+"deepCoadd/ids?";

        }
        else {
            return DAX_URL+"calexp/ids?";
        }
    }

}