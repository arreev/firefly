/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
package edu.caltech.ipac.util.download;

import edu.caltech.ipac.firefly.data.FileInfo;
import edu.caltech.ipac.firefly.server.util.Logger;
import edu.caltech.ipac.util.Base64;
import edu.caltech.ipac.util.FileUtil;
import edu.caltech.ipac.util.StringUtils;
import edu.caltech.ipac.util.UTCTimeUtil;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;
import java.util.zip.GZIPInputStream;
import java.util.zip.InflaterInputStream;


public class URLDownload {


    public static final String PATTERN_RFC1123 = "EEE, dd MMM yyyy HH:mm:ss zzz";
    /**
     * Date format pattern used to parse HTTP date headers in RFC 1036 format.
     */
    private static final String PATTERN_RFC1036 = "EEEE, dd-MMM-yy HH:mm:ss zzz";

    /**
     * Date format pattern used to parse HTTP date headers in ANSI C
     */
    private static final String PATTERN_ASCTIME = "EEE MMM d HH:mm:ss yyyy";

    public static final String DEFAULT_PATTERNS[] = {PATTERN_ASCTIME, PATTERN_RFC1036, PATTERN_RFC1123};
    private static final int BUFFER_SIZE = FileUtil.BUFFER_SIZE;

    private static final Logger.LoggerImpl _log = Logger.getLogger();


    private static final SimpleDateFormat _browserDateParser;

    static {
        _browserDateParser = new SimpleDateFormat();
        _browserDateParser.setTimeZone(TimeZone.getTimeZone("GMT"));
    }


//=====================================================================
//----------- Public static methods -----------------------------------
//=====================================================================

    /**
     * @param url      the url to get data from
     * @param postData ?
     * @param outfile  write the url data to this file
     * @param dl       listen for progress and cancel if necessary
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFileUsingPost(URL url,
                                                  String postData,
                                                  Map<String, String> cookies,
                                                  Map<String, String> requestHeader,
                                                  File outfile,
                                                  DownloadListener dl, int timeout) throws FailedRequestException {
        try {
            URLConnection c = makeConnection(url, cookies, requestHeader, false);
            //TODO test the timeout for sql query, added by LZ
             c.setConnectTimeout(timeout * 1000);//Sets a specified timeout value, in milliseconds
             c.setReadTimeout(timeout * 1000);
            if (c instanceof HttpURLConnection) {
                return getDataToFile(c, outfile, dl, false, true, true, 0L,postData);
            } else {
                FailedRequestException fe = new FailedRequestException("Can only do post with http: " + url.toString());
                logError(url, postData, fe);
                throw fe;
            }
        } catch (IOException e) {
            logError(url, postData, e);
            throw new FailedRequestException(ResponseMessage.getNetworkCallFailureMessage(e), e);
        }
    }

    public static String getUserFromUrl(String url) {
        try {
            String[] userInfo = getUserInfo(new URL(url));
            return userInfo == null ? null : userInfo[0];
        } catch (MalformedURLException e) {
            // don't need to handle
        }
        return null;
    }

    public static String getPasswordFromUrl(String url) {
        try {
            String[] userInfo = getUserInfo(new URL(url));
            return userInfo == null ? null : userInfo[1];
        } catch (MalformedURLException e) {
            // don't need to handle
        }
        return null;
    }

    private static String[] getUserInfo(URL url) {
        String auth = url.getAuthority();
        if (auth != null && auth.contains("@")) {
            String[] parts = auth.split("@", 2);
            String[] userInfo = parts[0].split(":", 2);
            if (userInfo.length == 2) {
                return userInfo;
            }
        }
        return null;
    }

    public static URLConnection makeConnection(URL url) throws IOException {
        return makeConnection(url, null, null, false);
    }

    public static URLConnection makeConnection(URL url, Map<String, String> cookies) throws IOException {
        return makeConnection(url, cookies, null, false);
    }

    public static URLConnection makeConnection(URL url,
                                               Map<String, String> cookies,
                                               Map<String, String> requestProperties,
                                               boolean setupCompression) throws IOException {
        try {

            URLConnection conn = url.openConnection();


            if (conn instanceof HttpURLConnection) addCookiesToConnection((HttpURLConnection) conn, cookies);
            String[] userInfo = getUserInfo(url);
            if (userInfo != null) {
                String authStringEnc = Base64.encode(userInfo[0] + ":" + userInfo[1]);
                conn.setRequestProperty("Authorization", "Basic " + authStringEnc);
            }
            if (requestProperties != null && requestProperties.size() > 0) {
                for (Map.Entry<String, String> entry : requestProperties.entrySet()) {
                    conn.setRequestProperty(entry.getKey(), entry.getValue());
                }
            }
            if (setupCompression) {
                conn.setRequestProperty("Accept-Encoding", "gzip, deflate");
            }
            return conn;
        } catch (IOException e) {
            logError(url,e);
            throw e;
        }
    }


    private static void addCookiesToConnection(HttpURLConnection conn, Map<String, String> cookies) {
        if (cookies != null) {
            StringBuilder sb = new StringBuilder(200);
            for (Map.Entry<String, String> entry : cookies.entrySet()) {
                if (sb.length() > 0) sb.append("; ");
                sb.append(entry.getKey());
                sb.append("=");
                sb.append(entry.getValue());
            }
            if (sb.length() > 0) {
                conn.setRequestProperty("Cookie", sb.toString());
            }
        }
    }

//======================================================================
//------------------ Public getDataToFile using a URL  -----------------
//======================================================================


    /**
     * Download data from the URL to the file. If this data appears to be compressed then uncompress it first
     *
     * @param url     the url to get data from
     * @param outfile The name of the file to write the data to. uncompress it first
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFile(URL url, File outfile) throws FailedRequestException {
        return getDataToFile(url, outfile, true);
    }

    /**
     * @param url        the url to get data from
     * @param outfile    The name of the file to write the data to.
     * @param uncompress if this data appears to be compressed then uncompress it first
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFile(URL url, File outfile, boolean uncompress) throws FailedRequestException {
        return getDataToFile(url, outfile, null, null, null, false, uncompress, 0);
    }


    /**
     * @param url                  the url to get data from
     * @param outfile              The name of the file to write the data to.  If useSuggestedFilename if false then the
     *                             data is written to this file.  If it is true then only the directory part of this
     *                             file is used and the filename come from the Content-disposition of the URL.
     * @param dl                   listen for progress and cancel if necessary
     * @param cookies              a map of cookies as name value pairs
     * @param useSuggestedFilename if true then use the name from the Content-disposition of the outfile parameter at
     *                             the directory. if false then the outfile parameter specifies the filename
     * @param uncompress           if this data appears to be compressed then uncompress it first
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFile(URL url,
                                           File outfile,
                                           Map<String, String> cookies,
                                           Map<String, String> requestProperties,
                                           DownloadListener dl,
                                           boolean useSuggestedFilename,
                                           boolean uncompress,
                                           long maxFileSize) throws FailedRequestException {
        try {
            URLConnection conn = makeConnection(url, cookies, requestProperties, true);
            return getDataToFile(conn, outfile, dl, useSuggestedFilename, uncompress, true, maxFileSize);
        } catch (IOException e) {
            throw new FailedRequestException(ResponseMessage.getNetworkCallFailureMessage(e), e);
        }
    }

//================================================================================
//------------------ Public getDataToFile using a URLConnection  -----------------
//================================================================================

    /**
     * Download from the URLConnection. If this data appears to be compressed then uncompress it first
     *
     * @param conn    the URLConnection
     * @param outfile write the url data to this file
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFile(URLConnection conn, File outfile) throws FailedRequestException{
        return getDataToFile(conn, outfile, true);
    }


    /**
     * @param conn       the URLConnection
     * @param outfile    write the url data to this file
     * @param uncompress if this data appears to be compressed then uncompress it first
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFile(URLConnection conn, File outfile, boolean uncompress)
            throws FailedRequestException {
        return getDataToFile(conn, outfile, null, false, uncompress, true, 0L,null);
    }

    /**
     * @param conn    the URLConnection
     * @param outfile write the url data to this file
     * @param dl      listen for progress and cancel if necessary
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFile(URLConnection conn, File outfile, DownloadListener dl)
            throws FailedRequestException {
        return getDataToFile(conn, outfile, dl, false, true, true, 0L,null);
    }

    public static FileInfo getDataToFile(URLConnection conn,
                                         File outfile,
                                         DownloadListener dl,
                                         boolean useSuggestedFilename,
                                         boolean uncompress,
                                         boolean onlyIfModified,
                                         long maxFileSize) throws FailedRequestException {
        return getDataToFile(conn, outfile, dl, useSuggestedFilename, uncompress, onlyIfModified, maxFileSize,null);
    }



    /**
     * @param conn                 the URLConnection
     * @param outfile              The name of the file to write the data to.  If useSuggestedFilename if false then the
     *                             data is written to this file.  If it is true then only the directory part of this
     *                             file is used and the filename come from the Content-disposition of the URL.
     * @param dl                   listen for progress and cancel if necessary
     * @param useSuggestedFilename if true then use the name from the Content-disposition of the outfile parameter at
     *                             the directory. if false then the outfile parameter specifies the filename
     * @param uncompress           if this data appears to be compressed then uncompress it first
     * @param onlyIfModified       get the file only if there is a newer version
     * @param maxFileSize          maximum that can be downloaded
     * @param postData             If non-null then send as post data
     * @return an array of FileInfo objects
     * @throws FailedRequestException Any Network Error with simple message, cause will probably be IOException
     */
    public static FileInfo getDataToFile(URLConnection conn,
                                         File outfile,
                                         DownloadListener dl,
                                         boolean useSuggestedFilename,
                                         boolean uncompress,
                                         boolean onlyIfModified,
                                         long maxFileSize,
                                         String postData) throws FailedRequestException {
        try {
            FileInfo outFileData;
            int responseCode= 200;
            Map<String,List<String>> sendHeaders= null;
            try {
                if (conn instanceof HttpURLConnection) {
                    HttpURLConnection httpConn= (HttpURLConnection) conn;
                    if (postData!=null) {
                        pushPostData(httpConn,postData);
                    }
                    if (uncompress) httpConn.setRequestProperty("Accept-Encoding", "gzip, deflate");
                    sendHeaders= httpConn.getRequestProperties();
                    if (onlyIfModified) {
                        outFileData = checkAlreadyDownloaded(httpConn, outfile);
                        if (outFileData != null) return outFileData;
                    }
                }
            } catch (IllegalStateException e) {
                // if I get this exception then the connection was already open and I can't set any more headers
                // If that happens then just go ahead with the download
            }
            //------
            //---From here on the server should be responding
            //------
            OutputStream out;
            DataInputStream in;
            logHeader(postData, conn, sendHeaders);
            if (conn instanceof HttpURLConnection) responseCode= ((HttpURLConnection)conn).getResponseCode();
            validFileSize(conn, maxFileSize);
            File f = useSuggestedFilename ? makeFile(conn, outfile) : outfile;
            String suggested = getSugestedFileName(conn);
            String encodeType = conn.getContentEncoding();
            String contentType = conn.getContentType();
            if (uncompress) {
                if (encodeType != null) {
                    if (encodeType.toLowerCase().endsWith("gzip")) {
                        in = makeGZipInStream(conn);
                    } else if (encodeType.toLowerCase().endsWith("deflate")) {
                        in = new DataInputStream(new InflaterInputStream(makeInStream(conn)));
                    } else {
                        in = makeDataInStream(conn);
                        _log.warn("unrecognized Content-encoding: " + encodeType, "cannot uncompress");
                    }
                } else if (contentType != null && contentType.toLowerCase().endsWith("gzip")) {
                    in = makeGZipInStream(conn);
                } else {
                    in = makeDataInStream(conn);
                }
            } else {
                if (encodeType != null && encodeType.endsWith("gzip")) { //can be x-gzip or gzip
                    if (suggested != null && !suggested.toLowerCase().endsWith(FileUtil.GZ)) {
                        _log.warn("content encoded without accepting it: " + encodeType,
                                          "added gz to file name");
                        suggested = suggested + "." + FileUtil.GZ;
                        if (useSuggestedFilename) f = new File(f.getParent(), suggested);
                    }
                }

                in = makeDataInStream(conn);
            }
            out = makeOutStream(f);


            long start = System.currentTimeMillis();
            netCopy(in, out, conn, maxFileSize, dl);

            long elapse = System.currentTimeMillis() - start;
            outFileData = new FileInfo(f, suggested, responseCode, ResponseMessage.getHttpResponseMessage(responseCode));
            logDownload(outFileData, conn.getURL().toString(), elapse );

            if (responseCode>=300 && responseCode<400) {
                throw new FailedRequestException(ResponseMessage.getHttpResponseMessage(responseCode),
                                                 "Response Code: "+responseCode);
            }

            return outFileData;
        } catch (IOException e) {
            logError(conn.getURL(), e);
            throw new FailedRequestException(ResponseMessage.getNetworkCallFailureMessage(e),e);
        }
    }

    private static void pushPostData(HttpURLConnection conn, String postData) throws IOException {
        conn.setRequestMethod("POST");
        conn.setRequestProperty( "Content-Type", "application/x-www-form-urlencoded" );
        conn.setRequestProperty( "Content-Length", String.valueOf(postData.length()));
        conn.setDoOutput(true);
        OutputStream os = conn.getOutputStream();
        OutputStreamWriter wr = new OutputStreamWriter(os);
        wr.write(postData);
        wr.flush();
        wr.close();
    }

    public static byte[] getDataFromURL(URL url, DownloadListener dl)
            throws FailedRequestException,
            IOException {
        URLConnection conn = makeConnection(url);
        return getDataFromOpenURL(conn, true, dl);
    }

    public static byte[] getDataFromURLUsingPost(URL url,
                                                 String postData,
                                                 DownloadListener dl) throws FailedRequestException,
            IOException {
        try {
            URLConnection conn = makeConnection(url);
            conn.setDoOutput(true);

            OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
            wr.write(postData);
            wr.flush();
            wr.close();
            return getDataFromOpenURL(conn, false, false, dl);
        } catch (IOException e) {
            logError(url, postData, e);
            throw e;
        }
    }

    private static void validFileSize(URLConnection conn, long maxFileSize) throws FailedRequestException {
        long contLen = conn.getContentLength();
        if (maxFileSize > 0 && contLen > 0 && contLen > maxFileSize) {
            throw new FailedRequestException(
                    "File too big to download, " + FileUtil.getSizeAsString(contLen) +
                            ", Max: " + FileUtil.getSizeAsString(maxFileSize),
                    "URL content length header reports content size greater then max size passed as parameter. " +
                            "Content length:  " + contLen + ", maxFileSize: " + maxFileSize);
        }
    }


    /**
     * Check if the file already has been downloaded. This method needs to be call at a very precise time in the life of
     * the URL.  If must be called before any data is retrieve. It will set a header and then retrieve the response. No
     * headers can be set after this method and no data or headers can be retrieved before this method.
     *
     * @param urlConn the connection
     * @param outfile target file download filename
     * @return the FileInfo if the file exist and is not out of date, otherwise null
     * @throws IOException if something goes wrong
     */
    private static FileInfo checkAlreadyDownloaded(HttpURLConnection urlConn, File outfile) throws IOException {
        FileInfo retval = null;
        try {
            if (outfile != null && outfile.canRead() && outfile.length() > 0) {
                urlConn.setIfModifiedSince(outfile.lastModified());
                if (urlConn.getResponseCode() == HttpURLConnection.HTTP_NOT_MODIFIED) {
                    _log.briefInfo(outfile.getName() + ": Not downloading, already have current version");
                    retval = new FileInfo(outfile, getSugestedFileName(urlConn), HttpURLConnection.HTTP_NOT_MODIFIED,
                                     ResponseMessage.getHttpResponseMessage(HttpURLConnection.HTTP_NOT_MODIFIED));
                    retval.putAttribute(FileInfo.FILE_DOWNLOADED,false+"");
                }
            }
        } catch (IllegalStateException e) {
            // if I get this exception then the connection was already open and I can't set any more headers
            // If that happens then just go ahead with the download
        }
        return retval;
    }

    public static byte[] getDataFromOpenURL(URLConnection conn,
                                            DownloadListener dl) throws FailedRequestException,
            IOException {
        return getDataFromOpenURL(conn, false, true, dl);
    }

    public static byte[] getDataFromOpenURL(URLConnection conn,
                                            boolean logHeader,
                                            DownloadListener dl)
            throws FailedRequestException,
            IOException {

        return getDataFromOpenURL(conn, logHeader, true, dl);
    }


    public static byte[] getDataFromOpenURL(URLConnection conn,
                                            boolean logHeader,
                                            boolean logError,
                                            DownloadListener dl)
            throws FailedRequestException,
            IOException {

        try {
            DataInputStream in;
            if (logHeader) logHeader(conn);
            try {
                in = new DataInputStream(new BufferedInputStream(
                        conn.getInputStream()));
            } catch (IOException e) {
                throw new FailedRequestException(
                        "URL not found",
                        "This file probably does not exist on the server", e);
            }
            ByteArrayOutputStream out = new ByteArrayOutputStream(4096);
            netCopy(in, out, conn, dl);
            byte[] retval = out.toByteArray();
            logCompletedDownload(conn.getURL(), retval.length); // always log this
            return retval;
        } catch (IOException e) {
            if (logError) logError(conn.getURL(), e);
            throw e;
        }
    }


    public static String getStringFromOpenURL(URLConnection conn, DownloadListener dl)
            throws FailedRequestException,
            IOException {
        return new String(getDataFromOpenURL(conn, dl));
    }

    public static String getStringFromURL(URL url, DownloadListener dl) throws FailedRequestException,
            IOException {
        return new String(getDataFromURL(url, dl));
    }

    public static String getStringFromURLUsingPost(URL url,
                                                   String postData,
                                                   DownloadListener dl) throws FailedRequestException,
            IOException {
        return new String(getDataFromURLUsingPost(url, postData, dl));
    }

    public static void netCopy(DataInputStream in, OutputStream out, URLConnection conn, DownloadListener dl)
            throws FailedRequestException, IOException {
        netCopy(in, out, conn, 0, dl);
    }

    private static void netCopy(DataInputStream in,
                               OutputStream out,
                               URLConnection conn,
                               long maxSize,
                               DownloadListener dl) throws FailedRequestException,
            IOException {
        Downloader downloader;

        downloader = new Downloader(in, out, conn);
        downloader.setMaxDownloadSize(maxSize);
        downloader.setDownloadListener(dl);
        try {
            downloader.download();
        } finally {
            FileUtil.silentClose(in);
            FileUtil.silentClose(out);
        }
    }


    public static String getSugestedFileName(URLConnection conn) {
        String retval = null;
        if (conn != null) {
            String disposition = conn.getHeaderField("Content-disposition");
            if (disposition != null) {
                String strs[] = disposition.split(";");
                if (strs.length == 2) {
                    String fname[] = strs[1].split("=");
                    if (fname[0].toLowerCase().contains("filename")) {
                        retval = fname[1];
                    }
                }
            }
        }
        return retval;
    }


    private static void logDownload(FileInfo retFile, String urlStr, long elapse) {
        String timeStr = (elapse>0) ? ", time: "+UTCTimeUtil.getHMSFromMills(elapse) : "";
        if (retFile != null) {
            List<String> outList = new ArrayList<>(2);
            outList.add(String.format("Download Complete: %s : %d bytes%s",
                          retFile.getFile().getName(), retFile.getFile().length(), timeStr));
            outList.add(urlStr);
            _log.info(outList.toArray(new String[outList.size()]));
        }
    }


    public static void logError(URL url, Exception e) {
        logError(url, null, e);
    }

    public static void logError(URL url, String postData, Exception e) {
        List<String> strList = new ArrayList<>(6);
        strList.add("----------Network Error-----------");
        if (url != null) {
            strList.add("----------Sending");
            strList.add(url.toString());
        }
        if (postData != null) {
            strList.add(StringUtils.pad(20, "Post Data ") + ": " + postData);
        }
        if (e != null) {
            strList.add(StringUtils.pad(20,"----------Exception "));
            strList.add(e.toString());
        }
        _log.warn(strList.toArray(new String[strList.size()]));
    }


    private static void logHeader(String postData, URLConnection conn, Map<String,List<String>> sendHeaders) {
        StringBuffer workBuff;
        try {
            Set hSet= null;
            List<String> outStr= new ArrayList<>(40);
            if (conn instanceof HttpURLConnection && ((HttpURLConnection)conn).getResponseCode()==-1) {
            }
            else {
                hSet = conn.getHeaderFields().entrySet();
            }
            Map.Entry entry;
            List values;
            int m;
            String key;
            Iterator k;
            if (conn.getURL() != null) {
                outStr.add("----------Sending");
                outStr.add( conn.getURL().toString());
                if (sendHeaders!=null) {
                    for(Map.Entry<String,List<String>> se: sendHeaders.entrySet()) {
                        workBuff = new StringBuffer(100);
                        if (se.getKey() == null) key = "<none>";
                        else key= se.getKey();
                        workBuff.append(StringUtils.pad(20,key));
                        workBuff.append(": ");
                        if (key.equalsIgnoreCase("cookie")) {
                            workBuff.append("not shown");
                        }
                        else {
                            workBuff.append(se.getValue());
                        }
                        outStr.add(workBuff.toString());
                    }

                }
            }
            if (postData != null) {
                outStr.add(StringUtils.pad(20,"Post Data ") + ": " + postData);
            }
            if (conn instanceof HttpURLConnection) {
                int responseCode= ((HttpURLConnection)conn).getResponseCode();
                outStr.add("----------Received Headers, response status code: " + responseCode);
            }
            else {
                outStr.add("----------Received Headers");
            }
            if (hSet!=null) {
                for (Iterator j = hSet.iterator(); (j.hasNext()); ) {
                    workBuff = new StringBuffer(100);
                    entry = (Map.Entry) j.next();
                    key = (String) entry.getKey();
                    if (key == null) key = "<none>";
                    workBuff.append(StringUtils.pad(20,key));
                    workBuff.append(": ");
                    values = (List) entry.getValue();
                    for (m = 0, k = values.iterator(); (k.hasNext()); m++) {
                        if (m > 0) workBuff.append("; ");
                        workBuff.append(k.next().toString());
                    }
                    outStr.add(workBuff.toString());
                }
            }
            else {
                outStr.add("No headers or status received, invalid http response, using work around");
            }
            _log.info(outStr.toArray(new String[outStr.size()]));
        } catch (Exception e) {
            _log.info(e.getMessage() + ":" + " url=" + (conn.getURL()!=null ? conn.getURL().toString() : "none"));
        }
    }

    public static void logHeader(URLConnection conn) {
        logHeader(null, conn, null);
    }

    private static void logCompletedDownload(URL url, long size) {
        logCompletedDownload(url, size, 0);
    }

    private static void logCompletedDownload(URL url, long size, long elapse) {

        String timeStr = (elapse>0) ? ", time: "+UTCTimeUtil.getHMSFromMills(elapse) : "";
        String s = String.format("Download Complete- %d bytes%s", size, timeStr);
        String urlString = null;
        if (url != null) urlString = url.toString();
        _log.info(s, urlString);
    }


//======================================================================
//------------------ Private / Protected / Methods -----------------------
//======================================================================


    private static OutputStream makeOutStream(File f) throws IOException {
        return new BufferedOutputStream(new FileOutputStream(f), BUFFER_SIZE);
    }

    private static DataInputStream makeGZipInStream(URLConnection conn) throws IOException {
        return new DataInputStream(new GZIPInputStream(conn.getInputStream(), BUFFER_SIZE));
    }

    private static DataInputStream makeDataInStream(URLConnection conn) throws IOException {
        if (conn instanceof HttpURLConnection && ((HttpURLConnection)conn).getResponseCode()==-1) {
            throw new IOException("Http Response Code is -1, invalid http protocol, " +
                                          "probably no status line in response headers");
        }
        else {
            try {
                return new DataInputStream(makeInStream(conn));
            }
            catch (IOException e) {
                if (conn instanceof HttpURLConnection) {
                    return new DataInputStream(makeErrStream((HttpURLConnection) conn));
                }
                else {
                    throw e;
                }
            }
        }
    }

    private static InputStream makeInStream(URLConnection conn) throws IOException {
        return new BufferedInputStream(conn.getInputStream(), BUFFER_SIZE);
    }

    private static InputStream makeErrStream(HttpURLConnection conn) {
        return new BufferedInputStream(conn.getErrorStream(), BUFFER_SIZE);
    }

    private static File makeFile(URLConnection conn, File outfile) {
        if (outfile == null) outfile = new File(".", "out.dat");
        File retval = outfile;
        String fStr = getSugestedFileName(conn);
        if (fStr != null) {
            retval = new File(outfile.getParentFile(), fStr);
        }
        return retval;
    }

}
