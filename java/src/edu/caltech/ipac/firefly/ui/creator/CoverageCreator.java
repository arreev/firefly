package edu.caltech.ipac.firefly.ui.creator;

import edu.caltech.ipac.firefly.ui.previews.CoveragePreview;
import edu.caltech.ipac.firefly.ui.previews.SimpleCoverageData;
import edu.caltech.ipac.firefly.ui.table.TablePreview;
import edu.caltech.ipac.firefly.visualize.WebPlotRequest;
import edu.caltech.ipac.firefly.visualize.ZoomType;
import edu.caltech.ipac.firefly.visualize.draw.DrawSymbol;

import java.util.List;
import java.util.Map;
/**
 * User: roby
 * Date: Apr 15, 2010
 * Time: 11:41:20 AM
 */



/**
 * @author Trey Roby
 */
public class CoverageCreator implements ObsResultCreator {

    public static final String QUERY_ID= "QUERY_ID";

    public TablePreview create(Map<String, String> params) {

        List<String> sources= DataViewCreator.getListParam(params,QUERY_ID);



        ZoomType hint= ZoomType.SMART;
        if (params.containsKey(WebPlotRequest.ZOOM_TYPE)) {
            try {
                hint= Enum.valueOf(ZoomType.class, params.get(WebPlotRequest.ZOOM_TYPE));
            } catch (Exception e) {
                hint= ZoomType.SMART;
            }
        }
        SimpleCoverageData covData= new SimpleCoverageData(sources,hint);
        covData.setEventWorkerList(DataViewCreator.getListParam(params,"EVENT_WORKER_ID"));

        boolean enableDefColumns=DataViewCreator.getBooleanParam(params,CommonParams.ENABLE_DEFAULT_COLUMNS);
        if (enableDefColumns) {
           covData.enableDefaultColumns();
        }
        else {
            List<String> cenC= DataViewCreator.getListParam(params,CommonParams.CENTER_COLUMNS);
            if (cenC!=null && cenC.size()==2) {
                List<String> corC= DataViewCreator.getListParam(params,CommonParams.CORNER_COLUMNS);
                if (corC!=null && corC.size()==8) {
                    covData.initFallbackCol(cenC.get(0), cenC.get(1),
                                            corC.get(0), corC.get(1),
                                            corC.get(2), corC.get(3),
                                            corC.get(4), corC.get(5),
                                            corC.get(6), corC.get(7) );
                }
                else {
                    covData.initFallbackCol(cenC.get(0),cenC.get(1));
                }
            }

        }

        String group= params.get(CommonParams.PLOT_GROUP);
        if (group!=null)  covData.setGroup(group);




        if (params.containsKey(CommonParams.MIN_SIZE)) {
            String s[]= params.get(CommonParams.MIN_SIZE).split("x",2);
            if (s.length==2) {
                try {
                    int minWidth= Integer.parseInt(s[0]);
                    int minHeight= Integer.parseInt(s[1]);
                    covData.setMinSize(minWidth,minHeight);
                } catch (NumberFormatException e) {
                    // do nothing
                }
            }
            params.remove(CommonParams.MIN_SIZE);
        }


        boolean multi= DataViewCreator.getBooleanParam(params,CommonParams.MULTI_COVERAGE);
        covData.setMultiCoverage(multi);

        boolean catalogsAsOverlays= DataViewCreator.getBooleanParam(params,CommonParams.CATALOGS_AS_OVERLAYS, true);
        covData.setTreatCatalogsAsOverlays(catalogsAsOverlays);



        if (params.containsKey(CommonParams.SHAPE)) {
            List<String> sList= DataViewCreator.getListParam(params,CommonParams.SHAPE);
            for(String s : sList) {
                String sAry[]= s.split("=");
                if (sAry.length==2) {
                   covData.setShape(sAry[0], DrawSymbol.getSymbol(sAry[1]));
                }
            }
        }

        if (params.containsKey(CommonParams.COLOR)) {
            List<String> sList= DataViewCreator.getListParam(params,CommonParams.COLOR);
            for(String s : sList) {
                String sAry[]= s.split("=");
                if (sAry.length==2) {
                    covData.setColor(sAry[0], sAry[1]);
                }
            }
        }

        boolean blank= DataViewCreator.getBooleanParam(params,CommonParams.BLANK);
        covData.setUseBlankPlot(blank);


        return new CoveragePreview(covData);



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
