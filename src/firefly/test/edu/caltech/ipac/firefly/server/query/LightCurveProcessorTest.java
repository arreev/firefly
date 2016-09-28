/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
package edu.caltech.ipac.firefly.server.query;

import edu.caltech.ipac.astro.IpacTableException;
import edu.caltech.ipac.astro.IpacTableReader;
import edu.caltech.ipac.firefly.server.query.lc.IrsaLightCurveHandler;
import edu.caltech.ipac.firefly.server.query.lc.LightCurveHandler;
import edu.caltech.ipac.firefly.server.query.lc.PeriodogramAPIRequest;
import edu.caltech.ipac.util.DataGroup;
import edu.caltech.ipac.util.DataObject;
import edu.caltech.ipac.util.DataType;
import edu.caltech.ipac.util.download.FailedRequestException;
import edu.caltech.ipac.util.download.URLDownload;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.List;

/**
 * Created by ejoliet on 8/23/16.
 */
public class LightCurveProcessorTest {


    private static PeriodogramAPIRequestTest req;
    private static File rawTable;


    @BeforeClass
    public static void setUp() {
        req = new PeriodogramAPIRequestTest();
        try {
            rawTable = File.createTempFile("phasefolded-temp-", ".tbl", new File("."));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testGetPeriodogram() {

        boolean deleteOnExit = true;
        LightCurveHandler t = new IrsaLightCurveHandler() {
            @Override
            protected URL buildUrl(PeriodogramAPIRequest req) throws MalformedURLException {
                return new URL(req.getResultTable());
            }

            @Override
            protected File makeResultTempFile(RESULT_TABLES_IDX resultTable) throws IOException {
                File tmpFile = File.createTempFile("period-test-", ".tbl", new File("."));
                if (deleteOnExit) {
                    tmpFile.deleteOnExit();
                }
                return tmpFile;
            }

            @Override
            protected File makeApiResultTempFile() throws IOException {
                File tmpFile = File.createTempFile("votable-test-", ".xml", new File("."));
                if (deleteOnExit) {
                    tmpFile.deleteOnExit();
                }
                return tmpFile;
            }
        };

        File p = t.getPeriodogramTable(req);

        try {
            DataGroup inDataGroup = IpacTableReader.readIpacTable(p, "periodogram");
            List<DataObject> dgjList = inDataGroup.values();
            DataType[] inColumns = inDataGroup.getDataDefinitions();
            Assert.assertTrue(inColumns.length + " is not 2", inColumns.length == 2);
            Assert.assertTrue("expected " + dgjList.size(), dgjList.size() == 390);
        } catch (IpacTableException e) {
            e.printStackTrace();
        }

    }

    @Test
    public void testGetPeaks() {

        boolean deleteOnExit = true;
        LightCurveHandler t = new IrsaLightCurveHandler() {
            @Override
            protected URL buildUrl(PeriodogramAPIRequest req) throws MalformedURLException {
                return new URL(req.getResultTable());
            }

            @Override
            protected File makeResultTempFile(RESULT_TABLES_IDX resultTable) throws IOException {
                File tmpFile = File.createTempFile("peaks-test-", ".tbl", new File("."));
                if (deleteOnExit) {
                    tmpFile.deleteOnExit();
                }
                return tmpFile;
            }

            @Override
            protected File makeApiResultTempFile() throws IOException {
                File tmpFile = File.createTempFile("votable-test-", ".xml", new File("."));
                if (deleteOnExit) {
                    tmpFile.deleteOnExit();
                }
                return tmpFile;
            }
        };

        File peaks = t.getPeaksTable(req);

        try {
            DataGroup inPeaksDataGroup = IpacTableReader.readIpacTable(peaks, "peaks");
            DataType[] inColumns = inPeaksDataGroup.getDataDefinitions();
            Assert.assertTrue(inColumns.length + " is not 5", inColumns.length == 5);
            List<DataObject> dgjList = inPeaksDataGroup.values();
            inColumns = inPeaksDataGroup.getDataDefinitions();
            Assert.assertTrue("expected " + dgjList.size(), dgjList.size() == req.getNumberPeaks());
        } catch (IpacTableException e) {
            e.printStackTrace();
        }

    }

    @Test
    public void testPhaseFoldedCurve() {

        boolean deleteOnExit = true;
        if(deleteOnExit){
            rawTable.deleteOnExit();
        }
        IrsaLightCurveHandler t = new IrsaLightCurveHandler(){
            @Override
            protected File createPhaseFoldedTempFile() throws IOException {
                File tmpFile = File.createTempFile("phase-folded-test-", ".tbl", new File("."));
                if (deleteOnExit) {
                    tmpFile.deleteOnExit();
                }
                return tmpFile;
            }
        };

        DataGroup inDataGroup = null;
        try {
            inDataGroup = IpacTableReader.readIpacTable(new File(req.getLcSource()), "lc_raw");
        } catch (IpacTableException e) {
            e.printStackTrace();
        }
        List<DataObject> dgjListOrigin = inDataGroup.values();
        DataType[] inColumns = inDataGroup.getDataDefinitions();

        File p = t.toPhaseFoldedTable(new File(req.getLcSource()), req.getPeriod(), req.getTimeColName());

        try {
            inDataGroup = IpacTableReader.readIpacTable(p, "phasefolded");
            List<DataObject> dgjList = inDataGroup.values();
            DataType[] inColumnsPhaseFolded = inDataGroup.getDataDefinitions();
            //should be one more extra column (Phase)
            Assert.assertTrue(inColumns.length + " is not correct", inColumns.length == inColumnsPhaseFolded.length-1);
            Assert.assertTrue("expected " + dgjList.size(), dgjList.size() == dgjListOrigin.size());
            double period = req.getPeriod();
            double tzero = (Double)dgjListOrigin.get(0).getDataElement("mjd");

            int iTest = 3;
            double mjdTest = (Double) dgjListOrigin.get(iTest).getDataElement("mjd");

            double phaseTested = (Double) inDataGroup.get(iTest).getDataElement("phase");
            double phaseExpected = (mjdTest-tzero)/period - Math.floor((mjdTest-tzero)/period);

            Assert.assertEquals(phaseExpected ,phaseTested,0.0001);
        } catch (IpacTableException e) {
            e.printStackTrace();
        }

    }

    /**
     * Class that will support a request test
     */
    static class PeriodogramAPIRequestTest extends PeriodogramAPIRequest {

        @Override
        public float getPeriod() {
            return 1.345f;
        }

        @Override
        public String getLcSource() {
            try {
                URL demo = new URL("http://web.ipac.caltech.edu/staff/ejoliet/demo/AllWISE-MEP-m82-2targets-10arsecs.tbl");
                URLConnection uc = URLDownload.makeConnection(demo);
                URLDownload.getDataToFile(uc, rawTable);
                return rawTable.getAbsolutePath();
            } catch (MalformedURLException e) {
                e.printStackTrace();
            } catch (FailedRequestException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return null;
        }

        @Override
        public int getNumberPeaks() {
            return 50;
        }

        @Override
        public String getTimeColName() {
            return "mjd";
        }

        /**
         * @return the built url api
         */
        @Override
        public String getUrl() {
            //As for the test, we return the result table
            return getResultTable();
        }

        @Override
        public String getResultTable() {
            return "http://web.ipac.caltech.edu/staff/ejoliet/demo/vo-nexsci-result-sample.xml";
        }
    }

    /**
     * Could be useful to define algorithm with default parameter and use them by mapping an enum from the request to ease the URL API building
     * Example of classes below
     */

    class LombScargle implements Periodogram {


        @Override
        public AlgorithmDefinition getAlgoDef() {
            return AlgorithmDefinition.LS;
        }

        @Override
        public int getNPeaks() {
            return 50;
        }

        @Override
        public Period getPeriod() {
            return new PeriodSample();
        }

        @Override
        public double[] getAlgoValues() {
            return new double[0];
        }

        @Override
        public StepMethod getStepMethod(StepMethod.STEPMETHOD_NAME sName) {
            return new FixedPeriodMethod(0.1f);
        }
    }

    class FixedPeriodMethod implements StepMethod {

        private final float val;

        FixedPeriodMethod(float step) {
            this.val = step;
        }

        @Override
        public String getName() {
            return STEPMETHOD_NAME.FIXED_PERIOD.name();
        }

        @Override
        public float getValue() {
            return val;
        }
    }

    protected class PeriodSample implements Period {

        @Override
        public float getMin() {
            return 0;
        }

        @Override
        public float getMax() {
            return 10;
        }

        @Override
        public float getPeakValue() {
            return 2;
        }
    }
}
