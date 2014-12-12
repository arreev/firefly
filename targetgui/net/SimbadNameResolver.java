package edu.caltech.ipac.targetgui.net;

import edu.caltech.ipac.astro.simbad.SimbadClient;
import edu.caltech.ipac.astro.simbad.SimbadException;
import edu.caltech.ipac.astro.simbad.SimbadObject;
import edu.caltech.ipac.client.net.FailedRequestException;
import edu.caltech.ipac.client.net.HostPort;
import edu.caltech.ipac.client.net.NetworkManager;
import edu.caltech.ipac.client.net.ThreadedService;
import edu.caltech.ipac.target.PositionJ2000;
import edu.caltech.ipac.target.ProperMotion;
import edu.caltech.ipac.target.SimbadAttribute;
import edu.caltech.ipac.util.action.ClassProperties;

import java.awt.Window;
import java.io.IOException;
import java.net.UnknownHostException;



/**
 * @author Xiuqin Wu, based on Trey Roby's CoordConvert
 */
public class SimbadNameResolver extends ThreadedService {


    private String          _in = null;
    private SimbadAttribute _out= null;
    private final static ClassProperties _prop= new ClassProperties(
                                           SimbadNameResolver.class);
    private static final String SERVER_SCRIPT= _prop.getName("script");
    private static final boolean USE_SIMBAD4_CLIENT = _prop.getSelected("use.simbad4");





    private SimbadNameResolver (String objname, Window w) {
        super(w);
        setOperationDesc(_prop.getName("desc"));
        setProcessingDesc(_prop.getName("searching") + objname);
        _in = objname;
//        if (USE_SIMBAD4_CLIENT) {
//            _simbad4Client = new Simbad4Client();
//        } else {
////            System.out.println("server: " + _server.getHost()
////                             + "  script: " + SERVER_SCRIPT);
//        }
    }

    protected void doService() throws Exception {
        _out= lowlevelNameResolver(_in);
    }

    public static SimbadAttribute getPosition(String objname, Window w)
                                         throws FailedRequestException {
        SimbadNameResolver action= new SimbadNameResolver(objname,w);
        action.execute();
        System.out.println("SimbadNameResolver: got results");
        return action._out;
    }

    public static SimbadAttribute lowlevelNameResolver(String objname)
                                         throws  FailedRequestException {
        SimbadObject simbadObject;
        PositionJ2000 positionOut = null;

        try {
            if (USE_SIMBAD4_CLIENT) {
                Simbad4Client simbad4Client = new Simbad4Client();
                simbadObject = simbad4Client.searchByName(objname);
            } else {
                HostPort server= NetworkManager.getInstance().getServer(
                        NetworkManager.SIMBAD_NAME_RESOLVER);
                SimbadClient simbadClient = new SimbadClient(server.getHost(), SERVER_SCRIPT);
                simbadObject = simbadClient.searchByName(objname);
            }
        }
        catch (UnknownHostException uhe) {
            throw new FailedRequestException(uhe.getMessage(), null, uhe);
        }
        catch (SimbadException ne) {
            throw new FailedRequestException(
                    "Simbad did not find the object: " + objname,
                    ne.getMessage(), ne);
        }
        catch (IOException ioe) {
            throw new FailedRequestException(
                    "Simbad Service Unavailable",  null, ioe);
        }


        ProperMotion pm= PositionJ2000.DEFAULT_PM;
        if (!Float.isNaN(simbadObject.getRaPM()) &&
                !Float.isNaN(simbadObject.getDecPM()) ) {
            pm=   new ProperMotion(simbadObject.getRaPM()/1000.0F,
                    simbadObject.getDecPM()/1000.0F);
        }

        positionOut = new PositionJ2000( simbadObject.getRa(),
                simbadObject.getDec(), pm);

        SimbadAttribute sa= new SimbadAttribute(positionOut);
        sa.setName(simbadObject.getName());
        sa.setType(simbadObject.getType());
        sa.setBMagnitude(simbadObject.getBMagnitude());
        if ( simbadObject.getMagBand().equalsIgnoreCase("V") ) {
            sa.setVMagnitude(simbadObject.getMagnitude());
        }
        sa.setSpectralType(simbadObject.getSpectralType());
        sa.setParallax(simbadObject.getParallax());
        return sa;
    }
}
/*
 * THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA 
 * INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. GOVERNMENT CONTRACT WITH 
 * THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION (NASA). THE SOFTWARE 
 * IS TECHNOLOGY AND SOFTWARE PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS 
 * AND IS PROVIDED AS-IS TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, 
 * INCLUDING ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR 
 * A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC 2312-2313) 
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
