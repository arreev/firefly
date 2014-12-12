package edu.caltech.ipac.visualize.plot;

import java.util.EventObject;

/**
 * The event that is passed when a PlotView status has a changed
 * @author Trey Roby
 */
public class PlotGroupStatusEvent extends EventObject {

    private Plot _plot;
    private int  _colorBand= -1;

    /**
     * Create a new PlotViewStatusEvent
     * @param plotGroup source of the event.
     * @param plot the plot added or removed
     */
    public PlotGroupStatusEvent (PlotGroup plotGroup, Plot plot) {
        this(plotGroup, plot, -1);
    }
    /**
     * Create a new PlotViewStatusEvent
     * @param plotGroup source of the event.
     * @param plot the plot added or removed
     */
    public PlotGroupStatusEvent (PlotGroup plotGroup, Plot plot, int band) {
        super(plotGroup);
        _plot= plot;
        _colorBand= band;
    }

    public Plot getPlot() {return _plot;}
    public int  getColorBand() {return _colorBand;}
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
