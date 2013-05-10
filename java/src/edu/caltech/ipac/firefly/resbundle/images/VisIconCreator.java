package edu.caltech.ipac.firefly.resbundle.images;

import com.google.gwt.core.client.GWT;
import com.google.gwt.resources.client.ClientBundle;
import com.google.gwt.resources.client.ImageResource;
/**
 * User: roby
 * Date: Sep 1, 2009
 * Time: 11:22:02 AM
 */


/**
 * @author Trey Roby
 */
public interface VisIconCreator extends ClientBundle {

    @Source("cyan_left_arrow_20x20.png")
    public ImageResource getSideLeftArrow();

    @Source("cyan_right_arrow_20x20.png")
    public ImageResource getSideRightArrow();

    @Source("cyan_down_arrow_20x20.png")
    public ImageResource getSideDownArrow();

    @Source("cyan_up_arrow_20x20.png")
    public ImageResource getSideUpArrow();

    @Source("zoom-out_v2-20x20.png")
    public ImageResource getZoomDown();

    @Source("zoom-in_v2-20x20.png")
    public ImageResource getZoomUp();

    @Source("zoom-original_v2-20x20.png")
    public ImageResource getZoomOriginal();

    @Source("zoom-fit_v2-20x20.png")
    public ImageResource getZoomFit();

    @Source("zoom-fill_v2-20x20.png")
    public ImageResource getZoomFill();

    @Source("restore-20x20.png")
    public ImageResource getRestore();

    @Source("lock_20x20.png")
    public ImageResource getLocked();

    @Source("unlock_20x20.png")
    public ImageResource getUnlocked();

    @Source("save_20x20.png")
    public ImageResource getSave();

    @Source("grid-20x20.png")
    public ImageResource getGridOff();

    @Source("grid-on-20x20.png")
    public ImageResource getGridOn();

    @Source("select-20x20.png")
    public ImageResource getSelectAreaOff();

    @Source("select-dark-20x20.png")
    public ImageResource getSelectAreaOn();

    @Source("settings-16x16.png")
    public ImageResource getSettings();

    @Source("crop_20x20.png")
    public ImageResource getCrop();

    @Source("crop_20x20.png")
    public ImageResource getCatalogWithSelection();

    @Source("palette_20x20.png")
    public ImageResource getColorTable();

    @Source("palette_stretch_20X20.png")
    public ImageResource getStretch();

    @Source("palette_stretch_quick_20X20.png")
    public ImageResource getStretchQuick();


    @Source("header-info_20x20.png")
    public ImageResource getFitsHeader();

    @Source("plot_layers_dim_20x20.png")
    public ImageResource getLayer();

    @Source("plot_layers_20x20.png")
    public ImageResource getLayerBright();

    @Source("rotate-north-on_20x20.png")
    public ImageResource getRotateNorthOn();

    @Source("rotate-north_20x20.png")
    public ImageResource getRotateNorth();

    @Source("flip_20x20.png")
    public ImageResource getFlip();

    @Source("rotate_20x20.png")
    public ImageResource getRotate();

    @Source("distance_20x20.png")
    public ImageResource getDistanceOn();


    @Source("marker-off-20x20.png")
    public ImageResource getMarkerOff();

    @Source("marker-on-20x20.png")
    public ImageResource getMarkerOn();



    @Source("distance-off_20x20.png")
    public ImageResource getDistanceOff();

    @Source("current-target_20x20.png")
    public ImageResource getCurrentTarget();

    @Source("catalog_20x20.png")
    public ImageResource getCatalog();

    @Source("image-working-background-24x24.png")
    public ImageResource getImageWorkingBackground();


    @Source("starry_sky_20x20.png")
    public ImageResource getStarrySky();

    @Source("statistics_20x20.png")
    public ImageResource getStatistics();


    @Source("compass-20x20.png")
    public ImageResource getCompass();

    @Source("compass-on-20x20.png")
    public ImageResource getCompassOn();


    @Source("sun-20x20.png")
    public ImageResource getDS9Symbol();

//    @Source("step-right-16x16.png")
//    public ImageResource getStepRight();
//
//    @Source("step-left-16x16.png")
//    public ImageResource getStepLeft();



    public static class Creator  {
        private final static VisIconCreator _instance=
                (VisIconCreator) GWT.create(VisIconCreator.class);
        public static VisIconCreator getInstance() {
            return _instance;
        }
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
