/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */

/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
package edu.caltech.ipac.firefly.visualize;

import com.google.gwt.core.client.Scheduler;
import com.google.gwt.event.dom.client.KeyCodes;
import com.google.gwt.event.dom.client.KeyDownEvent;
import com.google.gwt.event.dom.client.KeyDownHandler;
import com.google.gwt.user.client.ui.AbsolutePanel;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.FlowPanel;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.Widget;
import edu.caltech.ipac.firefly.ui.GwtUtil;
import edu.caltech.ipac.firefly.ui.input.InputField;
import edu.caltech.ipac.firefly.ui.input.SimpleInputField;
import edu.caltech.ipac.firefly.ui.input.TextBoxInputField;
import edu.caltech.ipac.firefly.ui.input.ValidationInputField;
import edu.caltech.ipac.firefly.visualize.draw.LayerDrawer;
import edu.caltech.ipac.firefly.visualize.draw.WebLayerItem;
import edu.caltech.ipac.firefly.visualize.ui.MaskAdjust;

/**
 * *
 */
public class OverlayPlotView extends Composite implements LayerDrawer {


    public static final String OVERLAY_ITEM_ID= "AnyMaskLayer";
    private WebPlotView pv;
    private AbsolutePanel rootPanel= new AbsolutePanel();
    private WebPlot      maskPlot= null;
    private String defaultColor= "ff0000";
    private final int maskValue;
    private final int imageNumber;
    private final WebLayerItem webLayerItem;
    private float opacity=.58F;
    private WebPlotRequest maskRequest;

    private static int idCnt=0;

    /**
     *
     */
    public OverlayPlotView(WebPlotView pv,
                           int maskValue,
                           int imageNumber,
                           String color,
                           String bitDesc,
                           WebPlotRequest maskRequest) {
        this.pv= pv;
        initWidget(rootPanel);
        this.maskValue= maskValue;
        this.imageNumber= imageNumber;
        this.defaultColor= color;
        this.maskRequest= maskRequest;
        rootPanel.setStyleName("OverlayPlotView");
        pv.addDrawingArea(this, false);

        String id= OverlayPlotView.OVERLAY_ITEM_ID+idCnt;
        pv.addOverlayPlotView(this);

        WebPlotView.MouseAll ma= new WebPlotView.DefMouseAll();
        WebPlotView.MouseInfo mi= new WebPlotView.MouseInfo(ma,"put more help here");
        webLayerItem= new WebLayerItem(id,null, "Mask: "+bitDesc,
                "Adjust opacity with up/down arrow",this,mi,null,null);

        if (!WebLayerItem.hasUICreator(id)) {
            WebLayerItem.addUICreator(id, new MaskUICreator());
        }

        idCnt++;
    }

    public WebLayerItem getWebLayerItem() { return webLayerItem; }

    public void setOpacity(float opacity) {
        this.opacity= opacity;
        if (pv!=null && maskPlot!=null) {
            maskPlot.setOpacity(opacity);
            pv.refreshDisplay();
        }
    }

    @Override
    public void setVisible(boolean visible) {
        super.setVisible(visible);
        if (visible && maskPlot!=null) maskPlot.refreshWidget();
    }

    @Override
    public boolean hasData() { return true; }

    @Override
    public String getDefaultColor() { return defaultColor; }

    @Override
    public void setDefaultColor(String c) {
        defaultColor= c;
        MaskAdjust.updateMask(this,pv,maskValue , imageNumber, c, maskRequest);
    }

    @Override
    public boolean getSupportsRegions() { return false; }

    public boolean isImageOverlay() { return true; }

    public void setPixelSize(int width, int height) {
        rootPanel.setPixelSize(width,height);
        super.setPixelSize(width, height);
    }

    public void clear() {
        if (this.maskPlot!=null) {
            this.maskPlot.freeResources();
            this.maskPlot= null;
            rootPanel.clear();
        }
    }

    public void dispose() {
        rootPanel.clear();
        pv.removeDrawingArea(this);
        pv.removeOverlayPlotView(this);
        maskPlot.freeResources();
        maskPlot= null;
    }

    public void setMaskPlot(WebPlot maskPlot) {
        if (maskPlot==null) return;
        if (this.maskPlot!=null) {
            this.maskPlot.freeResources();
        }
        maskPlot.getPlotGroup().setPlotView(pv);
        rootPanel.clear();
        rootPanel.add(maskPlot.getWidget(),0,0);
        maskPlot.setOpacity(opacity);
        this.maskPlot= maskPlot;
    }

    public WebPlot getMaskPlot() { return maskPlot; }

    private static class MaskUICreator extends WebLayerItem.UICreator {

        private MaskUICreator() { super(true,true); }


        public Widget makeExtraUI(final WebLayerItem item) {

            FlowPanel fp = new FlowPanel();
            final SimpleInputField field= SimpleInputField.createByProp("OverlayPlotView.opacity");
            field.setInternalCellSpacing(1);
            final InputField tIf= ((ValidationInputField)field.getField()).getIF();
            TextBox tb= ((TextBoxInputField)tIf).getTextBox();
            GwtUtil.setStyle(tb, "width", "2.5em");
            OverlayPlotView opv = (OverlayPlotView) item.getDrawer();
            field.setValue( ((int)(opv.opacity*100.0))+"");


            fp.add(field);


            tb.addKeyDownHandler(new KeyDownHandler() {
                public void onKeyDown(KeyDownEvent ev) {
                    final int keyCode= ev.getNativeKeyCode();
                    Scheduler.get().scheduleDeferred( new Scheduler.ScheduledCommand() {
                        public void execute() {
                            if (tIf.validate() && tIf.isNumber()) {
                                int v = tIf.getNumberValue().intValue();
                                OverlayPlotView opv = (OverlayPlotView) item.getDrawer();
                                opv.setOpacity((float) v / 100F);
                                if (keyCode== KeyCodes.KEY_DOWN&& v>0 && v<=100) {
                                    field.setValue((v-1)+"");
                                }
                                else if (keyCode== KeyCodes.KEY_UP && v>=0 && v<100) {
                                    field.setValue((v+1)+"");
                                }
                            }
                        }
                    });
                }
        });

            return fp;
        }

        public void delete(WebLayerItem item) {
            if (item.getDrawer() instanceof OverlayPlotView) {
                OverlayPlotView maskLayer= (OverlayPlotView)item.getDrawer();
                maskLayer.dispose();
                maskLayer.pv.removeWebLayerItem(item);
            }
        }
    }
}
