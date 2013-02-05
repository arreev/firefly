package edu.caltech.ipac.firefly.ui.table;

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.gen2.table.client.DefaultTableDefinition;
import com.google.gwt.gen2.table.client.FixedWidthFlexTable;
import com.google.gwt.gen2.table.client.FixedWidthGrid;
import com.google.gwt.gen2.table.client.ScrollTable;
import com.google.gwt.gen2.table.override.client.FlexTable;
import com.google.gwt.user.client.Event;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.CheckBox;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.HasVerticalAlignment;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import edu.caltech.ipac.firefly.core.Preferences;
import edu.caltech.ipac.firefly.data.table.TableData;
import edu.caltech.ipac.firefly.resbundle.images.TableImages;
import edu.caltech.ipac.firefly.ui.FormUtil;
import edu.caltech.ipac.firefly.ui.GwtUtil;
import edu.caltech.ipac.firefly.ui.PopupPane;
import edu.caltech.ipac.firefly.ui.input.SimpleInputField;
import edu.caltech.ipac.firefly.ui.table.filter.FilterPanel;
import edu.caltech.ipac.util.StringUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Date: May 21, 2009
 *
 * @author loi
 * @version $Id: TableOptions.java,v 1.18 2012/06/16 00:21:53 loi Exp $
 */
public class TableOptions extends Composite {

    private static final String VISI_COL_PREF = "-VisibleCols";
    private TablePanel table;
    private PopupPane main;
    private Map<ColDef, CheckBox> checkBoxes = new HashMap<ColDef, CheckBox>();
    private Widget popupContent;
    private String defVisibleCols;
    private CheckBox selectAllCheckBox = new CheckBox();
    private SimpleInputField pageSize;
    private FilterPanel filterPanel;


    public TableOptions(final TablePanel table) {
        this.table = table;

        defVisibleCols = getVisibleColStr(table.getTable());
        applyPrefVisibleColumns();
        Image options = new Image(TableImages.Creator.getInstance().getTableOptions());
//        options.setPixelSize(16, 16);
        options.setTitle("Edit Table Options");
        options.addClickHandler(new ClickHandler() {
            public void onClick(ClickEvent ev) {
                if (main.isVisible()) {
                    main.hide();
                } else {
                    if (popupContent == null) {
                        popupContent = makeContent();
                    }
                    // sync table's column with option's checkboxes
                    DefaultTableDefinition<TableData.Row> tdef = (DefaultTableDefinition<TableData.Row>) table.getTable().getTableDefinition();
                    for (ColDef col : checkBoxes.keySet()) {
                        CheckBox cb = checkBoxes.get(col);
                        cb.setValue(tdef.isColumnVisible(col));
                    }
                    
                    // sync table's filters with filterPanel
                    filterPanel.setFilters(table.getTable().getFilters());
                    
                    pageSize.setValue(String.valueOf(table.getLoader().getPageSize()));

                    main.alignTo(TableOptions.this.table, PopupPane.Align.TOP_RIGHT_OR_LEFT, 2, 0);
                    main.show();
                }
            }
        });

        main = new PopupPane("Table Options");
        initWidget(options);
    }

    public static void applyPrefVisibleColumns(TablePanel table) {

        String visibleCols = Preferences.get(table.getName() + VISI_COL_PREF);
        if (!StringUtils.isEmpty(visibleCols)) {
            List<String> vcols = Arrays.asList(visibleCols.split(";"));
            DefaultTableDefinition<TableData.Row> tdef = (DefaultTableDefinition<TableData.Row>) table.getTable().getTableDefinition();
            for (int i = 0; i < tdef.getColumnDefinitionCount(); i++) {
                ColDef col = (ColDef) tdef.getColumnDefinition(i);
                if (col.getColumn() != null) {
                    boolean isHidden = !vcols.contains(col.getName());
                    col.getColumn().setHidden(isHidden);
                    tdef.setColumnVisible(col, !isHidden);
                }
            }
        }
    }

    public static void setPrefVisibleColumns(String tableName, String[] vcols) {
        String visibleCols = "";
        for (String vcol : vcols) {
            visibleCols += (visibleCols.length() == 0 ? "" : ";") + vcol;
        }
        if (!StringUtils.isEmpty(visibleCols)) {
            Preferences.set(tableName + VISI_COL_PREF, visibleCols);
        }
    }

    public static String [] getPrefVisibleColumns(String tableName) {
        String visibleCols = Preferences.get(tableName + VISI_COL_PREF);
        if (!StringUtils.isEmpty(visibleCols)) {
            return visibleCols.split(";");
        } else {
            return new String[]{};
        }
    }


    private void applyPrefVisibleColumns() {
        applyPrefVisibleColumns(table);
    }

    private String getVisibleColStr(BasicPagingTable table) {
        String visibleCols = "";
        DefaultTableDefinition<TableData.Row> tdef = (DefaultTableDefinition<TableData.Row>) table.getTableDefinition();
        for (int i = 0; i < tdef.getColumnDefinitionCount(); i++) {
            ColDef col = (ColDef) tdef.getColumnDefinition(i);
            if (col.getColumn() != null) {
                if (!col.getColumn().isHidden()) {
                    visibleCols += (visibleCols.length() == 0 ? "" : ";") + col.getColumn().getName();
                }
            }
        }
        return visibleCols;
    }

    private void applyChanges() {
        DefaultTableDefinition<TableData.Row> tdef =
                (DefaultTableDefinition<TableData.Row>) table.getTable().getTableDefinition();
        boolean reloadNeeded = false;

        for (ColDef col : checkBoxes.keySet()) {
            CheckBox cb = checkBoxes.get(col);
            if (tdef.isColumnVisible(col) != cb.getValue()) {
                col.getColumn().setHidden(!cb.getValue());
                tdef.setColumnVisible(col, cb.getValue());
                reloadNeeded = true;
            }
        }

        int newPS = FormUtil.getIntValue(pageSize.getField());
        if ( newPS != table.getLoader().getPageSize()) {
            table.getPagingBar().reloadPageSize(newPS);
            return;
        }

        if (reloadNeeded) {
            String vcols = getVisibleColStr(table.getTable());
            if (vcols.equals(defVisibleCols)) {
                Preferences.set(table.getName() + VISI_COL_PREF, null);
            } else {
                Preferences.set(table.getName() + VISI_COL_PREF, vcols);
            }
            table.getTable().reloadPage();
        }
    }

    private Widget makeContent() {
        final ScrollTable colsTable = makeColsTable(table.getTable());
        colsTable.setSize("200px", "300px");

        Button reset = new Button("Reset", new ClickHandler() {
            public void onClick(ClickEvent ev) {
                if (!StringUtils.isEmpty(defVisibleCols)) {
                    List<String> vcols = Arrays.asList(defVisibleCols.split(";"));
                    for(Map.Entry<ColDef, CheckBox> entry : checkBoxes.entrySet()) {
                        entry.getValue().setValue(vcols.contains(entry.getKey().getName()));
                    }
                }
                ensureSelectAllCB();
            }
        });
        GwtUtil.setStyle(reset, "padding", "0 0");

        Button okBtn = new Button("Ok", new ClickHandler() {
            public void onClick(ClickEvent ev) {
                if (pageSize.validate()) {
                    applyChanges();
                    table.getTable().setFilters(filterPanel.getFilters());
                    table.doFilters();
                    main.hide();
                }
            }
        });

        Button cancelBtn = new Button("Cancel", new ClickHandler() {
            public void onClick(ClickEvent ev) {
                main.hide();
            }
        });

        GwtUtil.setStyles(okBtn, "paddingLeft", "15px", "paddingRight", "15px");
//
//        VerticalPanel vp = new VerticalPanel() {
//            protected void onLoad() {
//                colsTable.fillWidth();
//            }
//        };

        pageSize = makePageSizeField();
        filterPanel = new FilterPanel(table.getDataset().getColumns());

        FlexTable content = new FlexTable();
        content.setWidget(0, 0, pageSize);
        content.setWidget(1, 0, GwtUtil.makeHoriPanel(null, VerticalPanel.ALIGN_BOTTOM, new HTML("<b>Show/Hide column(s):</b> &nbsp;&nbsp;"), reset));
        content.setWidget(1, 1, new HTML("<b>Filters(s):</b>"));
        content.setWidget(2, 0, colsTable);
        content.setWidget(2, 1, filterPanel);
        content.setWidget(3, 0, GwtUtil.makeHoriPanel(null, null, cancelBtn, GwtUtil.getFiller(10, 0), okBtn));
        content.setCellSpacing(7);
        content.getRowFormatter().setVerticalAlign(1, VerticalPanel.ALIGN_BOTTOM);
        content.getFlexCellFormatter().setAlignment(2, 1, HorizontalPanel.ALIGN_LEFT, VerticalPanel.ALIGN_TOP);
        content.getFlexCellFormatter().setAlignment(3, 0, HorizontalPanel.ALIGN_RIGHT, VerticalPanel.ALIGN_MIDDLE);
        content.getFlexCellFormatter().setRowSpan(2, 0, 2);

        main.setWidget(content);
        ensureSelectAllCB();
        return content;
    }

    private SimpleInputField makePageSizeField() {
        final SimpleInputField pageSize = SimpleInputField.createByProp("TablePanel.pagesize");
        pageSize.setValue(table.getLoader().getPageSize()+"");
        return pageSize;
    }


    private ScrollTable makeColsTable(final BasicPagingTable table) {

        final FixedWidthFlexTable header = new FixedWidthFlexTable();
        header.setHTML(0, 0, "Column");
        header.setWidget(0, 1, selectAllCheckBox);
        selectAllCheckBox.addClickHandler(new ClickHandler() {
            public void onClick(ClickEvent ev) {
                boolean hasSel = false;
                for(Map.Entry<ColDef, CheckBox> entry : checkBoxes.entrySet()) {
                    if (entry.getValue().getValue()) {
                        hasSel = true;
                        break;
                    }
                }

                if (selectAllCheckBox.getValue() && !hasSel) {
                    for(Map.Entry<ColDef, CheckBox> entry : checkBoxes.entrySet()) {
                        entry.getValue().setValue(true);
                    }
                } else {
                    for(Map.Entry<ColDef, CheckBox> entry : checkBoxes.entrySet()) {
                        entry.getValue().setValue(false);
                    }
                    selectAllCheckBox.setValue(false);
                }
            }
        });

//        final SortableGrid.ColumnSorter[] origSorter = new SortableGrid.ColumnSorter[1];
        @SuppressWarnings("deprecation")
        final FixedWidthGrid data = new FixedWidthGrid(0, 2);
        data.unsinkEvents(Event.ONMOUSEOVER);
        data.setSelectionEnabled(false);

        final ScrollTable view = new ScrollTable(data, header, new BasicTable.Images());
        FlexTable.FlexCellFormatter formatter = header.getFlexCellFormatter();
        formatter.setHorizontalAlignment(0, 1,
                HasHorizontalAlignment.ALIGN_CENTER);
        view.setMaximumColumnWidth(1, 35);
        view.setMinimumColumnWidth(1, 35);
        view.setColumnSortable(1, false);

        final DefaultTableDefinition<TableData.Row> tdef = (DefaultTableDefinition<TableData.Row>) table.getTableDefinition();
        int cRowIdx = 0;
        for (int i = 0; i < tdef.getColumnDefinitionCount(); i++) {
            final ColDef col = (ColDef) tdef.getColumnDefinition(i);
            if (!col.isImmutable()) {
                data.insertRow(cRowIdx);
                data.setHTML(cRowIdx, 0, col.getTitle());

                CheckBox cb = new CheckBox();
                cb.setValue(tdef.isColumnVisible(col));
                checkBoxes.put(col, cb);
                data.setWidget(cRowIdx, 1, cb);
                data.getCellFormatter().setAlignment(cRowIdx, 1, HasHorizontalAlignment.ALIGN_CENTER, HasVerticalAlignment.ALIGN_MIDDLE);
                cRowIdx++;

                cb.addClickHandler(new ClickHandler(){
                    public void onClick(ClickEvent event) {
                        ensureSelectAllCB();
                    }
                });
            }
        }

        return view;
    }

    private void ensureSelectAllCB() {
        boolean selAll = true;

        for(Map.Entry<ColDef, CheckBox> entry : checkBoxes.entrySet()) {
            if (!entry.getValue().getValue()) {
                selAll = false;
                break;
            }
        }
        selectAllCheckBox.setValue(selAll);
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