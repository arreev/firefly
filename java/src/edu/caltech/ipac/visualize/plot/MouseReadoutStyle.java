/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
package edu.caltech.ipac.visualize.plot;

import java.awt.event.MouseEvent;

/**
 * Date: Sep 28, 2006
 *
 * @author Trey Roby
 * @version $id:$
 */
public class MouseReadoutStyle {


    private int _rows;
    private int _columns;
    private int _styleGrid[][];
//============================================================================
//---------------------------- Constructors ----------------------------------
//============================================================================
    public MouseReadoutStyle(int rows, int columns) {
        _styleGrid= new int[rows][columns];
    }

//============================================================================
//---------------------------- Public Methods --------------------------------
//============================================================================

    public void setCellStyle(int cellStyle, int row, int column) {
        _styleGrid[row][column]= cellStyle;
    }

//============================================================================
//---------------------------- Methods from xxx Interface --------------------
//============================================================================

//============================================================================
//---------------------------- Private / Protected Methods -------------------
//============================================================================

//============================================================================
//---------------------------- Factory Methods -------------------------------
//============================================================================

//============================================================================
//---------------------------- Inner Classes ---------------------------------
//============================================================================
    public interface ReadoutType {
        public String getCellValue(Plot plot, MouseEvent ev,  WorldPt wpt);
    }

}

