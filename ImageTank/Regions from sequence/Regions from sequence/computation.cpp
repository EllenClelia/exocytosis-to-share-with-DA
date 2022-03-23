#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

void Computation(const DTSet<DTImage> &everything,const DTTable &spots,double time,int timeback,
                 int timeforward,double pixels,DTMutableSet<DTImage> &output)
{
    // everything is a set
    // Path
    //   DTTableColumnPath2D paths = spots("Path");
    //   DTPath2D singlePath = paths(0);
    //   DTPath2D combined = paths.Path();
    //   DTIntArray starts = paths.StartsOfIntervals();
    // Center
    //   DTTableColumnPoint2D pointsColumn = spots("Center");
    //   DTPointCollection2D Center = pointsColumn.Points();
    // Area
    //   DTTableColumnNumber values = spots("Area");
    //   DTDoubleArray Area = values.DoubleVersion();
    // Length
    //   DTTableColumnNumber values = spots("Length");
    //   DTDoubleArray Length = values.DoubleVersion();
    // Orientation
    //   DTTableColumnNumber values = spots("Orientation");
    //   DTDoubleArray Orientation = values.DoubleVersion();
    // Closed
    //   DTTableColumnNumber values = spots("Closed");
    //   DTDoubleArray Closed = values.DoubleVersion();
    // ABox
    //   DTTableColumnPath2D paths = spots("ABox");
    //   DTPath2D singlePath = paths(0);
    //   DTPath2D combined = paths.Path();
    //   DTIntArray starts = paths.StartsOfIntervals();
    // ABox.width
    //   DTTableColumnNumber values = spots("ABox.width");
    //   DTDoubleArray ABoxwidth = values.DoubleVersion();
    // ABox.height
    //   DTTableColumnNumber values = spots("ABox.height");
    //   DTDoubleArray ABoxheight = values.DoubleVersion();
    // PBox
    //   DTTableColumnPath2D paths = spots("PBox");
    //   DTPath2D singlePath = paths(0);
    //   DTPath2D combined = paths.Path();
    //   DTIntArray starts = paths.StartsOfIntervals();
    // PBox.width
    //   DTTableColumnNumber values = spots("PBox.width");
    //   DTDoubleArray PBoxwidth = values.DoubleVersion();
    // PBox.height
    //   DTTableColumnNumber values = spots("PBox.height");
    //   DTDoubleArray PBoxheight = values.DoubleVersion();
    
    spots.pall();

    output.Add(everything(0));
    
    DTSet<DTImage> withCache = everything.WithCache(timeback+1+timeforward);
    DTTable parameters = everything.Parameters();
    DTTableColumnNumber timeValues = parameters("t");
    ssize_t where = timeValues.FindClosest(time);
    
    DTTableColumnPoint2D center = spots("Center");
    
    ssize_t row,howMany = center.NumberOfRows();
    for (row=0;row<howMany;row++) {
        
    }

    // Create the parameter table, typically fill along side the Add calls
    DTMutableList<DTTableColumn> columns(2);
    DTMutableDoubleArray values(1);
    values(0) = 5;
    columns(0) = CreateTableColumn("time",values);
    // At the end, wrap up everything. Can't add entries after that
    columns(1) = CreateTableColumn("intensity",values);
    output.Finish(DTTable(columns));
    // At the end, wrap up everything. Can't add entries after that
    output.Finish(DTTable(columns));
}
