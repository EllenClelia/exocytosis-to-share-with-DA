#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTSet<DTImage> &images)
{
    // images is a set
    ssize_t images_count = images.NumberOfItems();
    DTTable images_par = images.Parameters();
    // time
    //   DTTableColumnNumber values = images_par("time");
    //   DTDoubleArray time = values.DoubleVersion();
    // intensity
    //   DTTableColumnNumber values = images_par("intensity");
    //   DTDoubleArray intensity = values.DoubleVersion();
    // center
    //   DTTableColumnPoint2D pointsColumn = images_par("center");
    //   DTPointCollection2D center = pointsColumn.Points();
    // ptNumber
    //   DTTableColumnNumber values = images_par("ptNumber");
    //   DTDoubleArray ptNumber = values.DoubleVersion();
    // centerSpot
    //   DTTableColumnPoint2D pointsColumn = images_par("centerSpot");
    //   DTPointCollection2D centerSpot = pointsColumn.Points();
    // average
    //   DTTableColumnNumber values = images_par("average");
    //   DTDoubleArray average z= values.DoubleVersion();

    // Create the table
    return DTTable({
          CreateTableColumn("bckgrnd",DTDoubleArray()),
          CreateTableColumn("decay",DTDoubleArray()),
          CreateTableColumn("R2",DTDoubleArray()),
          CreateTableColumn("center",DTPointCollection2D())
    });
}
