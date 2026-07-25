#include "computation.h"

#include <math.h>
#include "DTDataFile.h"
#include "DTDoubleArrayOperators.h"
#include "DTImage.h"
#include "DTSet.h"
#include "DTTable.h"
#include "DTUtilities.h"

DTSet<DTImage> Computation(const DTSet<DTImage> &input,const DTTable &points,
                           double time,double before,double after,double width)
{
    DTSet<DTImage> withCache = input.WithCache(before+1+after);
    DTTable imageParameters = withCache.Parameters();
    DTTableColumnNumber timeValues = imageParameters("t");
    // int channel = parameters("channel");
    
    DTTableColumnNumber timeColumn;
    bool timeColumnExists = points.Contains("time");
    if (timeColumnExists) {
        timeColumn = points("time");
    }
    
    DTMutableSet<DTImage> output;
    
    // Flexible regarding the point to use.
    DTTableColumnPoint2D center;
    std::string tp = points("point").Type();
    if (points.Contains("point") && points("point").Type()=="Point2D") {
        center = points("point");
    }
    else {
        // Find the first point
        for (int cN=0;cN<points.NumberOfColumns();cN++) {
            if (points(cN).Type()=="Point2D") {
                center = points(cN);
                break;
            }
        }
        if (center.NumberOfRows()==0) {
            DTErrorMessage("Did not find the point column");
            return output;
        }
    }
    
    // Optional point number input
    DTTableColumnNumber pointNumberColumn;
    bool pointNumberSpecified = points.Contains("pointNumber");
    if (pointNumberSpecified) pointNumberColumn = points("pointNumber");

    int row,howMany = int(center.NumberOfRows());
    int posInOutput = 0;
    
    double w = width;
    
    DTMutableDoubleArray timeList(howMany*(before+1+after));
    DTMutableDoubleArray pointColumnNumberList(howMany*(before+1+after));
    DTMutableList<DTPoint2D> centerList(howMany*(before+1+after));
    
    DTProgress progress;
    for (row=0;row<howMany;row++) {
        progress.UpdatePercentage(row/double(howMany));
        
        if (timeColumnExists) {
            time = timeColumn(row);
        }
        int pointNumber = (pointNumberSpecified ? pointNumberColumn(row) : row+1);
        
        ssize_t where = timeValues.FindClosest(time);
        if (where<0) continue;
        DTPoint2D p = center(row);
        // The center value is not a good initial guess
        // Need to correct that by taking a window around the current center
        // find the maximum point of that window
        // and use that for the back/forward points - Crop(image,box) a few lines below
        
        // double w = width(row);
        DTRegion2D cropBox = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2); // The final image that is saved
        
        for (ssize_t index=where-before;index<=where+after;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {
                DTImage image = withCache(index);
                image = Crop(image,cropBox);
                output.Add(image);
                
                timeList(posInOutput) = index-where;
                centerList(posInOutput) = p;
                pointColumnNumberList(posInOutput) = pointNumber;
                
                posInOutput++;
            }
        }
    }
    
    if (posInOutput!=timeList.Length()) {
        timeList = TruncateSize(timeList,posInOutput);
        centerList = TruncateSize(centerList,posInOutput);
        pointColumnNumberList = TruncateSize(pointColumnNumberList,posInOutput);
    }
    
    // Create the parameter table, typically fill alongside the Add calls
    output.SetParameterTable(DTTable({
        CreateTableColumn("time",timeList),
        CreateTableColumn("center",DTPointCollection2D(centerList)),
        CreateTableColumn("pointNumber",pointColumnNumberList),
    }));
    
    return output;
}
