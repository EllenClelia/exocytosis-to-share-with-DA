#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

void Computation(const DTSet<DTImage> &everything,const DTTable &points,double time,
                 double timeback,double timeforward,double width,
                 DTMutableSet<DTImage> &output)
{
    DTSet<DTImage> withCache = everything.WithCache(timeback+1+timeforward);
    DTTable imageParameters = withCache.Parameters();
    DTTableColumnNumber timeValues = imageParameters("t");
    // int channel = parameters("channel");
    
    DTTableColumnNumber timeColumn;
    bool timeColumnExists = points.Contains("time");
    if (timeColumnExists) {
        timeColumn = points("time");
    }
    
    DTTableColumnPoint2D center = points("point");
    ssize_t row,howMany = center.NumberOfRows();
    int posInOutput = 0;
    
    double w = width;
    
    DTMutableDoubleArray timeList(howMany*(timeback+1+timeforward));
    DTMutableList<DTPoint2D> centerList(howMany*(timeback+1+timeforward));
    
    DTProgress progress;
    for (row=0;row<howMany;row++) {
        progress.UpdatePercentage(row/double(howMany));
        
        if (timeColumnExists) {
            time = timeColumn(row);
        }
        
        ssize_t where = timeValues.FindClosest(time);
        if (where<0) continue;
        DTPoint2D p = center(row);
        // The center value is not a good initial guess
        // Need to correct that by taking a window around the current center
        // find the maximum point of that window
        // and use that for the back/forward points - Crop(image,box) a few lines below
        
        // double w = width(row);
        DTRegion2D cropBox = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2); // The final image that is saved
        
        for (ssize_t index=where-timeback;index<=where+timeforward;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {
                DTImage image = withCache(index);
                image = Crop(image,cropBox);
                output.Add(image);
                
                timeList(posInOutput) = index-where;
                centerList(posInOutput) = p;
                
                posInOutput++;
            }
        }
    }
    
    if (posInOutput!=timeList.Length()) {
        timeList = TruncateSize(timeList,posInOutput);
        centerList = TruncateSize(centerList,posInOutput);
    }
    
    // Create the parameter table, typically fill alongside the Add calls
    output.Finish(DTTable({
        CreateTableColumn("time",timeList),
        CreateTableColumn("center",DTPointCollection2D(centerList)),
    }));
}
