#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

void Computation(const DTSet<DTImage> &everything,const DTTable &spots,double time,int timeback,
                 int timeforward,double pixels,DTMutableSet<DTImage> &output)
{
    DTSet<DTImage> withCache = everything.WithCache(timeback+1+timeforward);
    DTTable parameters = withCache.Parameters();
    DTTableColumnNumber timeValues = parameters("t");
    ssize_t where = timeValues.FindClosest(time);
    
    DTTableColumnPoint2D center = spots("center");
    DTTableColumnNumber width = spots("width");
    
    ssize_t row,howMany = center.NumberOfRows();
    ssize_t index;
    
    ssize_t posInOutput = 0;
    DTMutableDoubleArray timeList(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray intensityList(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray centerList(2,howMany*(timeback+1+timeforward));
    DTMutableDoubleArray pointNumber(howMany*(timeback+1+timeforward));

    for (row=0;row<howMany;row++) {
        DTPoint2D p = center(row);
        double w = width(row);
        DTRegion2D box = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2);
        for (index=where-timeback;index<=where+timeforward;index++) {
            if (index<0) continue;
            if (index>=withCache.NumberOfItems()) continue;
            DTImage image = withCache(index);
            image = Crop(image,box);
            output.Add(image);
            pointNumber(posInOutput) = row;
            timeList(posInOutput) = index-where;
            intensityList(posInOutput) = Maximum(image(0));
            centerList(0,posInOutput) = p.x;
            centerList(1,posInOutput) = p.y;
            posInOutput++;
        }
    }
    
    if (posInOutput!=timeList.Length()) {
        pointNumber = TruncateSize(pointNumber,posInOutput);
        timeList = TruncateSize(timeList,posInOutput);
        intensityList = TruncateSize(intensityList,posInOutput);
        centerList = TruncateSize(centerList,2*posInOutput);
    }

    // Create the parameter table, typically fill along side the Add calls
    output.Finish(DTTable({
        CreateTableColumn("time",timeList),
        CreateTableColumn("intensity",intensityList),
        CreateTableColumn("center",DTPointCollection2D(centerList)),
        CreateTableColumn("ptNumber",pointNumber)
    }));
    // At the end, wrap up everything. Can't add entries after that
    output.Finish(DTTable());
}
