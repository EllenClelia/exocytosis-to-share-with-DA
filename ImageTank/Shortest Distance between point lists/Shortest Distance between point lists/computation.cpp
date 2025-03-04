#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTTable &from,const DTTable &to,
                    const DTMask2D &inside)
{
    
    // center
    DTTableColumnPoint2D fromPoint = from("center");
    DTTableColumnNumber fromTime = from("time");
    DTTableColumnPoint2D toPoints = to("Center");
    
    ssize_t howManyFrom = from.NumberOfRows();
    ssize_t howManyTo = to.NumberOfRows();
    ssize_t startPtN,endPtN;

    DTMutableDoubleArray distanceList(howManyFrom);
    DTMutableList<DTPoint2D> closestPoint(howManyFrom);
    // DTMutableDoubleArray closestPoint(2,howManyFrom);
    DTMutableIntArray closestIndex(howManyFrom);
    
    for (startPtN=0;startPtN<howManyFrom;startPtN++) {
        DTPoint2D startPoint = fromPoint(startPtN);
        double minDistance = INFINITY;
        ssize_t minIndex = 0;
        for (endPtN=0;endPtN<howManyTo;endPtN++) {
            DTPoint2D endPoint = toPoints(endPtN);
            double distance = Distance(startPoint,endPoint);
            if (distance<minDistance) {
                minDistance = distance;
                minIndex = endPtN;
            }
        }
        
        distanceList(startPtN) = minDistance;
        closestPoint(startPtN) = toPoints(minIndex);
        closestIndex(startPtN) = int(minIndex)+1;
    }

    // Table is a list of columns
    return DTTable({
        CreateTableColumn("exocytosis",fromPoint.Points()),
        CreateTableColumn("time",fromTime.DoubleVersion()),
        CreateTableColumn("erpm",closestPoint),
        CreateTableColumn("index",closestIndex),
        CreateTableColumn("distance",distanceList)
    });
}
