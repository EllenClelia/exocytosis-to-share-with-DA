#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTImage &grid,const DTTable &points,int radius)
{
    DTMesh2DGrid gr = grid.Grid(); // The underlying grid
    
    if (radius<0) radius = 0;

    DTTableColumnPoint2D pointsColumn = points("point");
    DTPointCollection2D pointList = pointsColumn.Points();
    
    int m = int(gr.m());
    int n = int(gr.n());
    
    ssize_t howMany = pointList.NumberOfPoints();
    DTMutableIntArray startList(howMany+1);
    DTMutableIntArray intervals(2,(2*radius+1)*howMany);
    int posInIntervals = 0;
    
    int rSq = radius*radius;
    
    for (int ptN=0;ptN<howMany;ptN++) {
        startList(ptN) = posInIntervals;
        DTPoint2D P = gr.SpaceToGrid(pointList(ptN));
        int iC = int(round(P.x));
        int jC = int(round(P.y));
        
        int minJ = jC-radius;
        if (minJ<0) minJ = 0;
        int maxJ = jC+radius+1;
        if (maxJ>n) maxJ = n;
        for (int j=minJ;j<maxJ;j++) {
            // want (iC-minI)^2 + (jc-j)^2 = rSq
            int minI = round(iC - sqrt((rSq-(jC-j)*(jC-j))));
            int maxI = round(iC + sqrt((rSq-(jC-j)*(jC-j))));
            if (minI<0) minI = 0;
            if (maxI>m-1) maxI = m-1;
            if (minI<=maxI) {
                intervals(0,posInIntervals) = minI + j*m;
                intervals(1,posInIntervals) = maxI + j*m;
                posInIntervals++;
            }
        }
    }
    startList(howMany) = posInIntervals;
    if (posInIntervals<intervals.n()) {
        intervals = TruncateSize(intervals,2*posInIntervals);
    }

    // Table is a list of columns
    return DTTable({
        points("point"),
        CreateTableColumn("mask",gr,intervals,startList)
    });
}
