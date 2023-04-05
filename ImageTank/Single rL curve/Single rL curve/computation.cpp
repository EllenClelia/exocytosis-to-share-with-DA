#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTRandom.h"

DTPath2D ComputeL(const DTTable &points,
                 const DTMesh2DGrid &grid,
                 const DTCharArray &maskArray,
                 ssize_t howManyPoints,
                 const DTRegion1D &range,
                 int rCount)
{
    DTTableColumnPoint2D pointsColumn = points(0);
    DTPointCollection2D point = pointsColumn.Points();
    DTDoubleArray data = point.Data();
    const double *dataD = data.Pointer();
    
    ssize_t i,j,howMany = data.n();
    DTMutableDoubleArray distances((howMany*(howMany-1))/2);
    double *distancesD = distances.Pointer();
    ssize_t pos = 0;
    double x1,y1,r,dx,dy;
    for (i=1;i<howMany;i++) {
        x1 = data(0,i);
        y1 = data(1,i);
        for (j=0;j<i;j++) {
            dx = dataD[2*j]-x1;
            dy = dataD[1+2*j]-y1;
            r = dx*dx+dy*dy; // Don't compute the square root, that is expensive and will be handled later
            distancesD[pos++] = r;
        }
    }
    
    std::sort(distances.Pointer(),distances.Pointer()+distances.Length());
    
    DTMutableDoubleArray rList(rCount);
    DTMutableDoubleArray values(rCount);
    
    rList(0) = range.minV;
    values(0) = 0;
    double dr = (range.maxV-range.minV)/(rCount-1);
    pos = 0;
    ssize_t totalLength = distances.Length();
    // ssize_t maxCount = distances(totalLength-1)/dr;
    
    // Want to make sure that pos<totalLength when distances(pos)<until
    // that is I want distances(totalLength-1)>=until
    // that is distances(totalLength-1)>=range.minV + (maxCount-1)*dr
    // that is (distances(totalLength-1)-range.minV)/dr+1 >= maxCount
    
    ssize_t maxCount = floor(sqrt(distances(totalLength-1))-range.minV)/dr+1;
    if (maxCount>rCount) maxCount = rCount;
    // while (maxCount<totalLength && distances(maxCount)>=range.minV+maxCount) maxCount++;
    double until = range.minV + (maxCount-1)*dr;
    if (distances(totalLength-1)<until*until) {
        maxCount--;
    }
    for (i=1;i<maxCount;i++) {
        until = range.minV + i*dr;
        double untilSquared = until*until;
        while (distances(pos)<untilSquared) {
            pos++;
        }
        values(i) = pos;
        rList(i) = until;
    }
    for (i=maxCount;i<rCount;i++) {
        values(i) = pos;
        rList(i) = range.minV + i*dr;
    }
    
    //DTMesh2DGrid grid = mask.Grid();
    // ssize_t howManyPoints = grid.Mask().NumberOfPoints();
    double area = howManyPoints*grid.dx()*grid.dy();
    
    double scale = 1.0/totalLength/M_PI*area;
    for (i=0;i<rCount;i++) {
        values(i) *= scale;
    }
    

    DTMutableDoubleArray pathValues(2,rCount+1);
    pathValues(0,0) = 0;
    pathValues(1,0) = rCount;

    for (i=0;i<rCount;i++) {
        pathValues(0,i+1) = rList(i);
        pathValues(1,i+1) = sqrt(values(i));
        // lValues(i) = sqrt(values(i));
    }
    
    return DTPath2D(pathValues);
}

DTPath2D Computation(const DTMask2D &mask,const DTTable &points,int rCount,
                     DTRegion1D rRange)
{
    DTMesh2DGrid grid = mask.Grid();
    DTCharArray maskArray = mask.Mask().MaskArray();
    
    return ComputeL(points,grid,maskArray,mask.Mask().NumberOfPoints(),rRange,rCount);
}
