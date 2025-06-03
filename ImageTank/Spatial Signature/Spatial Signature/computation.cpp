#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTRandom.h"
#include "DTRegion1D.h"
#include "DTProgress.h"

DTTable Realization(const DTMesh2DGrid &grid,
                    const DTCharArray &mask,
                    DTRandom &randNumber,
                    int length)
{
    DTMutableDoubleArray timeList(length);
    DTMutableDoubleArray pointList(2,length);
    
    //DTMesh2DGrid grid = region.Grid();
    //DTCharArray mask = region.Mask().MaskArray();
    ssize_t m = grid.m();
    // ssize_t mn = m*n;
    double width = grid.m();
    double height = grid.n();
    double h = grid.dx();
    DTPoint2D origin = grid.Origin();
    double x0 = origin.x;
    double y0 = origin.y;

    // Assume that most of the entries are included
    // If that is not the case I can create a list of all of the interior points
    // and then draw a offset into that list at random
    
    ssize_t i,j,ij;
    double x,y;
    for (ssize_t pos=0;pos<length;pos++) {
        // Random in space
        while (1) {
            x = randNumber.UniformHalf()*width;
            y = randNumber.UniformHalf()*height;
            i = floor(x);
            j = floor(y);
            ij = i+j*m;
            if (mask(ij)) {
                pointList(0,pos) = x0 + x*h; // (i-0.5)*h;
                pointList(1,pos) = y0 + y*h; // (j-0.5)*h;
                break;
            }
        }
    }

    return DTTable({
        CreateTableColumn("xy",DTPointCollection2D(pointList))
    });
}

DTTable ComputeL(const DTMesh2DGrid &grid,
                 const DTCharArray &maskArray,
                 ssize_t howManyPoints,
                 DTRandom &randNumber,
                 int count,
                 const DTRegion1D &range,
                 int rCount)
{
    DTTable points = Realization(grid,maskArray,randNumber,count);
    if (points.IsEmpty()) return DTTable();
    
    DTTableColumnPoint2D pointsColumn = points("xy");
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
    
    // DTRegion1D range(0,1);
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

    DTMutableDoubleArray lValues(rCount);
    for (i=0;i<rCount;i++) {
        lValues(i) = sqrt(values(i));
    }

    return DTTable({
        CreateTableColumn("r",rList),
        CreateTableColumn("L",lValues)
    });
}

Group Computation(const DTMask2D &mask,int count,int seed,DTRegion1D rRange,
                  DTRegion1D yRange,int rCount,int yCount,int runs)
{
    DTRandom randNumber(seed);
    DTMesh2DGrid grid = mask.Grid();
    DTCharArray maskArray = mask.Mask().MaskArray();
    
    
    DTTable lValues = ComputeL(grid,maskArray,mask.Mask().NumberOfPoints(),randNumber,count,rRange,rCount);

    Group toReturn;

    toReturn.L = lValues;
    toReturn.L2 = ComputeL(grid,maskArray,mask.Mask().NumberOfPoints(),randNumber,count,rRange,rCount);
    
    DTMutableDoubleArray values(rCount,yCount);
    double dr = (rRange.maxV-rRange.minV)/(rCount-1);
    double dy = (yRange.maxV-yRange.minV)/(yCount-1);
    
    values = 0.0;
    
    int iteration;
    int i,j;
    double yScale = 1.0/dy;
    
    DTProgress progress;

    for (iteration=0;iteration<runs;iteration++) {
        progress.UpdatePercentage(iteration/(runs*1.0));
        lValues = ComputeL(grid,maskArray,mask.Mask().NumberOfPoints(),randNumber,count,rRange,rCount);
        
        DTTableColumnNumber listOfL = lValues("L");
        DTDoubleArray lList = listOfL.DoubleVersion();
        for (i=0;i<rCount;i++) {
            int yLoc = round((lList(i)-yRange.minV)*yScale);
            if (yLoc>=0 && yLoc<yCount) {
                values(i,yLoc)++;
            }
        }
    }
    
    // Compute the PDF
    double scaleCount = 1.0/runs;
    for (i=0;i<rCount;i++) {
        double sum = 0;
        for (j=0;j<yCount;j++) {
            sum += values(i,j);
            values(i,j) = sum*scaleCount;
        }
    }
    
    DTMesh2DGrid pdfGrid(DTPoint2D(rRange.minV,yRange.minV),dr,dy);
    toReturn.pdf = DTImage(pdfGrid,values,"value");

    toReturn.sample = Realization(grid,maskArray,randNumber,count);

    return toReturn;
}
