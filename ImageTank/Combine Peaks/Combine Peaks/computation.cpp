#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTImage Computation(DTPoint2D origin,double h,int m,int n,double radius,
                    const DTTable &points)
{
    // point
    DTTableColumnPoint2D pointsColumn = points("point");
    DTPointCollection2D point = pointsColumn.Points();
    
    DTMutableDoubleArray values(m,n);
    DTMesh2DGrid grid(origin,h);
    values = 0;
    
    double r = radius/h;
    
    ssize_t ptN,howManyPoints = point.NumberOfPoints();
    int i,j;
    double f,v;
    for (ptN=0;ptN<howManyPoints;ptN++) {
        DTPoint2D p = grid.SpaceToGrid(point(ptN));
        int iStart = ceil(p.x-r);
        int iEnd = floor(p.x+r)+1;
        int jStart = ceil(p.y-r);
        int jEnd = floor(p.y+r)+1;
        if (iStart<0) iStart = 0;
        if (jStart<0) jStart = 0;
        if (iEnd>m) iEnd = m;
        if (jEnd>n) jEnd = n;
        for (j=jStart;j<jEnd;j++) {
            for (i=iStart;i<iEnd;i++) {
                f = sqrt((i-p.x)*(i-p.x)+(j-p.y)*(j-p.y))/r; // Between 0 and 1
                if (f>1) f = 1;
                v = 1-f*f*f;
                v = v*v*v;
                values(i,j) += v;
            }
        }
    }

    return DTImage(grid,values,"intensity");
}
