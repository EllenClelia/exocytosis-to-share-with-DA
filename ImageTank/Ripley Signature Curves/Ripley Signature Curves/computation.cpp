#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTDoubleArray ComputeL(const DTTable &points,
                       double area,
                       const DTDoubleArray &rList)
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
    ssize_t rCount = rList.Length();
    DTMutableDoubleArray values(rCount);
    
    values(0) = 0;
    pos = 0;
    ssize_t totalLength = distances.Length();
    // ssize_t maxCount = distances(totalLength-1)/dr;
    
    // Want to make sure that pos<totalLength when distances(pos)<until
    // that is I want distances(totalLength-1)>=until
    // that is distances(totalLength-1)>=range.minV + (maxCount-1)*dr
    // that is (distances(totalLength-1)-range.minV)/dr+1 >= maxCount
    
    ssize_t distLength = distances.Length();
    for (i=1;i<rCount;i++) {
        double until = rList(i);
        double untilSquared = until*until;
        while (pos<distLength && distances(pos)<untilSquared) {
            pos++;
        }
        values(i) = pos;
    }

    double scale = 1.0/totalLength/M_PI*area;
    for (i=0;i<rCount;i++) {
        values(i) = sqrt(values(i)*scale);
    }
    
    return values;
}

DTDoubleArray FindP(const DTDoubleArray &pdfArray,double p,const DTMesh2DGrid &grid)
{
    ssize_t i,j,m = pdfArray.m(), n = pdfArray.n();
    DTMutableDoubleArray toReturn(m);
    double breakAt;
    double h = grid.dy();
    double origin = grid.Origin().y;
    for (i=0;i<m;i++) {
        // Go through the list pdfArray(i,:) and find where the value crosses p
        breakAt = NAN;
        if (pdfArray(i,0)>=p) {
            breakAt = 0;
        }
        else if (pdfArray(i,n-1)<=p) {
            breakAt = n-1;
        }
        else {
            for (j=0;j<n-1;j++) {
                if (pdfArray(i,j)<=p && p<pdfArray(i,j+1)) {
                    // Find where the value would be p exactly
                    breakAt = j + (p-pdfArray(i,j))/(pdfArray(i,j+1)-pdfArray(i,j));
                    break;
                }
            }
        }
        breakAt = origin + h*breakAt;
        toReturn(i) = breakAt;
    }
    
    return toReturn;
}

DTDoubleArray MAD(const DTImage &pdf)
{
    // norminv(0.75) = 0.674489750196082
    DTDoubleArray values = ConvertToDouble(pdf(0)).DoubleArray();
    ssize_t m = values.m();
    ssize_t n = values.n();
    ssize_t i,j,k,lower,upper;
    
    double dy = pdf.Grid().dy();
    
    DTMutableDoubleArray mad(m);

    for (i=0;i<m;i++) {
        // Median is where percentage is 0.5
        for (j=0;j<n;j++) {
            if (values(i,j)>0.5) {
                break;
            }
        }
        
        if (j==n) {
            // Didn't find the median
            mad(i) = NAN;
            continue;
        }
        if (j>0 && fabs(values(i,j-1)-0.5)<fabs(values(i,j))) {
            j--;
        }
        if (j>0) {
            // values(i,j) is closest to 0.5, so consider that the median.
            // Compute the median of the absolute distance from the median, that
            // means that I find k such that values(i,j+k)-values(i,j-k)==0.5
            for (k=1;k<n;k++) {
                lower = std::max(j-k,0L);
                upper = std::min(j+k,n);
                if (values(i,upper-1)-values(i,lower)>0.5) {
                    break;
                }
            }
            mad(i) = k*dy/0.674489750196082;
        }
        else {
            mad(i) = 0;
        }
    }

    return mad;
}

DTTable Computation(const DTImage &pdf,double p,const DTMask2D &mask,
                    const DTTable &points)
{
    // Extract channels from pdf
    // DTImageChannel channel = pdf("value") or pdf(0)
    // If you want to access this as a float even for 8,32 or 64 bit
    // channel = ConvertToFloat(channel); // no cost if already float
    // DTFloatArray values = channel.FloatArray(); // Complains if not stored as floats

    DTMesh2DGrid grid = pdf.Grid();
    DTCharArray maskArray = mask.Mask().MaskArray();
    
    
    // Compute the L for the point/mask combo
    ssize_t i,m = grid.m();
    double dr = grid.dx();
    double r0 = grid.Origin().x;
    DTMutableDoubleArray rList(m);
    for (i=0;i<m;i++) {
        rList(i) = r0 + dr*i;
    }

    ssize_t howManyPoints = mask.Mask().NumberOfPoints();
    DTMesh2DGrid maskGrid = mask.Grid();
    double area = howManyPoints*maskGrid.dx()*maskGrid.dy();
    DTDoubleArray lList = ComputeL(points,area,rList);
    
    
    // Compute the "level" curves for p, 0.5, and 1-p from the pdf
    DTDoubleArray pdfArray = pdf(0).DoubleArray();
    DTDoubleArray lowerList = FindP(pdfArray,p,grid);
    DTDoubleArray medianList = FindP(pdfArray,0.5,grid);
    DTDoubleArray upperList = FindP(pdfArray,1-p,grid);
    
    // Find the pdf at the y coordinates from lList
    DTMutableDoubleArray pdfAtL(m);
    double y0 = grid.Origin().y;
    double dy = grid.dy();
    ssize_t n = pdfArray.n();
    for (i=0;i<m;i++) {
        ssize_t j = round((lList(i)-y0)/dy);
        if (j<0) {
            pdfAtL(i) = 0;
        }
        else if (j>=n-1) {
            pdfAtL(i) = 1;
        }
        else {
            pdfAtL(i) = pdfArray(i,j);
        }
    }
    
    DTDoubleArray mad = MAD(pdf);

    // Table is a list of columns
    return DTTable({
        CreateTableColumn("r",rList),
        CreateTableColumn("lower",lowerList),
        CreateTableColumn("median",medianList),
        CreateTableColumn("mad",mad),
        CreateTableColumn("high",upperList),
        CreateTableColumn("Actual",lList),
        CreateTableColumn("pdfAtL",pdfAtL)
    });
}
