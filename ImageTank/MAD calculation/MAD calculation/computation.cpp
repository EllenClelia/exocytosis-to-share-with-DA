#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include <math.h>

#include <algorithm>

DTDoubleArray MAD(const DTImage &pdf)
{
    // norminv(0.75) = 0.674489750196082
    DTDoubleArray values = ConvertToDouble(pdf(0)).DoubleArray();
    ssize_t m = values.m();
    ssize_t n = values.n();
    ssize_t i,j,k;
    
    double dy = pdf.Grid().dy();
    
    DTMutableDoubleArray mad(m);

    for (i=0;i<m;i++) {
        // Median is where percentage is 0.5
        for (j=0;j<n;j++) {
            if (values(i,j)>0.5) {
                break;
            }
        }
        if (j>0 && fabs(values(i,j-1)-0.5)<fabs(values(i,j))) {
            j--;
        }
        if (j>0) {
            // values(i,j) is closest to 0.5, so consider that the median.
            // Compute the median of the absolute distance from the median, that
            // means that I find k such that values(i,j+k)-values(i,j-k)==0.5
            for (k=1;k<n;k++) {
                ssize_t lower = std::max(j-k,0L);
                ssize_t upper = std::min(j+k,n);
                if (values(i,upper)-values(i,lower)>0.5) {
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

DTTable Computation(const DTImage &pdf)
{
    DTMesh2DGrid grid = pdf.Grid();
    ssize_t i,m = pdf.m();
    double dr = grid.dx();
    double r0 = grid.Origin().x;
    DTMutableDoubleArray rList(m);
    for (i=0;i<m;i++) {
        rList(i) = r0 + dr*i;
    }
    
    DTDoubleArray mad = MAD(pdf);

    // Table is a list of columns
    return DTTable({
        CreateTableColumn("r",rList),
        CreateTableColumn("mad",mad)
    });
}
