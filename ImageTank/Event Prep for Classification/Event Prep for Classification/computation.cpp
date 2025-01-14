#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include "Utilities.h"

PeakGroup Computation(const DTSet<DTImage> &combined,int start)
{
    PeakGroup toReturn;

    int whichEntry = start;
    DTImage image = combined(whichEntry);
    
    DTTable setParameters = combined.Parameters();
    
    DTMesh2DGrid grid = image.Grid(); // The underlying grid
    // Extract channels from image
    DTImageChannel channel = image("difference");
    channel = ConvertToFloat(channel);
    DTFloatArray values = channel.FloatArray();
    
    DTPoint2D origin = grid.Origin();
    double h = grid.dx();
    ssize_t m = image.m();
    ssize_t n = image.m();
    
    DTTableColumnPoint2D centerColumn = setParameters("center");
    DTPoint2D centerPoint = centerColumn(whichEntry);
    
    {
        // The center point is where we think the center of the peak is.
        // Want to compute the
        DTRegion2D cropWith(centerPoint.x-10*h,centerPoint.x+10*h,
                            centerPoint.y-10*h,centerPoint.y+10*h);
        DTTableColumnNumber tColumn = setParameters("t");
        ssize_t howManyRows = setParameters.NumberOfRows();
        ssize_t rowN;
        DTMutableDoubleArray averageList(howManyRows);
        for (rowN=0;rowN<howManyRows;rowN++) {
            DTImage singleImage = combined(rowN);
            singleImage = Crop(singleImage,cropWith);
            DTImageChannel singleChannel = ConvertToFloat(image("difference"));
            DTFloatArray values = channel.FloatArray();
            averageList(rowN) = Mean(values);
        }
        toReturn.intensities = DTTable({
            setParameters("t"),
            CreateTableColumn("average",averageList)
         });
    }

    DTMutableDictionary parameters;
    parameters("channel") = 1; // "difference";
    parameters("R2 for Peak") = 0.5;

    LocalPeak peak = FindGaussianPeak(image,parameters);
    // base + scale*exp(-((x-x0)^2+(y-y0)^2)/(2radius^2));

    // Evaluate the fit function
    double base = peak.base;
    double scale = peak.height-base;
    double x0 = peak.center.x;
    double y0 = peak.center.y;
    double width = peak.width;
    
    // width = 2*sqrt(2*log(2))*fabs(returned("radius")); // FWHM - full width at half maximumum
    // so radius = width/sqrt(2*log(2))
    double radius = width/(2*sqrt(2*log(2)))*h;

    DTMutableDoubleArray peakValues(m,n);
    for (ssize_t j=0;j<n;j++) {
        double y = origin.y + j*h;
        for (ssize_t i=0;i<m;i++) {
            double x = origin.x + i*h;
            peakValues(i,j) = base + scale*exp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/(2*radius*radius));
        }
    }

    toReturn.peak = DTImage(grid,peakValues,"height");
    toReturn.base = base;
    toReturn.height = scale;
    toReturn.threshold = base + (scale-base)*0.1;
    toReturn.center = peak.center;

    return toReturn;
}
