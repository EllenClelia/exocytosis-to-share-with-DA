#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

#include "Utilities.h"

DTMutableDoubleArray EvaluateFit(const DTMesh2DGrid &grid,const LocalPeak &fit)
{
    // Base + scale*exp(-((x-xC)^2+(y-yC)^2)/2radius^2);
    int m = (int)grid.m();
    int n = (int)grid.n();
    
    DTMutableDoubleArray returnArray(m,n);
    
    DTPoint2D center = fit.center;
    center = grid.SpaceToGrid(center);
    double x0 = center.x;
    double y0 = center.y;
    
    // toReturn.width = 2*sqrt(2*log(2))*fabs(returned("radius")); // FWHM - full width at half maximumum
    
    double radius = fit.width/(2*sqrt(2*log(2)));
    double base = fit.base;
    double scale = fit.height-base;

    /*
    DTMutableDictionary guesses;
    guesses("x0") = maxI;
    guesses("y0") = maxJ;
    guesses("scale") = maxV-minV;
    guesses("base") = minV;
    guesses("radius") = 5;
    
    
    double base = constants("base");
    double scale = constants("scale");
    double x0 = constants("x0");
    double y0 = constants("y0");
    double radius = constants("radius");
     */

    int i,j,ij;
    double x,y,arg;
    double C = 1.0/(2*radius*radius);
    ij = 0;
    double *returnArrayD = returnArray.Pointer();
    for (j=0;j<n;j++) {
        y = (j-y0);
        for (i=0;i<m;i++) {
            x = (i-x0);
            arg = -(x*x+y*y)*C;
            returnArrayD[ij] = base + scale*exp(arg);
            ij++;
        }
    }
    
    return returnArray;
}

MyGroup Computation(const DTImage &image)
{
    MyGroup toReturn;
    toReturn.input = image;
    
    DTMutableDictionary parameters;
    parameters("channel") = 0;
    parameters("R2 for Peak") = 0;

    // Compute the main peak
    LocalPeak fit = FindGaussianPeak(image,parameters);
    DTMutableDoubleArray values = EvaluateFit(image.Grid(),fit);
    toReturn.first_peak = DTImage(image.Grid(),values,"height");

    // Compute the fit and subtract it from the channel
    DTDoubleArray leftAfterFirstFitRemoved = values-ConvertToDouble(image(0)).DoubleArray();
    DTImage differenceImage = DTImage(image.Grid(),leftAfterFirstFitRemoved,"height");
    toReturn.first_residual = differenceImage;

    // Compute the second peak
    LocalPeak fitInDifference = FindGaussianPeak(differenceImage,parameters);
    values = EvaluateFit(image.Grid(),fitInDifference);
    toReturn.second_peak = DTImage(image.Grid(),values,"height");
    
    DTDoubleArray leftAfterSecondFitRemoved = leftAfterFirstFitRemoved-values;
    toReturn.second_residual = DTImage(image.Grid(),leftAfterSecondFitRemoved,"height");
    
    toReturn.first_base = fit.base;
    toReturn.first_height = fit.height;
    toReturn.first_center = fit.center;
    
    toReturn.second_base = fitInDifference.base;
    toReturn.second_height = fitInDifference.height;
    toReturn.second_center = fitInDifference.center;
    
    toReturn.flagFirst = fit.failureMode;
    toReturn.flagSecond = fitInDifference.failureMode;
    toReturn.R2first = fit.R2;
    toReturn.RMSEfirst = fit.RMSE;
    toReturn.R2second = fitInDifference.R2;
    toReturn.RMSEsecond = fitInDifference.RMSE;
    
    toReturn.ratio = toReturn.second_height/toReturn.first_height;
    
    toReturn.distance = Distance(toReturn.first_center,toReturn.second_center);

    return toReturn;
}
