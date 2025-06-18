#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTTable.h"
#include "Adhesion.hpp"
#include "DTPath2DValues.h"

DTPath2DValues Computation(const DTImage &image)
{
    return DTPath2DValues();
}

double VariationPrev(const DTImageChannel &magnitude,DTPoint2D P,DTPoint2D Q)
{
    // 5/13 - disabled this.
    int N = 100;
    DTPoint2D delta = (Q-P)/N;
    DTDoubleArray data = magnitude.DoubleArray();
    
    // First and last value
    double firstValue = Interpolate(data,P);
    double lastValue = Interpolate(data,Q);
    
    double maxAboveEnds = std::max(firstValue,lastValue);
    double minBelowEnds = std::min(firstValue,lastValue);
    
    static bool display = false;
    
    for (int ptN=0;ptN<=N;ptN++) {
        DTPoint2D at = P + delta*ptN;
        double value = Interpolate(data,at);
        
        if (display) {
            std::cerr << value << std::endl;
        }

        if (value<minBelowEnds) minBelowEnds = value;
        if (value>maxAboveEnds) maxAboveEnds = value;
    }
    
    double risesAbove = maxAboveEnds-std::max(firstValue,lastValue);
    double goesBelow = std::min(firstValue,lastValue) - minBelowEnds;
    
    return std::max(risesAbove,goesBelow);
}

double Variation(const DTImageChannel &magnitude,const DTMesh2DGrid &grid,DTPoint2D P,DTPoint2D Q)
{
    DTTable interpolated = InterpolateSegment(magnitude,grid,P,Q,100);
    return Variation(interpolated);
}

struct PairingEntry
{
    int lowerIndex;
    int upperIndex;
    int rowIndex;
    
    bool operator<(const PairingEntry &other) const {return lowerIndex<other.lowerIndex;}
};



/*
double FindLocalMaxima(const DTImageChannel &magnitude,const DTMesh2DGrid &grid,const DTPoint2D &at,int r)
{
    DTDoubleArray data = magnitude.DoubleArray();
    
    int mMax = int(magnitude.m()-1);
    int nMax = int(magnitude.n()-1);
    int i,j;

    DTPoint2D atC = grid.SpaceToGrid(at);
    
    int startI = round(atC.x)-r;
    int endI = round(atC.x)+r;
    int startJ = round(atC.y)-r;
    int endJ = round(atC.y)+r;
    if (startI<0) startI = 0;
    if (startJ<0) startJ = 0;
    if (endI>mMax) endI = mMax;
    if (endJ>nMax) endJ = nMax;
    
    double maxV = data(startI,startJ);
    for (j=startJ;j<=endJ;j++) {
        for (i=startI;i<=endI;i++) {
            double v = data(i,j);
            if (v>maxV) maxV = v;
        }
    }

    return maxV;
}

double ComputeRatio(const DTImageChannel &magnitude,const DTMesh2DGrid &grid,const DTTable &interpolated,int r)
{
    DTTableColumnPoint2D points = interpolated("point");
    DTTableColumnNumber values = interpolated("value");
    ssize_t len = interpolated.NumberOfRows();

    DTPointCollection2D pointList = points.Points();

    int maxIndexAt = 0;
    double maxValue = values(0);
    for (int i=0;i<len;i++) {
        double v = values(i);
        if (v>maxValue) {
            maxValue = v;
            maxIndexAt = i;
        }
    }
    
    double ratio = FindLocalMaxima(magnitude,grid,points(maxIndexAt),2)/values(maxIndexAt);

    return ratio;
}
*/

MyGroup Computation(const DTImage &image,DTPoint2D from,DTPoint2D to,int N,
                    int Width)
{
    // Interpolate each channel
    DTImageChannel magnitudeChannel = image("magnitude");
    DTTable interpolated;
    
    if (Width==0) {
        interpolated = InterpolateSegment(magnitudeChannel,image.Grid(),from,to,N);
    }
    else {
        interpolated = InterpolateSegmentWide(magnitudeChannel,image.Grid(),from,to);
    }
    
    MyGroup toReturn;
    
    DTTableColumnPoint2D points = interpolated("point");
    DTTableColumnNumber values = interpolated("value");
    DTTableColumnNumber arcValues = interpolated("arc");
    DTDoubleArray arcList = arcValues.DoubleVersion();
    ssize_t len = interpolated.NumberOfRows();

    DTPointCollection2D pointList = points.Points();
    DTDoubleArray localMaximas = FindLocalMaxima(magnitudeChannel,image.Grid(),pointList,2);

    DTMutableDoubleArray pathArray(2,len+1);
    DTMutableDoubleArray valueArray(len+1);
    DTMutableDoubleArray maximaArray(len+1);
    DTMutableDoubleArray arcArray(len+1);
    pathArray(0,0) = 0;
    pathArray(1,0) = len;
    valueArray(0) = len;
    maximaArray(0) = len;
    arcArray(0) = len;
    
    double minRatioSoFar = INFINITY;
    int indexOfMinRatio = 0;
    for (int i=0;i<len;i++) {
        DTPoint2D p = points(i);
        pathArray(0,i+1) = p.x;
        pathArray(1,i+1) = p.y;
        valueArray(i+1) = values(i);
        maximaArray(i+1) = localMaximas(i);
        double currRatio = localMaximas(i)/values(i);
        if (currRatio<minRatioSoFar) {
            minRatioSoFar = currRatio;
            indexOfMinRatio = i;
        }
        arcArray(i+1) = arcList(i);
    }
    toReturn.ratio = minRatioSoFar;
    toReturn.arcForRatio = arcList(indexOfMinRatio);
    toReturn.ratio = ComputeRatio(magnitudeChannel,image.Grid(),interpolated,2);
    
    // ratioAtMax and arcForRatioAtMax
    // Find the maximum ratio. Start by computing the location and value of the maximum magnitude,
    int maxIndexAt = 0;
    double maxValue = values(0);
    for (int i=0;i<len;i++) {
        double v = values(i);
        if (v>maxValue) {
            maxValue = v;
            maxIndexAt = i;
        }
    }
    // double ratioAtMax = FindLocalMaxima(magnitudeChannel,image.Grid(),points(maxIndexAt),2)/values(maxIndexAt);
    toReturn.ratioAtMax = localMaximas(maxIndexAt)/maxValue;
    toReturn.arcForRatioAtMax = arcList(maxIndexAt);
    toReturn.ratioAtMax = ComputeRatioAtMax(magnitudeChannel,image.Grid(),interpolated,2);

    
    DTMutableList<DTPath2DValuesChannel> pathChannels(3);
    pathChannels(0) = DTPath2DValuesChannel("magnitude",valueArray);
    pathChannels(1) = DTPath2DValuesChannel("closeBy",maximaArray);
    pathChannels(2) = DTPath2DValuesChannel("arc",arcArray);
    toReturn.Values = DTPath2DValues(DTPath2D(pathArray),pathChannels);

    toReturn.Variation1 = Variation(interpolated);

    return toReturn;
                
}

