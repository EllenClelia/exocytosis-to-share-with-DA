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
    ssize_t len = interpolated.NumberOfRows();
    
    DTMutableDoubleArray pathArray(2,len+1);
    DTMutableDoubleArray valueArray(len+1);
    pathArray(0,0) = 0;
    pathArray(1,0) = len;
    valueArray(0) = len;
    for (int i=0;i<len;i++) {
        DTPoint2D p = points(i);
        pathArray(0,i+1) = p.x;
        pathArray(1,i+1) = p.y;
        valueArray(i+1) = values(i);
    }

    toReturn.Values = DTPath2DValues(DTPath2D(pathArray),valueArray,"magnitude");
    toReturn.Variation1 = Variation(interpolated);

    return toReturn;
                
}

