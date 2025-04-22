#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

double Interpolate(const DTDoubleArray &data,DTPoint2D at)
{
    int m = int(data.m());
    int n = int(data.n());

    double i = floor(at.x);
    double j = floor(at.y);
    if (i<0) i = 0;
    if (j<0) j = 0;
    if (i>=m-2) i = m-2;
    if (j>=n-2) j = n-2;
    
    double dx = at.x-i;
    double dy = at.y-j;

    double val = data(i,j)*(1-dx)*(1-dy) + data(i+1,j)*dx*(1-dy) + data(i,j+1)*(1-dx)*dy + data(i+1,j+1)*dx*dy;

    return val;
}

double Variation(const DTImageChannel &magnitude,DTPoint2D P,DTPoint2D Q)
{
    int N = 100;
    DTPoint2D delta = (Q-P)/N;
    DTDoubleArray data = magnitude.DoubleArray();
    int m = int(data.m());
    int n = int(data.n());
    double prev = NAN;
    double total = 0;
    double first = NAN, last = NAN;
    double minV = INFINITY;
    double maxV = -INFINITY;
    
    int i,j;
    
    
    // First and last value
    double firstValue = Interpolate(data,P);
    double lastValue = Interpolate(data,Q);
    double h = 1.0/N;
    double sumOfSquares = 0;
    
    double maxAboveEnds = std::max(firstValue,lastValue);
    double minBelowEnds = std::min(firstValue,lastValue);
    
    static bool display = false;
    
    for (int ptN=0;ptN<=N;ptN++) {
        DTPoint2D at = P + delta*ptN;
        double value = Interpolate(data,at);
        
        if (display) {
            std::cerr << value << std::endl;
        }
        /*
        double prediction = firstValue + double(ptN)/N*(lastValue-firstValue);
        double difference = value-prediction;
         sumOfSquares += difference*difference;
         */
        if (value<minBelowEnds) minBelowEnds = value;
        if (value>maxAboveEnds) maxAboveEnds = value;
        
        /*
        // Interpolate at the point
        double i = floor(at.x);
        double j = floor(at.y);
        if (i<0) i = 0;
        if (j<0) j = 0;
        if (i>=m-2) i = m-2;
        if (j>=n-2) j = n-2;
        double dx = at.x-i;
        double dy = at.y-j;
        double val = data(i,j)*(1-dx)*(1-dy) + data(1+1,j)*dx*(1-dy) + data(i,j+1)*(1-dx)*dy + data(i+1,j+1)*dx*dy;
        double add = fabs(val-prev);
        if (ptN==0) first = val;
        if (ptN==N) last = val;
        if (isfinite(add)) total+=add;
        if (val<minV) minV = val;
        if (val>maxV) maxV = val;
        prev = val;
         */
    }
    
    double risesAbove = maxAboveEnds-std::max(firstValue,lastValue);
    double goesBelow = std::min(firstValue,lastValue) - minBelowEnds;
    
    return std::max(risesAbove,goesBelow);
    
    // return sqrt(sumOfSquares/N);
}

DTTable Computation(const DTTable &extrema,const DTImage &gradient,
                    double distance)
{
    // point
    DTTableColumnPoint2D pointsColumn = extrema("point");
    DTPointCollection2D point = pointsColumn.Points();

    DTMesh2DGrid grid = gradient.Grid(); // The underlying grid
    // Extract channels from gradient
    // DTImageChannel channel = gradient("gx") or gradient(0)
    // DTImageChannel channel = gradient("gy") or gradient(1)
    DTImageChannel magnitudeChannel = gradient("magnitude");
    magnitudeChannel = ConvertToDouble(magnitudeChannel);
    DTDoubleArray magnitude = magnitudeChannel.DoubleArray();
    // If you want to access this as a float even for 8,32 or 64 bit
    // channel = ConvertToFloat(channel); // no cost if already float
    // DTFloatArray values = channel.FloatArray(); // Complains if not stored as floats
    
    ssize_t howManyPoints = point.NumberOfPoints();
    DTDoubleArray pointList = point.Data();
    ssize_t i,j;
    DTPoint2D Q;
    double dx,dy,d;
    double distSq = distance*distance;
    
    int lenOfTable = int(howManyPoints);
    DTMutableDoubleArray fromList(2,lenOfTable);
    DTMutableDoubleArray toList(2,lenOfTable);
    DTMutableDoubleArray variationList(lenOfTable);
    int posInTable = 0;
    DTMutableDoubleArray centerList(2,lenOfTable);
    DTMutableDoubleArray intervalList(lenOfTable);
    DTMutableDoubleArray startList(lenOfTable);
    DTMutableDoubleArray endList(lenOfTable);

    for (j=0;j<howManyPoints;j++) {
        DTPoint2D P = point(j);
        for (i=j+1;i<howManyPoints;i++) {
            Q.x = pointList(0,i);
            Q.y = pointList(1,i);
            
            dx = Q.x-P.x;
            dy = Q.y-P.y;
            d = dx*dx+dy*dy;
            if (d<distSq) {
                if (posInTable==lenOfTable) {
                    fromList = IncreaseSize(fromList);
                    toList = IncreaseSize(fromList);
                    variationList = IncreaseSize(variationList);
                    centerList = IncreaseSize(centerList);
                    intervalList = IncreaseSize(intervalList);
                    startList = IncreaseSize(startList);
                    endList = IncreaseSize(endList);
                }
                fromList(0,posInTable) = P.x;
                fromList(1,posInTable) = P.y;
                toList(0,posInTable) = Q.x;
                toList(1,posInTable) = Q.y;
                DTPoint2D c = (P+Q)/2;
                centerList(0,posInTable) = c.x;
                centerList(1,posInTable) = c.y;
                intervalList(posInTable) = posInTable+1;
                startList(posInTable) = i+1;
                endList(posInTable) = j+1;

                variationList(posInTable) = Variation(magnitudeChannel,grid.SpaceToGrid(P),grid.SpaceToGrid(Q));
                
                posInTable++;
            }
        }
    }
    
    // Table is a list of columns
    return DTTable({
        CreateTableColumn("from",DTPointCollection2D(TruncateSize(fromList,2*posInTable))),
        CreateTableColumn("to",DTPointCollection2D(TruncateSize(toList,2*posInTable))),
        CreateTableColumn("variation",TruncateSize(variationList,posInTable)),
        CreateTableColumn("center",DTPointCollection2D(TruncateSize(centerList,2*posInTable))),
        CreateTableColumn("interval",TruncateSize(intervalList,posInTable)),
        CreateTableColumn("start",TruncateSize(startList,posInTable)),
        CreateTableColumn("end",TruncateSize(endList,posInTable))
      });
}
