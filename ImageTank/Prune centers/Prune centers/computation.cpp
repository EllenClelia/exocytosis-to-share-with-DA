#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTTable.h"
#include "Adhesion.hpp"

double VariationPrev(const DTImageChannel &magnitude,DTPoint2D P,DTPoint2D Q)
{
    // 5/13 - disabled this.
    int N = 100;
    DTPoint2D delta = (Q-P)/N;
    DTDoubleArray data = magnitude.DoubleArray();
    /*
    int m = int(data.m());
    int n = int(data.n());
    double prev = NAN;
    double total = 0;
    double first = NAN, last = NAN;
    double minV = INFINITY;
    double maxV = -INFINITY;
    */

    
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

MyGroup Computation(const DTTable &extrema,const DTImage &magnitude,
                    double distance,double threshold)
{
    // point
    DTTableColumnPoint2D pointsColumn = extrema("point");
    DTPointCollection2D point = pointsColumn.Points();

    DTMesh2DGrid grid = magnitude.Grid(); // The underlying grid
    // Extract channels from gradient
    // DTImageChannel channel = gradient("gx") or gradient(0)
    // DTImageChannel channel = gradient("gy") or gradient(1)
    DTImageChannel magnitudeChannel = magnitude("magnitude");
    magnitudeChannel = ConvertToDouble(magnitudeChannel);
    // DTDoubleArray intensity = magnitudeChannel.DoubleArray();
    // If you want to access this as a float even for 8,32 or 64 bit
    // channel = ConvertToFloat(channel); // no cost if already float
    // DTFloatArray values = channel.FloatArray(); // Complains if not stored as floats
    
    ssize_t howManyPoints = point.NumberOfPoints();
    DTDoubleArray pointList = point.Data();
    int i,j;
    DTPoint2D Q;
    double dx,dy,d;
    double distSq = distance*distance;
    
    int lenOfTable = int(howManyPoints);
    DTMutableDoubleArray fromList(2,lenOfTable);
    DTMutableIntArray fromIndex(lenOfTable);
    DTMutableDoubleArray toList(2,lenOfTable);
    DTMutableIntArray toIndex(lenOfTable);
    DTMutableDoubleArray variationList(lenOfTable);
    DTMutableDoubleArray ratioList(lenOfTable);
    int posInTable = 0;
    DTMutableDoubleArray centerList(2,lenOfTable);
    DTMutableDoubleArray intervalList(lenOfTable);
    DTMutableDoubleArray startList(lenOfTable);
    DTMutableDoubleArray endList(lenOfTable);
    
    DTMutableList<DTPoint2D> finalListOfPoints(howManyPoints);
    int posInFinalList = 0;
    DTMutableIntArray addedThisPoint(howManyPoints);
    addedThisPoint = 0;

    bool addedPoint;
    DTPoint2D P;
    for (j=0;j<howManyPoints;j++) {
        DTPoint2D Q = point(j);
        addedPoint = false;
        for (i=j+1;i<howManyPoints;i++) {
            P.x = pointList(0,i);
            P.y = pointList(1,i);
            
            dx = Q.x-P.x;
            dy = Q.y-P.y;
            d = dx*dx+dy*dy;
            if (d<distSq) {
                if (posInTable==lenOfTable) {
                    fromList = IncreaseSize(fromList);
                    toList = IncreaseSize(toList);
                    variationList = IncreaseSize(variationList);
                    ratioList = IncreaseSize(ratioList);
                    centerList = IncreaseSize(centerList);
                    intervalList = IncreaseSize(intervalList);
                    startList = IncreaseSize(startList);
                    endList = IncreaseSize(endList);
                    fromIndex = IncreaseSize(fromIndex);
                    toIndex = IncreaseSize(toIndex);
                }
                fromList(0,posInTable) = Q.x;
                fromList(1,posInTable) = Q.y;
                toList(0,posInTable) = P.x;
                toList(1,posInTable) = P.y;
                DTPoint2D c = (P+Q)/2;
                centerList(0,posInTable) = c.x;
                centerList(1,posInTable) = c.y;
                intervalList(posInTable) = posInTable+1;
                startList(posInTable) = j;
                endList(posInTable) = i;
                
                fromIndex(posInTable) = i;
                toIndex(posInTable) = j;
                addedThisPoint(i)++;
                addedThisPoint(j)++;

                DTTable interpolated = InterpolateSegment(magnitudeChannel,grid,P,Q,100);
                variationList(posInTable) = Variation(interpolated);
                ratioList(posInTable) = ComputeRatio(magnitudeChannel,grid,interpolated,2);
                
                addedPoint = true;
                
                posInTable++;
            }
        }
        
        // if (addedPoint) addedThisPoint(j) = 1;
    }
    
    // Trim the result
    fromList = TruncateSize(fromList,2*posInTable);
    toList = TruncateSize(toList,2*posInTable);
    variationList = TruncateSize(variationList,posInTable);
    ratioList = TruncateSize(ratioList,posInTable);
    centerList = TruncateSize(centerList,2*posInTable);
    intervalList = TruncateSize(intervalList,posInTable);
    startList = TruncateSize(startList,posInTable);
    endList = TruncateSize(endList,posInTable);
    
    MyGroup toReturn;
    
    // Create the table process
    DTTable processTable = DTTable({
        CreateTableColumn("from",DTPointCollection2D(fromList)),
        CreateTableColumn("to",DTPointCollection2D(toList)),
        CreateTableColumn("variation",variationList),
        CreateTableColumn("ratio",ratioList),
        CreateTableColumn("center",DTPointCollection2D(centerList)),
        CreateTableColumn("interval",intervalList),
        CreateTableColumn("start",startList),
        CreateTableColumn("end",endList)
    });
    
    toReturn.process = processTable;
    
    // Prune out the rows where the variation is too large
    DTMutableIntArray rowsToUse(posInTable);
    int posInNewTable = 0;
    for (int rowNumber=0;rowNumber<posInTable;rowNumber++) {
        // if (variationList(rowNumber)<threshold) {
        if (ratioList(rowNumber)>threshold) {
            rowsToUse(posInNewTable++) = rowNumber;
        }
        else {
            addedThisPoint(fromIndex(rowNumber))--;
            addedThisPoint(toIndex(rowNumber))--;
        }
    }
    
    // Add the points that are not connected to any other points.
    // Leaves points that are a part of a group.
    for (i=0;i<howManyPoints;i++) {
        if (addedThisPoint(i)==0) {
            finalListOfPoints(posInFinalList++) = point(i);
        }
    }
    
    // rowsToUse is now the entries from the processing table where the variation is < threshold.
    // The next step is to remove from that list any row where the end points are singletons

    DTTable acceptedSegments = processTable.ExtractRows(TruncateSize(rowsToUse,posInNewTable));

    // Find clusters, i.e. points that are connected together in a loop.
    DTMutableCharArray includedInTheLoop(howManyPoints);
    includedInTheLoop = 0;

    DTMutableCharArray pointsThatHaveBeenHandled(howManyPoints);
    includedInTheLoop = 0;

    DTMutableIntArray segmentNumberForIndex(howManyPoints+1);
    segmentNumberForIndex = -1; // Not found yet

    DTTableColumnNumber start = acceptedSegments("start");
    DTTableColumnNumber end = acceptedSegments("end");
    
    DTMutableIntArray pendingPointsInLoop(howManyPoints);
    int posInPendingLoop = 0;
    

    DTMutableCharArray handledSegment(posInNewTable);
    handledSegment = 0;
    DTMutableCharArray pointsIncludedInLoop(howManyPoints);
    pointsIncludedInLoop = 0;
    DTMutableIntArray pointsFoundList(howManyPoints);
    int pointsFoundCount = 0;
    
    for (i=0;i<posInNewTable;i++) {
        if (handledSegment(i)) continue; // Already handled this segment
        
        // This segment hasn't been seen before, that means that the start and end points
        // don't belong to a segment.
        int startIndex = start(i);
        int endIndex = end(i);
        
        pointsFoundCount = 0;
        
        // These points belong to a new segment. Brute force go through the remaining segments.
        pendingPointsInLoop(posInPendingLoop++) = startIndex;
        pointsIncludedInLoop(startIndex) = 1;
        pointsFoundList(pointsFoundCount++) = startIndex;
        
        pendingPointsInLoop(posInPendingLoop++) = endIndex;
        pointsIncludedInLoop(endIndex) = 1;
        pointsFoundList(pointsFoundCount++) = endIndex;
        
        while (posInPendingLoop>0) {
            int findPoint = pendingPointsInLoop(posInPendingLoop-1);
            posInPendingLoop--;
            
            for (int search=i+1;search<posInNewTable;search++) {
                if (handledSegment(search)) continue;
                    
                if (start(search)==findPoint) {
                    // Include the end point in the segment, and add that to the list
                    endIndex = end(search);
                    if (pointsIncludedInLoop(endIndex)==0) {
                        pendingPointsInLoop(posInPendingLoop++) = endIndex;
                        pointsIncludedInLoop(endIndex) = 1;
                        pointsFoundList(pointsFoundCount++) = endIndex;
                    }
                    handledSegment(search) = 1;
                }
                if (end(search)==findPoint) {
                    startIndex = start(search);
                    // Include the end point in the segment, and add that to the list
                    if (pointsIncludedInLoop(startIndex)==0) {
                        pendingPointsInLoop(posInPendingLoop++) = startIndex;
                        pointsIncludedInLoop(startIndex) = 1;
                        pointsFoundList(pointsFoundCount++) = startIndex;
                    }
                    handledSegment(search) = 1;
                }
            }
        }
        
        // Marked all of the points that are in this segment, now compute the center of mass for those points
        DTPoint2D sumXY(0,0);
        for (int ptN=0;ptN<pointsFoundCount;ptN++) {
            int foundPoint = pointsFoundList(ptN);
            pointsFoundList(ptN) = -1; // Not needed, just for debugging
            pointsIncludedInLoop(foundPoint) = 0;
            sumXY += pointsColumn(foundPoint);
        }
        finalListOfPoints(posInFinalList++) = sumXY/pointsFoundCount;
    }
    
    // Now go through the connected components
    finalListOfPoints = TruncateSize(finalListOfPoints,posInFinalList);
    
    // Create the table centers
    toReturn.centers = DTTable({
        CreateTableColumn("center",finalListOfPoints)
    });

    return toReturn;
}

MyGroup ComputationOld(const DTTable &extrema,const DTImage &magnitude,
                    double distance,double threshold)
{
    // point
    DTTableColumnPoint2D pointsColumn = extrema("point");
    DTPointCollection2D point = pointsColumn.Points();

    DTMesh2DGrid grid = magnitude.Grid(); // The underlying grid
    // Extract channels from gradient
    // DTImageChannel channel = gradient("gx") or gradient(0)
    // DTImageChannel channel = gradient("gy") or gradient(1)
    DTImageChannel magnitudeChannel = magnitude("magnitude");
    magnitudeChannel = ConvertToDouble(magnitudeChannel);
    // DTDoubleArray intensity = magnitudeChannel.DoubleArray();
    // If you want to access this as a float even for 8,32 or 64 bit
    // channel = ConvertToFloat(channel); // no cost if already float
    // DTFloatArray values = channel.FloatArray(); // Complains if not stored as floats
    
    ssize_t howManyPoints = point.NumberOfPoints();
    DTDoubleArray pointList = point.Data();
    int i,j;
    DTPoint2D Q;
    double dx,dy,d;
    double distSq = distance*distance;
    
    int lenOfTable = int(howManyPoints);
    DTMutableDoubleArray fromList(2,lenOfTable);
    DTMutableIntArray fromIndex(lenOfTable);
    DTMutableDoubleArray toList(2,lenOfTable);
    DTMutableIntArray toIndex(lenOfTable);
    DTMutableDoubleArray variationList(lenOfTable);
    DTMutableDoubleArray ratioList(lenOfTable);
    int posInTable = 0;
    DTMutableDoubleArray centerList(2,lenOfTable);
    DTMutableDoubleArray intervalList(lenOfTable);
    DTMutableDoubleArray startList(lenOfTable);
    DTMutableDoubleArray endList(lenOfTable);
    
    DTMutableList<DTPoint2D> finalListOfPoints(howManyPoints);
    int posInFinalList = 0;
    DTMutableIntArray addedThisPoint(howManyPoints);
    addedThisPoint = 0;

    bool addedPoint;
    DTPoint2D P;
    for (j=0;j<howManyPoints;j++) {
        DTPoint2D Q = point(j);
        addedPoint = false;
        for (i=j+1;i<howManyPoints;i++) {
            P.x = pointList(0,i);
            P.y = pointList(1,i);
            
            dx = Q.x-P.x;
            dy = Q.y-P.y;
            d = dx*dx+dy*dy;
            if (d<distSq) {
                if (posInTable==lenOfTable) {
                    fromList = IncreaseSize(fromList);
                    toList = IncreaseSize(fromList);
                    variationList = IncreaseSize(variationList);
                    ratioList = IncreaseSize(ratioList);
                    centerList = IncreaseSize(centerList);
                    intervalList = IncreaseSize(intervalList);
                    startList = IncreaseSize(startList);
                    endList = IncreaseSize(endList);
                    fromIndex = IncreaseSize(fromIndex);
                    toIndex = IncreaseSize(toIndex);
                }
                fromList(0,posInTable) = Q.x;
                fromList(1,posInTable) = Q.y;
                toList(0,posInTable) = P.x;
                toList(1,posInTable) = P.y;
                DTPoint2D c = (P+Q)/2;
                centerList(0,posInTable) = c.x;
                centerList(1,posInTable) = c.y;
                intervalList(posInTable) = posInTable+1;
                startList(posInTable) = j;
                endList(posInTable) = i;
                
                fromIndex(posInTable) = i;
                toIndex(posInTable) = j;
                addedThisPoint(i)++;
                addedThisPoint(j)++;

                // variationList(posInTable) = Variation(magnitudeChannel,grid.SpaceToGrid(P),grid.SpaceToGrid(Q));
                
                addedPoint = true;
                
                posInTable++;
            }
        }
        
        // if (addedPoint) addedThisPoint(j) = 1;
    }
    
    // Trim the result
    fromList = TruncateSize(fromList,2*posInTable);
    toList = TruncateSize(toList,2*posInTable);
    variationList = TruncateSize(variationList,posInTable);
    ratioList = TruncateSize(ratioList,posInTable);
    centerList = TruncateSize(centerList,2*posInTable);
    intervalList = TruncateSize(intervalList,posInTable);
    startList = TruncateSize(startList,posInTable);
    endList = TruncateSize(endList,posInTable);
    
    MyGroup toReturn;
    
    // Create the table process
    DTTable processTable = DTTable({
        CreateTableColumn("from",DTPointCollection2D(fromList)),
        CreateTableColumn("to",DTPointCollection2D(toList)),
        CreateTableColumn("variation",variationList),
        CreateTableColumn("center",DTPointCollection2D(centerList)),
        CreateTableColumn("interval",intervalList),
        CreateTableColumn("start",startList),
        CreateTableColumn("end",endList)
    });
    
    toReturn.process = processTable;
    
    // Prune out the rows where the variation is too large
    DTMutableIntArray rowsToUse(posInTable);
    int posInNewTable = 0;
    for (int rowNumber=0;rowNumber<posInTable;rowNumber++) {
        if (variationList(rowNumber)<threshold) {
            rowsToUse(posInNewTable++) = rowNumber;
        }
        else {
            addedThisPoint(fromIndex(rowNumber))--;
            addedThisPoint(toIndex(rowNumber))--;
        }
    }
    
    // Add the points that are not connected to any other points.
    // Leaves points that are a part of a group.
    for (i=0;i<howManyPoints;i++) {
        if (addedThisPoint(i)==0) {
            finalListOfPoints(posInFinalList++) = point(i);
        }
    }
    
    // rowsToUse is now the entries from the processing table where the variation is < threshold.
    // The next step is to remove from that list any row where the end points are singletons

    DTTable acceptedSegments = processTable.ExtractRows(TruncateSize(rowsToUse,posInNewTable));

    // Find clusters, i.e. points that are connected together in a loop.
    // I know that the start list is sorted it increasing order and start<end for every row.
    DTMutableIntArray segmentNumberForIndex(howManyPoints+1);
    segmentNumberForIndex = -1; // Not found yet
    int thisSegmentIndex = -1;

    DTTableColumnNumber start = acceptedSegments("start");
    DTTableColumnNumber end = acceptedSegments("end");
    
    DTMutableIntArray remap(howManyPoints+1);
    remap = -1;

    for (i=0;i<posInNewTable;i++) {
        int startIndex = start(i);
        int endIndex = end(i);
        if (segmentNumberForIndex(startIndex)==-1) {
            // First time I see this number
            thisSegmentIndex++;
            segmentNumberForIndex(startIndex) = thisSegmentIndex;
            remap(thisSegmentIndex) = thisSegmentIndex;
            if (segmentNumberForIndex(endIndex)==-1) {
                segmentNumberForIndex(endIndex) = thisSegmentIndex;
            }
            else {
                DTErrorMessage("");
            }
        }
        else {
            if (segmentNumberForIndex(endIndex)==-1) {
                segmentNumberForIndex(endIndex) = segmentNumberForIndex(startIndex);
            }
            else {
                // Need to remap the higher number to the smaller number. Figure out what the number ends up being
                // Know that segmentNumberForIndex(endIndex) and segmentNumberForIndex(startIndex)
                // really should point to the same number. That number should be stored in remap(lower of the two)
                if (segmentNumberForIndex(endIndex)<segmentNumberForIndex(startIndex)) {
                    DTErrorMessage("deal with this");
                }
                else if (segmentNumberForIndex(endIndex)>segmentNumberForIndex(startIndex)) {
                    DTErrorMessage("deal with this");
                }
                else { // The same, nothing to do
                    
                }
            }
        }
    }

    // Count how many are in each segment
    DTMutableIntArray howManyForSegment(thisSegmentIndex+1);
    howManyForSegment = 0;
    DTMutableDoubleArray sumOfXY(2,thisSegmentIndex+1);
    sumOfXY = 0.0;
    DTTableColumnPoint2D centerColumn = acceptedSegments("center");
    for (i=0;i<howManyPoints;i++) {
        int segmentN = segmentNumberForIndex(i);
        if (segmentN>=0) {
            howManyForSegment(segmentN)++;
            DTPoint2D c = pointsColumn(i);
            sumOfXY(0,segmentN) += c.x;
            sumOfXY(1,segmentN) += c.y;
        }
    }
    
    for (i=0;i<howManyForSegment.Length();i++) {
        finalListOfPoints(posInFinalList++) = DTPoint2D(sumOfXY(0,i),sumOfXY(1,i))/howManyForSegment(i);
    }
    
    // Now go through the connected components
    finalListOfPoints = TruncateSize(finalListOfPoints,posInFinalList);
    
    // Create the table centers
    toReturn.centers = DTTable({
        CreateTableColumn("center",finalListOfPoints)
    });

    return toReturn;
}
