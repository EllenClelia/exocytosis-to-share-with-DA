#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTTable &events,double timeDelta,double spaceDelta)
{
    int pos = 0;
    ssize_t howMany = events.NumberOfRows();
    
    DTTableColumnPoint2D centerColumn = events("center");
    DTTableColumnNumber timeColumn = events("time");

    DTMutableCharArray includeRow(howMany);
    includeRow = 1;
    ssize_t secondIndex;
    ssize_t includeCount = 0;
    
    for (pos=0;pos<howMany;pos++) {
        double tv = timeColumn(pos);
        DTPoint2D pt = centerColumn(pos);
        if (includeRow(pos)) includeCount++;
        for (secondIndex=pos+1;secondIndex<howMany;secondIndex++) {
            double tComp = timeColumn(secondIndex);
            if (tComp>tv+timeDelta) break; // too far out
            DTPoint2D ptComp = centerColumn(secondIndex);
            if (Distance(ptComp,pt)<spaceDelta) {
                includeRow(secondIndex) = 0;
            }
        }
    }
    
    DTMutableIntArray selectRows(includeCount);
    ssize_t posInOutput = 0;
    for (pos=0;pos<howMany;pos++) {
        if (includeRow(pos)) selectRows(posInOutput++) = pos;
    }
    
    return events.ExtractRows(selectRows);
 }
