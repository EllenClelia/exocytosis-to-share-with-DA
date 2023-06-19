#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTFunctionFit.h"
#include "DTDictionary.h"
#include "DTFunction1D.h"
#include "Utilities.h"

Group Computation(const DTSet<DTImage> &images,int pt,
                  const DTDictionary &parameters)
{
    // images is a set
    ssize_t images_count = images.NumberOfItems();
    // int channel = parameters("channel");
    
    DTTable images_par = images.Parameters();
    DTTableColumnNumber ptNumber = images_par("ptNumber");
    
    // Just extract point == pt
    ssize_t startAt = 0;
    for (startAt=0;startAt<images_count;startAt++) {
        if (ptNumber(startAt)>=pt) break;
    }
    
    ssize_t i;
    for (i=startAt;i<images_count;i++) {
        if (ptNumber(i)!=pt) break;
    }
    ssize_t endAt = i;
    
    if (startAt==endAt) {
        // Empty
        return Group();
    }
    
    DTSet<DTImage> imagesToView = images.ExtractRows(DTRange(startAt,endAt-startAt));

    QuantifyEvent event = Quantify(imagesToView,parameters);
 
    Group toReturn;
    toReturn.average = event.average;
    toReturn.delay = event.delay;
    toReturn.decay = event.decay;
    toReturn.width = event.width;
    toReturn.histogram = event.histogram;
    DTFunction xv = DTFunction::Constant("x");
    double shift = event.shift;
    DTFunction1D fitFcn;
    
    double kinkAt = event.shift+event.delay;
    fitFcn = IfFunction(xv<kinkAt,DTFunction::Value(event.base+event.spike),event.base+event.spike*exp(-event.decay*(xv-kinkAt)));
    
    toReturn.fit = fitFcn;
    toReturn.R2 = event.R2;
    toReturn.shift = event.shift;
    
    // Create the table Drift
    DTTable eventParameters = imagesToView.Parameters();
    DTTableColumnNumber time = eventParameters("time");
    ssize_t startingIndex = time.FindClosest(event.shift);
    // The drift calculation can start before the startingIndex
    DTTableColumnNumber failure = eventParameters("failure");
    
    int howManyFailuresToAllowInDrift = parameters.GetNumber("Allow drift failures",0);
    bool checkDriftBefore = parameters.GetNumber("drift",0);
    
    DTTable driftPortion;
    if (howManyFailuresToAllowInDrift==0) {
        if (checkDriftBefore) {
            while (startingIndex-1>=0 && failure(startingIndex-1)==0) startingIndex--;
        }
        DTTable tail = eventParameters.ExtractRows(DTRange(startingIndex,eventParameters.NumberOfRows()-startingIndex));

        // The drift should only be the points until the first failure
        failure = tail("failure");
        int firstFailure = 1;
        while (firstFailure<tail.NumberOfRows() && failure(firstFailure)==0) firstFailure++;
        driftPortion = tail.ExtractRows(DTRange(0,firstFailure));
    }
    else {
        // Extract a subset of the eventParameters as follows:
        // if checkDriftBefore==true, then go back until I have a point that fails.
        // Go forward in time until I reach the end or have more than "howManyFailuresToAllowInDrift" failure points.
        DTMutableIntArray whichToTake(eventParameters.NumberOfRows());
        int posInWhichToTake = 0;
        int startAt = int(startingIndex);
        while (startAt-1>=0 && failure(startAt-1)==0) {
            startAt--;
        }
        while (startAt<startingIndex) {
            whichToTake(posInWhichToTake++) = startAt;
            startAt++;
        }
        
        int howManyFailedSoFar = 0;
        int check = int(startingIndex);
        while (check<failure.NumberOfRows() && howManyFailedSoFar<=howManyFailuresToAllowInDrift) {
            if (failure(check)==0) {
                whichToTake(posInWhichToTake++) = check;
            }
            else {
                howManyFailedSoFar++;
            }
            check++;
        }
        whichToTake = TruncateSize(whichToTake,posInWhichToTake);
        driftPortion = eventParameters.ExtractRows(whichToTake);
    }
        
    toReturn.Drift = DTTable({
        driftPortion("time"),
        driftPortion("centerSpot")
    });
    
    
    toReturn.Piecewise_Fit = event.piecewiseFitResults;

    return toReturn;
}

