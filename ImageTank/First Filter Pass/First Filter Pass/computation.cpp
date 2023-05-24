#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

#include "Utilities.h"

DTTable Computation(const DTSet<DTImage> &images,
                    double fromBkgrnd, // How much it has to rise from background
                    double tailThreshold, // When to stop tracking
                    double maxDrift, // Maximum Drift
                    int stepsForDrift,const DTDictionary &parameters)
{
    // images is a set
    ssize_t images_count = images.NumberOfItems();
    DTTable images_par = images.Parameters();
    
    DTTableColumnNumber T = images_par("T");
    DTTableColumnNumber ptNumber = images_par("ptNumber");
    // int channel = parameters("channel");
    bool useAverage = parameters("useAverage");
    int howManyRequiredAboveThresholdForDecay = parameters("Points Above Baseline for threshold");

    // Know that the T column is sorted in increasing order
    // inside that the ptNumber is increasing
    // and inside that time is increasing
    ssize_t i;
    
    ssize_t outputLineLength = images_count/20; // Don't know how many, but this is a safe upper bound
    DTMutableDoubleArray outputTime(outputLineLength);
    DTMutableDoubleArray outputShift(outputLineLength);
    DTMutableDoubleArray outputFlag(outputLineLength);
    DTMutableDoubleArray outputDelay(outputLineLength);
    DTMutableDoubleArray outputDecay(outputLineLength);
    DTMutableDoubleArray outputR2(outputLineLength);
    DTMutableDoubleArray outputBackground(outputLineLength);
    DTMutableDoubleArray outputWidth(outputLineLength);
    DTMutableDoubleArray outputCenter(2,outputLineLength);
    DTMutableDoubleArray outputDrift(outputLineLength);
    ssize_t posInOutput = 0;
    
    outputTime = NAN;
    outputShift = NAN;
    outputDecay = NAN;
    outputR2 = NAN;
    outputBackground = NAN;
    outputWidth = NAN;
    outputFlag = 0;
    outputCenter = NAN;
    outputDrift = NAN;

    DTProgress progress;
    
    DTMutableList<DTPoint2D> tempPoints(100);
    int posInTempPoints;
    
    // Loop through each event
    ssize_t startsAt = 0;
    while (startsAt<images_count) {
        progress.UpdatePercentage(startsAt/double(images_count));
        
        int ptN = ptNumber(startsAt);
        double Tval = T(startsAt);
        i = startsAt;
        while (i<images_count && ptNumber(i)==ptN && T(i)==Tval) {
            i++;
        }
        ssize_t endsAt = i;
        
        // startsAt<=i<endsAt is one event.
        DTSet<DTImage> event = images.ExtractRows(DTRange(startsAt,endsAt-startsAt));
        // DTSet<DTImage> subTable = images.ExtractRows(DTRange(80,100));

        QuantifyEvent info = Quantify(event,parameters);
        
        // See if the previous time value has a bigger intensity. In which case go back a step.
        outputTime(posInOutput) = Tval;
        
        outputBackground(posInOutput) = info.average;
        outputWidth(posInOutput) = info.width;
        outputShift(posInOutput) = info.shift;
        outputDelay(posInOutput) = info.delay;
        outputDecay(posInOutput) = info.decay;
        outputR2(posInOutput) = info.R2;

        DTTable eventParameters = event.Parameters();
        ssize_t lengthOfEvent = eventParameters.NumberOfRows();

        // Find time==0, i.e. the starting point
        DTTableColumnNumber time = eventParameters("time");
        DTTableColumnNumber average = eventParameters("average");
        DTTableColumnNumber failure = eventParameters("failure");
        DTTableColumnNumber intensityToUse;
        if (useAverage) {
            intensityToUse = eventParameters("average");
        }
        else {
            intensityToUse = eventParameters("intensity");
        }
        DTTableColumnNumber fittedIntensity = eventParameters("intensity");

        DTTableColumnPoint2D centerSpot = eventParameters("centerSpot");
        ssize_t startingIndex = time.FindClosest(info.shift); // Start where the event starts, not where the decay starts
        
        double valueAtStart = intensityToUse(startingIndex);
        // double valueAtNext = intensityToUse(startingIndex+1);
        //if (valueAtNext>valueAtStart) {
        //    // Still going up, need to look at this further
        //    outputFlag(posInOutput) += 1;
        //}
        
        // Make sure that the initial spike is high enough from the background
        if (valueAtStart<info.average + info.width*fromBkgrnd) {
            outputFlag(posInOutput) += 1; // Didn't rise high enough from the background
        }
        
        DTPoint2D startAt = centerSpot(startingIndex);
        outputCenter(0,posInOutput) = startAt.x;
        outputCenter(1,posInOutput) = startAt.y;
        

        // On 5/18/23, changed the drift to be the biggest distance between points
        // around startingIndex where failure(ptIndex)==0.
        // Do that by first computing that point list, and then find the maximum
        // pairwise distance beween the points.
        posInTempPoints = 0;
        
        DTMutableList<DTPoint2D> tempPoints(100);
        ssize_t lookAtPoint = startingIndex;
        if (failure(lookAtPoint)==0 && failure(lookAtPoint-1)==0) {
            lookAtPoint--;
            while (lookAtPoint>=0 && failure(lookAtPoint)==0) {
                tempPoints(posInTempPoints++) = centerSpot(lookAtPoint);
                lookAtPoint--;
            }
        }
        // Forwards
        ssize_t stopSearchingDrift = startingIndex+stepsForDrift;
        lookAtPoint = startingIndex;
        while (lookAtPoint<lengthOfEvent && lookAtPoint<stopSearchingDrift && failure(lookAtPoint)==0) {
            tempPoints(posInTempPoints++) = centerSpot(lookAtPoint);
            lookAtPoint++;
        }
        // Pairwise distance
        double drift = 0;
        for (int firstPt=0;firstPt<posInTempPoints;firstPt++) {
            for (int secondPt=firstPt+1;secondPt<posInTempPoints;secondPt++) {
                double dist = Distance(tempPoints(firstPt),tempPoints(secondPt));
                if (dist>drift) drift = dist;
            }
        }
        
        /*
         Before 5/18/23 we computed the distance from the centerSpot(startingIndex)
        // Find when the spike goes below the tail intensityToUse. Also check how much movement happened
        double drift = 0;
        // Find the drift
        ssize_t stopSearchingDrift = startingIndex+stepsForDrift;
        ssize_t lookAtPoint = startingIndex;
        
        // Compute the drift from previous values
        if (failure(lookAtPoint)==0 && failure(lookAtPoint-1)==0) {
            // The forward drift calculation fails in the first two steps of the next loop
            // Check to see if there was a valid point before it.
            ssize_t check = lookAtPoint-1;
            while (check>=0 && failure(check)==0) {
                double dist = Distance(startAt,centerSpot(check));
                if (dist>drift) drift = dist;
                check--;
            }
        }

        while (lookAtPoint<lengthOfEvent &&
               lookAtPoint<stopSearchingDrift &&
               failure(lookAtPoint)==0
               // The next line was commented out before April 5th 2023
               // it was put back in because when the intensity dips too low the drift is messed up.
               // On May 18th 23 we commented this out again
               // intensityToUse(lookAtPoint)>info.average+info.width*tailThreshold
               ) {
            // The center point is still valid, and intensity is still large enough, check the drift.
            double dist = Distance(startAt,centerSpot(lookAtPoint));
            if (dist>drift) drift = dist;
            lookAtPoint++;
        }
         
         */
        
        outputDrift(posInOutput) = drift;
        if (drift>maxDrift) {
            outputFlag(posInOutput) += 2; // Drifted too far
        }
        
        // See if the intensity decays too fast
        ssize_t tailEndsAt = startingIndex;
        while (tailEndsAt<lengthOfEvent &&
               tailEndsAt<stopSearchingDrift &&
               (useAverage || failure(tailEndsAt)==0) &&
               intensityToUse(tailEndsAt)>info.average+info.width*tailThreshold) {
            tailEndsAt++;
        }
        
        if (tailEndsAt<=startingIndex+howManyRequiredAboveThresholdForDecay) {
            // Only the first point was above the threshold
            outputFlag(posInOutput) += 4;
        }
        
        if (failure(startingIndex)!=0) {
            // Require at least the first point to look like a peak.
            outputFlag(posInOutput) += 8;
        }
                                
        // Ready for the next point
        startsAt = endsAt;
        
        posInOutput++;
    }
    
    outputTime = TruncateSize(outputTime,posInOutput);
    outputShift = TruncateSize(outputShift,posInOutput);
    outputDelay = TruncateSize(outputDelay,posInOutput);
    outputDecay = TruncateSize(outputDecay,posInOutput);
    outputR2 = TruncateSize(outputR2,posInOutput);
    outputBackground = TruncateSize(outputBackground,posInOutput);
    outputWidth = TruncateSize(outputWidth,posInOutput);
    outputFlag = TruncateSize(outputFlag,posInOutput);
    outputCenter = TruncateSize(outputCenter,2*posInOutput);
    outputDrift = TruncateSize(outputDrift,posInOutput);
    
    return DTTable({
        CreateTableColumn("time",outputTime),
        CreateTableColumn("shift",outputShift),
        CreateTableColumn("flag",outputFlag),
        CreateTableColumn("center",DTPointCollection2D(outputCenter)),
        CreateTableColumn("R2",outputR2),
        CreateTableColumn("drift",outputDrift),
        CreateTableColumn("background",outputBackground),
        CreateTableColumn("delay",outputDelay),
        CreateTableColumn("decay",outputDecay),
        CreateTableColumn("width",outputWidth)
    });
}
