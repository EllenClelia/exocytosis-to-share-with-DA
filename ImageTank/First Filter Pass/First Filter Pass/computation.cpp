#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

#include "Utilities.h"

DTTable Computation(const DTSet<DTImage> &images,
                    double fromBkgrnd, // How much it has to rise from background
                    double tailThreshold, // When to stop tracking
                    double maxDrift) // Maximum Drift
{
    // images is a set
    ssize_t images_count = images.NumberOfItems();
    DTTable images_par = images.Parameters();
    
    DTTableColumnNumber T = images_par("T");
    DTTableColumnNumber ptNumber = images_par("ptNumber");

    // Know that the T column is sorted in increasing order
    // inside that the ptNumber is increasing
    // and inside that time is increasing
    ssize_t i;
    
    ssize_t outputLineLength = images_count/20; // Don't know how many, but this is a safe upper bound
    DTMutableDoubleArray outputTime(outputLineLength);
    DTMutableDoubleArray outputShift(outputLineLength);
    DTMutableDoubleArray outputDecay(outputLineLength);
    DTMutableDoubleArray outputR2(outputLineLength);
    DTMutableDoubleArray outputBackground(outputLineLength);
    DTMutableDoubleArray outputWidth(outputLineLength);
    DTMutableDoubleArray outputSkip(outputLineLength);
    DTMutableDoubleArray outputFlag(outputLineLength);
    DTMutableDoubleArray outputCenter(2,outputLineLength);
    DTMutableDoubleArray outputDrift(outputLineLength);
    ssize_t posInOutput = 0;
    
    outputTime = NAN;
    outputShift = NAN;
    outputDecay = NAN;
    outputR2 = NAN;
    outputBackground = NAN;
    outputWidth = NAN;
    outputSkip = 0;
    outputFlag = 0;
    outputCenter = NAN;
    outputDrift = NAN;

    // Loop through each event
    ssize_t startsAt = 0;
    while (startsAt<images_count) {
        
        int ptN = ptNumber(startsAt);
        double Tval = T(startsAt);
        i = startsAt;
        while (i<images_count && ptNumber(i)==ptN && T(i)==Tval) {
            i++;
        }
        ssize_t endsAt = i;
        
        outputTime(posInOutput) = Tval;
        outputShift(posInOutput) = 0;
        
        // startsAt<=i<endsAt is one event.
        DTSet<DTImage> event = images.ExtractRows(DTRange(startsAt,endsAt-startsAt));
        // DTSet<DTImage> subTable = images.ExtractRows(DTRange(80,100));

        QuantifyEvent info = Quantify(event);
        
        // See if the previous time value has a bigger intensity. In which case go back a step.
        
        outputBackground(posInOutput) = info.average;
        outputWidth(posInOutput) = info.width;
        outputShift(posInOutput) = info.shift;
        outputDecay(posInOutput) = info.decay;
        outputR2(posInOutput) = info.R2;

        DTTable eventParameters = event.Parameters();
        ssize_t lengthOfEvent = eventParameters.NumberOfRows();

        // Find time==0, i.e. the starting point
        DTTableColumnNumber time = eventParameters("time");
        DTTableColumnNumber average = eventParameters("average");
        DTTableColumnNumber intensity = eventParameters("intensity");
        DTTableColumnPoint2D centerSpot = eventParameters("centerSpot");
        ssize_t startingIndex = time.FindClosest(0)+info.shift;
        
        double valueAtStart = intensity(startingIndex);
        double valueAtNext = intensity(startingIndex+1);
        
        if (valueAtNext>valueAtStart) {
            // Still going up, need to look at this further
            outputFlag(posInOutput) += 1;
        }
        
        // Make sure that the initial spike is high enough from the background
        if (valueAtStart<info.average + info.width*fromBkgrnd) {
            outputSkip(posInOutput) += 1; // Didn't rise high enough from the background
        }
        
        // Find when the spike goes below the tail intensity. Also check how much movement happened
        double drift = 0;
        ssize_t tailEndsAt = startingIndex;
        DTPoint2D startAt = centerSpot(startingIndex);
        outputCenter(0,posInOutput) = startAt.x;
        outputCenter(1,posInOutput) = startAt.y;
        while (tailEndsAt<lengthOfEvent && intensity(tailEndsAt)>info.average+info.width*tailThreshold) {
            double dist = Distance(startAt,centerSpot(tailEndsAt));
            if (dist>drift) drift = dist;
            tailEndsAt++;
        }
        outputDrift(posInOutput) = drift;
        if (drift>maxDrift) {
            outputSkip(posInOutput) += 2; // Drifted too far
        }
        if (startingIndex+1==tailEndsAt) {
            // Only the first point was above the threshold
            outputSkip(posInOutput) += 4;
        }
                                
        // Ready for the next point
        startsAt = endsAt;
        
        posInOutput++;
    }
    
    outputTime = TruncateSize(outputTime,posInOutput);
    outputShift = TruncateSize(outputShift,posInOutput);
    outputDecay = TruncateSize(outputDecay,posInOutput);
    outputR2 = TruncateSize(outputR2,posInOutput);
    outputBackground = TruncateSize(outputBackground,posInOutput);
    outputWidth = TruncateSize(outputWidth,posInOutput);
    outputSkip = TruncateSize(outputSkip,posInOutput);
    outputFlag = TruncateSize(outputFlag,posInOutput);
    outputCenter = TruncateSize(outputCenter,2*posInOutput);
    outputDrift = TruncateSize(outputDrift,posInOutput);
    
    return DTTable({
        CreateTableColumn("time",outputTime),
        CreateTableColumn("shift",outputShift),
        CreateTableColumn("skip",outputSkip),
        CreateTableColumn("center",DTPointCollection2D(outputCenter)),
        CreateTableColumn("R2",outputR2),
        CreateTableColumn("drift",outputDrift),
        CreateTableColumn("background",outputBackground),
        CreateTableColumn("decay",outputDecay),
        CreateTableColumn("width",outputWidth),
        CreateTableColumn("flag",outputFlag)
    });
}
