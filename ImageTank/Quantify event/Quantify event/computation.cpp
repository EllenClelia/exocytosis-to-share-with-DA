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
    int channel = parameters("channel");
    
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

    QuantifyEvent event = Quantify(imagesToView,channel);
 
    Group toReturn;
    toReturn.average = event.average;
    toReturn.decay = event.decay;
    toReturn.width = event.width;
    toReturn.histogram = event.histogram;
    DTFunction xv = DTFunction::Constant("x");
    double shift = event.shift;
    DTFunction1D fitFcn;
    if (shift==0) {
        fitFcn = event.base + event.spike*exp(-event.decay*xv);
    }
    else {
        fitFcn = IfFunction(xv<shift,DTFunction::Value(event.base+event.spike),event.base+event.spike*exp(-event.decay*(xv-shift)));
    }
    
    toReturn.fit = fitFcn;
    toReturn.R2 = event.R2;
    toReturn.shift = event.shift;
    
    // Create the table Drift
    DTTable eventParameters = imagesToView.Parameters();
    DTTableColumnNumber time = eventParameters("time");
    ssize_t startingIndex = time.FindClosest(0)+event.shift;
    DTTable tail = eventParameters.ExtractRows(DTRange(startingIndex,eventParameters.NumberOfRows()-startingIndex));
    
    toReturn.Drift = DTTable({
        tail("time"),
        tail("centerSpot")
    });
    
    
    toReturn.Piecewise_Fit = event.piecewiseFitResults;

    return toReturn;
}

