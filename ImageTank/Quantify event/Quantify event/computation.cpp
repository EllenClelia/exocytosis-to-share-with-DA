#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTFunctionFit.h"
#include "DTDictionary.h"
#include "DTFunction1D.h"
#include "Utilities.h"

Group Computation(const DTSet<DTImage> &images)
{
    // images is a set
    ssize_t images_count = images.NumberOfItems();
    
    DTTable images_par = images.Parameters();
    DTTableColumnNumber ptNumber = images_par("ptNumber");
    
    // Just extract point == 0
    ssize_t i;
    for (i=0;i<images_count;i++) {
        if (ptNumber(i)!=0) break;
    }
    
    QuantifyEvent event = Quantify(images.ExtractRows(DTRange(0,i)));
 
    Group toReturn;
    toReturn.average = event.average;
    toReturn.decay = event.decay;
    toReturn.width = event.width;
    toReturn.histogram = event.histogram;
    DTFunction1D xv = DTFunction1D::x();
    double shift = event.shift;
    DTFunction1D fitFcn;
    if (shift==0) {
        fitFcn = event.base + event.spike*exp(-event.decay*xv);
    }
    else if (shift>0) {
        fitFcn = event.base + event.spike*exp(-event.decay*(xv-shift));
    }
    else {
        fitFcn = event.base + event.spike*exp(-event.decay*(xv+(-shift)));
    }
    
    toReturn.fit = fitFcn;
    toReturn.R2 = event.R2;
    toReturn.shift = event.shift;
    
    return toReturn;
}

