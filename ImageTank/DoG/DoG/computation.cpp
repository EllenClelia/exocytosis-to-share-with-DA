#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTFloatArrayOperators.h"
#include "DTUtilities.h"
#include "Utilities.h"
#include "DTProgress.h"

DTImage Computation(const DTImage &image,double sigma,int inOctave,
                    int numOctaves)
{
    double mult = pow(2,1.0/inOctave);
    
    int octaveNumber;
    int stepNumber;
    
    DTImage imageToUse = ConvertToFloat(image);
    
    DTProgress progress;
    
    DTImage current,next;
    next = GaussianFilter(imageToUse,sigma);
    
    DTFloatArray difference;
    
    ssize_t m = image.m();
    ssize_t n = image.n();
    DTMutableFloatArray stack(m,n,inOctave*numOctaves);
    
    DTMutableList<DTImageChannel> outputImageChannels(3);

    int pos = 0;
    for (octaveNumber=0;octaveNumber<numOctaves;octaveNumber++) {
        for (stepNumber=0;stepNumber<inOctave;stepNumber++) {
            current = next;

            // Compute the next smoothing level and find the difference
            sigma *= mult;
            next = GaussianFilter(imageToUse,sigma);
            difference = current(0).FloatArray()-next(0).FloatArray();
            CopyIntoSlice(stack,difference,pos);
            
            pos++;
            
            progress.UpdatePercentage(pos/double(inOctave*numOctaves));
        }
    }

    DTFloatArray median = MedianOfStack(stack,pos);

    return DTImage(image.Grid(),median,"median");
}
