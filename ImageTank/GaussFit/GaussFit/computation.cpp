#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTProgress.h"
#include "DTFloatArrayOperators.h"

#include "Utilities.h"

void DoG(const DTImage &image,double sigma,int inOctave,int numOctaves,DTMutableSet<DTImage> &output)
{
    double mult = pow(2,1.0/inOctave);
    
    int octaveNumber;
    int stepNumber;
    DTMutableDoubleArray sigmaColumn(inOctave*numOctaves);
    DTMutableDoubleArray octaveColumn(inOctave*numOctaves);
    int pos = 0;
    
    DTImage imageToUse = ConvertToFloat(image);
    
    DTProgress progress;
    
    DTImage current,next;
    next = GaussianFilter(imageToUse,sigma);
    
    DTFloatArray difference;
    
    ssize_t m = image.m();
    ssize_t n = image.n();
    DTMutableFloatArray stack(m,n,inOctave*numOctaves);
    
    DTMutableList<DTImageChannel> outputImageChannels(3);

    for (octaveNumber=0;octaveNumber<numOctaves;octaveNumber++) {
        for (stepNumber=0;stepNumber<inOctave;stepNumber++) {
            // Save the images to the set.
            // Add this to the list (next entry)

            current = next;

            // Save the current smoothing level as the first channel
            outputImageChannels(0) = current(0); // Channel 1
            sigmaColumn(pos) = sigma;
            octaveColumn(pos) = (pos/inOctave)+1;
            
            // Compute the next smoothing level and find the difference
            sigma *= mult;
            next = GaussianFilter(imageToUse,sigma);
            difference = current(0).FloatArray()-next(0).FloatArray();
            outputImageChannels(1) = DTImageChannel("difference",difference); // Channel 2
            CopyIntoSlice(stack,difference,pos);
            
            // Third channel is the median up to this point
            DTFloatArray median = MedianOfStack(stack,pos+1);
            outputImageChannels(2) = DTImageChannel("median",median);

            // Create and save the image as part of the set
            DTImage imageToAdd(image.Grid(),outputImageChannels);
            output.Add(imageToAdd);

            pos++;
            
            progress.UpdatePercentage(pos/double(inOctave*numOctaves));
        }
    }
    
    // Add the meta data
    DTMutableList<DTTableColumn> columns(2);
    columns(0) = CreateTableColumn("sigma",sigmaColumn);
    columns(1) = CreateTableColumn("octave",octaveColumn);
    DTTable parameterTable(columns);
    output.Finish(parameterTable);
}
