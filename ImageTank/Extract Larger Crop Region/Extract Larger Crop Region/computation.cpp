#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

void Computation(const DTSet<DTImage> &images,const DTTable &parameters,
                 int width,double t,DTMutableSet<DTImage> &output)
{
    // images is a set
    ssize_t images_count = images.NumberOfItems();
    DTTableColumnNumber timeValuesOverall = images.Parameters()("t");
    int locationOfEntryInImageList = int(timeValuesOverall.FindClosest(t));
    
    // t
    //
    //   DTDoubleArray t = values.DoubleVersion();
    //   double single = values(3); // Higher overhead, simpler syntax.
    
    DTTableColumnPoint2D center = parameters("center");
    DTTableColumnNumber timeOffset = parameters("time");
    DTDoubleArray timeOff = timeOffset.Values();
    
    DTRegion2D croppingRegion;
    
    double radius = width/2+0.01;

    ssize_t numberOfEntries = parameters.NumberOfRows();
    for (ssize_t entryNumber=0;entryNumber<numberOfEntries;entryNumber++) {
        // Read the image
        int offset = timeOff(entryNumber);
        int locationInImages = locationOfEntryInImageList + offset;
        DTImage imageToCrop = images(locationInImages);
        
        // Find the region
        double h = imageToCrop.Grid().dx();
        DTPoint2D c = center(entryNumber);
        croppingRegion = DTRegion2D(c.x-radius*h,c.x+radius*h,
                                    c.y-radius*h,c.y+radius*h);
        
        DTImage croppedImage = Crop(imageToCrop,croppingRegion);
        output.Add(croppedImage);
    }
    
    // This call finishes the creation of the set variable.
    output.Finish(parameters);
}
