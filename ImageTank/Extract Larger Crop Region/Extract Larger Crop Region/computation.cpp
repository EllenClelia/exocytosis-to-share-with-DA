#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "Utilities.h"

void Computation(const DTSet<DTImage> &images,int width,double t,
                 const DTTable &parameters,DTMutableSet<DTImage> &output)
{
    // images is a set
    ssize_t images_count = images.NumberOfItems();
    DTTableColumnNumber timeValuesOverall = images.Parameters()("t");
    int locationOfEntryInImageList = int(timeValuesOverall.FindClosest(t));
    
    DTTableColumnPoint2D center = parameters("center");
    DTTableColumnNumber timeOffset = parameters("time");
    DTDoubleArray timeOff = timeOffset.Values();
    
    DTRegion2D croppingRegion;
    
    double radius = width/2+0.01;
    bool saveSmooth = true;
    
    // Extract point numbers until I don't have any more.
    DTTableColumnNumber pointNumber = parameters("ptNumber");
    int ptN = 0;
    while (true) {
        DTIntArray where = pointNumber.FindValue(ptN);
        if (where.Length()==0) break; // All of the points found
        
        DTTable forPointNumber = parameters.ExtractRows(where);
        
        DTTableColumnPoint2D centerForPoint = forPointNumber("center");
        DTIntArray timeOffsets = ConvertToInt(forPointNumber("time")).IntVersion();
        
        // Get the information about forward and backward times from the point information
        int timeback = -timeOffsets(0);
        int timeforward = timeOffsets(timeOffsets.Length()-1);
        DTMutableList<DTImage> smallImagesSmooth(timeback+timeforward+1);
        DTMutableList<DTImage> smallImagesRaw(timeback+timeforward+1);
        DTMutableList<DTImage> backgroundImages(timeback-1);
        
        DTPoint2D c = centerForPoint(0);

        ssize_t pos = 0;
        ssize_t posOfWhere = 0;
        for (int index=locationOfEntryInImageList-timeback;index<=locationOfEntryInImageList+timeforward;index++) {
            if (index>=0 && index<images_count) {
                
                if (index==locationOfEntryInImageList) { // Need to know where the starting image is in the sequence.
                    posOfWhere = pos;
                }
                
                DTImage image = images(index);
                double h = image.Grid().dx();
                croppingRegion = DTRegion2D(c.x-radius*h,c.x+radius*h,
                                            c.y-radius*h,c.y+radius*h);
                
                image = Crop(image,croppingRegion);
                image = ConvertToDouble(image);
                
                smallImagesRaw(pos) = image;
                smallImagesSmooth(pos++) = GaussianFilter(image,1);
            }
        }
        
        pos = 0;
        for (int index=locationOfEntryInImageList-timeback;index<locationOfEntryInImageList-1;index++) {
            if (index>=0 && index<images_count) {
                backgroundImages(pos) = smallImagesSmooth(pos);
                pos++;
            }
        }
        
        DTImage backgroundImage;
        if (pos<timeback-1) {
            backgroundImage = MedianOfImages(TruncateSize(backgroundImages,pos));
        }
        else {
            backgroundImage = MedianOfImages(backgroundImages);
        }

        // Now save the images, need to match the order for the original table
        for (pos=0;pos<timeOffsets.Length();pos++) {
            DTImage singleImageRaw = smallImagesRaw(pos);
            DTImage singleImageSmooth = smallImagesSmooth(pos);

            // Subtract the background and put the result into a two channel image
            DTImage diffSmooth = singleImageSmooth - backgroundImage;
            DTImage diffRaw = singleImageRaw - backgroundImage;

            // Combine the channels, use different names
            DTImage combinedSmooth(diffSmooth.Grid(),{ChangeName(singleImageSmooth(0),"intensity"),ChangeName(diffSmooth(0),"difference")});
            DTImage combinedRaw(diffRaw.Grid(),{ChangeName(singleImageRaw(0),"intensity"),ChangeName(diffRaw(0),"difference")});
            
            if (saveSmooth) {
                output.Add(combinedSmooth);
            }
            else {
                output.Add(combinedRaw);
            }
        }
                
        ptN++;
    }
    
    // This call finishes the creation of the set variable.
    output.Finish(parameters);
}
