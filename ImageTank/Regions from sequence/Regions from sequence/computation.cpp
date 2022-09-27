#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

#include "Utilities.h"

void Computation(const DTSet<DTImage> &everything,const DTTable &spots,
                 double time,int timeback,int timeforward,
                 const DTDictionary &parameters,DTMutableSet<DTImage> &output)
{
    DTSet<DTImage> withCache = everything.WithCache(timeback+1+timeforward);
    DTTable imageParameters = withCache.Parameters();
    DTTableColumnNumber timeValues = imageParameters("t");
    int channel = parameters("channel");
    
    DTTableColumnNumber timeColumn;
    bool timeColumnExists = spots.Contains("time");
    if (timeColumnExists) {
        timeColumn = spots("time");
    }
    
    DTTableColumnPoint2D center = spots("center");
    DTTableColumnNumber width = spots("width");
    
    ssize_t row,howMany = center.NumberOfRows();
    ssize_t index;
    
    ssize_t posInOutput = 0;
    DTMutableDoubleArray timeList(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray fitWorked(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray intensityList(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray centerList(2,howMany*(timeback+1+timeforward));
    DTMutableDoubleArray pointNumber(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray centerSpot(2,howMany*(timeback+1+timeforward));
    DTMutableDoubleArray averageValues(howMany*(timeback+1+timeforward));
    // double maxV = 0;
    // DTPoint2D maxP;
    
    DTMutableList<DTImage> smallImagesSmooth(timeback+timeforward+1);
    DTMutableList<DTImage> smallImagesRaw(timeback+timeforward+1);
    DTMutableList<DTImage> backgroundImages(timeback-1);
    
    DTImage backgroundImage;
    
    bool peakFit = ((parameters.GetNumber("peak",1))==1);
    bool analyzeSmooth = ((parameters.GetNumber("smooth",1))==1);
    bool saveSmooth = ((parameters.GetNumber("save",1))==1);
    LocalPeak peak;
    
    DTProgress progress;
    for (row=0;row<howMany;row++) {
        progress.UpdatePercentage(row/double(howMany));

        if (timeColumnExists) {
            time = timeColumn(row);
        }
        
        ssize_t where = timeValues.FindClosest(time);
        DTPoint2D p = center(row);
        // The center value is not a good initial guess
        // Need to correct that by taking a window around the current center
        // find the maximum point of that window
        // and use that for the back/forward points - Crop(image,box) a few lines below
        
        double w = width(row);
        DTRegion2D finalBox = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2); // The final image that is saved
        DTRegion2D backgroundBox; // The image where I compute the smoothing and background. The finalBox needs to live inside this box
        
        // New method. Take a bigger box
        backgroundBox = DTRegion2D(p.x-w,p.x+w,p.y-w,p.y+w);

        DTPoint2D startingPoint = p;
        
        // Compute all of the raw smooth images
        ssize_t pos = 0;
        ssize_t posOfWhere = 0;
        for (index=where-timeback;index<=where+timeforward;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {

                if (index==where) { // Need to know where the starting image is in the sequence.
                    posOfWhere = pos;
                }
                
                DTImage image = withCache(index);
                
                image = Crop(image,backgroundBox);
                image = ConvertToDouble(image);

                smallImagesRaw(pos) = image;
                smallImagesSmooth(pos++) = GaussianFilter(image,1);
            }
        }
        
        // Compute a background image using a median of the smooth images.
        pos = 0;
        for (index=where-timeback;index<where-1;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {
                backgroundImages(pos) = smallImagesSmooth(pos);
                pos++;
            }
        }
        
        if (pos<timeback-1) {
            backgroundImage = MedianOfImages(TruncateSize(backgroundImages,pos));
        }
        else {
            backgroundImage = MedianOfImages(backgroundImages);
        }
        
        // The center
        
        //Find the maxima of this image inside the middle of the big box
        //for either channel 0 or channel 1 depending on if I want the background removed.
        // Need to tweak the finalBox so that it is centered around the maxima at where
        DTImage startingImage = smallImagesSmooth(posOfWhere);
        DTImage diff = startingImage - backgroundImage;
        DTImage combined(diff.Grid(),{ChangeName(startingImage(0),"intensity"),ChangeName(diff(0),"difference")});
        
        //DTDataFile temp("/tmp/test.dtbin",DTFile::NewReadWrite);
        //WriteOne(temp, "Combined", combined);
        
        combined = Crop(combined,finalBox);
        //WriteOne(temp, "Cropped", combined);
        
        if (peakFit) {
            peak = FindGaussianPeak(combined,channel);
        }
        else {
            peak = FindMaximumPeak(combined,channel);
        }
        p = peak.center;
        
        startingPoint = p;
        finalBox = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2); // The final image that is saved
        
        pos = 0;
        for (index=where-timeback;index<=where+timeforward;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {
                DTImage singleImageRaw = smallImagesRaw(pos);
                DTImage singleImageSmooth = smallImagesSmooth(pos);

                // Subtract the background and put the result into a two channel image
                DTImage diffSmooth = singleImageSmooth - backgroundImage;
                DTImage diffRaw = singleImageRaw - backgroundImage;

                // Combine the channels, use different names
                DTImage combinedSmooth(diff.Grid(),{ChangeName(singleImageSmooth(0),"intensity"),ChangeName(diffSmooth(0),"difference")});
                combinedSmooth = Crop(combinedSmooth,finalBox);
                DTImage combinedRaw(diff.Grid(),{ChangeName(singleImageRaw(0),"intensity"),ChangeName(diffRaw(0),"difference")});
                combinedRaw = Crop(combinedRaw,finalBox);
                
                if (saveSmooth) {
                    output.Add(combinedSmooth);
                }
                else {
                    output.Add(combinedRaw);
                }

                if (peakFit) {
                    if (analyzeSmooth) {
                        peak = FindGaussianPeak(combinedSmooth,channel);
                        if (peak.worked==false) {
                            //DTDataFile temp("/tmp/test.dtbin",DTFile::NewReadWrite);
                            //WriteOne(temp, "Combined", combinedSmooth);
                            peak = FindMaximumPeak(combinedSmooth,channel);
                            peak.worked = false;
                        }
                    }
                    else {
                        peak = FindGaussianPeak(combinedRaw,channel);
                        if (peak.worked==false) {
                            peak = FindMaximumPeak(combinedSmooth,channel);
                            peak.worked = false;
                        }
                        if (peak.worked==false) {
                            peak = FindMaximumPeak(combinedSmooth,channel);
                            peak.worked = false;
                        }
                    }
                }
                else {
                    if (analyzeSmooth) {
                        peak = FindMaximumPeak(combinedSmooth,channel);
                    }
                    else {
                        peak = FindMaximumPeak(combinedRaw,channel);
                    }
                }
                p = peak.center;
                
                pointNumber(posInOutput) = row;
                timeList(posInOutput) = index-where;
                fitWorked(posInOutput) = peak.worked;
                intensityList(posInOutput) = peak.height;
                centerList(0,posInOutput) = startingPoint.x;
                centerList(1,posInOutput) = startingPoint.y;
                
                centerSpot(0,posInOutput) = p.x;
                centerSpot(1,posInOutput) = p.y;
                
                if (analyzeSmooth) {
                    averageValues(posInOutput) = Mean(combinedSmooth(channel));
                }
                else {
                    averageValues(posInOutput) = Mean(combinedRaw(channel));
                }
                
                // This new center spot should be used as the center of the cropping window
                
                // compute how to define the brightness of the spot
                
                posInOutput++;
                pos++;
            }
        }
    }
    
    if (posInOutput!=timeList.Length()) {
        pointNumber = TruncateSize(pointNumber,posInOutput);
        timeList = TruncateSize(timeList,posInOutput);
        fitWorked = TruncateSize(fitWorked,posInOutput);
        intensityList = TruncateSize(intensityList,posInOutput);
        centerList = TruncateSize(centerList,2*posInOutput);
        centerSpot = TruncateSize(centerList,2*posInOutput);
        averageValues = TruncateSize(averageValues,posInOutput);
    }

    // Create the parameter table, typically fill along side the Add calls
    output.Finish(DTTable({
        CreateTableColumn("time",timeList),
        CreateTableColumn("fitworked",fitWorked),
        CreateTableColumn("intensity",intensityList),
        CreateTableColumn("center",DTPointCollection2D(centerList)),
        CreateTableColumn("ptNumber",pointNumber),
        CreateTableColumn("centerSpot",DTPointCollection2D(centerSpot)),
        CreateTableColumn("average",averageValues)
    }));
}
