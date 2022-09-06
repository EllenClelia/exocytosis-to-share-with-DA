#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

#include "Utilities.h"

DTPoint2D FindLocalMaxima(const DTDoubleArray &values,double &maxV)
{
    // Find the maximum point in the center.
    // Need a buffer around it to find the optimal center using a least squares fit
    ssize_t maxI=-1,maxJ=-1;
    maxV = 0;
    ssize_t i,j;
    ssize_t m = values.m();
    ssize_t n = values.n();
    int bdry = 3;
    for (j=bdry;j<n-bdry;j++) {
        for (i=bdry;i<m-bdry;i++) {
            if (values(i,j)>maxV) {
                maxI = i;
                maxJ = j;
                maxV = values(i,j);
            }
        }
    }

    return DTPoint2D(maxI,maxJ);
}

void Computation(const DTSet<DTImage> &everything,const DTTable &spots,
                 double time,int timeback,int timeforward,
                 const DTDictionary &parameters,DTMutableSet<DTImage> &output)
{
    DTSet<DTImage> withCache = everything.WithCache(timeback+1+timeforward);
    DTTable imageParameters = withCache.Parameters();
    DTTableColumnNumber timeValues = imageParameters("t");
    ssize_t where = timeValues.FindClosest(time);
    int channel = parameters("channel");
    
    DTTableColumnPoint2D center = spots("center");
    DTTableColumnNumber width = spots("width");
    
    ssize_t row,howMany = center.NumberOfRows();
    ssize_t index;
    
    ssize_t posInOutput = 0;
    DTMutableDoubleArray timeList(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray intensityList(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray centerList(2,howMany*(timeback+1+timeforward));
    DTMutableDoubleArray pointNumber(howMany*(timeback+1+timeforward));
    DTMutableDoubleArray centerSpot(2,howMany*(timeback+1+timeforward));
    DTMutableDoubleArray averageValues(howMany*(timeback+1+timeforward));
    double maxV = 0;
    DTPoint2D maxP;
    
    DTMutableList<DTImage> smallImages(timeback+timeforward+1);
    DTMutableList<DTImage> backgroundImages(timeback-1);
    
    DTImage backgroundImage;

    for (row=0;row<howMany;row++) {
        DTPoint2D p = center(row);
        // The center value is not a good initial guess
        // Need to correct that by taking a window around the current center
        // find the maximum point of that window
        // and use that for the back/forward points - Crop(image,box) a few lines below
        
        double w = width(row);
        DTRegion2D box = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2);

        // Go to time==where, find the maximum of that frame, and use
        // that point as the center of the region instead of the center(row) point
        DTImage rawImage = withCache(where);
        double dx = rawImage.Grid().dx();
        double dy = rawImage.Grid().dy();
        DTRegion2D biggerBox = AddBorder(box,3*dx,3*dy);

        // Find the maxima close by the initial guess p
        DTImage image = Crop(rawImage,biggerBox);
        image = ConvertToDouble(image);
        DTImage smooth = GaussianFilter(image,1);
        maxP = FindLocalMaxima(smooth(0).DoubleArray(),maxV);
        p = image.Grid().GridToSpace(maxP);
        box = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2);
        biggerBox = AddBorder(box,3*dx,3*dy);

        // Do this again, might refine the point.
        image = Crop(rawImage,biggerBox);
        image = ConvertToDouble(image);
        smooth = GaussianFilter(image,1);
        maxP = FindLocalMaxima(smooth(0).DoubleArray(),maxV);
        p = image.Grid().GridToSpace(maxP);
        box = DTRegion2D(p.x-w/2,p.x+w/2,p.y-w/2,p.y+w/2);
        biggerBox = AddBorder(box,3*dx,3*dy);


        DTPoint2D startingPoint = p;
        
        // Compute all of the raw smooth images
        ssize_t pos = 0;
        for (index=where-timeback;index<=where+timeforward;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {
                
                DTImage image = withCache(index);
                
                image = Crop(image,biggerBox);
                image = ConvertToDouble(image);
                DTImage smooth = GaussianFilter(image,1);
                
                // Crop the image with the tight box, previously it was cropped with a slightly padded box.
                smooth = Crop(smooth,box);
                
                smallImages(pos++) = smooth;
            }
        }
        
        // Compute a background image using a median.
        pos = 0;
        for (index=where-timeback;index<where-1;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {
                backgroundImages(pos) = smallImages(pos);
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
        
        pos = 0;
        for (index=where-timeback;index<=where+timeforward;index++) {
            if (index>=0 && index<withCache.NumberOfItems()) {
                DTImage smooth = smallImages(pos);
                
                // Subtract the background and put the result into a two channel image
                DTImage diff = smooth - backgroundImage;
                
                // Combine the channels, use different names
                DTImage combined(diff.Grid(),{ChangeName(smooth(0),"intensity"),ChangeName(diff(0),"difference")});
                output.Add(combined);

                maxP = FindLocalMaxima(ConvertToDouble(combined(channel)).DoubleArray(),maxV);
                p = smooth.Grid().GridToSpace(maxP);
                
                pointNumber(posInOutput) = row;
                timeList(posInOutput) = index-where;
                intensityList(posInOutput) = maxV; // Maximum(image(0));
                centerList(0,posInOutput) = startingPoint.x;
                centerList(1,posInOutput) = startingPoint.y;
                
                // centerSpot is the brightest spot, using a Gaussian Peak Fit
                // just like the bead alignment. Uses the p.x as the initial guess
                // and a reasonable starting width, but finds that spike.
                centerSpot(0,posInOutput) = p.x;
                centerSpot(1,posInOutput) = p.y;
                
                averageValues(posInOutput) = Mean(combined(channel));
                
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
        intensityList = TruncateSize(intensityList,posInOutput);
        centerList = TruncateSize(centerList,2*posInOutput);
        centerSpot = TruncateSize(centerList,2*posInOutput);
        averageValues = TruncateSize(averageValues,posInOutput);
    }

    // Create the parameter table, typically fill along side the Add calls
    output.Finish(DTTable({
        CreateTableColumn("time",timeList),
        CreateTableColumn("intensity",intensityList),
        CreateTableColumn("center",DTPointCollection2D(centerList)),
        CreateTableColumn("ptNumber",pointNumber),
        CreateTableColumn("centerSpot",DTPointCollection2D(centerSpot)),
        CreateTableColumn("average",averageValues)
    }));
}
