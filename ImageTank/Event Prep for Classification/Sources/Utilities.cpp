//
//  Utilities.cpp
//  GaussFit
//
//  Created by David Adalsteinsson on 2/14/22.
//  Copyright © 2022 Visual Data Tools, Inc. All rights reserved.
//

#include "Utilities.h"

#include "DTFunction.h"
#include "DTDictionary.h"
#include "DTFunctionFit.h"
#include "DTFunction1D.h"
#include "DTDoubleArrayOperators.h"

DTImage GaussianFilter(const DTImage &image,double sigma)
{
    DTList<DTImageChannel> channels = image.Channels();
    DTMutableList<DTImageChannel> returnChannels(channels.Length());
    ssize_t i;
    for (i=0;i<channels.Length();i++) {
        returnChannels(i) = GaussianFilter(channels(i),sigma);
    }
    return DTImage(image.Grid(),returnChannels);
}

DTImageChannel GaussianFilter(const DTImageChannel &channel,double sigma)
{
    if (channel.IsEmpty()) {
        return channel;
    }
    else if (channel.isFloat()) {
        return DTImageChannel(channel.Name(),GaussianFilter(channel.FloatArray(),sigma));
    }
    else if (channel.isDouble()) {
        return DTImageChannel(channel.Name(),GaussianFilter(channel.DoubleArray(),sigma));
    }
    else {
        DTErrorMessage("GaussianFilter","Only supports double or float");
        return channel;
    }
}

DTFloatArray GaussianFilter(const DTFloatArray &arr,double sigma)
{
    if (arr.IsEmpty()) return arr;
    
    const ssize_t m = (int)arr.m();
    const ssize_t n = (int)arr.n();
    
    DTMutableFloatArray firstArray(m,n);
    
    ssize_t i;
    ssize_t padBy = int(floor(sigma*3.5));
    if (padBy==0) padBy = 1;

    if (arr.m()<=1+2*padBy || arr.n()<=1+2*padBy) {
        DTErrorMessage("GaussianFilter(image,sigma)","sigma too large");
        return DTFloatArray();
    }
    
        
    DTMutableFloatArray weights(1+padBy);
    double weight;
    double sum = 0;
    weight = erf(0.5/(sqrt(2.0)*sigma));
    weights(0) = weight;
    sum += weight;
    
    for (i=1;i<=padBy;i++) {
        weight = 0.5*(erfc((i-0.5)/(sqrt(2.0)*sigma)) - erfc((i+0.5)/(sqrt(2.0)*sigma)));
        weights(i) = weight;
        sum += 2*weight;
        // v = exp(-i*i/(2*sigma*sigma))/(sqrt(2*M_PI)*sigma);
        // cerr << weight << ", v = " << v << ", combined = " << sum << endl;
    }
    
    weights *= 1.0/sum;
    
    ssize_t j;
    
    float *retD = firstArray.Pointer();
    const float *D = arr.Pointer();
    const float *weightsD = weights.Pointer();
    
    // Now smooth each direction.  Use mirrored boundary conditions, but might want to
    // have an option of using periodic boundary conditions.
    
    ssize_t ij;
    ssize_t jm;
    ssize_t k;
    
    // Smooth (i-2:i+2,j).
    ij = 0;
    ssize_t irefl;
    for (j=0;j<n;j++) {
        // Left hand side
        jm = j*m;
        for (i=0;i<padBy;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                irefl = i-k;
                if (irefl<0) irefl = -irefl-1;
                sum += weightsD[k]*(D[ij+k] + D[irefl+jm]);
            }
            
            retD[ij] = sum;
            ij++;
        }
        
        // Bulk of the points
        ij = padBy + j*m;
        for (i=padBy;i<m-padBy;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                sum += weightsD[k]*(D[ij-k]+D[ij+k]);
            }
            retD[ij] = sum;
            ij++;
        }
        
        // Right side
        for (i=m-padBy;i<m;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                irefl = i+k;
                if (irefl>=m) irefl = m-irefl+m-1;
                sum += weightsD[k]*(D[irefl+jm] + D[ij-k]);
            }
            retD[ij] = sum;
            ij++;
        }
    }
    
    // Smooth (i,j-2:j+2)
    ij = 0;
    
    DTMutableFloatArray returnArray(m,n);
    retD = returnArray.Pointer();
    D = firstArray.Pointer();
    
    // Bottom
    j = 0;
    ssize_t jrefl;
    
    ij = 0;
    for (j=0;j<padBy;j++) {
        jm = j*m;
        for (i=0;i<m;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                jrefl = j-k;
                if (jrefl<0) jrefl = -jrefl-1;
                sum += weightsD[k]*(D[i+jrefl*m] + D[ij + k*m]);
            }
            retD[ij] = sum;
            ij++;
        }
    }
    
    // Bulk of the points
    ij = padBy*m;
    for (j=padBy;j<n-padBy;j++) {
        ij = 0 + j*m;
        for (i=0;i<m;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                sum += weightsD[k]*(D[ij-k*m]+D[ij+k*m]);
            }
            retD[ij] = sum;
            ij++;
        }
    }
    
    // Top
    for (j=n-padBy;j<n;j++) {
        for (i=0;i<m;i++) {
            sum = weightsD[0]*D[ij];
            
            for (k=1;k<=padBy;k++) {
                jrefl = j+k;
                if (jrefl>=n) jrefl = n-jrefl+n-1;
                sum += weightsD[k]*(D[i+jrefl*m] + D[ij-k*m]);
            }
            
            retD[ij] = sum;
            ij++;
        }
    }
    
    cerr << returnArray(300,300) << std::endl;
    
    return returnArray;
}

DTDoubleArray GaussianFilter(const DTDoubleArray &arr,double sigma)
{
    if (arr.IsEmpty()) return arr;
    
    const ssize_t m = (int)arr.m();
    const ssize_t n = (int)arr.n();
    
    DTMutableDoubleArray firstArray(m,n);
    
    ssize_t i;
    ssize_t padBy = int(floor(sigma*3.5));
    if (padBy==0) padBy = 1;

    if (arr.m()<=1+2*padBy || arr.n()<=1+2*padBy) {
        DTErrorMessage("GaussianFilter(image,sigma)","sigma too large");
        return DTDoubleArray();
    }
    
        
    DTMutableDoubleArray weights(1+padBy);
    double weight;
    double sum = 0;
    weight = erf(0.5/(sqrt(2.0)*sigma));
    weights(0) = weight;
    sum += weight;
    
    for (i=1;i<=padBy;i++) {
        weight = 0.5*(erfc((i-0.5)/(sqrt(2.0)*sigma)) - erfc((i+0.5)/(sqrt(2.0)*sigma)));
        weights(i) = weight;
        sum += 2*weight;
        // v = exp(-i*i/(2*sigma*sigma))/(sqrt(2*M_PI)*sigma);
        // cerr << weight << ", v = " << v << ", combined = " << sum << endl;
    }
    
    weights *= 1.0/sum;
    
    ssize_t j;
    
    double *retD = firstArray.Pointer();
    const double *D = arr.Pointer();
    const double *weightsD = weights.Pointer();
    
    // Now smooth each direction.  Use mirrored boundary conditions, but might want to
    // have an option of using periodic boundary conditions.
    
    ssize_t ij;
    ssize_t jm;
    ssize_t k;
    
    // Smooth (i-2:i+2,j).
    ij = 0;
    ssize_t irefl;
    for (j=0;j<n;j++) {
        // Left hand side
        jm = j*m;
        for (i=0;i<padBy;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                irefl = i-k;
                if (irefl<0) irefl = -irefl-1;
                sum += weightsD[k]*(D[ij+k] + D[irefl+jm]);
            }
            
            retD[ij] = sum;
            ij++;
        }
        
        // Bulk of the points
        ij = padBy + j*m;
        for (i=padBy;i<m-padBy;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                sum += weightsD[k]*(D[ij-k]+D[ij+k]);
            }
            retD[ij] = sum;
            ij++;
        }
        
        // Right side
        for (i=m-padBy;i<m;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                irefl = i+k;
                if (irefl>=m) irefl = m-irefl+m-1;
                sum += weightsD[k]*(D[irefl+jm] + D[ij-k]);
            }
            retD[ij] = sum;
            ij++;
        }
    }
    
    // Smooth (i,j-2:j+2)
    ij = 0;
    
    DTMutableDoubleArray returnArray(m,n);
    retD = returnArray.Pointer();
    D = firstArray.Pointer();
    
    // Bottom
    j = 0;
    ssize_t jrefl;
    
    ij = 0;
    for (j=0;j<padBy;j++) {
        jm = j*m;
        for (i=0;i<m;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                jrefl = j-k;
                if (jrefl<0) jrefl = -jrefl-1;
                sum += weightsD[k]*(D[i+jrefl*m] + D[ij + k*m]);
            }
            retD[ij] = sum;
            ij++;
        }
    }
    
    // Bulk of the points
    ij = padBy*m;
    for (j=padBy;j<n-padBy;j++) {
        ij = 0 + j*m;
        for (i=0;i<m;i++) {
            sum = weightsD[0]*D[ij];
            for (k=1;k<=padBy;k++) {
                sum += weightsD[k]*(D[ij-k*m]+D[ij+k*m]);
            }
            retD[ij] = sum;
            ij++;
        }
    }
    
    // Top
    for (j=n-padBy;j<n;j++) {
        for (i=0;i<m;i++) {
            sum = weightsD[0]*D[ij];
            
            for (k=1;k<=padBy;k++) {
                jrefl = j+k;
                if (jrefl>=n) jrefl = n-jrefl+n-1;
                sum += weightsD[k]*(D[i+jrefl*m] + D[ij-k*m]);
            }
            
            retD[ij] = sum;
            ij++;
        }
    }
    
    return returnArray;
}

DTFloatArray MedianOfStack(const DTFloatArray &stack,ssize_t slices)
{
    // Median so far
    ssize_t m = stack.m();
    ssize_t n = stack.n();
    DTMutableFloatArray toReturn(m,n);
    
    DTMutableFloatArray list(slices);
    ssize_t i,j,k;
    ssize_t half = (slices%2==0 ? slices/2-1 : slices/2);
    for (j=0;j<n;j++) {
        for (i=0;i<m;i++) {
            for (k=0;k<slices;k++) {
                list(k) = stack(i,j,k);
            }
            std::sort(list.Pointer(),list.Pointer()+k);
            if (slices%2==0) {
                toReturn(i,j) = (list(half)+list(half+1))*0.5;
            }
            else {
                toReturn(i,j) = list(half);
            }
        }
    }
    
    return toReturn;
}

DTImage MedianOfImages(const DTList<DTImage> &images)
{
    if (images.Length()==0) return DTImage();
    
    DTImage firstImage = images(0);
    ssize_t channelCount = firstImage.Channels().Length();
    ssize_t m = firstImage.m();
    ssize_t n = firstImage.n();

    DTMutableList<DTImageChannel> outputChannels(channelCount);
    
    DTMutableFloatArray stack(m,n,images.Length());
    for (ssize_t chNumber=0;chNumber<channelCount;chNumber++) {
        for (ssize_t imNumber=0;imNumber<images.Length();imNumber++) {
            CopyIntoSlice(stack,ConvertToFloat(images(imNumber)(chNumber)).FloatArray(),imNumber);
        }
        DTFloatArray median = MedianOfStack(stack,stack.o());
        outputChannels(chNumber) = DTImageChannel(firstImage(chNumber).Name(),median);
    }
    
    return DTImage(firstImage.Grid(),outputChannels);
}

double ComputeRMSE(const DTDoubleArray &yValuesList,const DTDoubleArray &fitValuesList)
{
    if (yValuesList.Length()!=fitValuesList.Length()) {
        DTErrorMessage("ComputeR2","Sizes don't match");
        return NAN;
    }
    
    double SSerr = 0.0;
    double temp;
    int residualCount = 0;
    ssize_t i;
    ssize_t howManyEntries = yValuesList.Length();
    for (i=0;i<howManyEntries;i++) {
        temp = yValuesList(i) - fitValuesList(i);
        if (isfinite(temp)) {
            residualCount++;
            SSerr += temp*temp;
        }
    }
    
    return sqrt(SSerr/howManyEntries);
}

double ComputeR2(const DTDoubleArray &yValuesList,const DTDoubleArray &fitValuesList)
{
    if (yValuesList.Length()!=fitValuesList.Length()) {
        DTErrorMessage("ComputeR2","Sizes don't match");
        return NAN;
    }
    
    double SSerr = 0.0;
    double temp;
    int residualCount = 0;
    ssize_t i;
    ssize_t howManyEntries = yValuesList.Length();
    for (i=0;i<howManyEntries;i++) {
        temp = yValuesList(i) - fitValuesList(i);
        if (isfinite(temp)) {
            residualCount++;
            SSerr += temp*temp;
        }
    }
    
    // Compute the average of y values
    double sumOfValues = 0.0;
    double sumW = 0.0;
    sumW = howManyEntries;
    for (i=0;i<howManyEntries;i++) {
        sumOfValues += yValuesList(i);
    }
    
    double yaverage = sumOfValues/sumW;
    
    // Compute SStot
    double SStot = 0.0;
    for (i=0;i<howManyEntries;i++) {
        temp = yValuesList(i)-yaverage;
        SStot += temp*temp;
    }
    
    double R2 = 1.0-SSerr/SStot;
    
    return R2;
}

DTTable Prune(const DTTable &table)
{
    double switchBelow = 1.0/2.0;
    double switchAbove = 2.0/3.0;
    
    if (table.NumberOfRows()<2) return table;

    // time
    DTTableColumnNumber timeValues = table("t");
    DTDoubleArray time = timeValues.DoubleVersion();
    DTTableColumnNumber intensityValues = table("y");
    DTDoubleArray value = intensityValues.DoubleVersion();
    
    ssize_t i,howMany = time.Length();
    
    DTMutableDoubleArray newValues(howMany);
    newValues(0) = std::max(value(0),value(1));
    for (i=1;i<howMany-1;i++) {
        newValues(i) = std::max(std::max(value(i-1),value(i)),value(i+1));
    }
    newValues(howMany-1) = std::max(value(howMany-2),value(howMany-1));
    
    // Look for time = 0
    ssize_t startPosition = timeValues.FindClosest(0);
    double vmax = newValues(startPosition);

    // Find the maximum
    for (i=startPosition;i<howMany;i++) {
        if (newValues(i)>vmax) {
            vmax = newValues(i);
        }
        else if (newValues(i)<vmax*switchBelow) {
            break;
        }
    }
    
    double minValue = NAN;
    double stopAbove = vmax*switchAbove;
    if (i<howMany) {
        minValue = newValues(i);
        while (i<howMany) {
            if (newValues(i)>stopAbove || (newValues(i)-minValue)>0.3*(vmax-minValue)) {
                break;
            }
            else if (newValues(i)<minValue) {
                minValue = newValues(i);
            }
            i++;
        }
    }
    
    return table.ExtractRows(DTRange(0,i));
}

QuantifyEvent Quantify(const DTSet<DTImage> &images,const DTDictionary &parameters)
{
    //DTSet<DTImage> imagesToView = images.ExtractRows("ptNumber",[](double v) {return (v==0);}); - To pick a sub-table with a particular value for ptNumber
    
    DTTable images_par = images.Parameters(); // In the debugger, use images_par.pinfo() and images_par.pall() to see the content.

    bool useAverage = parameters.GetNumber("useAverage",true);
    int channel = parameters("channel");
    
    DTTableColumnNumber time = images_par("time");
    DTTableColumnNumber average = images_par("average");
    DTTableColumnNumber intensity;
    if (useAverage) {
        intensity = images_par("average");
    }
    else {
        intensity = images_par("intensity");
    }
    DTTableColumnNumber failure = images_par("failure");

    ssize_t locOrigin = time.FindClosest(0);
    int shift = 0;
    
    // Check to see if the value before or after is larger. If so, use that
    if (locOrigin+1<images_par.NumberOfRows()) {
        // Might want to move the focus one point to the left or right.
        // Move it to the left if the initial screening was too late -> shift<0
        // move it to the right if the initial screening picked it up too early -> shift>0
        
        double sizeScale = std::max(intensity(locOrigin),intensity(locOrigin+1))-average(0);

        if (intensity(locOrigin-1)>sizeScale*0.8+average(0)) {
            // Previous value is big enough, back up one step.
            locOrigin--;
            shift--;
        }
        else if (intensity(locOrigin+1)-intensity(locOrigin)>0.6*sizeScale) {
            // Next value is much larger, consider this a false start.
            locOrigin++;
            shift++;
        }
    }

    // Compute the average for the background before the location. Also exclude the frame immediately before the starting frame.
    double sum = 0.0;
    ssize_t i;
    for (i=0;i<locOrigin-1;i++) {
        sum += average(i);
    }
    double mean = sum/(locOrigin-1);
    // Compute the standard deviation
    
    QuantifyEvent toReturn;

    toReturn.shift = shift; // The time when the event starts
    toReturn.average = mean; // The "size" of the signal before the start.
    
    // Two methods to compute the width
    double width;
    if (useAverage) {
        // Variability in the averages.
        double sumSq = 0;
        for (i=0;i<locOrigin-1;i++) {
            double val = average(i);
            sumSq = (val-mean)*(val-mean);
        }
        width = sqrt(sumSq/(locOrigin-1)); // standard deviation of the avarage
        toReturn.width = width;
    }
    else {
        // Compute the standard deviation of the image values (individual pixels, not just the averages).
        ssize_t firstBin = -3000;
        ssize_t lastBin = 3000+mean*10;
        ssize_t maxBin = lastBin-firstBin+1;
        DTMutableDoubleArray bins(maxBin);
        bins = 0;
        for (i=0;i<locOrigin-1;i++) {
            DTImage single = images(i);
            DTDoubleArray values = ConvertToDouble(single)(channel).DoubleArray();
            ssize_t ij,len = values.Length();
            for (ij=0;ij<len;ij++) {
                double v = values(ij);
                int binNumber = int(round(v));
                if (firstBin<=binNumber && binNumber<=lastBin) {
                    bins(binNumber-firstBin)++;
                }
            }
        }
        
        ssize_t minNonZero = 0;
        while (minNonZero<maxBin && bins(minNonZero)==0) minNonZero++;
        ssize_t maxNonZero = int(maxBin)-1;
        while (maxNonZero>=0 && bins(maxNonZero)==0) maxNonZero--;
        if (minNonZero+firstBin>0) minNonZero = -firstBin;
        
        // Compute the standard deviation
        sum = 0;
        double sumSq = 0.0;
        int totalCount = 0;
        // Average from before is going to work, in fact a little more accurate.
        // Compute the standard error
        ssize_t ival;
        for (i=minNonZero;i<=maxNonZero;i++) {
            // Should really do
            // sum (i-val)^2 bins(i) times
            ival = i+firstBin;
            sumSq += bins(i)*(ival-mean)*(ival-mean);
            totalCount += bins(i);
        }
        width = sqrt(sumSq/totalCount);
        
        // Save the histogram, in case you want to look at it in ImageTank.
        DTMutableDoubleArray intensityList(maxNonZero-minNonZero+1);
        DTMutableDoubleArray countList(maxNonZero-minNonZero+1);
        int pos = 0;
        for (i=minNonZero;i<=maxNonZero;i++) {
            intensityList(pos) = i+firstBin;
            countList(pos) = bins(i);
            pos++;
        }
        
        // The histogram of the intensity (combined from all of the images).
        toReturn.histogram = DTTable({
            CreateTableColumn("intensity",intensityList),
            CreateTableColumn("count",countList)
        });
        toReturn.width = width; // The variability in the signal before the event. The standard deviation from the mean.
    }

    // Create the fitting data table.
    ssize_t rowCount = images_par.NumberOfRows();
    
    int howFar = int(rowCount);
    double scaleForMin = parameters.GetNumber("minIntensityForFit",2.0);
    double minIntensityValueForFit = mean+scaleForMin*width;
    int allowIntensityFailures = parameters.GetNumber("Allow intensity failures",100);
    bool includeIntensityFailures = int(parameters.GetNumber("Include intensity failures",false));
    bool includePlateau = parameters.GetNumber("Include Plateau",0);

    DTMutableDoubleArray xvalTemp(rowCount), tvalTemp(rowCount), yvalTemp(rowCount);
    int pos = 0;
    int howManyIntensityFailures = 0;
    while (i<rowCount && time(i)<shift) i++; // Find the first index

    DTMutableIntArray skippedOver(allowIntensityFailures);
    int posInSkippedOver = 0;
    
    while (i<rowCount && time(i)<=howFar+shift) {
        double t = time(i);
        // Only include points that are coming from a valid fit
        if (useAverage || failure(i)==0) {
            // Don't include intensity values that are too low
            if (intensity(i)>minIntensityValueForFit) {
                if (posInSkippedOver>0 && includeIntensityFailures) {
                    // These entries were skipped over, but now need to be included
                    for (int tempI=0;tempI<posInSkippedOver;tempI++) {
                        int whichI = skippedOver(tempI);
                        xvalTemp(pos) = time(whichI)-shift;
                        tvalTemp(pos) = time(whichI);
                        yvalTemp(pos) = intensity(whichI);
                        pos++;
                    }
                    posInSkippedOver = 0;
                }
                
                xvalTemp(pos) = t-shift;
                tvalTemp(pos) = t;
                yvalTemp(pos) = intensity(i);
                pos++;
            }
            else {
                howManyIntensityFailures++;
                if (howManyIntensityFailures>allowIntensityFailures) break;
                skippedOver(posInSkippedOver++) = int(i);
            }
        }
        i++;
    }
    DTDoubleArray xval = TruncateSize(xvalTemp,pos);
    DTDoubleArray tval = TruncateSize(tvalTemp,pos);
    DTDoubleArray yval = TruncateSize(yvalTemp,pos);
    
    DTTable tableUsedForFit({CreateTableColumn("x",xval),CreateTableColumn("t",tval),CreateTableColumn("y",yval)});
    tableUsedForFit = Prune(tableUsedForFit);
    toReturn.pointsUsedForFit = tableUsedForFit;
    
    if (pos==0) {
        // Failed, don't want an error message to propagate to IT.
        toReturn.average = NAN;
        toReturn.width = NAN;
        toReturn.shift = 0;
        toReturn.delay = 0;
        toReturn.base = NAN;
        toReturn.spike = NAN;
        toReturn.decay = NAN;
        toReturn.R2 = NAN;

        return toReturn;
    }


    
    DTTableColumnNumber numCol = tableUsedForFit("x");
    xval = numCol.DoubleVersion();
    numCol = tableUsedForFit("t");
    tval = numCol.DoubleVersion();
    numCol = tableUsedForFit("y");
    yval = numCol.DoubleVersion();

    DTMutableDictionary knownConstants;
    DTMutableDictionary guesses;
    DTFunction1D fitFcn;
    DTDoubleArray fitValues;

    // Compute the piecewise fits
    
    // First a slight hack. The xval list before was initially shifted by "shift"
    // This allows the function fit to use "a + b*exp(-cx)" instead of "a + b*exp(-c*(x-shift))
    // However for the piecewise function we want to vary the shift, so the x values should be
    // changed back to the x values relative to the starting frame.
    xval = xval + shift; // Isn't this just tval?
    knownConstants("x") = xval; // Technically not needed since the xval is shared, but makes it more explicit.
    
    int maximumShift = std::min(20,howFar+shift-7);
    if (maximumShift<shift) maximumShift = shift;
    int howManyShifts = maximumShift-shift;
    
    if (howManyShifts<=1) {
        toReturn.average = NAN;
        toReturn.width = NAN;
        toReturn.shift = 0;
        toReturn.delay = 0;
        toReturn.base = NAN;
        toReturn.spike = NAN;
        toReturn.decay = NAN;
        toReturn.R2 = NAN;

        return toReturn;
    }
    
    DTMutableDoubleArray kinkList(howManyShifts);
    DTMutableDoubleArray baseList(howManyShifts);
    DTMutableDoubleArray spikeList(howManyShifts);
    DTMutableDoubleArray decayList(howManyShifts);
    DTMutableDoubleArray R2List(howManyShifts);
    DTMutableDoubleArray RMSEList(howManyShifts);

    double av = mean;
    double bv = yval(0)-mean;
    double cv = 1.0;
    
    // Components for the analytical function
    DTFunction a = DTFunction::Constant("a");
    DTFunction b = DTFunction::Constant("b");
    DTFunction c = DTFunction::Constant("c");
    DTFunction x = DTFunction::Constant("x");
    
    for (int kink=shift;kink<maximumShift;kink++) { // use delay, and then kink = shift + delay
        DTFunction combo = IfFunction(x<kink,a*a+b,a*a+b*exp(-c*(x-kink)));
        // if x<kink just a C, for x>=kink the function is Const + b*exp(-c*(distance from kink)).
        // Const = a^2+b to enforce that the shift is positive. b might be negative, but then the signal is not decaying.
        
        guesses("a") = sqrt(fabs(av));
        guesses("b") = bv;
        guesses("c") = cv;
        
        // This is a Levenberg–Marquardt method to solve the non-linear least squares problem
        fitFcn = FunctionFit(combo,yval,knownConstants,guesses);
        fitValues = fitFcn(xval);
        
        kinkList(kink-shift) = kink; // shift+delay
        baseList(kink-shift) = pow(guesses("a"),2.0);
        spikeList(kink-shift) = guesses("b");
        decayList(kink-shift) = guesses("c");
        R2List(kink-shift) = ComputeR2(yval,fitValues); // Quality of the fit. Not necessarily proper for a non-linear fit,
        RMSEList(kink-shift) = ComputeRMSE(yval,fitValues); // A better measurement of the quality of the fit.
    }
    
    // Find the best R2 value. A long kink will reduce the R2 and that is a reasonable
    // penalty, but a large flat portion after the decay portion will also lower R2
    // and that might be a problem since that fit is considered worse.
    double bestR = 0.0;
    ssize_t bestRindex = 0;
    for (i=0;i<R2List.Length();i++) {
        if (R2List(i)>bestR) {
            bestR = R2List(i);
            bestRindex = i;
        }
    }
    
    // To find the best fit, look at what minimizes the residual.
    double bestRMSE = INFINITY;
    // Find the best RMSE value, i.e. the smallest
    ssize_t bestRMSEindex = 0;
    for (i=0;i<RMSEList.Length();i++) {
        if (RMSEList(i)<bestRMSE) {
            bestRMSE = RMSEList(i);
            bestRMSEindex = i;
        }
    }
    
    ssize_t indexForKink = bestRMSEindex; // used to be bestRindex. That means that the R2List is ignored, but you can
    
    // Compute the R2 value for this fit. The R^2 for the fit should not include
    // the values before the kink.
    int kink = kinkList(indexForKink);
    DTFunction combo = IfFunction(x<kink,a*a+b,a*a+b*exp(-c*(x-kink)));
    guesses("a") = sqrt(fabs(av));
    guesses("b") = bv;
    guesses("c") = cv;
    fitFcn = FunctionFit(combo,yval,knownConstants,guesses);
    DTDoubleArray fitValuesForR2 = fitFcn(xval);
    // Trim away the flat part before the kink if there is any
    DTDoubleArray yvalForR2 = yval;

    // Compute the R^2 for values that happen at and after where the kink is.
    // Since the x values might have missing indices
    int startAt = 0;
    for (startAt=0;startAt<xval.Length();startAt++) {
        if (xval(startAt)>=kink) break;
    }
    // startAt = 0; // Include the flat part before the kink.
    if (startAt>0 && includePlateau==false) {
        // Strip out the part before the kink
        fitValuesForR2 = ExtractRows(fitValuesForR2,DTRange(startAt,yvalForR2.Length()-startAt));
        yvalForR2      = ExtractRows(yvalForR2,     DTRange(startAt,yvalForR2.Length()-startAt));
    }
    bestR = ComputeR2(yvalForR2,fitValuesForR2);
    
    //The returned fit, quality decay etc is the result of the best piecewise fit
    //not the original fit.
    // if (startAt<kinkList.Length()) {
    if (yvalForR2.Length()>0) {
        toReturn.delay = kinkList(indexForKink)-shift;
        toReturn.decay = decayList(indexForKink);
        toReturn.base = baseList(indexForKink);
        toReturn.spike = spikeList(indexForKink);
        toReturn.R2 = bestR;
    }
    else {
        toReturn.delay = NAN;
        toReturn.decay = NAN;
        toReturn.base = NAN;
        toReturn.spike = NAN;
        toReturn.R2 = NAN;
    }
    
    toReturn.piecewiseFitResults = DTTable({
        CreateTableColumn("kink",kinkList),
        CreateTableColumn("base",baseList),
        CreateTableColumn("spike",spikeList),
        CreateTableColumn("decay",decayList),
        CreateTableColumn("R2",R2List),
        CreateTableColumn("RMSE",RMSEList)});

    return toReturn;
}

bool EvaluateGaussianPeak(const DTDictionary &constants,DTMutableDoubleArray &returnArray);

LocalPeak FindGaussianPeak(const DTImage &image,const DTDictionary &parameters)
{
    int channel = parameters("channel");
    
    // Find the maximum point in the center.
    // Need a buffer around it to find the optimal center using a least squares fit
    DTDoubleArray values = ConvertToDouble(image(channel)).DoubleArray();
    
    ssize_t maxI=-1,maxJ=-1;
    ssize_t i,j;
    ssize_t m = values.m();
    ssize_t n = values.n();
    int bdry = 3;
    double maxV = -INFINITY;
    double minV = INFINITY;
    for (j=bdry;j<n-bdry;j++) {
        for (i=bdry;i<m-bdry;i++) {
            if (values(i,j)>maxV) {
                maxI = i;
                maxJ = j;
                maxV = values(i,j);
            }
            if (values(i,j)<minV) {
                minV = values(i,j);
            }
        }
    }
    
    DTMutableDictionary guesses;
    guesses("x0") = maxI;
    guesses("y0") = maxJ;
    guesses("scale") = maxV-minV;
    guesses("base") = minV;
    guesses("radius") = 5;

    DTMutableDictionary knownConstants;

    // Function form is
    // base + scale*exp(-((x-x0)^2+(y-y0)^2)/(2radius^2));
    DTDictionary returned = FunctionFit(EvaluateGaussianPeak,values,knownConstants,guesses);

    LocalPeak toReturn;
    
    toReturn.center = image.Grid().GridToSpace(DTPoint2D(returned("x0"),returned("y0")));
    toReturn.base = returned("base");
    toReturn.height = toReturn.base + double(returned("scale"));
    toReturn.width = 2*sqrt(2*log(2))*fabs(returned("radius")); // FWHM - full width at half maximumum
    toReturn.failureMode = 0;
    
    DTMutableDoubleArray fitValues(values.m(),values.n());
    EvaluateGaussianPeak(returned,fitValues);
    toReturn.R2 = ComputeR2(values,fitValues);
    
    double R2threshold = parameters("R2 for Peak");
    if (returned("LM::Status")==1) {
        // Might still fail if it is outside of the domain
        if (BoundingBox(image.Grid()).PointLiesInside(toReturn.center)==false) {
            // Went outside the domain
            toReturn.height = maxV;
            toReturn.failureMode = 2;
        }
        else if (toReturn.height<toReturn.base) {
            toReturn.failureMode = 5; // Inverted peak
        }
        else if (toReturn.R2<R2threshold) {
            toReturn.failureMode = 6; // Bad fit for a gaussian.
        }
        else if (toReturn.width>m) {
            // Too wide by far
            toReturn.failureMode = 4;
        }
        else if (toReturn.base<minV - 3*(maxV-minV)) {
            //
            toReturn.failureMode = 3;
        }
    }
    else {
        toReturn.failureMode = 1; // Convergence failed
    }
    
    return toReturn;
}

bool EvaluateGaussianPeak(const DTDictionary &constants,DTMutableDoubleArray &returnArray)
{
    // Function form is
    
    // Base + scale*exp(-((x-xC)^2+(y-yC)^2)/2radius^2);
    int m = (int)returnArray.m();
    int n = (int)returnArray.n();
    
    double base = constants("base");
    double scale = constants("scale");
    double x0 = constants("x0");
    double y0 = constants("y0");
    double radius = constants("radius");
    
    int i,j,ij;
    double x,y,arg;
    double C = 1.0/(2*radius*radius);
    ij = 0;
    double *returnArrayD = returnArray.Pointer();
    for (j=0;j<n;j++) {
        y = (j-y0);
        for (i=0;i<m;i++) {
            x = (i-x0);
            arg = -(x*x+y*y)*C;
            returnArrayD[ij] = base + scale*exp(arg);
            ij++;
        }
    }
    
    return true;
}

LocalPeak FindMaximumPeak(const DTImage &image,int channel)
{
    double maxV = 0.0;
    DTPoint2D maxP = FindLocalMaxima(ConvertToDouble(image(channel)).DoubleArray(),maxV);
    
    LocalPeak toReturn;
    toReturn.center = image.Grid().GridToSpace(maxP);
    toReturn.height = maxV;
    toReturn.base = 0.0;
    toReturn.width = 0.0;
    toReturn.failureMode = 0;
    toReturn.R2 = 0;

    return toReturn;
}

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

