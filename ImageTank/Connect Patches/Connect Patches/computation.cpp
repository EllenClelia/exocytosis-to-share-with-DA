#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTDoubleArray UniqueTimeValues(const DTTableColumnNumber &col)
{
    DTMutableDoubleArray sortV = Sort(col.DoubleVersion());
    ssize_t len = sortV.Length();
    int lengthOut = 1;
    double val = sortV(0);
    for (ssize_t posIn=1;posIn<len;posIn++) {
        if (sortV(posIn)!=val) {
            val = sortV(posIn);
            sortV(lengthOut++) = val;
        }
    }
    sortV = TruncateSize(sortV,lengthOut);
    return sortV;
}

DTIntArray TrackChanges(const DTTable &previousTime,const DTTable &currentTime)
{
    // Cross distances
    ssize_t previousCount = previousTime.NumberOfRows();
    ssize_t currentCount = currentTime.NumberOfRows();
    DTTableColumnMask2D previousMask = previousTime("Mask");
    DTTableColumnMask2D currentMask = currentTime("Mask");
    DTMutableDoubleArray crossDistance(previousCount,currentCount);
    for (int currentRow=0;currentRow<currentCount;currentRow++) {
        DTMask2D mask = currentMask(currentRow);
        for (int previousRow=0;previousRow<previousCount;previousRow++) {
            DTMask2D previous = previousMask(previousRow);
            DTMask2D intersection = Intersection(previousMask(previousRow),mask);
            // I have three masks here
            // mask = currentMask(currentRow)
            // previousMask(previousRow)
            // intersection - the overlap of the the above two masks.
            
            // The logic before 6/16/2025.
            // double fraction = intersection.Mask().NumberOfPoints()/(1.0*mask.Mask().NumberOfPoints());

            // this fraction is big if new overlaps highly with the previous time value.
            // double fraction = intersection.Mask().NumberOfPoints()/(1.0*previous.Mask().NumberOfPoints());
            double fraction = intersection.Mask().NumberOfPoints();
            crossDistance(previousRow,currentRow) = fraction;
        }
    }
    
    DTMutableIntArray matchedWithPreviousTimeValue(currentCount);
    matchedWithPreviousTimeValue = -1;
    
    double thresholdToAccept = 0.1; // Was 0.3
    
    // maximum of crossDistance(:,current entry) is the one at the previous time that matches the current one the best.
    DTMutableCharArray handled(previousCount,currentCount);
    handled = 0;
    // Find the best match that hasn't been handled
    int mn = int(previousCount*currentCount);
    int ij;
    while (1) {
        double bestSoFar = -1;
        int bestIJ = -1;
        for (ij=0;ij<mn;ij++) {
            if (handled(ij)==0 && crossDistance(ij)>bestSoFar) {
                bestSoFar = crossDistance(ij);
                bestIJ = ij;
            }
        }
        // Accept this as a good enough match, make sure that we don't try it again.
        if (bestSoFar<thresholdToAccept) {
            // Not enough of a match
            break;
        }
        int currentIndex = bestIJ/previousCount;
        int previousIndex = bestIJ%previousCount;
        
        // What this means is that 'currentIndex' at the current time value
        // is very close to the 'previousIndex' at the previous time value, and that
        // entry has not been matched with anything so far.
        matchedWithPreviousTimeValue(currentIndex) = previousIndex;
        
        // Ignore 'currentIndex' from consideration
        for (int i=0;i<previousCount;i++) {
            handled(i,currentIndex) = 1;
        }
        // The 'previousIndex' has been used, so it shouldn't be matched with any other at the current time value.
        for (int j=0;j<currentCount;j++) {
            handled(previousIndex,j) = 1;
        }
    }

    return matchedWithPreviousTimeValue;
}

std::pair<DTTable,DTTable> GetIDFromTable(const DTTable &currentPortion,const DTTable &previousTable,DTMutableIntArray &labels)
{
    // Input
    // currentPortion : a table (at the current time) that needs to be matched. Might be a sub-set of the original table.
    // previousTable : a sub-set of a previous table that needs to be matched, but has ID numbers
    // labels : the values for the starting table. If the currentPortion is not the full table, it has an "Index" column that specifies
    //          what the original location was for a given row.
    
    // Output
    // labels should be changed
    // return value
    //    .first  = is the sub-set of the current portion that is not matched (now with an "Index Column")
    //    .second = the sub set ofpreviousTable that was not matched.
    
    DTIntArray locationInPrevious = TrackChanges(previousTable,currentPortion);
    
    
    int countInCurrent = int(currentPortion.NumberOfRows());
    int countInPrevious = int(previousTable.NumberOfRows());
    
    DTMutableCharArray foundInPrevious(countInPrevious);
    foundInPrevious = 0;
    DTMutableCharArray foundInCurrent(countInCurrent);
    foundInCurrent = 0;
    
    DTMutableIntArray notFoundInCurrent(countInCurrent);
    int howManyNotFoundInCurrent = 0;
    DTIntArray trackNumberPrevious = DTTableColumnNumber(previousTable("ID")).IntVersion();
    
    bool hasIndex = currentPortion.Contains("Index");
    DTIntArray indexColumn;
    if (hasIndex) {
        indexColumn = DTTableColumnNumber(currentPortion("Index")).IntVersion();
    }
    DTMutableIntArray newIndexColumn(hasIndex ? 0 : countInCurrent);

    //int howManyFoundInPrevious = 0;
    int howManySet = 0;
    for (int i=0;i<countInCurrent;i++) {
        int loc = locationInPrevious(i);
        if (loc<0) {
            // not matched with the previous time, so keep it at -1
            notFoundInCurrent(howManyNotFoundInCurrent) = i;
            if (hasIndex) {
                //newIndexColumn(howManyNotFoundInCurrent) = indexColumn(i);
            }
            else {
                newIndexColumn(howManyNotFoundInCurrent) = i;
            }
            howManyNotFoundInCurrent++;
        }
        else {
            // found in the previous time
            foundInPrevious(loc) = 1;
            if (hasIndex) {
                labels(indexColumn(i)) = trackNumberPrevious(loc);
            }
            else {
                labels(i) = trackNumberPrevious(loc);
            }
            howManySet++;
        }
    }

    // What is not connected yet has an opportunity to connect to 2 steps back (and 3,4 etc in the future)
    notFoundInCurrent = TruncateSize(notFoundInCurrent,howManyNotFoundInCurrent);
    DTTable notConnectedInCurrent = currentPortion.ExtractRows(notFoundInCurrent);
    if (hasIndex) {
        //notConnectedInCurrent = notConnectedInCurrent.Replace(CreateTableColumn("Index",newIndexColumn));
    }
    else {
        newIndexColumn = TruncateSize(newIndexColumn,howManyNotFoundInCurrent);
        notConnectedInCurrent = notConnectedInCurrent.Append(CreateTableColumn("Index",newIndexColumn));
    }
    
    // The rows in previousTime that were not connected will be stopped2Back in the next iteration.
    DTMutableIntArray indexNotFoundInPrevious(countInPrevious-howManySet);
    int posInPrev = 0;
    for (int i=0;i<countInPrevious;i++) {
        if (foundInPrevious(i)==0) indexNotFoundInPrevious(posInPrev++) = i;
    }
    DTTable notUsedFromPrevious = previousTable.ExtractRows(TruncateSize(indexNotFoundInPrevious,posInPrev));
    
    return make_pair(notConnectedInCurrent,notUsedFromPrevious);
}

DTTable Computation(const DTTable &everything,int stepsBack)
{
    // time
    //
    //   DTDoubleArray time = values.DoubleVersion();
    //   double single = values(3); // Higher overhead, simpler syntax.
    // Value
    //   DTTableColumnNumber values = everything("Value");
    // Mask
    DTTableColumnNumber time = everything("time");
    DTTableColumnMask2D masks = everything("Mask");
    
    DTDoubleArray uniqueTimes = UniqueTimeValues(time);
    
    ssize_t totalRowCount = everything.NumberOfRows();
    
    int timeNumber = 0;
    DTTable previousTime = everything.ExtractRows(time.FindValue(uniqueTimes(timeNumber)));
    int howManyRowsPrevious = int(previousTime.NumberOfRows());
    DTMutableIntArray trackNumbers(howManyRowsPrevious);
    for (int i=0;i<howManyRowsPrevious;i++) trackNumbers(i) = i;
    previousTime = previousTime.Append(CreateTableColumn("ID",trackNumbers));
    
    DTMutableDoubleArray finalTrackNumbers(totalRowCount);
    finalTrackNumbers = -1;
    
    // Give them an initial track numbers
    for (int i = 0;i<howManyRowsPrevious;i++) {
        finalTrackNumbers(i) = trackNumbers(i);
    }
    int nextTrack = howManyRowsPrevious;
    int startOfEntries = howManyRowsPrevious;
    
    DTTable emptyTable = previousTime.ExtractRows(DTRange(0,0));
    DTMutableList<DTTable> stoppedInPreviousSteps(stepsBack);
    for (int i=0;i<stepsBack;i++) {
        stoppedInPreviousSteps(i) = emptyTable;
    }
    DTMutableList<DTTable> newHistory(stepsBack+1);
    
    DTProgress progress;

    for (timeNumber=1;timeNumber<uniqueTimes.Length();timeNumber++) {
        DTTable currentTime = everything.ExtractRows(time.FindValue(uniqueTimes(timeNumber)));
        
        progress.UpdatePercentage(timeNumber/(1.0*uniqueTimes.Length()));
        
        // Need to set labels at the current time based on the previous times.
        int countInCurrent = int(currentTime.NumberOfRows());
        DTMutableIntArray labelsAtCurrentTime(countInCurrent);
        labelsAtCurrentTime = -1;
        

        // Pick what we can from the previous time (has an ID column)
        auto information = GetIDFromTable(currentTime,previousTime,labelsAtCurrentTime);
        DTTable leftOverFromCurrent = information.first;
        newHistory(0) = information.second; // Left over from previous time, i.e. timeNumber-1
        
        // Go back and look for traces that stopped.
        for (int stepBack=0;stepBack<stepsBack;stepBack++) {
            information = GetIDFromTable(leftOverFromCurrent,stoppedInPreviousSteps(stepBack),labelsAtCurrentTime);
            leftOverFromCurrent = information.first;
            newHistory(stepBack+1) = information.second; // still left, just moves one step back
        }
        
        // Anything that is still not connected in labelsAtCurrentTime should get a new number
        for (int i=0;i<countInCurrent;i++) {
            if (labelsAtCurrentTime(i)<0) {
                labelsAtCurrentTime(i) = nextTrack++;
            }
        }

        // Get ready for the next step
        stoppedInPreviousSteps = Copy(newHistory);
        previousTime = currentTime.Append(CreateTableColumn("ID",labelsAtCurrentTime));

        // Put in the label IDs into the overall return table
        for (int i=0;i<countInCurrent;i++) {
            finalTrackNumbers(startOfEntries+i) = labelsAtCurrentTime(i);
        }
        startOfEntries += countInCurrent;
    }
    
    // Compress the finalTrackNumbers, might not be using a track number
    DTMutableCharArray usedTracks(nextTrack);
    usedTracks = 0;
    for (int i=0;i<totalRowCount;i++) {
        usedTracks(finalTrackNumbers(i)) = 1;
    }
    DTMutableIntArray enumerate(nextTrack);
    int newNumber = 0;
    for (int i=0;i<nextTrack;i++) {
        enumerate(i) = newNumber;
        if (usedTracks(i)) {
            newNumber++;
        }
    }
    if (newNumber!=nextTrack) {
        DTErrorMessage("Need to re-label");
    }
    
    // Add a column that counts how many entries have the same track number
    DTMutableIntArray countPerID(newNumber);
    countPerID = 0;
    for (int i=0;i<totalRowCount;i++) {
        countPerID(finalTrackNumbers(i))++;
    }
    DTMutableDoubleArray countColumn(totalRowCount);
    for (int i=0;i<totalRowCount;i++) {
        countColumn(i) = countPerID(finalTrackNumbers(i));
    }

    // Table is a list of columns
    return DTTable({
        everything("time"),
        everything("Mask"),
        CreateTableColumn("ID",finalTrackNumbers),
        CreateTableColumn("Count",countColumn)
    });
}
