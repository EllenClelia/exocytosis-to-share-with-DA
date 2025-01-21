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
            DTMask2D intersection = Intersection(previousMask(previousRow),mask);
            double fraction = intersection.Mask().NumberOfPoints()/(1.0*mask.Mask().NumberOfPoints());
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
    
    DTTable stopped2Back = previousTime.ExtractRows(DTRange(0,0));
    DTMutableList<DTTable> history(1);
    history(0) = stopped2Back;
    DTMutableList<DTTable> newHistory(1);

    for (timeNumber=1;timeNumber<uniqueTimes.Length();timeNumber++) {
        DTTable currentTime = everything.ExtractRows(time.FindValue(uniqueTimes(timeNumber)));
        
        // Need to set labels at the current time based on the previous times.
        int countInCurrent = int(currentTime.NumberOfRows());
        DTMutableIntArray labelsAtCurrentTime(countInCurrent);
        labelsAtCurrentTime = -1;

        // Pick what we can from the previous time (has an ID column)
        auto information = GetIDFromTable(currentTime,previousTime,labelsAtCurrentTime);
        DTTable leftOverFromCurrent = information.first;
        DTTable leftOverFromPrevious = information.second;
        newHistory(0) = leftOverFromPrevious;
        
        // One idea is to go to the accepted table from two steps back and see
        // find where the entries in leftOverFromPrevious were
                
        

        // Pick what we can from two steps back that was not found in the previous time
        information = GetIDFromTable(leftOverFromCurrent,history(0),labelsAtCurrentTime);
        DTTable leavingUnassigned = information.second;
        DTTable givingUpOn = information.second;
        
        // Anything that is still not connected in labelsAtCurrentTime should get a new number
        for (int i=0;i<countInCurrent;i++) {
            if (labelsAtCurrentTime(i)<0) {
                labelsAtCurrentTime(i) = nextTrack++;
            }
        }

        // stopped2Back = leftOverFromPrevious;
        
        // Get ready for the next step
        history = Copy(newHistory);
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
        
    // Table is a list of columns
    return DTTable({
        everything("time"),
        everything("Mask"),
        CreateTableColumn("ID",finalTrackNumbers)
    });
}

DTTable ComputationOld(const DTTable &everything,int stepsBack)
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
    
    DTMutableList<DTTable> history(1);
    DTTable stopped2Back = previousTime.ExtractRows(DTRange(0,0));
    
    for (timeNumber=1;timeNumber<uniqueTimes.Length();timeNumber++) {
        DTTable currentTime = everything.ExtractRows(time.FindValue(uniqueTimes(timeNumber)));
        
        // Need to set labels at the current time based on the previous times.
        int countInCurrent = int(currentTime.NumberOfRows());
        DTMutableIntArray labelsAtCurrentTime(countInCurrent);
        labelsAtCurrentTime = -1;

        // Pick what we can from the previous time (has an ID column)
        auto information = GetIDFromTable(currentTime,previousTime,labelsAtCurrentTime);
        DTTable leftOverFromCurrent = information.first;
        DTTable leftOverFromPrevious = information.second;

        // Pick what we can from two steps back
        information = GetIDFromTable(leftOverFromCurrent,stopped2Back,labelsAtCurrentTime);
        DTTable leavingUnassigned = information.second;
        DTTable givingUpOn = information.second;
        
        // Anything that is still not connected in labelsAtCurrentTime should get a new number
        for (int i=0;i<countInCurrent;i++) {
            if (labelsAtCurrentTime(i)<0) {
                labelsAtCurrentTime(i) = nextTrack++;
            }
        }

        // Prepare for the next iteration
        stopped2Back = leftOverFromPrevious;
        
        // Get ready for the next step
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
        
    // Table is a list of columns
    return DTTable({
        everything("time"),
        everything("Value"),
        everything("Mask"),
        CreateTableColumn("ID",finalTrackNumbers)
    });
}DTTable ComputationPrevious(const DTTable &everything,int stepsBack)
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
    
    DTTable stopped2Back;
    
    for (timeNumber=1;timeNumber<uniqueTimes.Length();timeNumber++) {
        DTTable currentTime = everything.ExtractRows(time.FindValue(uniqueTimes(timeNumber)));

        int countInPrevious = int(previousTime.NumberOfRows());
        int countInCurrent = int(currentTime.NumberOfRows());

        // find location of patches at the current time in the previous time value
        DTIntArray locationInPrevious = TrackChanges(previousTime,currentTime);
        
                
        // First pass is to connect what we can to the previous track numbers
        DTMutableIntArray labelsAtCurrentTime(countInCurrent);
        labelsAtCurrentTime = -1;
        
        DTMutableCharArray foundInPrevious(countInPrevious);
        foundInPrevious = 0;
        DTMutableCharArray foundInCurrent(countInCurrent);
        foundInCurrent = 0;
        
        DTMutableIntArray notFoundInCurrent(countInCurrent);
        int howManyNotFoundInCurrent = 0;
        DTIntArray trackNumberPrevious = DTTableColumnNumber(previousTime("ID")).IntVersion();

        //int howManyFoundInPrevious = 0;
        int howManySet = 0;
        for (int i=0;i<countInCurrent;i++) {
            int loc = locationInPrevious(i);
            if (loc<0) {
                // not matched with the previous time, so keep it at -1
                notFoundInCurrent(howManyNotFoundInCurrent++) = i;
            }
            else {
                // found in the previous time
                foundInPrevious(loc) = 1;
                labelsAtCurrentTime(i) = trackNumberPrevious(loc);
                howManySet++;
                //foundInCurrent(i) = 1;
                //foundInPrevious(loc) = 1;
                //howManyFoundInPrevious++;
            }
        }
        
        // What is not connected yet has an opportunity to connect to 2 steps back (and 3,4 etc in the future)
        notFoundInCurrent = TruncateSize(notFoundInCurrent,howManyNotFoundInCurrent);
        DTTable notConnectedInCurrent = currentTime.ExtractRows(notFoundInCurrent);
        
        // The rows in previousTime that were not connected will be stopped2Back in the next iteration.
        DTMutableIntArray indexNotFoundInPrevious(howManySet);
        int posInPrev = 0;
        for (int i=0;i<countInPrevious;i++) {
            if (foundInPrevious(i)==0) indexNotFoundInPrevious(posInPrev++) = i;
        }
        DTTable newStopped2Back = previousTime.ExtractRows(TruncateSize(indexNotFoundInPrevious,posInPrev));

        
        // Look to see if there is a match to entries 2 steps back (before the previous table)
        if (howManySet!=countInCurrent && stopped2Back.NumberOfRows()!=0) {
            // Missing some, look at indices two step back
            DTIntArray locationIn2Back = TrackChanges(stopped2Back,notConnectedInCurrent);
            
            DTMutableIntArray notFound2Back(howManyNotFoundInCurrent);
            int howManyNotFound2Back = 0;
            int countIn2Back = int(stopped2Back.NumberOfRows());
            DTMutableCharArray foundIn2Back(countIn2Back);
            foundIn2Back = 0;
            
            trackNumberPrevious = DTTableColumnNumber(stopped2Back("ID")).IntVersion();
            howManySet = 0;

            for (int i=0;i<howManyNotFoundInCurrent;i++) {
                int loc = locationIn2Back(i);
                if (loc<0) {
                    // Didn't match with two steps back, keep it at 1
                    notFoundInCurrent(howManyNotFound2Back++) = i;
                }
                else {
                    foundIn2Back(loc) = 1;
                    labelsAtCurrentTime(notFoundInCurrent(i)) = trackNumberPrevious(loc);
                    howManySet++;
                }
            }
        }
        
        // Anything that is still not connected in labelsAtCurrentTime should get a new number
        for (int i=0;i<countInCurrent;i++) {
            if (labelsAtCurrentTime(i)<0) {
                labelsAtCurrentTime(i) = nextTrack++;
            }
        }

        
        // Prepare for the next iteration
        stopped2Back = newStopped2Back;
        
        previousTime = currentTime.Append(CreateTableColumn("ID",labelsAtCurrentTime));

        
        // Put in the label IDs
        for (int i=0;i<countInCurrent;i++) {
            finalTrackNumbers(startOfEntries+i) = labelsAtCurrentTime(i);
        }
        startOfEntries += countInCurrent;
    }
    
    // Table is a list of columns
    return DTTable({
        everything("time"),
        everything("Value"),
        everything("Mask"),
        CreateTableColumn("ID",finalTrackNumbers)
    });
}
