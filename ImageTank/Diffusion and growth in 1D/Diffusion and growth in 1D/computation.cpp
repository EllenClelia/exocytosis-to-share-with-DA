#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTProgress.h"

#include "DTSeriesTable.h"
#include "DTSeriesOf.h"

#include "DTTable.h"
#include "localstructures.h"

double GrowthRate(double u,double v_max,double k)
{
    return v_max*(2/(1+exp(-k*u))-1);
}

void Computation(const DTTable &initial,double D,double dt,double until,
                 int stride,double maxVal,double v_max,double k,
                 double fluxOnLeft,
                 DTDataFile &output) // Write all output to this file
{
    // x
    DTTableColumnNumber xColumn = initial("x");
    DTMutableDoubleArray xValues = xColumn.DoubleVersion().Copy();
    DTTableColumnNumber values = initial("y");
    DTMutableDoubleArray yValues = values.DoubleVersion().Copy();
    
    DTSeriesTable result(output,"Output");
    DTSeriesOf<Values> valuesOut(output,"Values");

    ssize_t N = xValues.Length();
    double dx = xValues(1)-xValues(0);
    double x0 = xValues(0)-dx*0.5;

    DTProgress progress;
    
    DTMutableDoubleArray gradient(N+1);
    gradient = 0;

    double time = 0;
    result.Add(initial,time);
    
    // Overflow pixel
    double overflowValue = 0.0; // Total amount
    double overflowFraction = 0.0; // Fraction of a pixel that I've moved to the right (

    Values whatIsLeft;
    whatIsLeft.fraction = overflowFraction;
    whatIsLeft.overflow = overflowValue;
    whatIsLeft.rightHandSide = x0+(N+overflowFraction)*dx;
    whatIsLeft.u = yValues(N-1);
    valuesOut.Add(whatIsLeft,time);

    ssize_t i;
    int count = 0;
    
    DTMutableDoubleArray outputTime(1000);
    DTMutableDoubleArray outputVelocity(1000);
    DTMutableDoubleArray outputDT(1000);
    DTMutableDoubleArray outputU(1000);
    DTMutableDoubleArray outputLength(1000);
    int posInVelocityTable = 0;

    while (time<until) {
        
        // Compute the update flux
        for (i=1;i<N;i++) {
            gradient(i) = (yValues(i)-yValues(i-1))/dx;
        }
        
        // The assumption is that the overflow pixel has the same value as yValues(N-1) so there is
        // no flux over. That means that gradient(N) = 0.
        // But the overflow pixel has a width fraction that can increase (deal with decrease, i.e. gradient<0 later).
    // Might need to limit the time step, can't grow the front by more than a pixel
        double timestep = dt;

        // Compute the flow rate based on the value of yValues(N-1) (before update)
        double growthRate = GrowthRate(yValues(N-1),v_max,k);
        // cerr << growthRate << endl;
        if (growthRate*timestep>dx) {
            timestep = dx/growthRate;
        }
        whatIsLeft.speed = growthRate;
        
        // Record information for each time step
        if (outputTime.Length()==posInVelocityTable) {
            outputTime = IncreaseSize(outputTime);
            outputVelocity = IncreaseSize(outputVelocity);
            outputDT = IncreaseSize(outputDT);
            outputU = IncreaseSize(outputU);
            outputLength = IncreaseSize(outputLength);
        }
        outputTime(posInVelocityTable) = time;
        outputVelocity(posInVelocityTable) = growthRate;
        outputDT(posInVelocityTable) = timestep;
        outputU(posInVelocityTable) = yValues(N-1);
        outputLength(posInVelocityTable) = x0+(N+overflowFraction)*dx;
        posInVelocityTable++;

        double testSumBefore = overflowValue;
        for (i=0;i<N;i++) testSumBefore += dx*yValues(i);

        // Move between cells
        for (i=0;i<N;i++) {
            yValues(i) += D*(gradient(i+1)-gradient(i))*timestep/dx;
        }
        
        if (fluxOnLeft) {
            // Add material to the left side
            yValues(0) += fluxOnLeft*timestep/dx;
        }

        // Cell yValues(N-1) is getting bigger if growthRate>0
        // The total amount is now spread over a box with a greater width (1+overflowFraction)*dx
        // Know that growthRate*timestep/dx ≤ 1
        
        // Now handle the rightward movement, grows by growthRate to the right
        
        // Conservation of mass means y(N-1)*dx + overflowValue needs to be stay constant
        double newFraction = overflowFraction + growthRate*timestep/dx;
        double totalMass = yValues(N-1)*dx+overflowValue; // Should be conserved Conserved
        
        // Compute newAverage such that
        // newAverage*dx+newOverflowValue = yValues(N-1)*dx+overflowValue
        // newAverage is such that newAverage*(1+newFraction)*dx = yValues(N-1)*dx + overflowValue = totalMass
        
        double newAverage = totalMass/((1+newFraction)*dx);
        double newOverflowValue = newAverage*newFraction*dx;
        
        overflowValue = newOverflowValue;
        if (overflowValue<-1e-10) {
            DTErrorMessage("Should not be negative");
            overflowValue = 0;
        }
        else if (overflowValue<0) {
            overflowValue = 0;
        }
        
        overflowFraction = newFraction;
        yValues(N-1) = newAverage;
        
        double testSumAfter = overflowValue;
        for (i=0;i<N;i++) testSumAfter += dx*yValues(i);
        
        if (fabs(testSumBefore-testSumAfter)>1e-6) {
            DTErrorMessage("");
            // count = -1;
        }

        double testSumAfterChange = testSumAfter;
        if (overflowFraction>1) {
            // Add a grid point
            if (yValues.Length()==N) {
                yValues = IncreaseSize(yValues,N);
                xValues = IncreaseSize(xValues,N);
                gradient = IncreaseSize(gradient,N);
                gradient = 0;
            }
            xValues(N) = xValues(N-1) + dx;
            yValues(N) = yValues(N-1); // Copied right now, will change soon after.
            overflowValue -= newAverage*dx;
            overflowFraction -= 1.0;
            N++;
            
            testSumAfterChange = overflowValue;
            for (i=0;i<N;i++) testSumAfterChange += dx*yValues(i);
        }
        
        time += timestep;
        count++;
        
        if (count%stride==0) {
            DTTable saveThis({
                CreateTableColumn("x",TruncateSize(xValues,N)),
                CreateTableColumn("y",TruncateSize(yValues,N))
            });
            
            result.Add(saveThis,time);
            
            whatIsLeft.fraction = overflowFraction;
            whatIsLeft.overflow = overflowValue;
            whatIsLeft.rightHandSide = x0+(N+overflowFraction)*dx;
            whatIsLeft.u = yValues(N-1);
            valuesOut.Add(whatIsLeft,time);
        }
        
        progress.UpdatePercentage(time/until);
    }
    
    DTTable velocityTable({
        CreateTableColumn("time",TruncateSize(outputTime,posInVelocityTable)),
        CreateTableColumn("velocity",TruncateSize(outputVelocity,posInVelocityTable)),
        CreateTableColumn("dt",TruncateSize(outputDT,posInVelocityTable)),
        CreateTableColumn("u",TruncateSize(outputU,posInVelocityTable)),
        CreateTableColumn("length",TruncateSize(outputLength,posInVelocityTable))
    });
    WriteOne(output, "velocity", velocityTable);
}

void ComputationWithOverflow(const DTTable &initial,double D,double dt,double until,
                 int stride,double maxVal,
                 DTDataFile &output) // Write all output to this file
{
    // x
    DTTableColumnNumber xColumn = initial("x");
    DTDoubleArray xValues = xColumn.DoubleVersion();
    DTTableColumnNumber values = initial("y");
    DTMutableDoubleArray yValues = values.DoubleVersion().Copy();
    
    DTSeriesTable result(output,"Output");
    DTSeriesOf<Values> valuesOut(output,"Values");
    
    ssize_t N = xValues.Length();
    double dx = xValues(1)-xValues(0);
    double x0 = xValues(0)-dx*0.5;

    DTProgress progress;
    
    DTMutableDoubleArray gradient(N+1);
    gradient = 0;

    double time = 0;
    result.Add(initial,time);
    
    Values whatIsLeft;
    //whatIsLeft.left = 0.0;
    //whatIsLeft.added = 0.0;
    valuesOut.Add(whatIsLeft,time);

    ssize_t i;
    int count = 0;
    double absorbed = 0;
    while (time<until) {
        
        // Compute the update flux
        for (i=1;i<N;i++) {
            gradient(i) = (yValues(i)-yValues(i-1))/dx;
        }
        
        // Move between cells
        for (i=0;i<N;i++) {
            yValues(i) += D*(gradient(i+1)-gradient(i))*dt/dx;
        }
        
        // Last gradient that is computed is gradient(N-2)
        // That means that gradient(N-1) should have been computed, but isn't. That means that this flows out
        // double missingGradient = (yValues(N-1)-yValues(N-2))/dx;
        double valueAtEnd = yValues(N-1);
        double absorbThis = std::min(maxVal,valueAtEnd);
        absorbed += absorbThis*dx;
        //whatIsLeft.added = valueAtEnd;
        //whatIsLeft.left = absorbed;
        yValues(N-1) = yValues(N-1)-absorbThis;
        
        time += dt;
        count++;
        
        if (count%stride==0) {
            DTTable saveThis({
                CreateTableColumn("x",xValues),
                CreateTableColumn("y",yValues)
            });
            
            result.Add(saveThis,time);
            
            valuesOut.Add(whatIsLeft,time);
        }
        
        progress.UpdatePercentage(time/until);
    }
}

void ComputationNoFlux(const DTTable &initial,double D,double dt,double until,
                       int stride,DTDataFile &output) // Write all output to this file
{
    // x
    DTTableColumnNumber xColumn = initial("x");
    DTDoubleArray xValues = xColumn.DoubleVersion();
    DTTableColumnNumber values = initial("y");
    DTMutableDoubleArray yValues = values.DoubleVersion().Copy();
    
    DTSeriesTable result(output,"Output");
    DTSeriesOf<Values> valuesOut(output,"Values");
    
    ssize_t N = xValues.Length();
    double dx = xValues(1)-xValues(0);
    double x0 = xValues(0)-dx*0.5;

    DTProgress progress;
    
    DTMutableDoubleArray gradient(N+1);
    gradient = 0;

    double time = 0;
    result.Add(initial,time);
    
    Values whatIsLeft;
    //whatIsLeft.left = 0.0;
    //whatIsLeft.added = 0.0;
    valuesOut.Add(whatIsLeft,time);

    ssize_t i;
    int count = 0;
    double absorbed = 0;
    while (time<until) {
        
        // Compute the update flux
        for (i=1;i<N-1;i++) {
            gradient(i) = (yValues(i)-yValues(i-1))/dx;
        }
        
        // Move between cells
        for (i=0;i<N;i++) {
            yValues(i) += D*(gradient(i+1)-gradient(i))*dt/dx;
        }
        
        // Last gradient that is computed is gradient(N-2)
        // That means that gradient(N-1) should have been computed, but isn't. That means that this flows out
        // double missingGradient = (yValues(N-1)-yValues(N-2))/dx;
        absorbed += yValues(N-2);
        //whatIsLeft.added = yValues(N-2);
        //whatIsLeft.left = absorbed*dx;
        yValues(N-2) = 0;
        
        time += dt;
        count++;
        
        if (count%stride==0) {
            DTTable saveThis({
                CreateTableColumn("x",xValues),
                CreateTableColumn("y",yValues)
            });
            
            result.Add(saveThis,time);
            
            valuesOut.Add(whatIsLeft,time);
        }
        
        progress.UpdatePercentage(time/until);
    }
}

