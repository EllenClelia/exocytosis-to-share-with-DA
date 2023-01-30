#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTRandom.h"

DTTable Realization(const DTMask2D &region,
                    DTRandom &randNumber,
                    const DTDictionary &parameters)
{
    double Tmax = parameters("Tmax");
    int length = parameters("Number of events");

    DTMutableDoubleArray timeList(length);
    DTMutableDoubleArray pointList(2,length);
    
    randNumber.UniformClosed(timeList.Pointer(),0,Tmax,length);
    std::sort(timeList.Pointer(),timeList.Pointer()+length);
    
    DTMesh2DGrid grid = region.Grid();
    DTCharArray mask = region.Mask().MaskArray();
    ssize_t m = grid.m();
    ssize_t n = grid.n();
    ssize_t mn = m*n;
    double width = grid.m();
    double height = grid.n();
    double h = grid.dx();

    // Assume that most of the entries are included
    // If that is not the case I can create a list of all of the interior points
    // and then draw a offset into that list at random
    
    double t;
    ssize_t i,j,ij;
    double x,y;
    double val, product = m*n;
    for (ssize_t pos=0;pos<length;pos++) {
        // Random in space
        while (1) {
            x = randNumber.UniformHalf()*width;
            y = randNumber.UniformHalf()*height;
            i = floor(x);
            j = floor(y);
            ij = i+j*m;
            if (mask(ij)) {
                pointList(0,pos) = (i-0.5)*h;
                pointList(1,pos) = (j-0.5)*h;
                break;
            }
        }
    }

    //double p = parameters("probability");
    //double dt = parameters("dt");
    //double mean = parameters("Events per second");
    //double meanPerDT = mean*dt;
        
    //timeList = TruncateSize(timeList,pos);
    //pointList = TruncateSize(pointList,2*pos);

    return DTTable({
        CreateTableColumn("time",timeList),
        CreateTableColumn("point",DTPointCollection2D(pointList))
    });
}

void Computation(const DTMask2D &region,int seed,int count,
                 const DTDictionary &parameters,DTMutableSet<DTTable> &output)
{
    // DTCharArray onoff = region.MaskArray(); // array of 0 and 1 (included)
    // double a = parameters("a");
    // DTDoubleArray b = parameters("b")
    
    DTRandom randNumber(seed);
    
    DTMutableDoubleArray runN(count);
    // Every element should be a table. An example constructor is:
    for (int i=0;i<count;i++) {
        output.Add(Realization(region, randNumber, parameters));
        runN(i) = i+1;
    }
    
    output.Finish(DTTable({
        CreateTableColumn("run",runN)
     }));
}
