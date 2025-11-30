#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

// google "how to compute manders correlation coefficient"

Correlation Computation(const DTImage &image,const DTMask2D &mask)
{
    // DTMesh2DGrid grid = image.Grid(); // The underlying grid
    // Extract channels from image
    // DTImageChannel channel = image("red") or image(0)
    // DTImageChannel channel = image("green") or image(1)
    // If you want to access this as a float even for 8,32 or 64 bit
    // channel = ConvertToFloat(channel); // no cost if already float
    // DTFloatArray values = channel.FloatArray(); // Complains if not stored as floats

    // DTCharArray onoff = mask.MaskArray(); // array of 0 and 1 (included)
    
    DTImage asDouble = ConvertToDouble(image);
    
    DTDoubleArray red = asDouble(0).DoubleArray();
    DTDoubleArray green = asDouble(1).DoubleArray();
    
    // sum(R*G)/sqrt(sum(R^2)*sum(G^2))
    // if G = c*R, i.e. correlated
    // then this is 1.
    // i.e. sum(R*cR)/sqrt(sum(R^2)*sum(c^2R^2)) = 1.
    
    double sumRG = 0.0;
    double sumRSq = 0.0;
    double sumGSq = 0.0;
    
    ssize_t ij,mn = red.Length();
    
    if (mask.MaskArray().IsEmpty()) {
        for (ij=0;ij<mn;ij++) {
            double r = red(ij);
            double g = green(ij);
            sumRG += r*g;
            sumRSq += r*r;
            sumGSq += g*g;
        }
    }
    else {
        DTCharArray maskArray = mask.MaskArray();
        for (ij=0;ij<mn;ij++) {
            if (maskArray(ij)) {
                double r = red(ij);
                double g = green(ij);
                sumRG += r*g;
                sumRSq += r*r;
                sumGSq += g*g;
            }
        }
    }

    Correlation toReturn;

    toReturn.coeff = sumRG/sqrt(sumRSq*sumGSq);

    return toReturn;
}
