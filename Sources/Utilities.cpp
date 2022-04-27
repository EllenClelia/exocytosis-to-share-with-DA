//
//  Utilities.cpp
//  GaussFit
//
//  Created by David Adalsteinsson on 2/14/22.
//  Copyright © 2022 Visual Data Tools, Inc. All rights reserved.
//

#include "Utilities.h"

DTImage GaussianFilter(const DTImage &image,double sigma)
{
    DTFloatArray arr = image(0).FloatArray();
    
    if (arr.IsEmpty()) return image;
    
    const ssize_t m = (int)arr.m();
    const ssize_t n = (int)arr.n();
    
    DTMutableFloatArray firstArray(m,n);
    
    ssize_t i;
    ssize_t padBy = int(floor(sigma*3.5));
    if (padBy==0) padBy = 1;

    if (arr.m()<=1+2*padBy || arr.n()<=1+2*padBy) {
        DTErrorMessage("GaussianFilter(image,sigma)","sigma too large");
        return DTImage();
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
    
    DTMutableList<DTImageChannel> channels(1);
    channels(0) = DTImageChannel("field",returnArray);

    return DTImage(image.Grid(),channels);
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
                toReturn(i,j) = (list(half)+list(half)+1)*0.5;
            }
            else {
                toReturn(i,j) = list(half);
            }
        }
    }
    
    return toReturn;
}
