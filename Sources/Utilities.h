//
//  Utilities.h
//  GaussFit
//
//  Created by David Adalsteinsson on 2/14/22.
//  Copyright © 2022 Visual Data Tools, Inc. All rights reserved.
//

#ifndef Utilities_hpp
#define Utilities_hpp

#include "DTImage.h"
#include "DTSet.h"
#include "DTTable.h"
#include "DTList.h"

DTImage GaussianFilter(const DTImage &image,double sigma);
DTImageChannel GaussianFilter(const DTImageChannel &channel,double sigma);
DTFloatArray GaussianFilter(const DTFloatArray &image,double sigma);
DTDoubleArray GaussianFilter(const DTDoubleArray &image,double sigma);

DTFloatArray MedianOfStack(const DTFloatArray &stack,ssize_t slices);
DTImage MedianOfImages(const DTList<DTImage> &);

double ComputeR2(const DTDoubleArray &xValuesList,const DTDoubleArray &yValuesList,const DTDoubleArray &fitValuesList);

// Part of quantifying an event
struct QuantifyEvent {
    double average;
    double width;
    
    // fit with a + b*exp(-c*t)
    int shift;
    double base;
    double spike;
    double decay;
    
    double R2;
    DTTable histogram;
    
};
QuantifyEvent Quantify(const DTSet<DTImage> &,int channel);

//DTTable

#endif /* Utilities_hpp */
