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
#include "DTDictionary.h"

DTImage GaussianFilter(const DTImage &image,double sigma);
DTImageChannel GaussianFilter(const DTImageChannel &channel,double sigma);
DTFloatArray GaussianFilter(const DTFloatArray &image,double sigma);
DTDoubleArray GaussianFilter(const DTDoubleArray &image,double sigma);

DTFloatArray MedianOfStack(const DTFloatArray &stack,ssize_t slices);
DTImage MedianOfImages(const DTList<DTImage> &);

double ComputeR2(const DTDoubleArray &yValuesList,const DTDoubleArray &fitValuesList);
double ComputeRMSE(const DTDoubleArray &yValuesList,const DTDoubleArray &fitValuesList);

// Part of quantifying an event. Callers will extract different parts.
struct QuantifyEvent {
    double average;
    double width;
    
    // fit with a + b*exp(-c*(t-shift)) for x>shift, just the constant a+b before.
    int shift; // If the real event starts earlier this is <0. The shift from the function above
    
    // The function fit.
    int delay;
    double base;
    double spike;
    double decay;
    
    double R2;
    DTTable histogram;
  
    DTTable piecewiseFitResults;
    
    DTTable pointsUsedForFit;
};

QuantifyEvent Quantify(const DTSet<DTImage> &,const DTDictionary &parameters);

//DTTable

// Find a local maxima from an image using a peak fit.
struct LocalPeak
{
    DTPoint2D center;
    double height;
    double base;
    double width;
    int failureMode;
    double R2; // The R2 value for the fit, globally speaking
    double RMSE;
    // 0 means everything is ok,
    // 1 means optimization failed
    // 2 means went out of bounds
};

LocalPeak FindGaussianPeak(const DTImage &image,const DTDictionary &);
LocalPeak FindMaximumPeak(const DTImage &image,int channel);

DTPoint2D FindLocalMaxima(const DTDoubleArray &values,double &maxV);


#endif /* Utilities_hpp */
