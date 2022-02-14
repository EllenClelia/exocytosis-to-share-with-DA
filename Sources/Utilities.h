//
//  Utilities.hpp
//  GaussFit
//
//  Created by David Adalsteinsson on 2/14/22.
//  Copyright © 2022 Visual Data Tools, Inc. All rights reserved.
//

#ifndef Utilities_hpp
#define Utilities_hpp

#include "DTImage.h"

DTImage GaussianFilter(const DTImage &image,double sigma);
DTFloatArray MedianOfStack(const DTFloatArray &stack,ssize_t slices);

#endif /* Utilities_hpp */
