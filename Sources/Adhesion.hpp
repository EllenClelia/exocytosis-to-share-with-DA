//
//  Adhesion.hpp
//  Prune centers
//
//  Created by David Adalsteinsson on 5/29/25.
//  Copyright © 2025 Visual Data Tools, Inc. All rights reserved.
//

#ifndef Adhesion_hpp
#define Adhesion_hpp

#include "DTTable.h"
#include "DTImage.h"

double Interpolate(const DTDoubleArray &data,DTPoint2D at);

DTTable InterpolateSegment(const DTImageChannel &magnitude,const DTMesh2DGrid &grid,DTPoint2D Pc,DTPoint2D Qc,int N);
DTTable InterpolateSegmentWide(const DTImageChannel &magnitude,const DTMesh2DGrid &grid,DTPoint2D Pc,DTPoint2D Qc);

double Variation(const DTTable &interpolated);

#endif /* Adhesion_hpp */
