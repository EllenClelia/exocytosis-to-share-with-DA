//
//  Adhesion.cpp
//  Prune centers
//
//  Created by David Adalsteinsson on 5/29/25.
//  Copyright © 2025 Visual Data Tools, Inc. All rights reserved.
//

#include "Adhesion.hpp"

double Interpolate(const DTDoubleArray &data,DTPoint2D at)
{
    int m = int(data.m());
    int n = int(data.n());

    double i = floor(at.x);
    double j = floor(at.y);
    if (i<0) i = 0;
    if (j<0) j = 0;
    if (i>=m-2) i = m-2;
    if (j>=n-2) j = n-2;
    
    double dx = at.x-i;
    double dy = at.y-j;

    double val = data(i,j)*(1-dx)*(1-dy) + data(i+1,j)*dx*(1-dy) + data(i,j+1)*(1-dx)*dy + data(i+1,j+1)*dx*dy;

    return val;
}

DTTable InterpolateSegment(const DTImageChannel &magnitude,const DTMesh2DGrid &grid,DTPoint2D Pc,DTPoint2D Qc,int N)
{
    DTDoubleArray data = magnitude.DoubleArray();
    
    DTMutableDoubleArray pointList(2,N+1);
    DTMutableDoubleArray arclengthList(N+1);
    DTMutableDoubleArray valueList(N+1);

    // First and last value
    double h = 1.0/N;
    
    DTPoint2D P = grid.SpaceToGrid(Pc);
    DTPoint2D Q = grid.SpaceToGrid(Qc);

    DTPoint2D deltaC = (Qc-Pc)/N;
    DTPoint2D delta = (Q-P)/N;

    //double maxAboveEnds = std::max(firstValue,lastValue);
    //double minBelowEnds = std::min(firstValue,lastValue);
    
    for (int ptN=0;ptN<=N;ptN++) {
        DTPoint2D at = P + delta*ptN;
        double value = Interpolate(data,at);
        DTPoint2D atC = Pc + deltaC*ptN;
        pointList(0,ptN) = atC.x;
        pointList(1,ptN) = atC.y;
        valueList(ptN) = value;
        arclengthList(ptN) = ptN*h;
    }

    return DTTable({CreateTableColumn("point",DTPointCollection2D(pointList)),CreateTableColumn("arc",arclengthList),CreateTableColumn("value",valueList)});
}

DTTable InterpolateSegmentWide(const DTImageChannel &magnitude,const DTMesh2DGrid &grid,DTPoint2D Pc,DTPoint2D Qc)
{
    DTDoubleArray data = magnitude.DoubleArray();
    
    // First and last value
    
    DTPoint2D P = grid.SpaceToGrid(Pc);
    DTPoint2D Q = grid.SpaceToGrid(Qc);
    
    double len = Distance(P,Q);
    int N = int(std::max(len*6,20.0));
    double h = 1.0/N;

    DTMutableDoubleArray pointList(2,N+1);
    DTMutableDoubleArray arclengthList(N+1);
    DTMutableDoubleArray valueList(N+1);
    
    // Go in the normal direction
    DTPoint2D normal(Q.y-P.y,P.x-Q.x);
    normal = normal/Norm(normal);

    DTPoint2D delta = (Q-P)/N;

    //double maxAboveEnds = std::max(firstValue,lastValue);
    //double minBelowEnds = std::min(firstValue,lastValue);

    DTPoint2D minPoint(NAN,NAN);

    for (int ptN=0;ptN<=N;ptN++) {
        DTPoint2D at = P + delta*ptN;
        
        // Go in the normal direction
        double minValue = INFINITY;
        for (int offset=-15;offset<=15;offset++) {
            DTPoint2D test = at + normal*(offset/10.0);
            double value = Interpolate(data,test);
            if (value<minValue) {
                minValue = value;
                minPoint = test;
            }
        }

        DTPoint2D atC = grid.GridToSpace(minPoint);
        pointList(0,ptN) = atC.x;
        pointList(1,ptN) = atC.y;
        valueList(ptN) = minValue;
        
        arclengthList(ptN) = ptN*h;
    }

    return DTTable({CreateTableColumn("point",DTPointCollection2D(pointList)),CreateTableColumn("arc",arclengthList),CreateTableColumn("value",valueList)});
}

double Variation(const DTTable &interpolated)
{
    DTTableColumnNumber interpolatedValues = interpolated("value");
    
    ssize_t N = interpolatedValues.NumberOfRows()-1;
    
    double h = 1.0/N;
    
    double firstValue = interpolatedValues(0);
    double lastValue = interpolatedValues(N);

    double maxDeviation = 0;
    double minDeviation = 0;

    for (int ptN=0;ptN<=N;ptN++) {
        // DTPoint2D at = P + delta*ptN;
        // double value = Interpolate(data,at);
        double value = interpolatedValues(ptN);
        
        // Do a linear interpolation between the values at the start&end and see how much the
        // the intensity differst from that
        double prediction = firstValue + ptN*h*(lastValue-firstValue);
        double difference = value-prediction;
        
        if (difference<minDeviation) minDeviation = difference;
        if (difference>maxDeviation) maxDeviation = difference;
    }
    
    return (maxDeviation-minDeviation);
}

