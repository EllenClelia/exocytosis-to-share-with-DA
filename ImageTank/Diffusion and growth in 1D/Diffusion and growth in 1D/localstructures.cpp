//
//  localstructures.cpp
//  Diffusion and growth in 1D
//
//  Created by David Adalsteinsson on 10/7/25.
//  Copyright © 2025 Visual Data Tools, Inc. All rights reserved.
//

#include "localstructures.h"

//////////////////////////////////////////////////////////////////////////////
//    Values
//////////////////////////////////////////////////////////////////////////////

void Values::pinfo(void) const
{
    pinfoIndent("");
}

void Values::pinfoIndent(std::string pad) const
{
    std::cerr << pad << "overflow = " << overflow << std::endl;
    std::cerr << pad << "fraction = " << fraction << std::endl;
    std::cerr << pad << "rightHandSide = " << rightHandSide << std::endl;
    std::cerr << pad << "u = " << u << std::endl;
    std::cerr << pad << "speed = " << speed << std::endl;
}

void Values::WriteStructure(DTDataStorage &output,std::string name)
{
    // Structure for "overflow"
    output.Save("overflow",name+"_1N");
    output.Save("Number",name+"_1T");

    // Structure for "fraction"
    output.Save("fraction",name+"_2N");
    output.Save("Number",name+"_2T");

    // Structure for "rightHandSide"
    output.Save("rightHandSide",name+"_3N");
    output.Save("Number",name+"_3T");

    // Structure for "u"
    output.Save("u",name+"_4N");
    output.Save("Number",name+"_4T");

    // Structure for "speed"
    output.Save("speed",name+"_5N");
    output.Save("Number",name+"_5T");

    output.Save(5,name+"_N");
    output.Save("Values",name+"_Name");
    output.Save("Group",name);
}

void Write(DTDataStorage &output,std::string name,const Values &var)
{
    output.Save(var.overflow,name+"_overflow");
    output.Save(var.fraction,name+"_fraction");
    output.Save(var.rightHandSide,name+"_rightHandSide");
    output.Save(var.u,name+"_u");
    output.Save(var.speed,name+"_speed");
    Write(output,name,DTDoubleArray());
}

void WriteOne(DTDataStorage &output,std::string name,const Values &var)
{
    Write(output,name,var);
    output.Save("Group","Seq_"+name);
    Values::WriteStructure(output,"SeqInfo_"+name);
}

void Read(DTDataStorage &input,std::string name,Values &var)
{
    var.overflow = input.ReadNumber(name+"_overflow");
    var.fraction = input.ReadNumber(name+"_fraction");
    var.rightHandSide = input.ReadNumber(name+"_rightHandSide");
    var.u = input.ReadNumber(name+"_u");
    var.speed = input.ReadNumber(name+"_speed");
}
