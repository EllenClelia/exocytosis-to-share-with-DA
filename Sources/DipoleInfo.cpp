#include "DipoleInfo.h"

#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

#include "Utilities.h"

//////////////////////////////////////////////////////////////////////////////
//    DipoleInfo
//////////////////////////////////////////////////////////////////////////////

void DipoleInfo::pinfo(void) const
{
    pinfoIndent("");
}

void DipoleInfo::pinfoIndent(std::string pad) const
{
    std::cerr << pad << "input = "; input.pinfo();
    std::cerr << pad << "first_peak = "; first_peak.pinfo();
    std::cerr << pad << "first_residual = "; first_residual.pinfo();
    std::cerr << pad << "second_peak = "; second_peak.pinfo();
    std::cerr << pad << "second_residual = "; second_residual.pinfo();
    std::cerr << pad << "first_center = "; first_center.pinfo();
    std::cerr << pad << "first_base = " << first_base << std::endl;
    std::cerr << pad << "first_height = " << first_height << std::endl;
    std::cerr << pad << "second_center = "; second_center.pinfo();
    std::cerr << pad << "second_base = " << second_base << std::endl;
    std::cerr << pad << "second_height = " << second_height << std::endl;
    std::cerr << pad << "distance = " << distance << std::endl;
    std::cerr << pad << "R2first = " << R2first << std::endl;
    std::cerr << pad << "RMSEfirst = " << RMSEfirst << std::endl;
    std::cerr << pad << "R2second = " << R2second << std::endl;
    std::cerr << pad << "RMSEsecond = " << RMSEsecond << std::endl;
    std::cerr << pad << "flagFirst = " << flagFirst << std::endl;
    std::cerr << pad << "flagSecond = " << flagSecond << std::endl;
    std::cerr << pad << "ratio = " << ratio << std::endl;
}

void DipoleInfo::WriteStructure(DTDataStorage &output,std::string name)
{
    // Structure for "input"
    output.Save("input",name+"_1N");
    output.Save("height",name+"_1T_1N");
    output.Save(1,name+"_1T_N");
    output.Save("Image",name+"_1T");

    // Structure for "first peak"
    output.Save("first peak",name+"_2N");
    output.Save("height",name+"_2T_1N");
    output.Save(1,name+"_2T_N");
    output.Save("Image",name+"_2T");

    // Structure for "first residual"
    output.Save("first residual",name+"_3N");
    output.Save("height",name+"_3T_1N");
    output.Save(1,name+"_3T_N");
    output.Save("Image",name+"_3T");

    // Structure for "second peak"
    output.Save("second peak",name+"_4N");
    output.Save("height",name+"_4T_1N");
    output.Save(1,name+"_4T_N");
    output.Save("Image",name+"_4T");

    // Structure for "second residual"
    output.Save("second residual",name+"_5N");
    output.Save("height",name+"_5T_1N");
    output.Save(1,name+"_5T_N");
    output.Save("Image",name+"_5T");

    // Structure for "first center"
    output.Save("first center",name+"_6N");
    output.Save("Point2D",name+"_6T");

    // Structure for "first base"
    output.Save("first base",name+"_7N");
    output.Save("Number",name+"_7T");

    // Structure for "first height"
    output.Save("first height",name+"_8N");
    output.Save("Number",name+"_8T");

    // Structure for "second center"
    output.Save("second center",name+"_9N");
    output.Save("Point2D",name+"_9T");

    // Structure for "second base"
    output.Save("second base",name+"_10N");
    output.Save("Number",name+"_10T");

    // Structure for "second height"
    output.Save("second height",name+"_11N");
    output.Save("Number",name+"_11T");

    // Structure for "distance"
    output.Save("distance",name+"_12N");
    output.Save("Number",name+"_12T");

    // Structure for "R2first"
    output.Save("R2first",name+"_13N");
    output.Save("Number",name+"_13T");

    // Structure for "RMSEfirst"
    output.Save("RMSEfirst",name+"_14N");
    output.Save("Number",name+"_14T");

    // Structure for "R2second"
    output.Save("R2second",name+"_15N");
    output.Save("Number",name+"_15T");

    // Structure for "RMSEsecond"
    output.Save("RMSEsecond",name+"_16N");
    output.Save("Number",name+"_16T");

    // Structure for "flagFirst"
    output.Save("flagFirst",name+"_17N");
    output.Save("Number",name+"_17T");

    // Structure for "flagSecond"
    output.Save("flagSecond",name+"_18N");
    output.Save("Number",name+"_18T");

    // Structure for "ratio"
    output.Save("ratio",name+"_19N");
    output.Save("Number",name+"_19T");

    output.Save(19,name+"_N");
    output.Save("DipoleInfo",name+"_Name");
    output.Save("Group",name);
}

void Write(DTDataStorage &output,std::string name,const DipoleInfo &var)
{
    Write(output,name+"_input",var.input);
    Write(output,name+"_first peak",var.first_peak);
    Write(output,name+"_first residual",var.first_residual);
    Write(output,name+"_second peak",var.second_peak);
    Write(output,name+"_second residual",var.second_residual);
    Write(output,name+"_first center",var.first_center);
    output.Save(var.first_base,name+"_first base");
    output.Save(var.first_height,name+"_first height");
    Write(output,name+"_second center",var.second_center);
    output.Save(var.second_base,name+"_second base");
    output.Save(var.second_height,name+"_second height");
    output.Save(var.distance,name+"_distance");
    output.Save(var.R2first,name+"_R2first");
    output.Save(var.RMSEfirst,name+"_RMSEfirst");
    output.Save(var.R2second,name+"_R2second");
    output.Save(var.RMSEsecond,name+"_RMSEsecond");
    output.Save(var.flagFirst,name+"_flagFirst");
    output.Save(var.flagSecond,name+"_flagSecond");
    output.Save(var.ratio,name+"_ratio");
    Write(output,name,DTDoubleArray());
}

void WriteOne(DTDataStorage &output,std::string name,const DipoleInfo &var)
{
    Write(output,name,var);
    output.Save("Group","Seq_"+name);
    DipoleInfo::WriteStructure(output,"SeqInfo_"+name);
}

void Read(DTDataStorage &input,std::string name,DipoleInfo &var)
{
    Read(input,name+"_input",var.input);
    Read(input,name+"_first peak",var.first_peak);
    Read(input,name+"_first residual",var.first_residual);
    Read(input,name+"_second peak",var.second_peak);
    Read(input,name+"_second residual",var.second_residual);
    Read(input,name+"_first center",var.first_center);
    var.first_base = input.ReadNumber(name+"_first base");
    var.first_height = input.ReadNumber(name+"_first height");
    Read(input,name+"_second center",var.second_center);
    var.second_base = input.ReadNumber(name+"_second base");
    var.second_height = input.ReadNumber(name+"_second height");
    var.distance = input.ReadNumber(name+"_distance");
    var.R2first = input.ReadNumber(name+"_R2first");
    var.RMSEfirst = input.ReadNumber(name+"_RMSEfirst");
    var.R2second = input.ReadNumber(name+"_R2second");
    var.RMSEsecond = input.ReadNumber(name+"_RMSEsecond");
    var.flagFirst = input.ReadNumber(name+"_flagFirst");
    var.flagSecond = input.ReadNumber(name+"_flagSecond");
    var.ratio = input.ReadNumber(name+"_ratio");
}

DTMutableDoubleArray EvaluateDipoleFit(const DTMesh2DGrid &grid,const LocalPeak &fit)
{
    // Base + scale*exp(-((x-xC)^2+(y-yC)^2)/2radius^2);
    int m = (int)grid.m();
    int n = (int)grid.n();
    
    DTMutableDoubleArray returnArray(m,n);
    
    DTPoint2D center = fit.center;
    center = grid.SpaceToGrid(center);
    double x0 = center.x;
    double y0 = center.y;
    
    // toReturn.width = 2*sqrt(2*log(2))*fabs(returned("radius")); // FWHM - full width at half maximumum
    
    double radius = fit.width/(2*sqrt(2*log(2)));
    double base = fit.base;
    double scale = fit.height-base;

    /*
    DTMutableDictionary guesses;
    guesses("x0") = maxI;
    guesses("y0") = maxJ;
    guesses("scale") = maxV-minV;
    guesses("base") = minV;
    guesses("radius") = 5;
    
    
    double base = constants("base");
    double scale = constants("scale");
    double x0 = constants("x0");
    double y0 = constants("y0");
    double radius = constants("radius");
     */

    int i,j,ij;
    double x,y,arg;
    double C = 1.0/(2*radius*radius);
    ij = 0;
    double *returnArrayD = returnArray.Pointer();
    for (j=0;j<n;j++) {
        y = (j-y0);
        for (i=0;i<m;i++) {
            x = (i-x0);
            arg = -(x*x+y*y)*C;
            returnArrayD[ij] = base + scale*exp(arg);
            ij++;
        }
    }
    
    return returnArray;
}

DipoleInfo ComputeDipole(const DTImage &image)
{
    if (image.NumberOfChannels()!=1) {
        DTErrorMessage("ComputeDipole(Image)","Requires just a single channel");
    }
    
    DipoleInfo toReturn;
    toReturn.input = image;
    
    DTMutableDictionary parameters;
    parameters("channel") = 0;
    parameters("R2 for Peak") = 0;

    // Compute the main peak
    LocalPeak fit = FindGaussianPeak(image,parameters);
    DTMutableDoubleArray values = EvaluateDipoleFit(image.Grid(),fit);
    toReturn.first_peak = DTImage(image.Grid(),values,"height");

    // Compute the fit and subtract it from the channel
    DTDoubleArray leftAfterFirstFitRemoved = values-ConvertToDouble(image(0)).DoubleArray();
    DTImage differenceImage = DTImage(image.Grid(),leftAfterFirstFitRemoved,"height");
    toReturn.first_residual = differenceImage;

    // Compute the second peak
    LocalPeak fitInDifference = FindGaussianPeak(differenceImage,parameters);
    values = EvaluateDipoleFit(image.Grid(),fitInDifference);
    toReturn.second_peak = DTImage(image.Grid(),values,"height");
    
    DTDoubleArray leftAfterSecondFitRemoved = leftAfterFirstFitRemoved-values;
    toReturn.second_residual = DTImage(image.Grid(),leftAfterSecondFitRemoved,"height");
    
    toReturn.first_base = fit.base;
    toReturn.first_height = fit.height;
    toReturn.first_center = fit.center;
    
    toReturn.second_base = fitInDifference.base;
    toReturn.second_height = fitInDifference.height;
    toReturn.second_center = fitInDifference.center;
    
    toReturn.flagFirst = fit.failureMode;
    toReturn.flagSecond = fitInDifference.failureMode;
    toReturn.R2first = fit.R2;
    toReturn.RMSEfirst = fit.RMSE;
    toReturn.R2second = fitInDifference.R2;
    toReturn.RMSEsecond = fitInDifference.RMSE;
    
    toReturn.ratio = toReturn.second_height/toReturn.first_height;
    
    toReturn.distance = Distance(toReturn.first_center,toReturn.second_center);

    return toReturn;
}
