#ifndef Dipole_hpp
#define Dipole_hpp

#include "DTDataStorage.h"
#include "DTDoubleArray.h"
#include "DTImage.h"
#include "DTIntArray.h"
#include "DTPoint2D.h"


//////////////////////////////////////////////////////////////////////////////
//    DipoleInfo
//////////////////////////////////////////////////////////////////////////////

struct DipoleInfo {
    DTImage input;
    DTImage first_peak;
    DTImage first_residual;
    DTImage second_peak;
    DTImage second_residual;
    DTPoint2D first_center;
    double first_base;
    double first_height;
    DTPoint2D second_center;
    double second_base;
    double second_height;
    double distance;
    double R2first;
    double RMSEfirst;
    double R2second;
    double RMSEsecond;
    double flagFirst;
    double flagSecond;
    double ratio;

    void pinfo(void) const;
    void pinfoIndent(std::string) const;
    static void WriteStructure(DTDataStorage &,std::string);
};

extern void Write(DTDataStorage &,std::string name,const DipoleInfo &);
extern void Read(DTDataStorage &,std::string name,DipoleInfo &);
extern void WriteOne(DTDataStorage &,std::string name,const DipoleInfo &);

#endif /* Dipole_hpp */
