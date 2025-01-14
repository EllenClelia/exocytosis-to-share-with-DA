// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#ifndef IT_structure_h
#define IT_structure_h

#include "DTDataStorage.h"
#include "DTDoubleArray.h"
#include "DTImage.h"
#include "DTIntArray.h"
#include "DTPoint2D.h"
#include "DTSet.h"
#include "DTTable.h"


//////////////////////////////////////////////////////////////////////////////
//    PeakGroup
//////////////////////////////////////////////////////////////////////////////

struct PeakGroup {
    double height;
    double base;
    double threshold;
    DTImage peak;
    DTPoint2D center;
    DTTable intensities;

    void pinfo(void) const;
    void pinfoIndent(std::string) const;
    static void WriteStructure(DTDataStorage &,std::string);
};

extern void Write(DTDataStorage &,std::string name,const PeakGroup &);
extern void Read(DTDataStorage &,std::string name,PeakGroup &);
extern void WriteOne(DTDataStorage &,std::string name,const PeakGroup &);

#endif /* IT_structure_h */ 
