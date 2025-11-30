#include "DTDoubleArray.h"
#include "DTIntArray.h"
#include "DTDataStorage.h"

#ifndef IT_localstructure_h
#define IT_localstructure_h

//////////////////////////////////////////////////////////////////////////////
//    Values
//////////////////////////////////////////////////////////////////////////////

struct Values {
    double overflow;
    double fraction;
    double rightHandSide;
    double u;
    double speed;

    void pinfo(void) const;
    void pinfoIndent(std::string) const;
    static void WriteStructure(DTDataStorage &,std::string);
};

extern void Write(DTDataStorage &,std::string name,const Values &);
extern void Read(DTDataStorage &,std::string name,Values &);
extern void WriteOne(DTDataStorage &,std::string name,const Values &);

#endif /* IT_structure_h */
