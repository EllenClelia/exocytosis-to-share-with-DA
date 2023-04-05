// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#ifndef IT_structure_h
#define IT_structure_h

#include "DTDataStorage.h"
#include "DTDoubleArray.h"
#include "DTImage.h"
#include "DTIntArray.h"
#include "DTMask2D.h"
#include "DTRegion1D.h"
#include "DTTable.h"


//////////////////////////////////////////////////////////////////////////////
//    Group
//////////////////////////////////////////////////////////////////////////////

struct Group {
    DTTable values;
    DTTable L;
    DTTable L2;
    DTImage pdf;
    DTTable sample;

    void pinfo(void) const;
    void pinfoIndent(std::string) const;
    static void WriteStructure(DTDataStorage &,std::string);
};

extern void Write(DTDataStorage &,std::string name,const Group &);
extern void Read(DTDataStorage &,std::string name,Group &);
extern void WriteOne(DTDataStorage &,std::string name,const Group &);

#endif /* IT_structure_h */ 
