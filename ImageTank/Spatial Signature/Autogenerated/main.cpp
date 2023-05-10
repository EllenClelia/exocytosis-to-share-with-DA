// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#include "computation.h"

#include "DTArguments.h"
#include "DTTimer.h"
#include "DTDataFile.h"
#include "DTError.h"

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTMask2D mask;
    DTRegion1D rRange, yRange;
    int count, seed, rCount, yCount, runs;

    {
        // Inside a scope so that the data files will be closed before the computation starts.
        DTDataFile inputDataFile("Input.dtbin",DTFile::ReadOnly);
        if (inputDataFile.IsOpen()==false) {
            std::cerr << "No input file found. Might have to save input for debugging." << std::endl;
        }

        Read(inputDataFile,"mask",mask);
        count = inputDataFile.ReadNumber("count");
        seed = inputDataFile.ReadNumber("seed");
        Read(inputDataFile,"rRange",rRange);
        Read(inputDataFile,"yRange",yRange);
        rCount = inputDataFile.ReadNumber("rCount");
        yCount = inputDataFile.ReadNumber("yCount");
        runs = inputDataFile.ReadNumber("runs");
    }

    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    //DTTimer timer;
    //timer.Start();
    Group output = Computation(mask,count,seed,rRange,yRange,rCount,yCount,
                               runs);

    //timer.Stop(); // Use timer.Time() to get the elapsed time
    if (DTHowManyErrors()>0) outputFile.Save(DTHowManyErrors(),"ErrorCount"); // For error logging

    WriteOne(outputFile,"Var",output);
    // The structure, to make it easy to open the output file
    Group::WriteStructure(outputFile,"SeqInfo_Var");
    outputFile.Save("Group","Seq_Var");

    // To speed up reading.
    outputFile.SaveIndex();

    return 0;
}