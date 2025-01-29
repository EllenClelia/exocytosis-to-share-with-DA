// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#include "computation.h"

#include "DTArguments.h"
#include "DTTimer.h"
#include "DTDataFile.h"
#include "DTError.h"

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTImage grid;
    DTTable points;
    int radius;

    {
        // Inside a scope so that the data files will be closed before the computation starts.
        DTDataFile inputDataFile("Input.dtbin",DTFile::ReadOnly);
        if (inputDataFile.IsOpen()==false) {
            std::cerr << "No input file found. Might have to save input for debugging." << std::endl;
        }
        DTDataFile variableDataFile;

        variableDataFile = DTDataFile("grid.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"grid",grid);
        variableDataFile = DTDataFile("points.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"points",points);
        radius = inputDataFile.ReadNumber("radius");
    }

    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    //DTTimer timer;
    //timer.Start();
    DTTable output = Computation(grid,points,radius);

    //timer.Stop(); // Use timer.Time() to get the elapsed time
    if (DTHowManyErrors()>0) outputFile.Save(DTHowManyErrors(),"ErrorCount"); // For error logging

    WriteOne(outputFile,"Var",output);
    // To speed up reading.
    outputFile.SaveIndex();

    return 0;
}
