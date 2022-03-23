// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#include "computation.h"

#include "DTArguments.h"
#include "DTTimer.h"
#include "DTDataFile.h"
#include "DTError.h"

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTSet<DTImage> everything;
    DTTable spots;
    double time, pixels;
    int timeback, timeforward;
    DTMutableList<std::string> _channelNames;

    {
        // Inside a scope so that the data files will be closed before the computation starts.
        DTDataFile inputDataFile("Input.dtbin",DTFile::ReadOnly);
        if (inputDataFile.IsOpen()==false) {
            std::cerr << "No input file found. Might have to save input for debugging." << std::endl;
        }
        DTDataFile variableDataFile;

        variableDataFile = DTDataFile("everything.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"everything",everything);
        // Need the channel names for the output structure
        std::string _cName = "SeqInfo_everything_E_";
        int _howM = variableDataFile.ReadNumber(_cName+"N");
        _channelNames = DTMutableList<std::string>(_howM);
        for (int _count=0;_count<_howM;_count++) {
            _channelNames(_count) = variableDataFile.ReadString(_cName+DTInt2String(_count+1)+"N");
        }

        variableDataFile = DTDataFile("spots.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"spots",spots);
        time = inputDataFile.ReadNumber("time");
        timeback = inputDataFile.ReadNumber("timeback");
        timeforward = inputDataFile.ReadNumber("timeforward");
        pixels = inputDataFile.ReadNumber("pixels");
    }

    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    DTMutableSet<DTImage> output(outputFile,"Var");
    Computation(everything,spots,time,timeback,timeforward,pixels,output);

    if (DTHowManyErrors()>0) outputFile.Save(DTHowManyErrors(),"ErrorCount"); // For error logging

    {
        // Structure information for the set
        std::string baseName = "SeqInfo_Var";
        std::string eName = baseName+"_E";
        std::string pName = baseName+"_P";

        // Structure for parameters
        outputFile.Save("time",pName+"_1N");
        outputFile.Save("Number",pName+"_1T");
        outputFile.Save("intensity",pName+"_2N");
        outputFile.Save("Number",pName+"_2T");
        outputFile.Save(2,pName+"_N");

        // Structure for element
        for (int _count=0;_count<_channelNames.Length();_count++) {
            outputFile.Save(_channelNames(_count),eName+"_"+DTInt2String(_count+1)+"N");
        }
        outputFile.Save((int)_channelNames.Length(),eName+"_N");
        outputFile.Save("Image",eName);

        outputFile.Save("Image Set","Seq_Var");
    }
    // To speed up reading.
    outputFile.SaveIndex();

    return 0;
}
