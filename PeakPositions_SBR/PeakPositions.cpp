//------------------------------------------------------------------------------------------------------------
//
//  PeakPositions: a generalised code to calculate mapping parameters (per run) such that the peaks match up
//
//  Author: K.C.W. Li
//  Institutions at the time of writing this code: Stellenbosch University and iThemba LABS (2017)
//
//------------------------------------------------------------------------------------------------------------

#include "PeakPositions.hh"

void PeakPositions() {
    
    
    //------------------------
    //      Run Numbers
    //------------------------
    
    //------------------------
    //      PR194: WE1
    //      12C target, 0.3 mg/cm^2
    //std::vector<int> runNumbers{1090, 1114, 1115, 1120, 1121, 1125, 1131, 1132, 1138, 1156, 1157};
    std::vector<int> runNumbers{1090, 1114, 1115, 1120, 1121, 1125, 1131, 1138, 1156, 1157};
    
    
    //std::vector<int> runNumbers{1090, 1120};

    //      Removed run 1132 (insufficient data)
    
    SetRunNumbers(runNumbers);
        
    //SetDirectoryPath("../../../sortedNoXMapping/");
    SetDirectoryPath("../../../sortedNew/");
    
    //  This is something simple I just implemented so you can place your basic cut here instead of hardcoding in the header file.
    //  Haven't tested it yet but it should be working.
    SetPrimitiveCut("cut_Pad1vsTOF && goodSbrEvent && U1flag==0 && X2flag==0 && Y1>-18 && Y1<18 && thetaSCAT>-2.0 && thetaSCAT<2.0");
    
    //--------------------------------
    //      Initialise Variables
    InitialiseVariables();
    

    //--------------------------------------------------------------------
    //      Peak Position Mapping
    //--------------------------------------------------------------------
    //      This is to provide offsets between different runs
    //      Unneccesary to run if only a single run is being used.
    //--------------------------------------------------------------------
    //      DefineMinimizationPeak(xMin, xMax, *referencePosition)
    //      xMin: Minimum peak fit range for the mapping peak
    //      xMax: Maximum peak fit range for the mapping peak
    //      *referencePosition: optional reference position
    
    //----------------------------------------------------
    //      1090-1157
    DefineMappingPeak(430.0, 470.0);
    DefineMappingPeak(625.0, 633.0);
    DefineMappingPeak(704.5, 708.5);

    
    MapPeakPositions();

    

    c1->Close();
    
}


