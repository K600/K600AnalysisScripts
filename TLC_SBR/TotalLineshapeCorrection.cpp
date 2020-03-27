//----------------------------------------------------------------------------------------------------
//
//  TotalLineshapeCorrection: a generalised code for iterative lineshape corrections
//
//  Author: K.C.W. Li
//  Institutions at the time of writing this code: Stellenbosch University and iThemba LABS (2017)
//
//----------------------------------------------------------------------------------------------------

#include "TotalLineshapeCorrection.hh"

void TotalLineshapeCorrection() {
    
    
    //------------------------
    //      Run Numbers
    //------------------------
    
    //------------------------
    //      PR194: WE1
    //      12C target, 0.3 mg/cm^2
    //std::vector<int> runNumbers{1090, 1114, 1115, 1120, 1121, 1125, 1131, 1132, 1138, 1156, 1157};
    std::vector<int> runNumbers{1138, 1156, 1157};
    
    SetRunNumbers(runNumbers);
        
    SetDirectoryPath("../../../sortedNew/");
    
    //--------------------------------
    //      Initialise Variables
    InitialiseVariables();
    
    //==================================================//
    //              LINESHAPE CORRECTION                //
    //==================================================//
    
    //--------------------------------------------------------------
    //      Define the peaks used for the lineshape correction
    //      DefineLineshapeCorrectionPeak(double a, double b, int c)
    //--------------------------------------------------------------
    //      a: Lower limit for X1pos fit range
    //      b: Upper limit for X1pos fit range
    //      *c: Polynomial background order: -1 => no polynomial.
    //      *Optional argument: ommitance results in default 1st order poynomial background.
    
    //--------------------------------------------------------------
    //      Experimental resolution Llimits (optional)
    //      SetExperimentalResolutionLimits(double expResMin, double expResMax)
    //--------------------------------------------------------------
    //      expResMin: Lower limit for the Gaussian standard deviation
    //      expResMax: Lower limit for the Gaussian standard deviation
    //      By default: the Gaussian standard deviation parameter is not limited
    
    //--------------------------------------------------------------
    //      Set the number of rebin attempts (optional)
    //      SetNRebinAttempts(int nRebinAttempts)
    //--------------------------------------------------------------
    //      nRebinAttempts: the number of rebin attempts
    //      The largest bin height (of a given peak) is a chosen metric to test whether a peak is valid for fitting.
    //      The nRebinAttempts parameter determines the amount of times the relevant histogram can
    //      be rebinned (by a factor two) in order to pass this minimum bin height metric.
    //      By default: nRebinAttempts = 0.

    //--------------------------------------------------------------
    //      Set the maximum reduced chi square (optional)
    //      SetMaximumReducedChiSquare(double maxRedChiSquare)
    //--------------------------------------------------------------
    //      maxRedChiSquare: the maximum reduced chi square that is accepted as a valid result
    //      The reduced chi square is a chosen metric to test whether a fit of a peak is valid.
    //      By default: maxRedChiSquare = 5.0
    //      Whilst the default maxRedChiSquare value may appear gratuitous, maxRedChiSquare>>1 is often neccesary
    //      since the intrinsic lineshape of some excited states/resonances are not always Gaussian.
    //      The assymmetry from target-related energy loss can also cause a non-negligible deviation from a Gaussian lineshape

    //------------------------------------------------------------------------------------
    //      Total Lineshape Correction
    //      TLC(int correctionType, int nDivisions_TLC, double min_TLC, double max_TLC)
    //------------------------------------------------------------------------------------
    //      correctionType: the parameter with thich the correction is performed. 0 => thetSCAT correction, 1 => Y1 correction
    //      nDivisions_TLC: N divisions for the amount of cuts for the correction parameter
    //      min_TLC: lower limit for the correction parameter
    //      max_TLC: upper limit for the correction parameter

    //------------------------------------------------------------------------------------
    //      Define a peak for mapping (post Total Lineshape Correction)
    //      DefineMappingPeak(double min_X1posC, double max_X1posC)
    //------------------------------------------------------------------------------------
    //      nDivisions_TLC: N divisions for the amount of cuts for the correction parameter

    
    //==================================================//
    //              PEAK POSITION MAPPING               //
    //==================================================//

    //------------------------------------------------------------------------------------
    //      Define a peak for mapping (post Total Lineshape Correction)
    //      DefineMappingPeak(double a, double b)
    //------------------------------------------------------------------------------------
    //      a: Lower limit for X1posC fit range
    //      b: Upper limit for X1posC fit range
    
    //------------------------------------------------------------------------------------
    //      Define a peak to test the peak mean after the peak position mapping (optional)
    //      PostCorrectionPeakFit(double a, double b)
    //------------------------------------------------------------------------------------
    //      These test fits are applied to a histogram with all the analysed data runs which
    //      have undergone their respective peak mapping corrections.
    //      a: Lower limit for X1posC fit range
    //      b: Upper limit for X1posC fit range

    //------------------------------------------------------------------------------------
    //      Set the data merging method
    //      SetDataMergingMethod(int dataMergeMethod) (optional)
    //      dataMergeMethod: the data merging method
    //      By default: dataMergeMethod = 0.
    //------------------------------------------------------------------------------------
    //      dataMergeMethod = 0: Implicit data merging. If one first merges runs {0, 1} and then subsequently
    //      merges runs {1, 2}, then runs {0, 1, 2} will all share the same data and will have
    //      the peak mapping correction factors.
    //      Consequently, runs {0, 1, 2} will all have identical peak mapping correction factors.
    //------------------------------------------------------------------------------------
    //      dataMergeMethod = 1: Explicit data merging. If one first merges runs {0, 1} and then subsequently
    //      merges runs {1, 2}, then run {0} will exclusively share data from run {1}, run {1}
    //      will share data with runs {0, 2} and run 2 will exclusively share data from run {1}.
    //      Consequently, runs {0, 1, 2} will all have different peak mapping correction factors.
    
    //------------------------------------------------------------------------------------
    //      Set the maximum polynomial order for the peak mapping
    //      SetMaxMappingOrder(int a)
    //------------------------------------------------------------------------------------
    //      a: Maximum polynomial order for the peak mapping

    //------------------------------------------------------------------------------------
    //      Map Peak Positions
    //      MapPeakPositions()
    //------------------------------------------------------------------------------------
    //      The function to calculate the peak mapping parameters for the pre-defined peaks
    //      with the DefineMappingPeak() function.
    //      After the mapping parameters have been calculated, the corrected data is then
    //      summated and each peak defined with PostCorrectionPeakFit() is fitted.
    
    //------------------------
    //      The typical ID for the type of lineshape correction
    //      For example:
    //      SetNumberOfExtrapolationPoints(X, 3, 2);
    //      SetMinimumPeakHeight(X, 5);
    //      SetNExtensionPoints(X, 1, 0, 2.0);
    //
    //      X == 0: theta_SCAT
    //      X == 1: Y1
    //      X == 2: tofCal, this is a calibration time-of-flight where the tof of the locus is centered over 0. A drawback of polynomials is their divergence towards -infty/infty
    //      X == 3: X1thCal, raytrace angle in X1 wireplane
    //      X == 4: U1thCal, raytrace angle in U1 wireplane
    //
    //      In principle, any correction variable can be added with any ID, as long as you are consistent.
    //      It should be noted that this is not always the same between all the experiments as some experiments required new correction variables to be used (due to issues).
    
    //----------------------------------------------------
    //      Runs 1138-1157
    
    SetNRebinAttempts(1);
    SetMaximumReducedChiSquare(40.0);
    ExtendToLimits(true);
    SetVisuals_thetaSCAT(400, -2.5, 2.5);
    SetVisuals_Y1(400, -50.0, 50.0);
    SetVisuals_tofCal(50, -5.0, 5.0);
    SetVisuals_X1thCal(200, -10.0, 5.0);
    SetVisuals_U1thCal(200, -15.0, 15.0);
    
    SetNumberOfExtrapolationPoints(0, 3, 2);
    SetNumberOfExtrapolationPoints(1, 3, 3);
    SetNumberOfExtrapolationPoints(2, 3, 3);
    SetNumberOfExtrapolationPoints(3, 3, 3);
    
    SetMinimumPeakHeight(0, 5);
    SetMinimumPeakHeight(1, 5);
    SetMinimumPeakHeight(2, 5);
    SetMinimumPeakHeight(3, 5);
    
    SetNExtensionPoints(0, 1, 0, 2.0);
    SetNExtensionPoints(1, 5, 1, 5.0);
    SetNExtensionPoints(2, 0, 0, 2.0);
    

    //----------------------------------------------------
    SetNumberOfExtrapolationPoints(0, 3, 3);
    DefineLineshapeCorrectionPeak(615.0, 640.0);
    DefineLineshapeCorrectionPeak(695.0, 715.0);
    TLC(0, 8, -1.7, 1.8);

    
    //--------------------------------------------------------
    //      PrintTLCResults()
    //--------------------------------------------------------
    //      Saves the output of the lineshape correction
    PrintTLCResults();
    
    c1->Close();
    
}


