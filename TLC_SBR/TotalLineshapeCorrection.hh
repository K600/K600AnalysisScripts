
//======================================//
//      TotalLineshapeCorrection.hh     //
//======================================//

#ifndef TOTALLINESHAPECORRECTION_H
#define TOTALLINESHAPECORRECTION_H

#include "time.h"
#include "iostream"
#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TTree.h"
#include "TCutG.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TMath.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

#include "../../../Analysis/cuts/CUTS.hh"

char name[512];
char condition[512];
char command[512];
std::string directoryPath;

std::vector<int> runNumbers;
int nrOfRuns;

std::string primitiveCut;

//----------------------------------
int nBins_X1pos;
double min_X1pos, max_X1pos;

int nBins_X1posC;
double min_X1posC, max_X1posC;

int nBins_thetaSCAT;
double min_thetaSCAT, max_thetaSCAT;

int nBins_Y1;
double min_Y1, max_Y1;

int nBins_tofCal;
double min_tofCal, max_tofCal;

int nBins_X1thCal;
double min_X1thCal, max_X1thCal;

int nBins_U1thCal;
double min_U1thCal, max_U1thCal;

int nBins_thetaRaytraceCal;
double min_thetaRaytraceCal, max_thetaRaytraceCal;

//----------------------------------
std::unique_ptr<TCanvas> c1;

std::vector<std::unique_ptr<TChain>> DATAChains;
std::vector<TEventList*> eventLists;

//----------------------------------
//      Peak position mapping
int nMappingPeaks;
TF1 *fMapping;
std::vector<double> referencePeakPositions;
std::vector<std::tuple<double, double>> defaultX1posCMappingRanges;
std::vector<std::tuple<double, double>> postCorrectionPeakFits;
std::vector<std::tuple<int, int, double, double>> X1posCRanges;
std::vector<std::vector<double>> peakPositions;
std::vector<std::vector<double>> peakPositionErrors;
std::vector<std::vector<double>> mappingParameters;
std::vector<double> initialCorrectionPeakMeans;
std::vector<std::vector<int>> mergedDataRunsStorage;
std::vector<std::vector<int>> mergedDataRuns;

bool setMaxMappingOrder;
int maxMappingOrder;
int dataMergingMethod;

//--------------------------------
//      Lineshape correction
std::vector<std::tuple<int, int, int, double>> nExtensionPointsPerCorrectionType;
bool extendToLimits;
std::vector<std::tuple<int, int, int>> extrapolationPoints;

int nRebinAttempts;

bool setMaximumRedChiSquare;
double maximumRedChiSquare;

int nLineshapeCorrectionPeaks;
std::vector<std::tuple<double, double, int>> lineshapeCorrectionPeaks;
std::tuple<double, double, bool> expRes;
std::vector<std::vector<std::shared_ptr<TGraphErrors>>> grTLCCoefficientsStorage;
std::vector<std::vector<std::shared_ptr<TGraphErrors>>> grTLCExtendedCoefficientsStorage;
std::vector<std::tuple<int, double>> minimumPeakHeights;

//----------------------------------------
//  finalLineshapeCorrectionParameters[a][b][c]
//  a: Corresponds to the correction iteration (thetaSCAT/Y1)
//  b: Corresponds to the order of the polynomial coefficient for the correction itself
//  c: Corresponds to the order of the polynomial coefficient used to calculate b
std::vector<int> finalLineshapeCorrectionTypes;
std::vector<std::vector<std::vector<double>>> finalLineshapeCorrectionParameters;

void GenerateGoodEventLists();
double PreviousLineshapeCorrections(double fX1pos, std::vector<double> correctionVariables);
double fPeakPosition_voigt(double *x, double *par);

//--------------------------------------------
void SetPrimitiveCut(std::string s)
{
    primitiveCut = s;
}

//--------------------------------------------
void InitialiseVariables()
{
    //------------------------------------------------
    nBins_X1pos = 2000;
    min_X1pos = 0.0;
    max_X1pos = 1000.0;
    
    nBins_X1posC = 2000;
    min_X1posC = 0.0;
    max_X1posC = 1000.0;
    
    nBins_thetaSCAT = 400;
    min_thetaSCAT = -3.0;
    max_thetaSCAT = 3.0;
    
    nBins_Y1 = 400;
    min_Y1 = -40.0;
    max_Y1 = 40.0;
    
    nBins_tofCal = 400;
    min_tofCal = -3.0;
    max_tofCal = 3.0;
    
    nBins_X1thCal = 400;
    min_X1thCal = -10.0;
    max_X1thCal = 10.0;
    
    nBins_U1thCal = 400;
    min_U1thCal = -10.0;
    max_U1thCal = 10.0;

    nBins_thetaRaytraceCal = 400;
    min_thetaRaytraceCal = -10.0;
    max_thetaRaytraceCal = 10.0;

    //------------------------------------------------
    c1 = std::unique_ptr<TCanvas>(new TCanvas("TotalLineshapeCorrection","",10,10,700,700));
    
    nLineshapeCorrectionPeaks = 0;
    nMappingPeaks = 0;
    extendToLimits = false;
    nRebinAttempts = 0;
    setMaximumRedChiSquare = false;
    maximumRedChiSquare = 0.0;
    expRes = std::make_tuple(0.0, 0.0, false);
    
    setMaxMappingOrder = false;
    maxMappingOrder = 0;
    dataMergingMethod = 0;
    
    initialiseCuts();
    
    GenerateGoodEventLists();
}

//--------------------------------------------
void SetExperimentalResolutionLimits(double min, double max)
{
    std::get<0>(expRes) = min;
    std::get<1>(expRes) = max;
    std::get<2>(expRes) = true;
}

//--------------------------------------------
void SetVisuals_thetaSCAT(int nBins, double min, double max)
{
    nBins_thetaSCAT = nBins;
    min_thetaSCAT = min;
    max_thetaSCAT = max;
}

//--------------------------------------------
void SetVisuals_Y1(int nBins, double min, double max)
{
    nBins_Y1 = nBins;
    min_Y1 = min;
    max_Y1 = max;
}

//--------------------------------------------
void SetVisuals_tofCal(int nBins, double min, double max)
{
    nBins_tofCal = nBins;
    min_tofCal = min;
    max_tofCal = max;
}

//--------------------------------------------
void SetVisuals_X1thCal(int nBins, double min, double max)
{
    nBins_X1thCal = nBins;
    min_X1thCal = min;
    max_X1thCal = max;
}

//--------------------------------------------
void SetVisuals_U1thCal(int nBins, double min, double max)
{
    nBins_U1thCal = nBins;
    min_U1thCal = min;
    max_U1thCal = max;
}

//--------------------------------------------
void SetVisuals_ThetaRaytraceCal(int nBins, double min, double max)
{
    nBins_thetaRaytraceCal = nBins;
    min_thetaRaytraceCal = min;
    max_thetaRaytraceCal = max;
}

//--------------------------------------------
void DefineMappingPeak(double min_X1posC, double max_X1posC)
{
    defaultX1posCMappingRanges.push_back(std::make_tuple(min_X1posC, max_X1posC));
    nMappingPeaks++;
}

//--------------------------------------------
void DefineMappingPeak(double min_X1posC, double max_X1posC, double referenceValue)
{
    defaultX1posCMappingRanges.push_back(std::make_tuple(min_X1posC, max_X1posC));
    referencePeakPositions.push_back(referenceValue);
    nMappingPeaks++;
}

//--------------------------------------------
void SetNExtensionPoints(int correctionType, int nPoints, int nSkippedPoints, double errorWeighting)
{
    bool found = false;
    
    for(int i=0; i<(int) nExtensionPointsPerCorrectionType.size(); i++)
    {
        if(std::get<0>(nExtensionPointsPerCorrectionType[i])==correctionType)
        {
            found = true;
            
            std::get<1>(nExtensionPointsPerCorrectionType[i]) = nPoints;
            std::get<2>(nExtensionPointsPerCorrectionType[i]) = nSkippedPoints;
            std::get<3>(nExtensionPointsPerCorrectionType[i]) = errorWeighting;
        }
    }
    
    if(!found)
    {
        nExtensionPointsPerCorrectionType.push_back(std::make_tuple(correctionType, nPoints, nSkippedPoints, errorWeighting));
    }
}

//--------------------------------------------
void PostCorrectionPeakFit(double min_X1posC, double max_X1posC)
{
    postCorrectionPeakFits.push_back(std::make_tuple(min_X1posC, max_X1posC));
}

//--------------------------------------------
void DefineLineshapeCorrectionPeak(double min_X1pos, double max_X1pos, int polBkgOrder)
{
    lineshapeCorrectionPeaks.push_back(std::make_tuple(min_X1pos, max_X1pos, polBkgOrder));
    nLineshapeCorrectionPeaks++;
}

//--------------------------------------------
void DefineLineshapeCorrectionPeak(double min_X1pos, double max_X1pos)
{
    DefineLineshapeCorrectionPeak(min_X1pos, max_X1pos, 1);
}

//--------------------------------------------
void RemoveAllLineshapeCorrectionPeaks()
{
    lineshapeCorrectionPeaks.clear();
    nLineshapeCorrectionPeaks = 0;
}

//--------------------------------------------
void RemoveLatestLineshapeCorrectionPeak()
{
    lineshapeCorrectionPeaks.pop_back();
    nLineshapeCorrectionPeaks--;
}

//--------------------------------------------
void RemoveLineshapeCorrectionPeak(int peakN)
{
    if(peakN<(int) lineshapeCorrectionPeaks.size())
    {
        lineshapeCorrectionPeaks.erase(lineshapeCorrectionPeaks.begin() + peakN);
        nLineshapeCorrectionPeaks--;
    }
}

//--------------------------------------------
void RemoveLineshapeCorrectionPeak_ForHeathensWhoNeglectCNumbering(int peakN)
{
    if(peakN<=(int) lineshapeCorrectionPeaks.size())
    {
        lineshapeCorrectionPeaks.erase(lineshapeCorrectionPeaks.begin() + (peakN-1));
        nLineshapeCorrectionPeaks--;
    }
}

//--------------------------------------------
void SetRunNumbers(std::vector<int> givenRunNumbers)
{
    runNumbers = givenRunNumbers;
    nrOfRuns = (int) givenRunNumbers.size();
}

//------------------------------------------------------------------------
void SetDirectoryPath(std::string path)
{
    directoryPath = path;
}

//------------------------------------------------------------------------
void ExtendToLimits(bool a)
{
    extendToLimits = a;
}

//------------------------------------------------------------------------
void SetNRebinAttempts(int a)
{
    nRebinAttempts = a;
}

//------------------------------------------------------------------------
void SetMaximumReducedChiSquare(double a)
{
    setMaximumRedChiSquare = true;
    maximumRedChiSquare = a;
}

//------------------------------------------------------------------------
void SetNumberOfExtrapolationPoints(int correctionType, int a, int b)
{
    bool found = false;
    
    for(int i=0; i<(int) extrapolationPoints.size(); i++)
    {
        if(std::get<0>(extrapolationPoints[i])==correctionType)
        {
            found = true;
            
            if(a>=0)
            {
                std::get<1>(extrapolationPoints[i]) = a;
            }
            else
            {
                std::get<1>(extrapolationPoints[i]) = 0;
            }
            
            if(b>=0)
            {
                std::get<2>(extrapolationPoints[i]) = b;
            }
            else
            {
                std::get<2>(extrapolationPoints[i]) = 0;
            }
        }
    }
    
    if(!found)
    {
        extrapolationPoints.push_back(std::make_tuple(correctionType, a, b));
    }
}

//------------------------------------------------------------------------
void SetMinimumPeakHeight(int correctionType, double minimumPeakHeight)
{
    bool found = false;
    
    for(int i=0; i<(int) minimumPeakHeights.size(); i++)
    {
        if(std::get<0>(minimumPeakHeights[i])==correctionType)
        {
            found = true;
            std::get<1>(minimumPeakHeights[i]) = minimumPeakHeight;
        }
    }
    
    if(!found)
    {
        minimumPeakHeights.push_back(std::make_tuple(correctionType, minimumPeakHeight));
    }
}

//------------------------------------------------------------------------
void SetMaxMappingOrder(int a)
{
    if(a>=0)
    {
        setMaxMappingOrder = true;
        maxMappingOrder = a;
    }
}

//------------------------------------------------------------------------
void SetDataMergingMethod(int a)
{
    dataMergingMethod = a;
}

//------------------------------------------------------------------------
void MergeDataRuns(std::vector<int> dataRuns)
{
    for(auto i=dataRuns.begin(); i!=dataRuns.end(); i++)
    {
        for(auto j=i+1; j!=dataRuns.end();)
        {
            if((*i)==(*j))
            {
                j = dataRuns.erase(j);
            }
            else
            {
                ++j;
            }
        }
    }
    
    /*
     for(auto i=dataRuns.begin(); i!=dataRuns.end(); i++)
     {
     std::cout << "dataRuns: " << (*i) << std::endl;
     }
     */
    
    mergedDataRunsStorage.push_back(dataRuns);
}

//------------------------------------------------------------------------
void ConsolidateMergedDataRuns()
{
    if(dataMergingMethod==0)
    {
        //--------------------------------------------------------
        for(int g=0; g<(int) mergedDataRunsStorage.size(); g++)
        {
            bool mergeWithPreExistingSet = false;
            int preExistingSetNumber = 0;
            
            for(int h=0; h<(int) mergedDataRunsStorage[g].size(); h++)
            {
                //--------------------------------------------------------
                for(int i=0; i<(int) mergedDataRuns.size(); i++)
                {
                    for(int j=0; j<(int) mergedDataRuns[i].size(); j++)
                    {
                        if(mergedDataRuns[i][j] == mergedDataRunsStorage[g][h])
                        {
                            mergeWithPreExistingSet = true;
                            preExistingSetNumber = h;
                        }
                    }
                }
            }
            
            if(mergeWithPreExistingSet)
            {
                for(int k=0; k<(int) mergedDataRunsStorage[g].size(); k++)
                {
                    bool repeated = false;
                    
                    for(int l=0; l<(int) mergedDataRuns[preExistingSetNumber].size(); l++)
                    {
                        if(mergedDataRunsStorage[g][k]==mergedDataRuns[preExistingSetNumber][l])
                        {
                            repeated = true;
                        }
                    }
                    
                    if(!repeated)
                    {
                        mergedDataRuns[preExistingSetNumber].push_back(mergedDataRunsStorage[g][k]);
                    }
                }
            }
            else
            {
                mergedDataRuns.push_back(mergedDataRunsStorage[g]);
            }
        }
    }
    
    if(dataMergingMethod==1)
    {
        mergedDataRuns = mergedDataRunsStorage;
    }
    
    for(int g=0; g<(int) mergedDataRuns.size(); g++)
    {
        std::cout << "Merged data runs:" ;
        
        for(int h=0; h<(int) mergedDataRuns[g].size(); h++)
        {
            std::cout << " " << mergedDataRuns[g][h];
        }
        
        std::cout << std::endl;
    }
    
    
    /*
     for(int k=0; k<(int) dataRuns.size(); k++)
     {
     std::cout << " " << dataRuns[k];
     }
     std::cout << std::endl;
     */
}

//--------------------------------------------
void GenerateGoodEventLists()
{
    std::cout << "//=====================================================//" << std::endl;
    std::cout << "      Generating Event Lists (with specified cuts)" << std::endl;
    
    for(int i=0; i<nrOfRuns; i++)
    {
        char filePath[1028];
        
        std::cout << "\rRun (" << i+1 << "/" << nrOfRuns << ")" << std::flush;
        
        //----------------------------------------------------------------
        
        DATAChains.push_back(std::unique_ptr<TChain>(new TChain("DATA")));
        
        
        if(directoryPath.find_last_of("/") != (directoryPath.length()-1))
        {
            sprintf(filePath, "%s/sorted", directoryPath.c_str());
        }
        else
        {
            sprintf(filePath, "%ssorted", directoryPath.c_str());
        }
        
        int nMissingZeros = 0;
        if(runNumbers[i]<10000) {nMissingZeros++;}
        if(runNumbers[i]<1000) {nMissingZeros++;}
        if(runNumbers[i]<100) {nMissingZeros++;}
        if(runNumbers[i]<10) {nMissingZeros++;}
        
        for(int j=0; j<nMissingZeros; j++)
        {
            sprintf(filePath, "%s0", filePath);
        }
        
        if(runNumbers[i]<10) {sprintf(filePath, "%s%1d.root", filePath, runNumbers[i]);}
        else if(runNumbers[i]<100) {sprintf(filePath, "%s%2d.root", filePath, runNumbers[i]);}
        else if(runNumbers[i]<1000) {sprintf(filePath, "%s%3d.root", filePath, runNumbers[i]);}
        else if(runNumbers[i]<10000) {sprintf(filePath, "%s%4d.root", filePath, runNumbers[i]);}
        else if(runNumbers[i]<100000) {sprintf(filePath, "%s%5d.root", filePath, runNumbers[i]);}
        
        DATAChains[i]->Add(filePath);
        
        //----------------------------------------------------------------
        
        sprintf(name, "goodEventList_%d", i);
        sprintf(command, ">>%s", name);

//        if(runNumbers[i]>=23153 && runNumbers[i]<=23191)
//        {
//            DATAChains[i]->Draw(command, "cut_Pad1vsTOF && X1flag==0 && X2flag==0 && U2flag==0");
//        }
//        else
//        {
//            DATAChains[i]->Draw(command, "cut_Pad1vsTOF && X1flag==0");
//        }
        
        DATAChains[i]->Draw(command, primitiveCut.c_str());
        
        //DATAChains[i]->SetEventList((TEventList*)gDirectory->Get(name));
        eventLists.push_back((TEventList*)gDirectory->Get(name));
    }
    
    std::cout << std::endl;
}

//--------------------------------------------
std::vector<int> GenerateMergedDataRunIndices(int runNumber)
{
    std::vector<int> mergedDataRunIndices;
    
    for(int j=0; j<(int) mergedDataRuns.size(); j++)
    {
        bool found = false;
        
        for(int k=0; k<(int) mergedDataRuns[j].size(); k++)
        {
            if(runNumber==mergedDataRuns[j][k])
            {
                found = true;
            }
        }
        
        if(found)
        {
            for(int k=0; k<(int) mergedDataRuns[j].size(); k++)
            {
                bool validRun = false;
                int runIndex = 0;
                
                for(int l=0; l<nrOfRuns; l++)
                {
                    if(runNumbers[l] == mergedDataRuns[j][k] && mergedDataRuns[j][k]!=runNumber)
                    {
                        validRun = true;
                        runIndex = l;
                    }
                }
                
                if(validRun)
                {
                    mergedDataRunIndices.push_back(runIndex);
                }
            }
        }
    }
    
    //--------------------------------------
    //      Removing repeated indices
    for(auto j=mergedDataRunIndices.begin(); j!=mergedDataRunIndices.end(); j++)
    {
        for(auto k=j+1; k!=mergedDataRunIndices.end();)
        {
            if((*j)==(*k))
            {
                k = mergedDataRunIndices.erase(k);
            }
            else
            {
                ++k;
            }
        }
    }
    
    return mergedDataRunIndices;
}


//--------------------------------------------
//void MapPeakPositions()
//{
//    ConsolidateMergedDataRuns();
//    
//    TFile f("Output/MapPeakPositions.root", "recreate");
//    
//    std::vector<std::unique_ptr<TGraphErrors>> grMappings;
//    
//    std::vector<std::unique_ptr<TH1F>> hX1posC_thetaSCATgated;
//    std::vector<std::unique_ptr<TH1F>> hX1posC_thetaSCATgated_corrected;
//    
//    //  Default
//    nBins_X1posC = 2000;
//    min_X1posC = 0.0;
//    max_X1posC = 1000.0;
//    
//    std::unique_ptr<TH1F> hX1posC_thetaSCATgated_corrected_total(new TH1F("hX1posC_thetaSCATgated_corrected_total", "", nBins_X1posC, min_X1posC, max_X1posC));
//    std::unique_ptr<TH1F> hX1posC_thetaSCATgated_total(new TH1F("hX1posC_thetaSCATgated_total", "", nBins_X1posC, min_X1posC, max_X1posC));
//    
//    //------------------------------------------------
//    nBins_X1posC = 2000;
//    min_X1posC = 0.0;
//    max_X1posC = 1000.0;
//    
//    std::unique_ptr<TH2F> hThetaSCATvsX1pos(new TH2F("hThetaSCATvsX1pos", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_thetaSCAT, min_thetaSCAT, max_thetaSCAT));
//    std::unique_ptr<TH2F> hThetaSCATvsX1pos_TLC(new TH2F("hThetaSCATvsX1pos_TLC", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_thetaSCAT, min_thetaSCAT, max_thetaSCAT));
//    std::unique_ptr<TH2F> hThetaSCATvsX1pos_TLC_mapped(new TH2F("hThetaSCATvsX1pos_TLC_mapped", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_thetaSCAT, min_thetaSCAT, max_thetaSCAT));
//
//    //------------------------------------------------
//    nBins_X1posC = 2000;
//    min_X1posC = 0.0;
//    max_X1posC = 1000.0;
//    
//    std::unique_ptr<TH2F> hY1vsX1pos(new TH2F("hY1vsX1pos", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_Y1, min_Y1, max_Y1));
//    std::unique_ptr<TH2F> hY1vsX1pos_TLC(new TH2F("hY1vsX1pos_TLC", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_Y1, min_Y1, max_Y1));
//    std::unique_ptr<TH2F> hY1vsX1pos_TLC_mapped(new TH2F("hY1vsX1pos_TLC_mapped", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_Y1, min_Y1, max_Y1));
//    
//    //------------------------------------------------
//    nBins_X1posC = 2000;
//    min_X1posC = 0.0;
//    max_X1posC = 1000.0;
//    
//    std::unique_ptr<TH2F> hTOFCalvsX1pos(new TH2F("hTOFCalvsX1pos", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_tofCal, min_tofCal, max_tofCal));
//    std::unique_ptr<TH2F> hTOFCalvsX1pos_TLC(new TH2F("hTOFCalvsX1pos_TLC", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_tofCal, min_tofCal, max_tofCal));
//    std::unique_ptr<TH2F> hTOFCalvsX1pos_TLC_mapped(new TH2F("hTOFCalvsX1pos_TLC_mapped", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_tofCal, min_tofCal, max_tofCal));
//    
//    //------------------------------------------------
//    nBins_X1posC = 2000;
//    min_X1posC = 0.0;
//    max_X1posC = 1000.0;
//    
//    std::unique_ptr<TH2F> hX1thCalvsX1pos(new TH2F("hX1thCalvsX1pos", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_X1thCal, min_X1thCal, max_X1thCal));
//    std::unique_ptr<TH2F> hX1thCalvsX1pos_TLC(new TH2F("hX1thCalvsX1pos_TLC", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_X1thCal, min_X1thCal, max_X1thCal));
//    std::unique_ptr<TH2F> hX1thCalvsX1pos_TLC_mapped(new TH2F("hX1thCalvsX1pos_TLC_mapped", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_X1thCal, min_X1thCal, max_X1thCal));
//
//    //------------------------------------------------
//    nBins_X1posC = 2000;
//    min_X1posC = 0.0;
//    max_X1posC = 1000.0;
//    
//    std::unique_ptr<TH2F> hU1thCalvsX1pos(new TH2F("hU1thCalvsX1pos", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_U1thCal, min_U1thCal, max_U1thCal));
//    std::unique_ptr<TH2F> hU1thCalvsX1pos_TLC(new TH2F("hU1thCalvsX1pos_TLC", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_U1thCal, min_U1thCal, max_U1thCal));
//    std::unique_ptr<TH2F> hU1thCalvsX1pos_TLC_mapped(new TH2F("hU1thCalvsX1pos_TLC_mapped", "", nBins_X1posC, min_X1posC, max_X1posC, nBins_U1thCal, min_U1thCal, max_U1thCal));
//
//    //------------------------------------------------
//    //std::unique_ptr<TF1> fitPeakPosition(new TF1("fitPeakPosition", "gaus(0) + pol1(3)", 0., 1000.));
//    std::unique_ptr<TF1> fitPeakPosition(new TF1("fitPeakPosition", fPeakPosition_voigt, 0., 1000., 6));
//    fitPeakPosition->SetParameters(1.0e+01, 7.20e+02, 1.0e+00, 1.0e+00, 0.0e+00, 0.0e+00);
//    fitPeakPosition->SetLineColor(kRed);
//    fitPeakPosition->SetNpx(4*nBins_X1posC);
//    fitPeakPosition->SetParLimits(0, 0.0, 1.0e+08);
//    
//    if((int) referencePeakPositions.size() != nMappingPeaks)
//    {
//        referencePeakPositions.clear();
//    }
//    
//    double par[5];
//    
//    //================================
//    //      Mapping Fit Function
//    
//    int mappingOrder;
//    
//    if(setMaxMappingOrder && (nMappingPeaks-1)>maxMappingOrder)
//    {
//        mappingOrder = maxMappingOrder;
//    }
//    else if(!setMaxMappingOrder && nMappingPeaks>=1)
//    {
//        mappingOrder = nMappingPeaks-1;
//    }
//    
//    if(mappingOrder==1)
//    {
//        sprintf(condition, "x + [0]");
//    }
//    else
//    {
//        sprintf(condition, "pol%d(0)", mappingOrder);
//    }
//    
//    fMapping = new TF1("fMapping", condition, 0.0, 1000.0);
//    
//    for(int k=0; k<mappingOrder+1; k++)
//    {
//        fMapping->SetParameter(k, 0.0e+00);
//    }
//    
//    //========================================================================
//    std::cout << "//==============================================//" << std::endl;
//    std::cout << "      Calculating X1posC mapping parameters" << std::endl;
//    
//    for(int i=0; i<nrOfRuns; i++)
//    {
//        sprintf(name, "hX1posC_thetaSCATgated_run%d", runNumbers[i]);
//        hX1posC_thetaSCATgated.push_back(std::unique_ptr<TH1F>(new TH1F(name, "", nBins_X1posC, min_X1posC, max_X1posC)));
//        
//        sprintf(name, "hX1posC_thetaSCATgated_run%d_corrected", runNumbers[i]);
//        hX1posC_thetaSCATgated_corrected.push_back(std::unique_ptr<TH1F>(new TH1F(name, "", nBins_X1posC, min_X1posC, max_X1posC)));
//    }
//    
//    for(int i=0; i<nrOfRuns; i++)
//    {
//        std::cout << "\rRun (" << i+1 << "/" << nrOfRuns << ")" << std::flush;
//        
//        TTreeReader                 reader("DATA", DATAChains[i]->GetFile());
//        TTreeReaderValue<double>    fX1pos(reader, "X1pos.t_X1pos");
//        TTreeReaderValue<double>    fThetaSCAT(reader, "thetaSCAT.t_thetaSCAT");
//        TTreeReaderValue<double>    fY1(reader, "Y1.t_Y1");
//        
//        TTreeReaderValue<double>    fTofCal(reader, "tofCal.t_tofCal");
//        TTreeReaderValue<double>    fX1thCal(reader, "X1thCal.t_X1thCal");
//        TTreeReaderValue<double>    fU1thCal(reader, "U1thCal.t_U1thCal");
//        
//        TTreeReaderValue<double>    fX1chisq(reader, "X1chisq.t_X1chisq");
//        
//        TTreeReaderArray<double>    fxPositions_TLC(reader, "xPositions_TLC");
//        TTreeReaderArray<double>    fthetaRaytraceCals(reader, "thetaRaytraceCals");
//
//        
//        //----------------------------------------------------
//        std::vector<int> mergedDataRunIndices = GenerateMergedDataRunIndices(runNumbers[i]);
//        
//        int eN = 0;
//        while(eventLists[i]->GetEntry(eN)>0)
//        {
//            reader.SetEntry(eventLists[i]->GetEntry(eN));
//         
//            //----------------------------------------------------
//            std::vector<double> correctionVariables;
//            
//            for(auto j=finalLineshapeCorrectionTypes.begin(); j!=finalLineshapeCorrectionTypes.end(); j++)
//            {
//                if((*j)==0)
//                {
//                    correctionVariables.push_back(*fThetaSCAT);
//                }
//                else if((*j)==1)
//                {
//                    correctionVariables.push_back(*fY1);
//                }
//                else if((*j)==2)
//                {
//                    correctionVariables.push_back(*fTofCal);
//                }
//                else if((*j)==3)
//                {
//                    correctionVariables.push_back(*fX1thCal);
//                }
//                else if((*j)==4)
//                {
//                    correctionVariables.push_back(*fU1thCal);
//                }
//                else if((*j)==5)
//                {
//                    correctionVariables.push_back(fthetaRaytraceCals[0]);
//                }
//            }
//            
//            //----------------------------------------------------
//            //double originalX1pos = (*fX1pos);
//            double originalX1pos = fxPositions_TLC[0];
//            double originalX1posC = PreviousLineshapeCorrections(originalX1pos, correctionVariables);
//            
//            hX1posC_thetaSCATgated[i]->Fill(originalX1posC);
//            
//            for(int j=0; j<(int) mergedDataRunIndices.size(); j++)
//            {
//                hX1posC_thetaSCATgated[mergedDataRunIndices[j]]->Fill(originalX1posC);
//            }
//            
//            //------------------------------------------------------------
//            hX1posC_thetaSCATgated_total->Fill(originalX1posC);
//            
//            hThetaSCATvsX1pos->Fill(originalX1pos, (*fThetaSCAT));
//            hY1vsX1pos->Fill(originalX1pos, (*fY1));
//            
//            hTOFCalvsX1pos->Fill(originalX1pos, (*fTofCal));
//            hX1thCalvsX1pos->Fill(originalX1pos, (*fX1thCal));
//            hU1thCalvsX1pos->Fill(originalX1pos, (*fU1thCal));
//            
//            
//            hThetaSCATvsX1pos_TLC->Fill(originalX1posC, (*fThetaSCAT));
//            hY1vsX1pos_TLC->Fill(originalX1posC, (*fY1));
//            
//            hTOFCalvsX1pos_TLC->Fill(originalX1posC, (*fTofCal));
//            hX1thCalvsX1pos_TLC->Fill(originalX1posC, (*fX1thCal));
//            hU1thCalvsX1pos_TLC->Fill(originalX1posC, (*fU1thCal));
//            
//            
//            eN++;
//        }
//    }
//    
//    for(int i=0; i<nrOfRuns; i++)
//    {
//        TTreeReader                 reader("DATA", DATAChains[i]->GetFile());
//        TTreeReaderValue<double>    fX1pos(reader, "X1pos.t_X1pos");
//        TTreeReaderValue<double>    fThetaSCAT(reader, "thetaSCAT.t_thetaSCAT");
//        TTreeReaderValue<double>    fY1(reader, "Y1.t_Y1");
//        
//        TTreeReaderValue<double>    fTofCal(reader, "tofCal.t_tofCal");
//        TTreeReaderValue<double>    fX1thCal(reader, "X1thCal.t_X1thCal");
//        TTreeReaderValue<double>    fU1thCal(reader, "U1thCal.t_U1thCal");
//        
//        TTreeReaderValue<double>    fX1chisq(reader, "X1chisq.t_X1chisq");
//        
//        TTreeReaderArray<double>    fxPositions_TLC(reader, "xPositions_TLC");
//        TTreeReaderArray<double>    fthetaRaytraceCals(reader, "thetaRaytraceCals");
//
//        std::vector<int> mergedDataRunIndices = GenerateMergedDataRunIndices(runNumbers[i]);
//        
//        //========================================================================
//        std::vector<double> tempPeakPositions;
//        std::vector<double> tempPeakPositionErrors;
//        
//        for(int j=0; j<nMappingPeaks; j++)
//        {
//            double tempmin_X1posC, tempmax_X1posC;
//            
//            bool found = false;
//            for(auto k = X1posCRanges.begin(); k != X1posCRanges.end(); ++k)
//            {
//                if(std::get<0>((*k))==j && std::get<1>((*k))==runNumbers[i])
//                {
//                    found = true;
//                    tempmin_X1posC = std::get<2>((*k));
//                    tempmax_X1posC = std::get<3>((*k));
//                }
//            }
//            
//            if(!found)
//            {
//                tempmin_X1posC = std::get<0>(defaultX1posCMappingRanges[j]);
//                tempmax_X1posC = std::get<1>(defaultX1posCMappingRanges[j]);
//            }
//            
//            //========================================================================
//            hX1posC_thetaSCATgated[i]->GetXaxis()->SetRangeUser(tempmin_X1posC, tempmax_X1posC);
//            
//            double initialPeakHeight = hX1posC_thetaSCATgated[i]->GetBinContent(hX1posC_thetaSCATgated[i]->GetMaximumBin());
//            double initialMean = hX1posC_thetaSCATgated[i]->GetBinCenter(hX1posC_thetaSCATgated[i]->GetMaximumBin());
//            double initialStdDev = hX1posC_thetaSCATgated[i]->GetStdDev();
//            
//            fitPeakPosition->SetParameter(0, initialPeakHeight);
//            fitPeakPosition->SetParameter(1, initialMean);
//            fitPeakPosition->SetParameter(2, 0.7*initialStdDev);
//            fitPeakPosition->SetParameter(3, 0.1);
//            fitPeakPosition->SetParameter(4, 0.0e+00);
//            fitPeakPosition->SetParameter(5, 0.0e+00);
//            
//            fitPeakPosition->SetParLimits(0, 0.0e+00, 3.0*initialPeakHeight);
//            fitPeakPosition->SetParLimits(1, tempmin_X1posC, tempmax_X1posC);
//            fitPeakPosition->SetParLimits(2, 0.02*initialStdDev, 2.0*initialStdDev);
//            fitPeakPosition->SetParLimits(3, 0.0e+00, 5.0e+01);
//            
//            //========================================================================
//            hX1posC_thetaSCATgated[i]->Fit("fitPeakPosition", "Q+", "", tempmin_X1posC, tempmax_X1posC);
//            
//            tempPeakPositions.push_back(fitPeakPosition->GetParameter(1));
//            tempPeakPositionErrors.push_back(fitPeakPosition->GetParError(1));
//        }
//        
//        hX1posC_thetaSCATgated[i]->GetXaxis()->SetRangeUser(min_X1posC, max_X1posC);
//        
//        peakPositions.push_back(tempPeakPositions);
//        peakPositionErrors.push_back(tempPeakPositionErrors);
//        
//        //========================================================================
//        sprintf(name, "grMapping_run%d", runNumbers[i]);
//        grMappings.push_back(std::unique_ptr<TGraphErrors>(new TGraphErrors()));
//        grMappings[i]->SetName(name);
//        
//        for(int j=0; j<nMappingPeaks; j++)
//        {
//            if(!referencePeakPositions.empty())
//            {
//                grMappings[i]->SetPoint(j, peakPositions[i][j], referencePeakPositions[j]);
//            }
//            else
//            {
//                grMappings[i]->SetPoint(j, peakPositions[i][j], peakPositions[0][j]);
//            }
//            
//            grMappings[i]->SetPointError(j, peakPositionErrors[i][j], 0.0);
//        }
//        
//        for(int k=0; k<mappingOrder+1; k++)
//        {
//            fMapping->SetParameter(k, 0.0e+00);
//        }
//        
//        if(nMappingPeaks>1)
//        {
//            double x1, x2, y1, y2;
//            grMappings[i]->GetPoint(0, x1, y1);
//            grMappings[i]->GetPoint(grMappings[i]->GetN()-1, x2, y2);
//            
//            double a, b;
//            
//            a = (y2-y1)/(x2-x1);
//            b = (-a*x1) + y1;
//            
//            fMapping->SetParameter(0, b);
//            fMapping->SetParameter(1, a);
//        }
//        
//        //----------------------------------------------------------------------------------------------------------------
//        //      This is to mitigate a strange issue where the fit of the first run exhibits strange behaviour:
//        //      For the very first fit, the >1st order polynomial parameters appear to not play a role in optomization.
//        if(!referencePeakPositions.empty() && i==0)
//        {
//            grMappings[i]->Fit("fMapping", "Q", "", 0.0, 1000.0);
//        }
//        
//        grMappings[i]->Fit("fMapping", "Q+", "", 0.0, 1000.0);
//        
//        /*
//         double chiSquare = fMapping->GetChisquare();
//         int nDF = fMapping->GetNDF();
//         double reducedChiSquare;
//         
//         if(nDF>=1)
//         {
//         reducedChiSquare = (chiSquare/nDF);
//         std::cout << "Run: " << runNumbers[i] << ", reducedChiSquare: " << reducedChiSquare << std::endl;
//         }
//         else
//         {
//         std::cout << "Run: " << runNumbers[i] << ", chiSquare: " << chiSquare << std::endl;
//         }
//         */
//        
//        grMappings[i]->SetMarkerStyle(kCircle);
//        grMappings[i]->SetMarkerSize(3);
//        grMappings[i]->SetMarkerColor(kBlack);
//        
//        grMappings[i]->Write();
//        
//        //------------------------------------------------------------------------
//        std::vector<double> tempMappingParameters;
//        
//        for(int j=0; j<mappingOrder+1; j++)
//        {
//            tempMappingParameters.push_back(fMapping->GetParameter(j));
//        }
//        
//        mappingParameters.push_back(tempMappingParameters);
//        
//        //========================================================================
//        
//        reader.Restart();
//        
//        int eN = 0;
//        while(eventLists[i]->GetEntry(eN)>0)
//        {
//            reader.SetEntry(eventLists[i]->GetEntry(eN));
//            
//            //----------------------------------------------------
//            std::vector<double> correctionVariables;
//            
//            for(auto j=finalLineshapeCorrectionTypes.begin(); j!=finalLineshapeCorrectionTypes.end(); j++)
//            {
//                if((*j)==0)
//                {
//                    correctionVariables.push_back(*fThetaSCAT);
//                }
//                else if((*j)==1)
//                {
//                    correctionVariables.push_back(*fY1);
//                }
//                else if((*j)==2)
//                {
//                    correctionVariables.push_back(*fTofCal);
//                }
//                else if((*j)==3)
//                {
//                    correctionVariables.push_back(*fX1thCal);
//                }
//                else if((*j)==4)
//                {
//                    correctionVariables.push_back(*fU1thCal);
//                }
//                else if((*j)==5)
//                {
//                    correctionVariables.push_back(fthetaRaytraceCals[0]);
//                }
//            }
//
//            //----------------------------------------------------
//            //double originalX1pos = (*fX1pos);
//            double originalX1pos = fxPositions_TLC[0];
//            double originalX1posC = PreviousLineshapeCorrections(originalX1pos, correctionVariables);
//            
//            double correctedX1pos_TLC_mapped = 0.0;
//            
//            for(int j=0; j<(int) mappingParameters[i].size(); j++)
//            {
//                correctedX1pos_TLC_mapped += mappingParameters[i][j]*pow(originalX1posC, j);
//            }
//            
//            hX1posC_thetaSCATgated_corrected[i]->Fill(correctedX1pos_TLC_mapped);
//            
//            /*
//             for(int j=0; j<(int) mergedDataRunIndices.size(); j++)
//             {
//             hX1posC_thetaSCATgated_corrected[mergedDataRunIndices[j]]->Fill(correctedX1pos_TLC_mapped);
//             }
//             */
//            
//            hX1posC_thetaSCATgated_corrected_total->Fill(correctedX1pos_TLC_mapped);
//            
//            hY1vsX1pos_TLC_mapped->Fill(correctedX1pos_TLC_mapped, (*fY1));
//            hThetaSCATvsX1pos_TLC_mapped->Fill(correctedX1pos_TLC_mapped, (*fThetaSCAT));
//            
//            eN++;
//        }
//    }
//    std::cout << std::endl;
//    
//    /*
//     if(referencePeakPositions.empty())
//     {
//     for(int j=0; j<(int) peakPositions[0].size(); j++)
//     {
//     std::cout << "Peak positions of first run: ";
//     std::cout << std::setprecision(12) << peakPositions[0][j] << " " << std::endl;
//     }
//     }
//     */
//    
//    //--------------------------------------------------------------------------------
//    std::cout << "Peak positions of the corrected spectra: " << std::endl;
//    
//    for(int j=0; j<(int) postCorrectionPeakFits.size(); j++)
//    {
//        double tempmin_X1posC, tempmax_X1posC;
//        
//        tempmin_X1posC = std::get<0>(postCorrectionPeakFits[j]);
//        tempmax_X1posC = std::get<1>(postCorrectionPeakFits[j]);
//        
//        hX1posC_thetaSCATgated_corrected_total->GetXaxis()->SetRangeUser(tempmin_X1posC, tempmax_X1posC);
//        
//        double initialPeakHeight = hX1posC_thetaSCATgated_corrected_total->GetBinContent(hX1posC_thetaSCATgated_corrected_total->GetMaximumBin());
//        double initialMean = hX1posC_thetaSCATgated_corrected_total->GetBinCenter(hX1posC_thetaSCATgated_corrected_total->GetMaximumBin());
//        double initialStdDev = hX1posC_thetaSCATgated_corrected_total->GetStdDev();
//        
//        fitPeakPosition->SetParameter(0, initialPeakHeight);
//        fitPeakPosition->SetParameter(1, initialMean);
//        fitPeakPosition->SetParameter(2, 0.8*initialStdDev);
//        fitPeakPosition->SetParameter(3, 0.1);
//        fitPeakPosition->SetParameter(4, 0.0e+00);
//        fitPeakPosition->SetParameter(5, 0.0e+00);
//        
//        fitPeakPosition->SetParLimits(0, 0.0e+00, 3.0*initialPeakHeight);
//        fitPeakPosition->SetParLimits(1, tempmin_X1posC, tempmax_X1posC);
//        fitPeakPosition->SetParLimits(2, 0.02*initialStdDev, 2.0*initialStdDev);
//        fitPeakPosition->SetParLimits(3, 0.0e+00, 5.0e+01);
//        
//        //========================================================================
//        hX1posC_thetaSCATgated_corrected_total->Fit("fitPeakPosition", "Q+", "", tempmin_X1posC, tempmax_X1posC);
//        
//        std::cout << "Peak mean (" << j << ") = " << std::setprecision(12) << fitPeakPosition->GetParameter(1) << " " << std::endl;
//    }
//    
//    hX1posC_thetaSCATgated_corrected_total->GetXaxis()->SetRangeUser(0.0, 1000.0);
//    
//    //--------------------------------------------------------------------------------
//    std::ofstream fileX1MappingParameters("Output/PR166_2011_XMappingParameters.dat", std::ofstream::out);
//    
//    for(int i=0; i<nrOfRuns; i++)
//    {
//        fileX1MappingParameters << runNumbers[i];
//        
//        fileX1MappingParameters << " " << mappingOrder+1;
//        
//        for(int j=0; j<mappingOrder+1; j++)
//        {
//            fileX1MappingParameters << " " << mappingParameters[i][j];
//        }
//        
//        fileX1MappingParameters << "\n";
//    }
//    
//    fileX1MappingParameters << "eof\n";
//    
//    fileX1MappingParameters.close();
//    f.Write();
//}

double PreLineshapeCorrectionMap(double origX1pos, std::vector<double> mapPars)
{
    double mappedX1pos = 0.0;
    
    if((int) mapPars.size()==1)
    {
        mappedX1pos = origX1pos + mapPars[0];
    }
    else
    {
        for(int i=0; i<(int) mapPars.size(); i++)
        {
            mappedX1pos += mapPars[i]*pow(origX1pos, i);
        }
    }
    
    return mappedX1pos;
}

//--------------------------------------------
double singleLineshapeCorrection(double fX1pos, double correctionParValue, std::vector<std::vector<double>> singleLineshapeCorrectionPars)
{
    double correctedPosition = 0.0;
    
    double previousCorrectedPosition;
    std::vector<double> correctionPars_PolCoefficients;
    
    previousCorrectedPosition = fX1pos;
    
    //--------------------------------------------------------------------
    for(int j=0; j<(int) singleLineshapeCorrectionPars.size(); j++)
    {
        double correctionPar_PolCoefficient = 0.0;
        
        //--------------------------------------------------------------------
        for(int k=0; k<(int)singleLineshapeCorrectionPars[j].size(); k++)
        {
            correctionPar_PolCoefficient += singleLineshapeCorrectionPars[j][k]*pow(correctionParValue, k);

            /*
            if(correctionType==0)
            {
                correctionPar_PolCoefficient += singleLineshapeCorrectionPars[j][k]*pow(fThetaSCAT, k);
            }
            
            if(correctionType==1)
            {
                correctionPar_PolCoefficient += singleLineshapeCorrectionPars[j][k]*pow(fY1, k);
            }
            */
            
        }
        
        correctionPars_PolCoefficients.push_back(correctionPar_PolCoefficient);
    }
    
    correctedPosition = 0.0;
    
    for(int j=0; j<(int) correctionPars_PolCoefficients.size(); j++)
    {
        if((int) correctionPars_PolCoefficients.size()==1)
        {
            correctedPosition = previousCorrectedPosition + correctionPars_PolCoefficients[j];
        }
        else
        {
            correctedPosition += correctionPars_PolCoefficients[j]*pow(previousCorrectedPosition, j);
        }
    }
    
    return correctedPosition;
}

//--------------------------------------------
void PrintTLCResults()
{
    //----------------------------------------
    TFile f("Output/TLC_results.root", "recreate");
    
    std::vector<std::unique_ptr<TH1F>> hX1_TLCIteration;
    std::vector<std::unique_ptr<TH2F>> hThetaSCATvsX1pos_TLCIteration;
    std::vector<std::unique_ptr<TH2F>> hY1vsX1pos_TLCIteration;
    std::vector<std::unique_ptr<TH2F>> hTOFCalvsX1pos_TLCIteration;
    std::vector<std::unique_ptr<TH2F>> hX1thCalvsX1pos_TLCIteration;
    std::vector<std::unique_ptr<TH2F>> hU1thCalvsX1pos_TLCIteration;
    std::vector<std::unique_ptr<TH2F>> hThetaRaytraceCalvsX1pos_TLCIteration;

    for(int j=0; j<(int) finalLineshapeCorrectionParameters.size() + 1; j++)
    {
        nBins_X1pos = 5000;
        min_X1pos = 0.0;
        max_X1pos = 1000.0;
        
        sprintf(name, "hX1_TLCIteration%d", j);
        hX1_TLCIteration.push_back(std::unique_ptr<TH1F>(new TH1F(name, "", nBins_X1pos, min_X1pos, max_X1pos)));
        
        //------------------------------------------------
        nBins_X1pos = 2000;
        min_X1pos = 0.0;
        max_X1pos = 1000.0;
        nBins_thetaSCAT = 400;
        
        sprintf(name, "hThetaSCATvsX1pos_TLCIteration_%d", j);
        hThetaSCATvsX1pos_TLCIteration.push_back(std::unique_ptr<TH2F>(new TH2F(name, "", nBins_X1pos, min_X1pos, max_X1pos, nBins_thetaSCAT, min_thetaSCAT, max_thetaSCAT)));
        
        //------------------------------------------------
        nBins_X1pos = 2000;
        min_X1pos = 0.0;
        max_X1pos = 1000.0;
        nBins_Y1 = 400;
        
        sprintf(name, "hY1vsX1pos_TLCIteration_%d", j);
        hY1vsX1pos_TLCIteration.push_back(std::unique_ptr<TH2F>(new TH2F(name, "", nBins_X1pos, min_X1pos, max_X1pos, nBins_Y1, min_Y1, max_Y1)));
        
        //------------------------------------------------
        nBins_X1pos = 2000;
        min_X1pos = 0.0;
        max_X1pos = 1000.0;
        
        sprintf(name, "hTOFCalvsX1pos_TLCIteration_%d", j);
        hTOFCalvsX1pos_TLCIteration.push_back(std::unique_ptr<TH2F>(new TH2F(name, "", nBins_X1pos, min_X1pos, max_X1pos, nBins_tofCal, min_tofCal, max_tofCal)));

        //------------------------------------------------
        nBins_X1pos = 2000;
        min_X1pos = 0.0;
        max_X1pos = 1000.0;
        
        sprintf(name, "hX1thCalvsX1pos_TLCIteration_%d", j);
        hX1thCalvsX1pos_TLCIteration.push_back(std::unique_ptr<TH2F>(new TH2F(name, "", nBins_X1pos, min_X1pos, max_X1pos, nBins_X1thCal, min_X1thCal, max_X1thCal)));

        //------------------------------------------------
        nBins_X1pos = 2000;
        min_X1pos = 0.0;
        max_X1pos = 1000.0;
        
        sprintf(name, "hU1thCalvsX1pos_TLCIteration_%d", j);
        hU1thCalvsX1pos_TLCIteration.push_back(std::unique_ptr<TH2F>(new TH2F(name, "", nBins_X1pos, min_X1pos, max_X1pos, nBins_U1thCal, min_U1thCal, max_U1thCal)));

        //------------------------------------------------
        nBins_X1pos = 4000;
        min_X1pos = 0.0;
        max_X1pos = 1000.0;
        
        sprintf(name, "hThetaRaytraceCalvsX1pos_TLCIteration_%d", j);
        hThetaRaytraceCalvsX1pos_TLCIteration.push_back(std::unique_ptr<TH2F>(new TH2F(name, "", nBins_X1pos, min_X1pos, max_X1pos, nBins_thetaRaytraceCal, min_thetaRaytraceCal, max_thetaRaytraceCal)));

        //------------------------------------------------
        hX1_TLCIteration[j]->SetStats(0);
        hThetaSCATvsX1pos_TLCIteration[j]->SetStats(0);
        hY1vsX1pos_TLCIteration[j]->SetStats(0);
        hTOFCalvsX1pos_TLCIteration[j]->SetStats(0);
        hX1thCalvsX1pos_TLCIteration[j]->SetStats(0);
        hU1thCalvsX1pos_TLCIteration[j]->SetStats(0);
        hThetaRaytraceCalvsX1pos_TLCIteration[j]->SetStats(0);
        
        //------------------------------------------------
        hX1_TLCIteration[j]->GetXaxis()->CenterTitle(true);
        hThetaSCATvsX1pos_TLCIteration[j]->GetXaxis()->CenterTitle(true);
        hY1vsX1pos_TLCIteration[j]->GetXaxis()->CenterTitle(true);
        hTOFCalvsX1pos_TLCIteration[j]->GetXaxis()->CenterTitle(true);
        hX1thCalvsX1pos_TLCIteration[j]->GetXaxis()->CenterTitle(true);
        hU1thCalvsX1pos_TLCIteration[j]->GetXaxis()->CenterTitle(true);
        hThetaRaytraceCalvsX1pos_TLCIteration[j]->GetXaxis()->CenterTitle(true);
        
        //------------------------------------------------
        hX1_TLCIteration[j]->GetXaxis()->SetTitleSize(0.05);
        hThetaSCATvsX1pos_TLCIteration[j]->GetXaxis()->SetTitleSize(0.05);
        hY1vsX1pos_TLCIteration[j]->GetXaxis()->SetTitleSize(0.05);
        hTOFCalvsX1pos_TLCIteration[j]->GetXaxis()->SetTitleSize(0.05);
        hX1thCalvsX1pos_TLCIteration[j]->GetXaxis()->SetTitleSize(0.05);
        hU1thCalvsX1pos_TLCIteration[j]->GetXaxis()->SetTitleSize(0.05);
        hThetaRaytraceCalvsX1pos_TLCIteration[j]->GetXaxis()->SetTitleSize(0.05);
        
        //------------------------------------------------
        hX1_TLCIteration[j]->GetYaxis()->CenterTitle(true);
        hThetaSCATvsX1pos_TLCIteration[j]->GetYaxis()->CenterTitle(true);
        hY1vsX1pos_TLCIteration[j]->GetYaxis()->CenterTitle(true);
        hTOFCalvsX1pos_TLCIteration[j]->GetYaxis()->CenterTitle(true);
        hX1thCalvsX1pos_TLCIteration[j]->GetYaxis()->CenterTitle(true);
        hU1thCalvsX1pos_TLCIteration[j]->GetYaxis()->CenterTitle(true);
        hThetaRaytraceCalvsX1pos_TLCIteration[j]->GetYaxis()->CenterTitle(true);
        
        //------------------------------------------------
        hX1_TLCIteration[j]->GetYaxis()->SetTitleSize(0.05);
        hThetaSCATvsX1pos_TLCIteration[j]->GetYaxis()->SetTitleSize(0.05);
        hY1vsX1pos_TLCIteration[j]->GetYaxis()->SetTitleSize(0.05);
        hTOFCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitleSize(0.05);
        hX1thCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitleSize(0.05);
        hU1thCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitleSize(0.05);
        hThetaRaytraceCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitleSize(0.05);
        
        if(j==0)
        {
            hX1_TLCIteration[j]->GetXaxis()->SetTitle("X1pos");
            hThetaSCATvsX1pos_TLCIteration[j]->GetXaxis()->SetTitle("X1pos");
            hY1vsX1pos_TLCIteration[j]->GetXaxis()->SetTitle("X1pos");
        }
        else
        {
            hX1_TLCIteration[j]->GetXaxis()->SetTitle("X1posC");
            hThetaSCATvsX1pos_TLCIteration[j]->GetXaxis()->SetTitle("X1posC");
            hY1vsX1pos_TLCIteration[j]->GetXaxis()->SetTitle("X1posC");
        }
        
        hX1_TLCIteration[j]->GetYaxis()->SetTitle("Counts");
        hThetaSCATvsX1pos_TLCIteration[j]->GetYaxis()->SetTitle("#theta_{SCAT}");
        hY1vsX1pos_TLCIteration[j]->GetYaxis()->SetTitle("Y1");
        hTOFCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitle("tofCal (Calibrated TOF)");
        hX1thCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitle("X1thCal (Calibrated X1 raytrace angle)");
        hU1thCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitle("U1thCal (Calibrated u1 raytrace angle)");
        hThetaRaytraceCalvsX1pos_TLCIteration[j]->GetYaxis()->SetTitle("thetaRaytraceCak (Calibrated SBR theta raytrace angle)");
    }
    
    for(int i=0; i<nrOfRuns; i++)
    {
        TTreeReader                 reader("DATA", DATAChains[i]->GetFile());
        TTreeReaderValue<double>    fX1pos(reader, "X1pos.t_X1pos");
        TTreeReaderValue<double>    fThetaSCAT(reader, "thetaSCAT.t_thetaSCAT");
        TTreeReaderValue<double>    fY1(reader, "Y1.t_Y1");
        
        TTreeReaderValue<double>    fTofCal(reader, "tofCal.t_tofCal");
        TTreeReaderValue<double>    fX1thCal(reader, "X1thCal.t_X1thCal");
        TTreeReaderValue<double>    fU1thCal(reader, "U1thCal.t_U1thCal");
        
        TTreeReaderValue<double>    fX1chisq(reader, "X1chisq.t_X1chisq");
        
        TTreeReaderArray<double>    fxPositions_TLC(reader, "xPositions_TLC");
        //TTreeReaderArray<double>    fthetaRaytraceCals(reader, "thetaRaytraceCals");

        
        double mfX1pos; // mfX1pos: mapped fX1pos using parameters from InitialPeakPositionMapping()
        double correctionParValue;
        
        int eN = 0;
        while(eventLists[i]->GetEntry(eN)>0)
        {
            reader.SetEntry(eventLists[i]->GetEntry(eN));
            //mfX1pos = (*fX1pos);
            mfX1pos = fxPositions_TLC[0];
            
            //------------------------------------------------
            //      Iterating the lineshape corrections
            
            //      Uncorrected spectra
            hX1_TLCIteration[0]->Fill(mfX1pos);
            hThetaSCATvsX1pos_TLCIteration[0]->Fill(mfX1pos, (*fThetaSCAT));
            hY1vsX1pos_TLCIteration[0]->Fill(mfX1pos, (*fY1));
            
            hTOFCalvsX1pos_TLCIteration[0]->Fill(mfX1pos, (*fTofCal));
            hX1thCalvsX1pos_TLCIteration[0]->Fill(mfX1pos, (*fX1thCal));
            hU1thCalvsX1pos_TLCIteration[0]->Fill(mfX1pos, (*fU1thCal));
            //hThetaRaytraceCalvsX1pos_TLCIteration[0]->Fill(mfX1pos, (fthetaRaytraceCals[0]));
            
            
            for(int j=0; j<(int) finalLineshapeCorrectionParameters.size(); j++)
            {
                //mfX1pos = singleLineshapeCorrection(mfX1pos, (*fThetaSCAT), (*fY1), finalLineshapeCorrectionTypes[j], finalLineshapeCorrectionParameters[j]);
                
                //------------------------------------------------
                double correctionParValue = 0.0;
                
                //--------------------------------------------
                if(finalLineshapeCorrectionTypes[j]==0)
                {
                    correctionParValue = (*fThetaSCAT);
                }
                else if(finalLineshapeCorrectionTypes[j]==1)
                {
                    correctionParValue = (*fY1);
                }
                else if(finalLineshapeCorrectionTypes[j]==2)
                {
                    correctionParValue = (*fTofCal);
                }
                else if(finalLineshapeCorrectionTypes[j]==3)
                {
                    correctionParValue = (*fX1thCal);
                }
                else if(finalLineshapeCorrectionTypes[j]==4)
                {
                    correctionParValue = (*fU1thCal);
                }
                else if(finalLineshapeCorrectionTypes[j]==5)
                {
                    //correctionParValue = (fthetaRaytraceCals[0]);
                }

                
                //------------------------------------------------
                mfX1pos = singleLineshapeCorrection(mfX1pos, correctionParValue, finalLineshapeCorrectionParameters[j]);
                
                hX1_TLCIteration[j+1]->Fill(mfX1pos);
                hThetaSCATvsX1pos_TLCIteration[j+1]->Fill(mfX1pos, (*fThetaSCAT));
                hY1vsX1pos_TLCIteration[j+1]->Fill(mfX1pos, (*fY1));
                
                hTOFCalvsX1pos_TLCIteration[j+1]->Fill(mfX1pos, (*fTofCal));
                hX1thCalvsX1pos_TLCIteration[j+1]->Fill(mfX1pos, (*fX1thCal));
                hU1thCalvsX1pos_TLCIteration[j+1]->Fill(mfX1pos, (*fU1thCal));
                //hThetaRaytraceCalvsX1pos_TLCIteration[j+1]->Fill(mfX1pos, (fthetaRaytraceCals[0]));
            }
            
            eN++;
        }
        
        reader.Restart();
    }
    
    std::unique_ptr<TF1> fitPeakRes(new TF1("fitPeakRes", "gaus(0) + pol1(3)", 0., 1000.));
    fitPeakRes->SetParameters(1.0e+01, 7.20e+02, 1.0e+00, 0.0e+00, 0.0e+00);
    fitPeakRes->SetLineColor(kRed);
    fitPeakRes->SetNpx(1000);
    fitPeakRes->SetParLimits(0, 0.0, 1.0e+08);
    
    
    f.Write();
    
    //--------------------------------------------------------------------------------
    bool producedSummarySpectra = false;
    
    if(producedSummarySpectra)
    {
        int nRows = (int) finalLineshapeCorrectionParameters.size() + 1;
        int nColumns = 3 + (int) lineshapeCorrectionPeaks.size();
        
        std::unique_ptr<TCanvas> c2(new TCanvas("c2","",10,10,500*nColumns,350*nRows));
        c2->Divide(nColumns, nRows, 0.0001, 0.0001);
        
        for(int h=0; h<(int) lineshapeCorrectionPeaks.size(); h++)
        {
            for(int i=0; i<(int) finalLineshapeCorrectionParameters.size() + 1; i++)
            {
                hX1_TLCIteration[i]->GetXaxis()->SetRangeUser(std::get<0>(lineshapeCorrectionPeaks[h]), std::get<1>(lineshapeCorrectionPeaks[h]));
                hThetaSCATvsX1pos_TLCIteration[i]->GetXaxis()->SetRangeUser(std::get<0>(lineshapeCorrectionPeaks[h]), std::get<1>(lineshapeCorrectionPeaks[h]));
                hY1vsX1pos_TLCIteration[i]->GetXaxis()->SetRangeUser(std::get<0>(lineshapeCorrectionPeaks[h]), std::get<1>(lineshapeCorrectionPeaks[h]));
                
                c2->cd(i*nColumns + 1);
                hX1_TLCIteration[i]->Draw("");
                
                c2->cd(i*nColumns + 2);
                hThetaSCATvsX1pos_TLCIteration[i]->Draw("COL");
                
                c2->cd(i*nColumns + 3);
                hY1vsX1pos_TLCIteration[i]->Draw("COL");
                
                if(i<(int) finalLineshapeCorrectionParameters.size())
                {
                    for(int j=0; j<(int) grTLCCoefficientsStorage[i].size(); j++)
                    {
                        grTLCCoefficientsStorage[i][j]->SetMarkerStyle(1);
                        grTLCCoefficientsStorage[i][j]->SetMarkerSize(3);
                        grTLCCoefficientsStorage[i][j]->SetMarkerColor(kBlack);
                        
                        grTLCCoefficientsStorage[i][j]->GetXaxis()->CenterTitle(true);
                        
                        grTLCCoefficientsStorage[i][j]->GetXaxis()->SetTitleSize(0.05);
                        
                        grTLCCoefficientsStorage[i][j]->GetYaxis()->CenterTitle(true);
                        grTLCCoefficientsStorage[i][j]->GetYaxis()->SetTitleSize(0.05);
                        
                        if(finalLineshapeCorrectionTypes[i]==0)
                        {
                            grTLCCoefficientsStorage[i][j]->GetXaxis()->SetTitle("#theta_{SCAT}");
                        }
                        else if(finalLineshapeCorrectionTypes[i]==1)
                        {
                            grTLCCoefficientsStorage[i][j]->GetXaxis()->SetTitle("Y1");
                        }
                        
                        sprintf(name, "Polynomial coefficient (order %d)", j);
                        grTLCCoefficientsStorage[i][j]->GetYaxis()->SetTitle(name);
                        
                        c2->cd(i*nColumns + 4 + j);
                        grTLCCoefficientsStorage[i][j]->Draw("AP");
                        
                        //------------------------------------------------------------------------
                        grTLCExtendedCoefficientsStorage[i][j]->SetMarkerStyle(kFullCircle);
                        grTLCExtendedCoefficientsStorage[i][j]->SetMarkerSize(1);
                        grTLCExtendedCoefficientsStorage[i][j]->SetMarkerColor(kOrange-2);
                        grTLCExtendedCoefficientsStorage[i][j]->SetLineColor(kBlack);
                        
                        grTLCExtendedCoefficientsStorage[i][j]->Draw("P");
                    }
                }
            }
            
            sprintf(name, "Output/TLC_summary_peak%d.root", h);
            c2->Print(name);
            
            sprintf(name, "Output/TLC_summary_peak%d.png", h);
            c2->Print(name);
        }
    }
    
    //--------------------------------------------------------------------------------
    std::ofstream fileTLCParameters("Output/PR166_2011_TLCParameters.dat", std::ofstream::out);
    
    //--------------------------------
    //      Correction Types
    fileTLCParameters << (int) finalLineshapeCorrectionTypes.size();
    for(int i=0; i<(int) finalLineshapeCorrectionTypes.size(); i++)
    {
        fileTLCParameters << " " << finalLineshapeCorrectionTypes[i];
    }
    fileTLCParameters << "\n";
    
    for(int i=0; i<(int) finalLineshapeCorrectionParameters.size(); i++)
    {
        fileTLCParameters << (int) finalLineshapeCorrectionParameters[i].size();
        
        for(int j=0; j<(int) finalLineshapeCorrectionParameters[i].size(); j++)
        {
            fileTLCParameters << " " << (int) finalLineshapeCorrectionParameters[i][j].size();
            
            for(int k=0; k<(int) finalLineshapeCorrectionParameters[i][j].size(); k++)
            {
                fileTLCParameters << " " << finalLineshapeCorrectionParameters[i][j][k];
            }
        }
        
        fileTLCParameters << "\n";
    }
    
    fileTLCParameters.close();
}

//--------------------------------------------
//double PreviousLineshapeCorrections(double fX1pos, double fThetaSCAT, double fY1)
double PreviousLineshapeCorrections(double fX1pos, std::vector<double> correctionVariables)
{
    double correctedPosition = 0.0;
    
    for(int i=0; i<(int) finalLineshapeCorrectionParameters.size(); i++)
    {
        double previousCorrectedPosition;
        std::vector<double> correctionPars_PolCoefficients;
        
        if(i==0)
        {
            previousCorrectedPosition = fX1pos;
        }
        else
        {
            previousCorrectedPosition = correctedPosition;
        }
        
        //--------------------------------------------------------------------
        for(int j=0; j<(int) finalLineshapeCorrectionParameters[i].size(); j++)
        {
            double correctionPar_PolCoefficient = 0.0;
            
            //--------------------------------------------------------------------
            for(int k=0; k<(int) finalLineshapeCorrectionParameters[i][j].size(); k++)
            {
                correctionPar_PolCoefficient += finalLineshapeCorrectionParameters[i][j][k]*pow(correctionVariables[i], k);

                /*
                if(finalLineshapeCorrectionTypes[i]==0)
                {
                    correctionPar_PolCoefficient += finalLineshapeCorrectionParameters[i][j][k]*pow(fThetaSCAT, k);
                }
                
                if(finalLineshapeCorrectionTypes[i]==1)
                {
                    correctionPar_PolCoefficient += finalLineshapeCorrectionParameters[i][j][k]*pow(fY1, k);
                }
                */
            }
            
            correctionPars_PolCoefficients.push_back(correctionPar_PolCoefficient);
        }
        
        correctedPosition = 0.0;
        
        for(int j=0; j<(int) correctionPars_PolCoefficients.size(); j++)
        {
            if((int) correctionPars_PolCoefficients.size()==1)
            {
                correctedPosition = previousCorrectedPosition + correctionPars_PolCoefficients[j];
            }
            else
            {
                correctedPosition += correctionPars_PolCoefficients[j]*pow(previousCorrectedPosition, j);
            }
        }
    }
    
    if(finalLineshapeCorrectionParameters.empty())
    {
        correctedPosition = fX1pos;
    }
    
    return correctedPosition;
}

//--------------------------------------------
void TLC(int correctionType, int nDivisions_TLC, double min_TLC, double max_TLC) // correctionType: 0 => thetaSCAT correction, 1 => Y1 correction
{
    //--------------------------------------------------------
    //      First check the amount of extension points
    int nExtensionPoints = 0;
    int nSkippedPoints = 0;
    double errorWeighting = 1.0;
    
    for(int i=0; i<(int) nExtensionPointsPerCorrectionType.size(); i++)
    {
        if(std::get<0>(nExtensionPointsPerCorrectionType[i]) == correctionType)
        {
            //------------------------------------------------
            nExtensionPoints = std::get<1>(nExtensionPointsPerCorrectionType[i]);
            nSkippedPoints = std::get<2>(nExtensionPointsPerCorrectionType[i]);
            errorWeighting = std::get<3>(nExtensionPointsPerCorrectionType[i]);
            
            //------------------------------------------------
            double divisionSize = (double) abs(max_TLC-min_TLC)/nDivisions_TLC;
            nDivisions_TLC += 2*nExtensionPoints;
            
            min_TLC -= nExtensionPoints*divisionSize;
            max_TLC += nExtensionPoints*divisionSize;
        }
    }
    
    //------------------------------------------------
    static int iterationTLC = 0;
    if(iterationTLC==0)
    {
        std::cout << "//==========================================================//" << std::endl;
    }
    
    std::cout << "      Commencing TLC iteration: " << iterationTLC;
    
    if(correctionType==0)
    {
        std::cout << ", thetaSCAT correction" << std::endl;
    }
    else if(correctionType==1)
    {
        std::cout << ", Y1 correction" << std::endl;
    }
    else if(correctionType==2)
    {
        std::cout << ", tofCal correction" << std::endl;
    }
    else if(correctionType==3)
    {
        std::cout << ", X1thCal correction" << std::endl;
    }
    else if(correctionType==1)
    {
        std::cout << ", U1thCal correction" << std::endl;
    }

    iterationTLC++;
    
    double division_TLC = (max_TLC-min_TLC)/nDivisions_TLC;
    
    //----------------------------------------
    bool forceExtrapolationPoints_lower = false;
    bool forceExtrapolationPoints_higher = false;
    
    int nForcedExtrapolationPoints_lower = 0;
    int nForcedExtrapolationPoints_higher = 0;
    
    for(int i=0; i<(int) extrapolationPoints.size(); i++)
    {
        if(std::get<0>(extrapolationPoints[i])==correctionType)
        {
            if(std::get<1>(extrapolationPoints[i])>0)
            {
                forceExtrapolationPoints_lower = true;
            }
            if(std::get<2>(extrapolationPoints[i])>0)
            {
                forceExtrapolationPoints_higher = true;
            }
            
            nForcedExtrapolationPoints_lower = std::get<1>(extrapolationPoints[i]);
            nForcedExtrapolationPoints_higher = std::get<2>(extrapolationPoints[i]);
        }
    }
    
    if(nForcedExtrapolationPoints_lower>nDivisions_TLC)
    {
        nForcedExtrapolationPoints_lower = nDivisions_TLC;
    }
    if(nForcedExtrapolationPoints_higher>nDivisions_TLC)
    {
        nForcedExtrapolationPoints_higher = nDivisions_TLC;
    }
    
    //----------------------------------------
    static int correctionNumber = 0;
    sprintf(name, "Output/TLC_fitting_iteration%d.root", correctionNumber);
    TFile f(name, "recreate");
    correctionNumber++;
    
    nBins_X1pos = 2000;
    
    min_X1pos = 0.0;
    max_X1pos = 1000.0;
    
    std::unique_ptr<TH1F> hX1pos_centralCut(new TH1F("hX1pos_centralCut", "", nBins_X1pos, min_X1pos, max_X1pos));
    std::vector<std::unique_ptr<TH1F>> hX1pos_TLCcuts;
    
    std::vector<std::unique_ptr<TF1>> fitTLC_peakPositions;
    
    std::vector<std::vector<std::pair<double, double>>> lineshapeCorrectionPositions;
    std::vector<bool> validTLCDataPoints_perCut;
    std::vector<std::vector<bool>> validTLCDataPoints_perPeak;
    std::vector<double> referenceX1posPositions_centralCut;
    std::vector<std::shared_ptr<TGraphErrors>> grLineshapeCorrections;
    
    //  These vectors contain the correction parameters for each thetaSCAT/Y1 cut
    std::vector<std::vector<std::pair<double, double>>> lineshapeCorrectionPars;
    
    //  These vectors contain the polynomial coefficients (0th order, 1st order...) for the thetaSCAT/Y1 correction factors
    std::vector<std::vector<double>> lineshapeCorrectionPolPars;
    
    //--------------------------------------------
    std::unique_ptr<TF1> fLineshapeCorrection;
    
    if(nLineshapeCorrectionPeaks==1)
    {
        fLineshapeCorrection = std::unique_ptr<TF1>(new TF1("fLineshapeCorrection", "x[0] + pol0(0)", 0., 1000.));
        fLineshapeCorrection->SetParameter(0, 0.0);
    }
    else
    {
        sprintf(condition, "pol%d(0)", nLineshapeCorrectionPeaks-1);
        fLineshapeCorrection = std::unique_ptr<TF1>(new TF1("fLineshapeCorrection", condition, 0., 1000.));
        
        for(int i=0; i<nLineshapeCorrectionPeaks; i++)
        {
            fLineshapeCorrection->SetParameter(i, 0.0);
        }
    }
    
    for(int i=0; i<nDivisions_TLC; i++)
    {
        sprintf(name, "hX1pos_TLCcuts%d", i);
        hX1pos_TLCcuts.push_back(std::unique_ptr<TH1F>(new TH1F(name, "", nBins_X1pos, min_X1pos, max_X1pos)));
    }
    
    for(int i=0; i<nLineshapeCorrectionPeaks; i++)
    {
        sprintf(name, "fitTLC_peakPositions%d", i);
        if(std::get<2>(lineshapeCorrectionPeaks[i])>=0)
        {
            sprintf(condition, "gaus(0) + pol%d(3)", std::get<2>(lineshapeCorrectionPeaks[i]));
            fitTLC_peakPositions.push_back(std::unique_ptr<TF1>(new TF1(name, "gaus(0) + pol1(3)", 0., 1000.)));
        }
        else
        {
            fitTLC_peakPositions.push_back(std::unique_ptr<TF1>(new TF1(name, "gaus(0)", 0., 1000.)));
        }
        
        fitTLC_peakPositions[i]->SetNpx(2000);
        fitTLC_peakPositions[i]->SetParameters(1.0e+01, 100.0, 1.0, 0.0, 0.0);
    }
    
    //------------------------------------------------
    nBins_X1pos = 2000;
    min_X1pos = 0.0;
    max_X1pos = 1000.0;
    nBins_thetaSCAT = 400;
    
    std::unique_ptr<TH2F> hThetaSCATvsX1pos(new TH2F("hThetaSCATvsX1pos", "", nBins_X1pos, min_X1pos, max_X1pos, nBins_thetaSCAT, min_thetaSCAT, max_thetaSCAT));
    
    //------------------------------------------------
    nBins_X1pos = 2000;
    min_X1pos = 0.0;
    max_X1pos = 1000.0;
    nBins_Y1 = 400;
    
    std::unique_ptr<TH2F> hY1vsX1pos(new TH2F("hY1vsX1pos", "", nBins_X1pos, min_X1pos, max_X1pos, nBins_Y1, min_Y1, max_Y1));
    
    //------------------------------------------------
    nBins_X1pos = 2000;
    min_X1pos = 0.0;
    max_X1pos = 1000.0;
    
    std::unique_ptr<TH2F> hTofCalvsX1pos(new TH2F("hTofCalvsX1pos", "", nBins_X1pos, min_X1pos, max_X1pos, nBins_tofCal, min_tofCal, max_tofCal));

    //------------------------------------------------
    nBins_X1pos = 2000;
    min_X1pos = 0.0;
    max_X1pos = 1000.0;
    
    std::unique_ptr<TH2F> hX1thCalvsX1pos(new TH2F("hX1thCalvsX1pos", "", nBins_X1pos, min_X1pos, max_X1pos, nBins_X1thCal, min_X1thCal, max_X1thCal));
    
    //------------------------------------------------
    nBins_X1pos = 2000;
    min_X1pos = 0.0;
    max_X1pos = 1000.0;
    
    std::unique_ptr<TH2F> hU1thCalvsX1pos(new TH2F("hU1thCalvsX1pos", "", nBins_X1pos, min_X1pos, max_X1pos, nBins_U1thCal, min_U1thCal, max_U1thCal));
    
    //------------------------------------------------
    std::cout << std::endl;
    
    for(int i=0; i<nrOfRuns; i++)
    {
        TTreeReader                 reader("DATA", DATAChains[i]->GetFile());
        TTreeReaderValue<double>    fX1pos(reader, "X1pos.t_X1pos");
        TTreeReaderValue<double>    fThetaSCAT(reader, "thetaSCAT.t_thetaSCAT");
        TTreeReaderValue<double>    fY1(reader, "Y1.t_Y1");
        
        TTreeReaderValue<double>    fTofCal(reader, "tofCal.t_tofCal");
        TTreeReaderValue<double>    fX1thCal(reader, "X1thCal.t_X1thCal");
        TTreeReaderValue<double>    fU1thCal(reader, "U1thCal.t_U1thCal");
        
        TTreeReaderValue<double>    fX1chisq(reader, "X1chisq.t_X1chisq");
        
        TTreeReaderArray<double>    fxPositions_TLC(reader, "xPositions_TLC");
        //TTreeReaderArray<double>    fthetaRaytraceCals(reader, "thetaRaytraceCals");

        
        double mfX1pos; // mfX1pos: mapped fX1pos using parameters from InitialPeakPositionMapping()
        double correctionParValue;
        
        int eN = 0;
        while(eventLists[i]->GetEntry(eN)>0)
        {
            reader.SetEntry(eventLists[i]->GetEntry(eN));
            
            //----------------------------------------------------
            std::vector<double> correctionVariables;
            
            for(auto j=finalLineshapeCorrectionTypes.begin(); j!=finalLineshapeCorrectionTypes.end(); j++)
            {
                if((*j)==0)
                {
                    correctionVariables.push_back(*fThetaSCAT);
                }
                else if((*j)==1)
                {
                    correctionVariables.push_back(*fY1);
                }
                else if((*j)==2)
                {
                    correctionVariables.push_back(*fTofCal);
                }
                else if((*j)==3)
                {
                    correctionVariables.push_back(*fX1thCal);
                }
                else if((*j)==4)
                {
                    correctionVariables.push_back(*fU1thCal);
                }
                else if((*j)==5)
                {
                    //correctionVariables.push_back(fthetaRaytraceCals[0]);
                }
            }
            //--------------------------------------------
            //      Pre-lineshape Correction Mapping
            //mfX1pos = (*fX1pos);
            mfX1pos = fxPositions_TLC[0];
            mfX1pos = PreviousLineshapeCorrections(mfX1pos, correctionVariables);
            
            //--------------------------------------------
            if(correctionType==0)
            {
                correctionParValue = (*fThetaSCAT);
            }
            else if(correctionType==1)
            {
                correctionParValue = (*fY1);
            }
            
            else if(correctionType==2)
            {
                correctionParValue = (*fTofCal);
            }
            else if(correctionType==3)
            {
                correctionParValue = (*fX1thCal);
            }
            else if(correctionType==4)
            {
                correctionParValue = (*fU1thCal);
            }
            else if(correctionType==5)
            {
                //correctionParValue = (fthetaRaytraceCals[0]);
            }

            
            //--------------------------------------------
            for(int j=0; j<nDivisions_TLC; j++)
            {
                if((correctionParValue > min_TLC + (j*division_TLC)) && (correctionParValue < min_TLC + ((j+1)*division_TLC)))
                {
                    hX1pos_TLCcuts[j]->Fill(mfX1pos);
                }
            }
            
            if(correctionParValue>-(0.5*division_TLC) && correctionParValue<(0.5*division_TLC))
            {
                hX1pos_centralCut->Fill(mfX1pos);
            }
            
            //--------------------------------------------
            hThetaSCATvsX1pos->Fill(mfX1pos, (*fThetaSCAT));
            hY1vsX1pos->Fill(mfX1pos, (*fY1));
            
            hTofCalvsX1pos->Fill(mfX1pos, (*fTofCal));
            hX1thCalvsX1pos->Fill(mfX1pos, (*fX1thCal));
            hU1thCalvsX1pos->Fill(mfX1pos, (*fU1thCal));
            
            
            eN++;
        }
        
        reader.Restart();
    }
    
    //----------------------------------------------------------------------------------------
    //  Central cut for thetaSCAT and Y1 to provide which reference positions to map to.
    double initialMean;
    double initialPeakHeight;
    double initialStdDev;
    double minimumRequiredHeight = 10.0;
    bool setMinimumRequiredHeight = false;
    
    for(int j=0; j<(int) minimumPeakHeights.size(); j++)
    {
        if(std::get<0>(minimumPeakHeights[j])==correctionType)
        {
            setMinimumRequiredHeight = true;
            minimumRequiredHeight = std::get<1>(minimumPeakHeights[j]);
        }
    }
    
    for(int j=0; j<nLineshapeCorrectionPeaks; j++)
    {
        hX1pos_centralCut->GetXaxis()->SetRangeUser(std::get<0>(lineshapeCorrectionPeaks[j]), std::get<1>(lineshapeCorrectionPeaks[j]));
        initialMean = hX1pos_centralCut->GetBinCenter(hX1pos_centralCut->GetMaximumBin());
        initialPeakHeight = hX1pos_centralCut->GetBinContent(hX1pos_centralCut->GetMaximumBin());
        initialStdDev = hX1pos_centralCut->GetStdDev();
        
        fitTLC_peakPositions[j]->SetParameter(0, initialPeakHeight);
        fitTLC_peakPositions[j]->SetParLimits(0, 0.0e+00, 3.0*initialPeakHeight);
        
        fitTLC_peakPositions[j]->SetParameter(1, initialMean);
        fitTLC_peakPositions[j]->SetParLimits(1, initialMean - 3.0*initialStdDev, initialMean + 3.0*initialStdDev);
        
        sprintf(name, "fitTLC_peakPositions%d", j);
        hX1pos_centralCut->Fit(name, "Q+", "", std::get<0>(lineshapeCorrectionPeaks[j]), std::get<1>(lineshapeCorrectionPeaks[j]));
        
        referenceX1posPositions_centralCut.push_back(fitTLC_peakPositions[j]->GetParameter(1));
    }
    
    //----------------------------------------------------------------------------------------
    //      Rebinning for some histograms which do not meet the minimum peak height test
    for(int i=0; i<nDivisions_TLC; i++)
    {
        int nAttempts = 0;
        
        for(int j=0; j<nLineshapeCorrectionPeaks; j++)
        {
            hX1pos_TLCcuts[i]->GetXaxis()->SetRangeUser(std::get<0>(lineshapeCorrectionPeaks[j]), std::get<1>(lineshapeCorrectionPeaks[j]));
            initialMean = hX1pos_TLCcuts[i]->GetBinCenter(hX1pos_TLCcuts[i]->GetMaximumBin());
            initialPeakHeight = hX1pos_TLCcuts[i]->GetBinContent(hX1pos_TLCcuts[i]->GetMaximumBin());
            initialStdDev = hX1pos_TLCcuts[i]->GetStdDev();
            
            while(setMinimumRequiredHeight && initialPeakHeight<minimumRequiredHeight && nAttempts<nRebinAttempts)
            {
                hX1pos_TLCcuts[i]->Rebin(2);
                
                initialMean = hX1pos_TLCcuts[i]->GetBinCenter(hX1pos_TLCcuts[i]->GetMaximumBin());
                initialPeakHeight = hX1pos_TLCcuts[i]->GetBinContent(hX1pos_TLCcuts[i]->GetMaximumBin());
                initialStdDev = hX1pos_TLCcuts[i]->GetStdDev();
                
                nAttempts++;
            }
        }
    }
    
    //----------------------------------------------------------------------------------------
    for(int i=0; i<nDivisions_TLC; i++)
    {
        bool validDataPoint = true;
        std::vector<std::pair<double, double>> tempLineshapeCorrectionPositions;
        
        std::vector<bool> validTLCDataPoints;
        
        for(int j=0; j<nLineshapeCorrectionPeaks; j++)
        {
            hX1pos_TLCcuts[i]->GetXaxis()->SetRangeUser(std::get<0>(lineshapeCorrectionPeaks[j]), std::get<1>(lineshapeCorrectionPeaks[j]));
            initialMean = hX1pos_TLCcuts[i]->GetBinCenter(hX1pos_TLCcuts[i]->GetMaximumBin());
            //initialMean = hX1pos_TLCcuts[i]->GetMean();
            initialPeakHeight = hX1pos_TLCcuts[i]->GetBinContent(hX1pos_TLCcuts[i]->GetMaximumBin());
            initialStdDev = hX1pos_TLCcuts[i]->GetStdDev();
            
            //------------------------------------------------
            //      Amplitude
            fitTLC_peakPositions[j]->SetParameter(0, initialPeakHeight);
            fitTLC_peakPositions[j]->SetParLimits(0, 0.0e+00, 3.0*initialPeakHeight);
            
            //------------------------------------------------
            //      Mean
            double lowerLimit, upperLimit;
            
            if((initialMean - 3.0*initialStdDev)>std::get<0>(lineshapeCorrectionPeaks[j]))
            {
                lowerLimit = initialMean - 3.0*initialStdDev;
            }
            else
            {
                lowerLimit = std::get<0>(lineshapeCorrectionPeaks[j]);
            }
            
            if((initialMean + 3.0*initialStdDev)<std::get<1>(lineshapeCorrectionPeaks[j]))
            {
                upperLimit = initialMean + 3.0*initialStdDev;
            }
            else
            {
                upperLimit = std::get<1>(lineshapeCorrectionPeaks[j]);
            }
            
            fitTLC_peakPositions[j]->SetParameter(1, initialMean);
            //  TEMP OLD
            //fitTLC_peakPositions[j]->SetParLimits(1, initialMean - 3.0*initialStdDev, initialMean + 3.0*initialStdDev);
            fitTLC_peakPositions[j]->SetParLimits(1, lowerLimit, upperLimit);
            
            //------------------------------------------------
            //      Standard Deviation
            fitTLC_peakPositions[j]->SetParameter(2, initialStdDev);
            
            if((std::get<2>(expRes) && initialStdDev>std::get<0>(expRes) && initialStdDev<std::get<1>(expRes)) || !std::get<2>(expRes))
            {
                fitTLC_peakPositions[j]->SetParLimits(2, 0.2*initialStdDev, 2.0*initialStdDev);
            }
            else if(std::get<2>(expRes))
            {
                fitTLC_peakPositions[j]->SetParLimits(2, std::get<0>(expRes), std::get<1>(expRes));
            }
            
            fitTLC_peakPositions[j]->SetParameter(3, 0.0);
            fitTLC_peakPositions[j]->SetParameter(4, 0.0);
            
            
            //------------------------------------------------
            bool useDataPoint = false;
            
            //------------------------
            //      Inclusions
            if(setMinimumRequiredHeight && initialPeakHeight>minimumRequiredHeight)
            {
                useDataPoint = true;
            }
            else if(!setMinimumRequiredHeight)
            {
                useDataPoint = true;
            }
            
            //------------------------
            //      Exclusions
            if(i<nExtensionPoints || i>((nDivisions_TLC-1)-nExtensionPoints))
            {
                useDataPoint = false;
            }
            
            //------------------------
            if(useDataPoint)
            {
                sprintf(name, "fitTLC_peakPositions%d", j);
                hX1pos_TLCcuts[i]->Fit(name, "Q+", "", std::get<0>(lineshapeCorrectionPeaks[j]), std::get<1>(lineshapeCorrectionPeaks[j]));
            }
            else
            {
                validDataPoint = false;
            }
            
            //------------------------------------------------------------------------
            double chiSquare = fitTLC_peakPositions[j]->GetChisquare();
            int nDF = fitTLC_peakPositions[j]->GetNDF();
            double reducedChiSquare = (chiSquare/nDF);
            
            double defaultReducedChiSquareLimit = 5.0;
            if((setMaximumRedChiSquare && reducedChiSquare>maximumRedChiSquare) || (!setMaximumRedChiSquare && reducedChiSquare>defaultReducedChiSquareLimit))
            {
                validDataPoint = false;
            }
            
            double amplitude = fitTLC_peakPositions[j]->GetParameter(0);
            double position = fitTLC_peakPositions[j]->GetParameter(1);
            double positionError = fitTLC_peakPositions[j]->GetParError(1);
            double stdDev = fitTLC_peakPositions[j]->GetParameter(2);
            
            //----------------------------------------------------------------------
            //      Diagnostics to determine whether a fitted position is valid
            if(setMinimumRequiredHeight && amplitude<amplitude)
            {
                validDataPoint = false;
            }
            if(position<std::get<0>(lineshapeCorrectionPeaks[j]) || position>std::get<1>(lineshapeCorrectionPeaks[j]))
            {
                validDataPoint = false;
            }
            if(std::get<2>(expRes) && (stdDev<std::get<0>(expRes) || stdDev>std::get<1>(expRes)))
            {
                validDataPoint = false;
            }
            
            tempLineshapeCorrectionPositions.push_back(std::make_pair(position, positionError));
            validTLCDataPoints.push_back(validDataPoint);
        }
        
        validTLCDataPoints_perPeak.push_back(validTLCDataPoints);
        lineshapeCorrectionPositions.push_back(tempLineshapeCorrectionPositions);
        validTLCDataPoints_perCut.push_back(validDataPoint);
    }
    
    //------------------------------------------------
    //      Resetting the Histogram X-axis range
    hX1pos_centralCut->GetXaxis()->SetRangeUser(min_X1pos, max_X1pos);
    
    for(int i=0; i<nDivisions_TLC; i++)
    {
        hX1pos_TLCcuts[i]->GetXaxis()->SetRangeUser(min_X1pos, max_X1pos);
    }
    
    //--------------------------------------------------------------
    //      Setting the maximum amount of extrapolation points
    std::vector<bool> originalValidityTLCDataPoints_perCut;
    
    int nInitialValidDataPoints = 0;
    for(int i=0; i<nDivisions_TLC; i++)
    {
        if(validTLCDataPoints_perCut[i])
        {
            originalValidityTLCDataPoints_perCut.push_back(true);
            nInitialValidDataPoints++;
        }
        else
        {
            originalValidityTLCDataPoints_perCut.push_back(false);
        }
    }
    
    if(nForcedExtrapolationPoints_lower>nInitialValidDataPoints)
    {
        nForcedExtrapolationPoints_lower = nInitialValidDataPoints;
    }
    if(nForcedExtrapolationPoints_higher>nInitialValidDataPoints)
    {
        nForcedExtrapolationPoints_higher = nInitialValidDataPoints;
    }
    
    //----------------------------------------------------
    //      Extension of the positions to the limits
    std::vector<std::vector<double>> bestExtrapolationParameters;
    std::vector<std::vector<double>> bestExtrapolationParameters_lower;
    std::vector<std::vector<double>> bestExtrapolationParameters_higher;
    std::vector<double> averageErrors;
    std::vector<double> averageErrors_lower;
    std::vector<double> averageErrors_higher;
    
    for(int j=0; j<nLineshapeCorrectionPeaks; j++)
    {
        //--------------------------------
        //      Lower Extrapolation
        std::unique_ptr<TGraphErrors> grPositionExtrapolation_lower(new TGraphErrors());
        sprintf(name, "grPositionExtrapolation_lower_peak%d", j);
        grPositionExtrapolation_lower->SetName(name);
        grPositionExtrapolation_lower->SetMarkerStyle(1);
        grPositionExtrapolation_lower->SetMarkerSize(3);
        grPositionExtrapolation_lower->SetMarkerColor(kBlack);
        
        int nBestFitPoints_lower = 0;
        double averageError_lower = 0.0;
        
        double bestReducedChiSquare_lower = 0.0;
        std::vector<double> bestParameters_lower;
        
        int nExtrapolationPoints_lower = 0;
        for(int i=0; i<nDivisions_TLC; i++)
        {
            double cutCenterValue = min_TLC + ((i+0.5)*division_TLC);
            
            if(validTLCDataPoints_perCut[i])
            {
                grPositionExtrapolation_lower->SetPoint(nExtrapolationPoints_lower, cutCenterValue, lineshapeCorrectionPositions[i][j].first);
                grPositionExtrapolation_lower->SetPointError(nExtrapolationPoints_lower, division_TLC/2.0, lineshapeCorrectionPositions[i][j].second);
                
                //--------------------------------------------
                std::unique_ptr<TF1> linearFit_lower(new TF1("linearFit_lower", "pol1(0)", -1000., 1000.));
                linearFit_lower->SetParameters(0, 0.0e+00);
                linearFit_lower->SetParameters(1, 0.0e+00);
                
                if(nExtrapolationPoints_lower==0)
                {
                    linearFit_lower->FixParameter(1, 0.0e+00);
                }
                linearFit_lower->SetParameters(0, 0.0e+00);
                linearFit_lower->SetParameters(1, 0.0e+00);
                
                if(nExtrapolationPoints_lower>=1)
                {
                    double x1, x2, y1, y2;
                    grPositionExtrapolation_lower->GetPoint(0, x1, y1);
                    grPositionExtrapolation_lower->GetPoint(grPositionExtrapolation_lower->GetN()-1, x2, y2);
                    
                    double a, b;
                    
                    a = (y2-y1)/(x2-x1);
                    b = (-a*x1) + y1;
                    
                    linearFit_lower->SetParameter(0, b);
                    linearFit_lower->SetParameter(1, a);
                }
                
                grPositionExtrapolation_lower->Fit("linearFit_lower", "Q+", "", min_TLC, max_TLC);
                
                if(nExtrapolationPoints_lower==0)
                {
                    linearFit_lower->ReleaseParameter(1);
                }
                
                //--------------------------------------------
                double chiSquare = linearFit_lower->GetChisquare();
                int nDF = linearFit_lower->GetNDF();
                double reducedChiSquare = (chiSquare/nDF);
                
                nExtrapolationPoints_lower++;
                
                //--------------------------------------------
                //  If a better Reduced Chi Square is found
                if((forceExtrapolationPoints_lower && nExtrapolationPoints_lower==nForcedExtrapolationPoints_lower) || (!forceExtrapolationPoints_lower && (nExtrapolationPoints_lower==3 || (nExtrapolationPoints_lower>3 && reducedChiSquare<bestReducedChiSquare_lower))))
                {
                    bestReducedChiSquare_lower = reducedChiSquare;
                    bestParameters_lower.clear();
                    
                    bestParameters_lower.push_back(linearFit_lower->GetParameter(0));
                    bestParameters_lower.push_back(linearFit_lower->GetParameter(1));
                    
                    nBestFitPoints_lower = nExtrapolationPoints_lower;
                    
                    averageError_lower += lineshapeCorrectionPositions[i][j].second;
                }
            }
        }
        
        grPositionExtrapolation_lower->Write();
        averageError_lower *= (1.0/nBestFitPoints_lower);
        averageErrors_lower.push_back(averageError_lower);
        bestExtrapolationParameters_lower.push_back(bestParameters_lower);
        
        //--------------------------------
        //      Higher Extrapolation
        std::unique_ptr<TGraphErrors> grPositionExtrapolation_higher(new TGraphErrors());
        sprintf(name, "grPositionExtrapolation_higher_peak%d", j);
        grPositionExtrapolation_higher->SetName(name);
        grPositionExtrapolation_higher->SetMarkerStyle(1);
        grPositionExtrapolation_higher->SetMarkerSize(3);
        grPositionExtrapolation_higher->SetMarkerColor(kBlack);
        
        int nBestFitPoints_higher = 0;
        double averageError_higher = 0.0;
        
        double bestReducedChiSquare_higher = 0.0;
        std::vector<double> bestParameters_higher;
        
        int nExtrapolationPoints_higher = 0;
        for(int i=0; i<nDivisions_TLC; i++)
        {
            int h = nDivisions_TLC - 1 - i;
            
            double cutCenterValue = min_TLC + ((h+0.5)*division_TLC);
            
            if(validTLCDataPoints_perCut[h])
            {
                grPositionExtrapolation_higher->SetPoint(nExtrapolationPoints_higher, cutCenterValue, lineshapeCorrectionPositions[h][j].first);
                grPositionExtrapolation_higher->SetPointError(nExtrapolationPoints_higher, division_TLC/2.0, lineshapeCorrectionPositions[h][j].second);
                
                //--------------------------------------------
                std::unique_ptr<TF1> linearFit_higher(new TF1("linearFit_higher", "pol1(0)", -1000., 1000.));
                linearFit_higher->SetParameters(0, 0.0e+00);
                linearFit_higher->SetParameters(1, 0.0e+00);
                
                if(nExtrapolationPoints_higher==0)
                {
                    linearFit_higher->FixParameter(1, 0.0e+00);
                }
                linearFit_higher->SetParameters(0, 0.0e+00);
                linearFit_higher->SetParameters(1, 0.0e+00);
                
                if(nExtrapolationPoints_higher>=1)
                {
                    double x1, x2, y1, y2;
                    grPositionExtrapolation_higher->GetPoint(0, x1, y1);
                    grPositionExtrapolation_higher->GetPoint(grPositionExtrapolation_higher->GetN()-1, x2, y2);
                    
                    double a, b;
                    
                    a = (y2-y1)/(x2-x1);
                    b = (-a*x1) + y1;
                    
                    linearFit_higher->SetParameter(0, b);
                    linearFit_higher->SetParameter(1, a);
                }
                
                grPositionExtrapolation_higher->Fit("linearFit_higher", "Q+", "", min_TLC, max_TLC);
                
                if(nExtrapolationPoints_higher==0)
                {
                    linearFit_higher->ReleaseParameter(1);
                }
                
                //--------------------------------------------
                double chiSquare = linearFit_higher->GetChisquare();
                int nDF = linearFit_higher->GetNDF();
                double reducedChiSquare = (chiSquare/nDF);
                
                nExtrapolationPoints_higher++;
                
                //--------------------------------------------
                //  If a better Reduced Chi Square is found
                if((forceExtrapolationPoints_higher && nExtrapolationPoints_higher==nForcedExtrapolationPoints_higher) || (!forceExtrapolationPoints_higher && (nExtrapolationPoints_higher==3 || (nExtrapolationPoints_higher>3 && reducedChiSquare<bestReducedChiSquare_higher))))
                {
                    bestReducedChiSquare_higher = reducedChiSquare;
                    bestParameters_higher.clear();
                    
                    bestParameters_higher.push_back(linearFit_higher->GetParameter(0));
                    bestParameters_higher.push_back(linearFit_higher->GetParameter(1));
                    
                    nBestFitPoints_higher = nExtrapolationPoints_higher;
                    
                    averageError_higher += lineshapeCorrectionPositions[h][j].second;
                }
            }
        }
        
        if(bestParameters_higher.empty())
        {
            std::cout << "bestParameters_higher: EMPTY!" << std::endl;
        }
        
        if(bestParameters_lower.empty())
        {
            std::cout << "bestParameters_lower: EMPTY!" << std::endl;
        }
        
        grPositionExtrapolation_higher->Write();
        averageError_higher *= (1.0/nBestFitPoints_higher);
        averageErrors_higher.push_back(averageError_higher);
        bestExtrapolationParameters_higher.push_back(bestParameters_higher);
    }
    
    if(extendToLimits)
    {
        //      This process forces the extrapolation of all the peak positions for a given cut if the fit for any individual peak is invalid
        for(int i=0; i<nDivisions_TLC; i++)
        {
            bool foundInvalid = false;
            
            for(int j=0; j<nLineshapeCorrectionPeaks; j++)
            {
                if(!validTLCDataPoints_perPeak[i][j])
                {
                    foundInvalid = true;
                }
            }
            
            if(foundInvalid)
            {
                for(int j=0; j<nLineshapeCorrectionPeaks; j++)
                {
                    validTLCDataPoints_perPeak[i][j] = false;
                }
            }
        }
        
        for(int j=0; j<nLineshapeCorrectionPeaks; j++)
        {
            for(int i=0; i<nDivisions_TLC; i++)
            {
                double cutCenterValue = min_TLC + ((i+0.5)*division_TLC);
                double extrapolatedValue = 0.0;
                double extrapolatedValueError = 0.0;
                
                if(!validTLCDataPoints_perPeak[i][j])
                {
                    bool validExtrapolationParameters = false;
                    
                    if(cutCenterValue<0.0 && !bestExtrapolationParameters_lower.empty())
                    {
                        validExtrapolationParameters = true;
                        
                        extrapolatedValue = bestExtrapolationParameters_lower[j][0] + (bestExtrapolationParameters_lower[j][1]*cutCenterValue);
                        extrapolatedValueError = averageErrors_lower[j];
                    }
                    else if(cutCenterValue>=0.0 && !bestExtrapolationParameters_higher.empty())
                    {
                        validExtrapolationParameters = true;
                        
                        extrapolatedValue = bestExtrapolationParameters_higher[j][0] + (bestExtrapolationParameters_higher[j][1]*cutCenterValue);
                        extrapolatedValueError = averageErrors_higher[j];
                    }

                    if(validExtrapolationParameters)
                    {
                        lineshapeCorrectionPositions[i][j].first = extrapolatedValue;
                        lineshapeCorrectionPositions[i][j].second = extrapolatedValueError;
                        
                        validTLCDataPoints_perPeak[i][j] = true;
                    }
                }
            }
        }
        
        //----------------------------------------------------
        //      Setting each data point (per cut) as valid
        for(int i=0; i<nDivisions_TLC; i++)
        {
            bool valid = true;
            
            for(int j=0; j<nLineshapeCorrectionPeaks; j++)
            {
                if(!validTLCDataPoints_perPeak[i][j])
                {
                    valid = false;
                }
            }
            
            if(valid)
            {
                validTLCDataPoints_perCut[i] = true;
            }
        }
    }
    
    //--------------------------------------------------------------------
    for(int i=0; i<(int) lineshapeCorrectionPositions.size(); i++)
    {
        std::vector<std::pair<double, double>> tempLineshapeCorrectionPars;
        
        grLineshapeCorrections.push_back(std::shared_ptr<TGraphErrors>(new TGraphErrors()));
        sprintf(name, "grPositions_peak%d", i);
        grLineshapeCorrections[i]->SetName(name);
        
        int grPointN = 0;
        for(int j=0; j<(int) lineshapeCorrectionPositions[i].size(); j++)
        {
            grLineshapeCorrections[i]->SetPoint(grPointN, lineshapeCorrectionPositions[i][j].first, referenceX1posPositions_centralCut[j]);
            grPointN++;
        }
        
        //--------------------------------------------------------------
        //  Set all initial parameters of fLineshapeCorrection to 0.0
        for(int i=0; i<nLineshapeCorrectionPeaks; i++)
        {
            fLineshapeCorrection->SetParameter(i, 0.0);
        }
        
        grLineshapeCorrections[i]->Fit("fLineshapeCorrection", "Q+", "", 0.0, 1000.0);
        
        for(int k=0; k<nLineshapeCorrectionPeaks; k++)
        {
            tempLineshapeCorrectionPars.push_back(std::make_pair(fLineshapeCorrection->GetParameter(k), fLineshapeCorrection->GetParError(k)));
        }
        
        lineshapeCorrectionPars.push_back(tempLineshapeCorrectionPars);
        grLineshapeCorrections[i]->Write();
    }
    
    //----------------------------------------------------------------------------
    std::vector<std::shared_ptr<TGraphErrors>> grTLCCoefficients_perIteration;
    std::vector<std::shared_ptr<TGraphErrors>> grTLCExtendedCoefficients_perIteration;
    
    for(int j=0; j<nLineshapeCorrectionPeaks; j++)
    {
        std::shared_ptr<TGraphErrors> grTLCCoefficients(new TGraphErrors());
        sprintf(name, "grTLCCoefficients_peak%d", j);
        grTLCCoefficients->SetName(name);
        
        grTLCCoefficients->GetXaxis()->CenterTitle(true);
        grTLCCoefficients->GetYaxis()->CenterTitle(true);
        
        if(correctionType==0)
        {
            grTLCCoefficients->GetXaxis()->SetTitle("#theta_{SCAT}");
        }
        else if(correctionType==1)
        {
            grTLCCoefficients->GetXaxis()->SetTitle("Y1");
        }
        else if(correctionType==2)
        {
            grTLCCoefficients->GetXaxis()->SetTitle("tofCal");
        }
        else if(correctionType==3)
        {
            grTLCCoefficients->GetXaxis()->SetTitle("X1thCal");
        }
        else if(correctionType==4)
        {
            grTLCCoefficients->GetXaxis()->SetTitle("U1thCal");
        }

        sprintf(name, "Polynomial coefficient (order %d)", j);
        grTLCCoefficients->GetYaxis()->SetTitle(name);
        
        //------------------------------------------------------------------------------------
        //  Filling grTLCCoefficients to extract the polynomial coefficients
        double smallestValidPoint;
        double largestValidPoint;
        double summedError = 0.0;
        double averageError = 0.0;
        
        int grPointN = 0;
        for(int i=0; i<nDivisions_TLC; i++)
        {
            if(validTLCDataPoints_perCut[i])
            {
                if(grPointN==0)
                {
                    smallestValidPoint = min_TLC + ((i+0.5)*division_TLC);
                    largestValidPoint = min_TLC + ((i+0.5)*division_TLC);
                }
                if(min_TLC + ((i+0.5)*division_TLC) < smallestValidPoint)
                {
                    smallestValidPoint = min_TLC + ((i+0.5)*division_TLC);
                }
                if(min_TLC + ((i+0.5)*division_TLC) > largestValidPoint)
                {
                    largestValidPoint = min_TLC + ((i+0.5)*division_TLC);
                }
                
                grTLCCoefficients->SetPoint(grPointN, min_TLC + ((i+0.5)*division_TLC), lineshapeCorrectionPars[i][j].first);
                
                if(originalValidityTLCDataPoints_perCut[i])
                {
                    grTLCCoefficients->SetPointError(grPointN, division_TLC/2.0, lineshapeCorrectionPars[i][j].second);
                    summedError += lineshapeCorrectionPars[i][j].second;
                }
                else if(i%(nSkippedPoints+1)==0)
                {
                    double scaledError = errorWeighting*lineshapeCorrectionPars[i][j].second;
                    grTLCCoefficients->SetPointError(grPointN, division_TLC/2.0, scaledError);
                    summedError += scaledError;
                }
                
                grPointN++;
            }
        }
        
        averageError = summedError/grPointN;
        
        //--------------------------------------------------------------------------------------------
        std::vector<double> linearFitPars;
        std::vector<double> linearFitPars_lower;
        std::vector<double> linearFitPars_higher;
        
        std::unique_ptr<TF1> linearFit(new TF1("linearFit", "pol1(0)", -1000., 1000.));
        linearFit->SetParameters(0, 0.0e+00);
        linearFit->SetParameters(1, 0.0e+00);
        
        //------------------------------------------------------------------------------------------------
        grTLCCoefficients->Fit("linearFit", "Q", "", -100.0, 100.0);
        linearFitPars.push_back(linearFit->GetParameter(0));
        linearFitPars.push_back(linearFit->GetParameter(1));
        grTLCCoefficients->GetListOfFunctions()->Clear();
        
        //------------------------------------------------------------------------
        std::shared_ptr<TGraphErrors> grTLCExtendedCoefficients(new TGraphErrors());
        sprintf(name, "grTLCExtendedCoefficients_peak%d", j);
        grTLCExtendedCoefficients->SetName(name);
        
        int grExtendedPointN = 0;
        int grPointNCounter = 0;
        
        if(extendToLimits)
        {
            for(int i=0; i<nDivisions_TLC; i++)
            {
                if(validTLCDataPoints_perCut[i])
                {
                    if(!originalValidityTLCDataPoints_perCut[i] && i%(nSkippedPoints+1)==0)
                    {
                        double x, y;
                        double eX, eY;
                        
                        grTLCCoefficients->GetPoint(grPointNCounter, x, y);
                        eX = grTLCCoefficients->GetErrorX(grPointNCounter);
                        eY = grTLCCoefficients->GetErrorY(grPointNCounter);
                        
                        grTLCExtendedCoefficients->SetPoint(grExtendedPointN, x, y);
                        grTLCExtendedCoefficients->SetPointError(grExtendedPointN, eX, eY);
                        
                        grExtendedPointN++;
                    }
                    
                    grPointNCounter++;
                }
            }
        }
        
        grTLCExtendedCoefficients->Write();
        
        //------------------------------------------------------------------------
        std::vector<double> bestParameters;
        double bestReducedChiSquare = 0.0;
        std::vector<std::unique_ptr<TF1>> fitLineshapeCorrectionCoefficients;
        std::string bestFitName;
        
        for(int i=0; i<grPointN-1; i++)
        {
            sprintf(name, "fitLineshapeCorrectionCoefficients%d", i);
            sprintf(condition, "pol%d(0)", i);
            fitLineshapeCorrectionCoefficients.push_back(std::unique_ptr<TF1>(new TF1(name, condition, -1000., 1000.)));
            
            for(int k=0; k<i; k++)
            {
                fitLineshapeCorrectionCoefficients[i]->SetParameter(k, 0.0);
            }
            
            fitLineshapeCorrectionCoefficients[i]->SetNpx(5*grPointN);
            fitLineshapeCorrectionCoefficients[i]->SetLineStyle(7);
            fitLineshapeCorrectionCoefficients[i]->SetLineWidth(1);
            fitLineshapeCorrectionCoefficients[i]->SetLineColor(kBlue);
            
            grTLCCoefficients->Fit(name, "Q+", "", -100.0, 100.0);
            
            double chiSquare = fitLineshapeCorrectionCoefficients[i]->GetChisquare();
            int nDF = fitLineshapeCorrectionCoefficients[i]->GetNDF();
            double reducedChiSquare = (chiSquare/nDF);
            
            //--------------------------------------------
            //  If a better Reduced Chi Square is found
            if(i==0 || reducedChiSquare<bestReducedChiSquare)
            {
                bestReducedChiSquare = reducedChiSquare;
                bestParameters.clear();
                
                for(int k=0; k<i+1; k++)
                {
                    bestParameters.push_back(fitLineshapeCorrectionCoefficients[i]->GetParameter(k));
                }
                
                bestFitName = name;
            }
        }
        
        //--------------------------------------------------------------------------
        //  Highlighting the best fitting polynomial (smallest reduced chi square)
        (((TF1*) grTLCCoefficients->GetListOfFunctions()->FindObject(bestFitName.c_str())))->SetLineColor(kRed);
        (((TF1*) grTLCCoefficients->GetListOfFunctions()->FindObject(bestFitName.c_str())))->SetLineStyle(1);
        (((TF1*) grTLCCoefficients->GetListOfFunctions()->FindObject(bestFitName.c_str())))->SetLineWidth(4);
        
        (((TF1*) grTLCCoefficients->GetListOfFunctions()->FindObject(bestFitName.c_str())))->SetNpx(5*grPointN);
        
        grTLCCoefficients->Write();
        
        lineshapeCorrectionPolPars.push_back(bestParameters);
        grTLCCoefficients_perIteration.push_back(grTLCCoefficients);
        grTLCExtendedCoefficients_perIteration.push_back(grTLCExtendedCoefficients);
    }
    
    finalLineshapeCorrectionTypes.push_back(correctionType);
    finalLineshapeCorrectionParameters.push_back(lineshapeCorrectionPolPars);
    
    grTLCCoefficientsStorage.push_back(grTLCCoefficients_perIteration);
    grTLCExtendedCoefficientsStorage.push_back(grTLCExtendedCoefficients_perIteration);
    
    f.Write();
}

double fPeakPosition_voigt(double *x, double *par)
{
    double fitFunctionValue = 0.0;
    
    fitFunctionValue += par[0]*TMath::Voigt(x[0] - par[1], par[2], par[3], 4) + par[4] + pow(par[5], 2.0);
    
    return fitFunctionValue;
}


#endif


