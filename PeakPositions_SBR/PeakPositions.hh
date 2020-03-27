
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
#include "TF1Convolution.h"
#include "TRandom.h"
#include "TMath.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

#include "../../../Analysis/cuts/CUTS.hh"

double gaussian(double *x, double *par);
double landau(double *x, double *par);
double landauMirror(double *x, double *par);
double cauchy(double *x, double *par);

char name[512];
char condition[512];
char command[512];
std::string directoryPath;

std::string primitiveCut;

std::vector<int> runNumbers;
int nrOfRuns;

int nBins_X1posC;
double min_X1posC, max_X1posC;

std::unique_ptr<TCanvas> c1;

std::vector<std::unique_ptr<TChain>> DATAChains;
std::vector<TEventList*> eventLists;

//----------------------------------
//      Peak position mapping
bool usePeakPositionMapping;
int nMappingPeaks;
TF1 *fMapping;
std::vector<double> referencePeakPositions;
std::vector<std::tuple<double, double>> defaultX1posCMappingRanges;
std::vector<std::tuple<int, int, double, double>> X1posCRanges;
std::vector<std::vector<double>> peakPositions;
std::vector<std::vector<double>> peakPositionErrors;
std::vector<std::vector<double>> mappingParameters;
std::vector<double> initialCorrectionPeakMeans;


void GenerateGoodEventLists();

void SetPrimitiveCut(std::string s)
{
    primitiveCut = s;
}

//--------------------------------------------
void InitialiseVariables()
{
    nMappingPeaks = 0;
    
    c1 = std::unique_ptr<TCanvas>(new TCanvas("TotalLineshapeCorrection","",10,10,700,700));
    
    usePeakPositionMapping = false;
    
    initialiseCuts();

    GenerateGoodEventLists();
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
        
        sprintf(filePath, "%s%3d.root", filePath, runNumbers[i]);
        DATAChains[i]->Add(filePath);
        
        //----------------------------------------------------------------
        
        sprintf(name, "goodEventList_%d", i);
        sprintf(command, ">>%s", name);
//        DATAChains[i]->Draw(command, "cut_Pad1vsTOF && goodSbrEvent && U1flag==0 && X2flag==0 && Y1>-18 && Y1<18 && thetaSCAT>-2.0 && thetaSCAT<2.0");
        DATAChains[i]->Draw(command, primitiveCut.c_str());
        DATAChains[i]->SetEventList((TEventList*)gDirectory->Get(name));
        eventLists.push_back((TEventList*)gDirectory->Get(name));
    }
    
    std::cout << std::endl;
}

//--------------------------------------------
void MapPeakPositions()
{
    usePeakPositionMapping = true;
    
    TFile f("Output/PeakPositions.root", "recreate");
    
    std::vector<std::unique_ptr<TH1F>> hX1posCThetSCATgated;
    std::vector<std::unique_ptr<TH1F>> hX1posCThetSCATgated_corrected;
    
    nBins_X1posC = 4000;
    min_X1posC = 0.0;
    max_X1posC = 1000.0;

    std::unique_ptr<TH1F> hX1posCThetSCATgated_corrected_total(new TH1F("hX1posCThetSCATgated_corrected_total", "", nBins_X1posC, min_X1posC, max_X1posC));
    std::unique_ptr<TH1F> hX1posCThetSCATgated_total(new TH1F("hX1posCThetSCATgated_total", "", nBins_X1posC, min_X1posC, max_X1posC));
    
    std::unique_ptr<TF1> fitPeakPosition(new TF1("fitPeakPosition", "gaus(0) + pol1(3)", 0., 1000.));
    fitPeakPosition->SetParameters(1.0e+01, 7.20e+02, 1.0e+00, 0.0e+00, 0.0e+00);
    fitPeakPosition->SetLineColor(kRed);
    fitPeakPosition->SetNpx(4*nBins_X1posC);
    fitPeakPosition->SetParLimits(0, 0.0, 1.0e+08);
    
    if((int) referencePeakPositions.size() != nMappingPeaks)
    {
        referencePeakPositions.clear();
    }
    
    double par[5];
    
    //================================
    //      Mapping Fit Function
    
    if(nMappingPeaks==1)
    {
        sprintf(condition, "x + [0]");
    }
    else
    {
        sprintf(condition, "pol%d(0)", nMappingPeaks-1);
    }
    
    fMapping = new TF1("fMapping", condition, 0.0, 1000.0);
    
    int mappingOrder;
    int maxMappingOrder = 7;
    
    if(nMappingPeaks<=maxMappingOrder && nMappingPeaks>=1)
    {
        mappingOrder = nMappingPeaks-1;
    }
    else
    {
        mappingOrder = maxMappingOrder-1;
    }
    //==========================================
    //      Alternate mapping Fit Function
    
    auto fGaus = std::make_shared<TF1>("fGaus", "[0] + [1]*x + [2]*TMath::Gaus(x, [3], [4])", -50.0, 1000.0);
    fGaus->SetNpx(1000);
    fGaus->SetParameter(0, 1.0);
    fGaus->SetParameter(1, 0.0);
    fGaus->SetParameter(2, 1.0);
    
    auto fLandau = std::make_shared<TF1>("fLandau", landau, -50., 1000., 3);
    fLandau->SetNpx(1000);
    
    auto fLandauMirror = std::make_shared<TF1>("fLandauMirror", "[0]*TMath::Landau((2.0*[1]) - x, [1], [2])", -50.0, 1000.0);
    fLandauMirror->SetNpx(1000);
    fLandauMirror->SetParameter(0, 1.0);
    fLandauMirror->SetParameter(1, 0.0);
    fLandauMirror->SetParameter(2, 1.0);

    //--------------------------------------------------------------------
    TF1Convolution *f_conv = new TF1Convolution(fGaus.get(),fLandauMirror.get(),-50.0,1000.0,true);
    //TF1Convolution *f_conv = new TF1Convolution("[0]*TMath::Gaus(x, [1], [2])", "[0]*TMath::Landau((2.0*[1]) - x, [1], [2])",-50.0,1000.0,true);
    f_conv->SetRange(-100.0,1000.0);
    f_conv->SetNofPointsFFT(100000);
    
    auto fConvolution = std::make_shared<TF1>("fConvolution",*f_conv, -50., 1000.0, f_conv->GetNpar());
    fConvolution->SetNpx(100000);

    //========================================================================
    
    std::cout << "//==================================//" << std::endl;
    std::cout << "      Calculating X1posC offsets" << std::endl;
    
    std::vector<std::unique_ptr<TGraphErrors>> grMappings;
    
    for(int i=0; i<nrOfRuns; i++)
    {
        std::cout << "\rRun (" << i+1 << "/" << nrOfRuns << ")" << std::flush;
        
        //========================================================================
        
        sprintf(name, "hX1posCThetSCATgated_run%d", runNumbers[i]);
        hX1posCThetSCATgated.push_back(std::unique_ptr<TH1F>(new TH1F(name, "", nBins_X1posC, min_X1posC, max_X1posC)));
        
        sprintf(name, "hX1posCThetSCATgated_run%d_corrected", runNumbers[i]);
        hX1posCThetSCATgated_corrected.push_back(std::unique_ptr<TH1F>(new TH1F(name, "", nBins_X1posC, min_X1posC, max_X1posC)));
        
        sprintf(name, "hX1posCThetSCATgated_run%d", runNumbers[i]);
        //sprintf(command, "X1posC>>%s", name);
        //sprintf(command, "xPositions_TLC[0]>>%s", name);
        sprintf(command, "sbrXPosition_TLC>>%s", name);
        
        DATAChains[i]->Draw(command,"cut_Pad1vsTOF && goodSbrEvent && Y1>-18 && Y1<18 && thetaSCAT>-2.0 && thetaSCAT<2.0","");
        
        hX1posCThetSCATgated_total->Add(hX1posCThetSCATgated[i].get(), 1.0);
        
        //========================================================================
        
        std::vector<double> tempPeakPositions;
        std::vector<double> tempPeakPositionErrors;
        
        for(int j=0; j<nMappingPeaks; j++)
        {
            double tempmin_X1posC, tempmax_X1posC;
            
            bool found = false;
            for(auto k = X1posCRanges.begin(); k != X1posCRanges.end(); ++k)
            {
                if(std::get<0>((*k))==j && std::get<1>((*k))==runNumbers[i])
                {
                    found = true;
                    tempmin_X1posC = std::get<2>((*k));
                    tempmax_X1posC = std::get<3>((*k));
                }
            }
            
            if(!found)
            {
                tempmin_X1posC = std::get<0>(defaultX1posCMappingRanges[j]);
                tempmax_X1posC = std::get<1>(defaultX1posCMappingRanges[j]);
            }
            
            //========================================================================
            hX1posCThetSCATgated[i]->GetXaxis()->SetRangeUser(tempmin_X1posC, tempmax_X1posC);
            double initialPeakHeight = hX1posCThetSCATgated[i]->GetBinContent(hX1posCThetSCATgated[i]->GetMaximumBin());
            double initialMean = hX1posCThetSCATgated[i]->GetBinCenter(hX1posCThetSCATgated[i]->GetMaximumBin());
            double initialStdDev = hX1posCThetSCATgated[i]->GetStdDev();

            //----------------------------
            if(j==1)
            {
                tempmin_X1posC = initialMean-9.0;
                tempmax_X1posC = initialMean+6.0;
            }
            if(j==2)
            {
                tempmin_X1posC = initialMean-1.0;
                tempmax_X1posC = initialMean+1.5;
            }

            //----------------------------
            fitPeakPosition->SetParameter(0, initialPeakHeight);
            fitPeakPosition->SetParameter(1, initialMean);
            fitPeakPosition->SetParameter(2, initialStdDev);
            
            fitPeakPosition->SetParLimits(0, 0.0e+00, 3.0*initialPeakHeight);
            fitPeakPosition->SetParLimits(1, tempmin_X1posC, tempmax_X1posC);
            fitPeakPosition->SetParLimits(2, 0.2*initialStdDev, 2.0*initialStdDev);
            
            //========================================================================
            /*
            fConvolution->SetParameter(0, 0.0);
            fConvolution->SetParameter(1, 0.0);
            fConvolution->SetParameter(2, initialPeakHeight);
            fConvolution->SetParameter(3, initialMean);
            fConvolution->SetParameter(4, initialStdDev);
            fConvolution->FixParameter(5, 1.0e+00);
            fConvolution->FixParameter(6, 0.0e+00);
            fConvolution->SetParameter(7, 3.0e-01);
            
            fConvolution->SetParLimits(0, 0.0, 2.0e+03);
            fConvolution->SetParLimits(1, -1.0e-04, 1.0e-04);
            fConvolution->SetParLimits(2, 0.0, 1.0e+01*initialPeakHeight);
            */
            
            fConvolution->SetParameter(0, 0.0);
            fConvolution->SetParameter(1, 0.0);
            fConvolution->SetParameter(2, initialPeakHeight);
            fConvolution->SetParameter(3, initialMean);
            //fConvolution->SetParameter(4, initialStdDev);
            fConvolution->SetParameter(4, 1.5e+00);
            fConvolution->FixParameter(5, 1.0e+00);
            fConvolution->FixParameter(6, 0.0e+00);
            fConvolution->SetParameter(7, 1.0e-01);
            
            fConvolution->SetParLimits(0, 0.0, 1.0e+01);
            fConvolution->SetParLimits(1, -1.0e-05, 1.0e-05);
            fConvolution->SetParLimits(2, 0.0, 1.0e+01*initialPeakHeight);
            fConvolution->SetParLimits(4, 2.0e-01, 3.5e+00);

            if(runNumbers[i]==1090)
            {
                if(j==1)
                {
                    fConvolution->SetParameter(0, 2.0e+01);
                    fConvolution->SetParLimits(0, 0.0, 4.0e+01);
                }
                else if(j==2)
                {
                    fConvolution->SetParameter(0, 0.0e+00);
                    fConvolution->SetParLimits(0, 0.0, 1.0e+01);
                }
            }
            else if(runNumbers[i]==1114)
            {
                if(j==1)
                {
                    fConvolution->SetParameter(0, 1.1e+02);
                    fConvolution->SetParLimits(0, 5.0e+01, 1.5e+02);
                }
                else if(j==2)
                {
                    fConvolution->SetParameter(0, 1.1e+02);
                    fConvolution->SetParLimits(0, 5.0e+01, 1.5e+02);
                }
            }
            else if(runNumbers[i]==1115 || runNumbers[i]==1120)
            {
                if(j==1)
                {
                    fConvolution->SetParameter(0, 3.0e+01);
                    fConvolution->SetParLimits(0, 0.0, 4.0e+01);
                }
                else if(j==2)
                {
                    fConvolution->SetParameter(0, 0.0e+00);
                    fConvolution->SetParLimits(0, 0.0, 1.0e+01);
                }
            }
            else if(runNumbers[i]==1121 || runNumbers[i]==1131 || runNumbers[i]==1138)
            {
                if(j==1)
                {
                    fConvolution->SetParameter(0, 3.0e+01);
                    fConvolution->SetParLimits(0, 0.0, 5.0e+01);
                }
                else if(j==2)
                {
                    fConvolution->SetParameter(0, 1.0e+01);
                    fConvolution->SetParLimits(0, 0.0, 1.5e+01);
                }
            }
            else if(runNumbers[i]==1125)
            {
                if(j==1)
                {
                    fConvolution->SetParameter(0, 3.0e+01);
                    fConvolution->SetParLimits(0, 0.0, 5.0e+01);
                }
                else if(j==2)
                {
                    fConvolution->SetParameter(0, 2.0e+01);
                    fConvolution->SetParLimits(0, 0.0, 3.0e+01);
                }
            }
            else if(runNumbers[i]==1156)
            {
                if(j==1)
                {
                    fConvolution->SetParameter(0, 3.5e+01);
                    fConvolution->SetParLimits(0, 0.0, 5.0e+01);
                }
                else if(j==2)
                {
                    fConvolution->SetParameter(0, 2.5e+01);
                    fConvolution->SetParLimits(0, 0.0, 3.0e+01);
                }
            }
            else if(runNumbers[i]==1157)
            {
                if(j==1)
                {
                    fConvolution->SetParameter(0, 6.0e+01);
                    fConvolution->SetParLimits(0, 2.0e+01, 1.2e+02);
                }
                else if(j==2)
                {
                    fConvolution->SetParameter(0, 0.0e+00);
                    fConvolution->SetParLimits(0, 0.0, 2.0e+01);
                }
            }

            //========================================================================
            
            if(j>0)
            {
                hX1posCThetSCATgated[i]->Fit("fConvolution", "Q+", "", tempmin_X1posC, tempmax_X1posC);
            }
            else if(j==0)
            {
                //tempmin_X1posC -= 22.0;
                
                if(runNumbers[i]==1119)
                {
                    tempmin_X1posC = 431.0;
                    tempmax_X1posC = 470.0;
                }

                if(runNumbers[i]==1120)
                {
                    tempmin_X1posC = 440.0;
                    tempmax_X1posC = 470.0;
                }
                if(runNumbers[i]==1131 || runNumbers[i]==1156)
                {
                    tempmin_X1posC = 437.0;
                    tempmax_X1posC = 470.0;
                }

                sprintf(name, "hTemp_run%d", runNumbers[i]);
                TH1F *hTemp = (TH1F*) hX1posCThetSCATgated[i]->Rebin(16,name);
                hTemp->GetXaxis()->SetRangeUser(tempmin_X1posC, tempmax_X1posC);
                hTemp->Smooth(2);
                
                double initialPeakHeight = hTemp->GetBinContent(hTemp->GetMaximumBin());
                double initialMean = hTemp->GetBinCenter(hTemp->GetMaximumBin());
                double initialStdDev = hTemp->GetStdDev();

                fConvolution->SetParameter(0, 1.5*initialPeakHeight);
                fConvolution->SetParameter(1, -1.0e+00);
                fConvolution->SetParameter(2, initialPeakHeight);
                fConvolution->SetParameter(3, initialMean);
                fConvolution->SetParameter(4, 6.0);
                fConvolution->FixParameter(5, 1.0e+00);
                fConvolution->FixParameter(6, 0.0e+00);
                fConvolution->SetParameter(7, 2.0e-01);
                
                fConvolution->SetParLimits(0, 0.0, 6.0*initialPeakHeight);
                fConvolution->SetParLimits(1, -7.0e+00, 0.0e+00);
                fConvolution->SetParLimits(2, 0.0, 1.0e+01*initialPeakHeight);
                fConvolution->SetParLimits(4, 4.0, 9.0);
                fConvolution->SetParLimits(7, 5.0e-02, 4.0e-01);

                /*
                tempmin_X1posC = initialMean - 30.0;
                tempmax_X1posC = initialMean + 25.0;
                */
                
                if(runNumbers[i]==1090)
                {
                    fConvolution->SetParameter(2, initialPeakHeight);
                    fConvolution->SetParLimits(2, 0.0, 2.0e+01*initialPeakHeight);

                    tempmin_X1posC = initialMean - 30.0;
                    tempmax_X1posC = initialMean + 25.0;
                }
                else if(runNumbers[i]==1114)
                {
                    tempmin_X1posC = initialMean - 22.0;
                    tempmax_X1posC = initialMean + 18.0;
                }
                else if(runNumbers[i]==1115)
                {
                    tempmin_X1posC = initialMean - 25.0;
                    tempmax_X1posC = initialMean + 16.0;
                }
                else if(runNumbers[i]==1120)
                {
                    fConvolution->SetParameter(2, initialPeakHeight);
                    fConvolution->SetParLimits(2, 0.0, 2.0e+01*initialPeakHeight);

                    tempmin_X1posC = initialMean - 30.0;
                    tempmax_X1posC = initialMean + 25.0;
                }
                else if(runNumbers[i]==1121)
                {
                    tempmin_X1posC = initialMean - 30.0;
                    tempmax_X1posC = initialMean + 15.0;
                }
                else if(runNumbers[i]==1125)
                {
                    tempmin_X1posC = initialMean - 30.0;
                    tempmax_X1posC = initialMean + 20.0;
                }
                else if(runNumbers[i]==1131)
                {
                    tempmin_X1posC = initialMean - 30.0;
                    tempmax_X1posC = initialMean + 16.0;
                }
                else if(runNumbers[i]==1138)
                {
                    tempmin_X1posC = initialMean - 20.0;
                    tempmax_X1posC = initialMean + 20.0;
                }
                else if(runNumbers[i]==1156)
                {
                    tempmin_X1posC = initialMean - 30.0;
                    tempmax_X1posC = initialMean + 18.0;
                }
                else if(runNumbers[i]==1157)
                {
                    tempmin_X1posC = initialMean - 25.0;
                    tempmax_X1posC = initialMean + 25.0;
                }

                
                hTemp->Fit("fConvolution", "Q+", "", tempmin_X1posC, tempmax_X1posC);
                
                hTemp->GetXaxis()->SetRangeUser(tempmin_X1posC - 20.0, tempmax_X1posC + 20.0);

                std::cout << "Run Number: " << runNumbers[i] << ", Position: " << fConvolution->GetParameter(3) << std::endl;
            }

            //hX1posCThetSCATgated[i]->Fit("fitPeakPosition", "Q+", "", tempmin_X1posC, tempmax_X1posC);
            //hX1posCThetSCATgated[i]->Fit("fConvolution", "Q+", "", tempmin_X1posC, tempmax_X1posC);
            
            tempPeakPositions.push_back(fConvolution->GetParameter(3));
            tempPeakPositionErrors.push_back(fConvolution->GetParError(3));
        }
        
        hX1posCThetSCATgated[i]->GetXaxis()->SetRangeUser(min_X1posC, max_X1posC);
        
        peakPositions.push_back(tempPeakPositions);
        peakPositionErrors.push_back(tempPeakPositionErrors);
        
        //========================================================================
        
        sprintf(name, "grMapping_run%d", runNumbers[i]);
        grMappings.push_back(std::unique_ptr<TGraphErrors>(new TGraphErrors()));
        grMappings[i]->SetName(name);

        for(int j=0; j<nMappingPeaks; j++)
        {
            if(!referencePeakPositions.empty())
            {
                grMappings[i]->SetPoint(j, peakPositions[i][j], referencePeakPositions[j]);
            }
            else
            {
                grMappings[i]->SetPoint(j, peakPositions[i][j], peakPositions[0][j]);
            }
            
            //grMappings[i]->SetPointError(j, 0.0, peakPositionErrors[i][j]);
        }
        
        grMappings[i]->SetMarkerStyle(1);
        grMappings[i]->SetMarkerSize(7);
        grMappings[i]->SetMarkerColor(kBlack);

        for(int k=0; k<mappingOrder; k++)
        {
            fMapping->SetParameter(k, 0.0e+00);
        }
        
        if(nMappingPeaks>=1)
        {
            double x1, x2, y1, y2;
            grMappings[i]->GetPoint(0, x1, y1);
            grMappings[i]->GetPoint(grMappings[i]->GetN()-1, x2, y2);
            
            double a, b;
            
            a = (y2-y1)/(x2-x1);
            b = (-a*x1) + y1;
            
            fMapping->SetParameter(0, b);
            fMapping->SetParameter(1, a);
        }
        
        grMappings[i]->Fit("fMapping", "Q+", "", 0.0, 1000.0);
        
        double chiSquare = fMapping->GetChisquare();
        int nDF = fMapping->GetNDF();
        double reducedChiSquare;
        
        /*
        if(nDF>=1)
        {
            reducedChiSquare = (chiSquare/nDF);
            std::cout << "Run: " << runNumbers[i] << ", reducedChiSquare: " << reducedChiSquare << std::endl;
        }
        else
        {
            std::cout << "Run: " << runNumbers[i] << ", chiSquare: " << chiSquare << std::endl;
        }
        */

        grMappings[i]->Write();
        
        std::vector<double> tempMappingParameters;
        
        for(int j=0; j<nMappingPeaks; j++)
        {
            tempMappingParameters.push_back(fMapping->GetParameter(j));
            //std::cout << "fMapping->GetParameter(j): " << fMapping->GetParameter(j) << std::endl;
        }
        
        
        mappingParameters.push_back(tempMappingParameters);
        
        //========================================================================
        /*
        char drawCondition[1028];
        
        if(nMappingPeaks==1)
        {
            sprintf(condition, "X1posC + (%f)", mappingParameters[i][0]);
        }
        else
        {
            for(int j=0; j<nMappingPeaks; j++)
            {
                if(j==0)
                {
                    sprintf(condition, "(%f)", mappingParameters[i][0]);
                }
                else
                {
                    sprintf(condition, "%s + (%f)*pow(X1posC, %d)", condition, mappingParameters[i][j], j);
                }
            }
        }
        
        sprintf(name, "hX1posCThetSCATgated_run%d_corrected", runNumbers[i]);
        sprintf(drawCondition, "%s >> %s", condition, name);
        
        DATAChains[i]->Draw(drawCondition, "cut_Pad1vsTOF && X1flag==0 && U1flag==0 && X2flag==0 && Y1>-18 && Y1<18 && thetaSCAT>-2.0 && thetaSCAT<2.0", "");
        hX1posCThetSCATgated_corrected_total->Add(hX1posCThetSCATgated_corrected[i].get(), 1.0);
        */
        
        TTreeReader reader("DATA", DATAChains[i]->GetFile());
        
        TTreeReaderValue<double> fX1posC(reader, "X1posC.t_X1posC");
        //TTreeReaderValue<double> fThetaSCAT(reader, "thetaSCAT.t_thetaSCAT");
        //TTreeReaderValue<double> fY1(reader, "Y1.t_Y1");
        
        TTreeReaderValue<bool>      fgoodSbrEvent(reader, "goodSbrEvent");
        TTreeReaderValue<int>       fgoodSbrEventNr(reader, "goodSbrEventNr");
        TTreeReaderValue<double>    fsbrXPosition(reader, "sbrXPosition");
        TTreeReaderValue<double>    fsbrXPosition_TLC(reader, "sbrXPosition_TLC");
        
        int eN = 0;
        while(eventLists[i]->GetEntry(eN)>0)
        {
            reader.SetEntry(eventLists[i]->GetEntry(eN));
            
            double originalX1posC = (*fX1posC);
            double correctedX1posC = 0.0;
            
            bool    goodSbrEvent = *fgoodSbrEvent;
            int     goodSbrEventNr = *fgoodSbrEventNr;
            
            double  sbrXPosition = *fsbrXPosition;
            double  sbrXPosition_TLC = *fsbrXPosition_TLC;
            double  sbrXPosition_TLC_mapped = 0.0;
            
            for(int j=0; j<(int) mappingParameters[i].size(); j++)
            {
                sbrXPosition_TLC_mapped += mappingParameters[i][j]*pow(sbrXPosition_TLC, j);
            }
            
            //  TEMP
            //correctedX1posC = (*fX1posC);
            
            //hX1posCThetSCATgated[i]->Fill(originalX1posC);
            hX1posCThetSCATgated_corrected[i]->Fill(sbrXPosition_TLC_mapped);
            hX1posCThetSCATgated_corrected_total->Fill(sbrXPosition_TLC_mapped);
            
            eN++;
        }
        
        reader.Restart();
    }
    
    std::cout << std::endl;

    for(int j=0; j<(int) peakPositions[0].size(); j++)
    {
        std::cout << "Peak positions of first run: " << std::endl;
        std::cout << std::setprecision(12) << peakPositions[0][j] << " " << std::endl;
    }
    
    //--------------------------------------------------------------------------------
    std::ofstream fileX1offsets("Output/PR194_XMappings.dat", std::ofstream::out);
    
    for(int i=0; i<nrOfRuns; i++)
    {
        fileX1offsets << runNumbers[i];
        
        fileX1offsets << " " << nMappingPeaks;
        
        for(int j=0; j<nMappingPeaks; j++)
        {
            fileX1offsets << " " << mappingParameters[i][j];
        }
        
        fileX1offsets << "\n";
    }
    
    fileX1offsets << "eof\n";
    
    fileX1offsets.close();
    f.Write();
}

double PreLineshapeCorrectionMap(double origX1posC, std::vector<double> mapPars)
{
    double mappedX1posC = 0.0;
    
    if((int) mapPars.size()==1)
    {
        mappedX1posC = origX1posC + mapPars[0];
    }
    else
    {
        for(int i=0; i<(int) mapPars.size(); i++)
        {
            mappedX1posC += mapPars[i]*pow(origX1posC, i);
        }
    }
    
    return mappedX1posC;
}

double gaussian(double *x, double *par) {
    
    double fitFunctionValue = par[0]*TMath::Gaus(x[0], par[1], par[2]);
    
    return fitFunctionValue;
}

double landauMirror(double *x, double *par) {
    
    double fitFunctionValue = par[0]*TMath::Landau((2.0*par[1]) - x[0], par[1], par[2]);
    
    return fitFunctionValue;
}

double landau(double *x, double *par) {
    
    double fitFunctionValue = par[0]*TMath::Landau(x[0], par[1], par[2]);
    
    return fitFunctionValue;
}

double cauchy(double *x, double *par) {
    
    double fitFunctionValue = par[0]*TMath::CauchyDist(x[0], par[1], par[2]);
    
    return fitFunctionValue;
}


#endif


