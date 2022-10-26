/**
  \file
  Analysis script for calculating y-dependent TPC sector drit loss. These gains 
  are stored in a text file in the database and applied by the TPCDEDXCalculatorSG
  module. Total cluster charges are accumulated and averaged as a function of binned
  y-position. Asymmetry around the beam Y (y=0) is the calculated, and an asymmetry
  factor is stored in the text file.

  \author B. Rumberger
  \version $Id:    $
  \date 24 Oct 2019
*/

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <set>

#include <det/TPCConst.h>
#include <modutils/PeakFinder.h>

#include "TH1D.h"

//Typedefs and containers for holding histograms.
typedef std::unordered_map<unsigned int, TH1D*> PadHistograms;
typedef std::unordered_map<unsigned int, PadHistograms> PadrowHistograms;
typedef std::unordered_map<unsigned int, PadrowHistograms> SectorHistograms;
typedef std::unordered_map<unsigned int, SectorHistograms> DetectorHistograms;
//The container itself.
DetectorHistograms fSpectraHistograms;

//Typedefs and containers for peak finders.
typedef std::unordered_map<unsigned int, modutils::PeakFinder> PadPeakFinders;
typedef std::unordered_map<unsigned int, PadPeakFinders> PadrowPeakFinders;
typedef std::unordered_map<unsigned int, PadrowPeakFinders> SectorPeakFinders;
typedef std::unordered_map<unsigned int, SectorPeakFinders> DetectorPeakFinders;
DetectorPeakFinders fPeakFinders;

//Typedefs and containers for average sector peaks.
typedef std::unordered_map<unsigned int, double> SectorPeaks;
typedef std::unordered_map<unsigned int, SectorPeaks> DetectorPeaks;
DetectorPeaks fAverageSectorPeaks;

//Config parameters.
std::set<det::TPCConst::EId> fTPCIdList;
std::string fFitFunction;
double fMinAcceptableGain;
double fMaxAcceptableGain;
unsigned int fMinHistogramEntries;
unsigned int fHistogramBins;
double fHistogramPadding;
unsigned int fMinPads;
unsigned int fMaxPads;
unsigned int fMinTimeSliceNumber;
unsigned int fMinTimeSlices;
unsigned int fMaxTimeSlices;
double fMaxADCCut;
double fChargeCut;
double fMinADCPeakSearch;
double fMinADCPeakSearchVTPC1Upstream;


/// Main function.
int main(int argc, char* argv[]);

/// Configuration file parsing function.
void ParseConfigFile(const std::string& configFile);

/// Function for replaing default pad gain XML with user-defined XML.
void ReplacePadGainPath(const std::string& bootstrap,
                        const std::string& newPadGainXML);

/// Function for replaing default pad gain XML with user-defined XML.
std::map<int,int> GetPadrowColorMap(const int maxPadrows);

// Gauss function.
double Gauss(const double x,
             const double mean,
             const double sigma,
             const double amplitude = 1.)
{
  return amplitude*std::exp(-(x-mean)*(x-mean)/(2.0*sigma*sigma));
}


// Display usage.
void DisplayUsage()
{
  std::cerr << "\nUsage:\n\tKryptonAnalyzer -o outputPrefix "
    "[ (-c / --config) configFilePath] [ (-u / --updateGains) previousPadGainsFile] "
    "-i rootFiles \n"
            << std::endl;
  exit(-1);
}

