/**
  \file Analysis script for analyzing Krypton decay data and
  calculating pad-by-pad gains. These gains are used during
  reconstruction, as the cluster position (calculated using the
  weighted mean) can be heavily influenced by large gain
  variations. The analysis contains options for peak detection or
  upper edge detection. In the case that the main Krypton decay peak
  (41.6 keV) is not visible, the edge detection option should be used.

  \author B. Rumberger
  \version $Id: 2   $
  \date 20 April 2022
*/

#include "KryptonAnalyzer.h"

#include <fwk/CentralConfig.h>
#include <det/Detector.h>
#include <det/TPC.h>
#include <det/TPCSector.h>
#include <modutils/DEDXTools.h>
#include <utl/ShineUnits.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TKey.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TPaveText.h>
#include <TMultiGraph.h>
#include <TTree.h>
#include <TStyle.h>

#include <cmath>
#include <fstream>
#include <string>

#include <boost/filesystem.hpp>

using namespace std;
using namespace modutils;

/// Main function.
int main(int argc, char* argv[])
{
  int exitCode = 0;

  const vector<string> argumentsVector(argv + 1, argv + argc);
  vector<string> filenamesVector;

  string configFilename = "Config.txt";
  string outputPrefix;
  string previousGainsFilename;
  bool updateGains = false;
  for (auto it = argumentsVector.begin(), itEnd = argumentsVector.end(); 
       it != itEnd; ++it) {
    if (*it == string("-h") || *it == string("--help")) {
      DisplayUsage();
    }
    else if (*it == string("-o")) {
      if (next(it) == itEnd || next(it)->rfind("-",0) != string::npos) {
	cout << "[ERROR] No output prefix provided with argument -o!" << endl;
	DisplayUsage();
      }
      advance(it,1);
      outputPrefix = *it + "-KryptonAnalysis";
    }
    else if (*it == string("-c") || *it == string("--config")) {
      if (next(it) == itEnd || next(it)->rfind("-",0) != string::npos) {
	cout << "[ERROR] No config filename provided with argument -c!" << endl;
	DisplayUsage();
      }
      advance(it,1);
      configFilename = *it;
      cout << "[INFO] User-provided config file: " << configFilename << endl;
    }
    else if (*it == string("-u") || *it == string("--updateGains")) {
      if (next(it) == itEnd || next(it)->rfind("-",0) != string::npos) {
	cout << "[ERROR] No pad gains XML path provided with argument -u!" << endl;
	DisplayUsage();
      }
      advance(it,1);
      previousGainsFilename = *it;
      updateGains = true;
      cout << "[INFO] User-provided gains file: " << previousGainsFilename << endl;
    }
    else if (*it == string("-i") || *it == string("--inputFiles")) {
      filenamesVector.assign(next(it),itEnd);
      break;
    }
    else {
      cout << "[ERROR] Invalid argument " << *it << "!" << endl;
      DisplayUsage();
    }
  }

  if (filenamesVector.size() == 0) {
    cout << "[ERROR] No input filenames provided!" << endl;
    DisplayUsage();
  }
  if (outputPrefix.size() == 0) {
    cout << "[ERROR] No output prefix provided!" << endl;
    DisplayUsage();
  }
  cout << "[INFO] Number of input files: " << filenamesVector.size()
       << ". Config file: " << configFilename 
       << ". Update previously-calculated gains? " << updateGains << endl;

  //Manage our own object ownership. Grrr.
  TH1::AddDirectory(kFALSE);

  //Parse configuration file.
  ParseConfigFile(configFilename);

  //Bootstrap XML path is by default in this directory.
  string bootstrapPath = "bootstrap.xml";

  if (updateGains) {
    ReplacePadGainPath(bootstrapPath,string(previousGainsFilename));
  }

  //Name and create output file. Use full path.
  const string& firstInputFile = filenamesVector.front();
  boost::filesystem::path currentPath( boost::filesystem::current_path() );
  const string& currentWorkingDirectory = currentPath.string() + "/";
  boost::filesystem::path inputNameAndPath(firstInputFile);

  //First histogram: No cuts. Second histogram: With cuts.
  unordered_map<int,unordered_map<int,pair<TH1D*,TH1D*> > > sectorSpectraHistograms;
  unordered_map<int,unordered_map<int,pair<TH2D*,TH2D*> > > sectorPadEntries;
  unordered_map<int,unordered_map<int,pair<TH1D*,TH1D*> > > sectorTimeSlices;
  unordered_map<int,unordered_map<int,pair<TH2D*,TH2D*> > > sectorChargeVsMaxADC;
  unordered_map<int,unordered_map<int,pair<TH2D*,TH2D*> > > sectorNPadsVsNTimeSlices;
    
  //Prepare PDF file.
  TString gainsPDFName = currentWorkingDirectory.c_str() + outputPrefix + ".pdf";
  TString gainsOpenString = gainsPDFName + "[";
  TString gainsCloseString = gainsPDFName + "]";
  TCanvas dummy;
  dummy.SaveAs(gainsOpenString);

  //Get parameters from XML file.
  fwk::CentralConfig::GetInstance(bootstrapPath);

  //Get detector and event interfaces.
  det::Detector& detector  = det::Detector::GetInstance();
  const unsigned int dummyRun = 1;
  const utl::TimeStamp dummyTime = utl::TimeStamp(1);
  detector.Update(dummyTime,dummyRun);
  const det::TPC& tpc = detector.GetTPC();

  //Create output file.
  TString outputFilename = currentWorkingDirectory + outputPrefix + ".root";
  cout << "[INFO] Output filename: " << outputFilename.Data() << endl;  
  TFile* outputFile = new TFile(outputFilename,"RECREATE");  


  //Create one histogram per active pad.
  for (auto chamberIt = tpc.ChambersBegin(), chamberEnd = tpc.ChambersEnd();
       chamberIt != chamberEnd; ++chamberIt) {
    const det::TPCChamber& chamber = *chamberIt;
    const unsigned int tpcId = (unsigned int)chamber.GetId();
    if (fTPCIdList.find(chamber.GetId()) == fTPCIdList.end())
      continue;
    for (auto sectorIt = chamber.SectorsBegin(), sectorEnd = chamber.SectorsEnd();
         sectorIt != sectorEnd; ++sectorIt) {
      const det::TPCSector& sector = *sectorIt;
      const unsigned int sectorId = sector.GetId();
      for (auto padrowIt = sector.PadrowsBegin(), padrowEnd = sector.PadrowsEnd();
           padrowIt != padrowEnd; ++padrowIt) {
        const det::TPCPadrow& padrow = *padrowIt;
        const unsigned int padrowId = padrow.GetId();
        for (unsigned int padId = 1; padId <= padrow.GetNPads(); ++padId) {
          //Create name and title.
          TString nameString =
            det::TPCConst::GetName(chamber.GetId()) +
            TString("Sector") + Form("%i",(unsigned int)sector.GetId()) +
            TString("Padrow") + Form("%i",(unsigned int)padrow.GetId()) +
            TString("Pad") + Form("%i",padId);
          TString titleString = "Krypton decay cluster charges, " + 
            det::TPCConst::GetName(chamber.GetId()) +
            TString(" Sector ") + Form("%i",(unsigned int)sector.GetId()) +
            TString(" Padrow ") + Form("%i",(unsigned int)padrow.GetId()) +
            TString(" Pad ") + Form("%i",padId) +
            TString(";Cluster Charge [ADC];Entries");

	  const double minADCPeakSearch = 
	    (tpcId == det::TPCConst::eVTPC1 && (sectorId == 1 || sectorId == 4)) ? 
	    fMinADCPeakSearchVTPC1Upstream : fMinADCPeakSearch;

          const double histogramMax = minADCPeakSearch*fHistogramPadding;
          TH1D* histogram = new TH1D(nameString,titleString,fHistogramBins,0,histogramMax);
          fSpectraHistograms[tpcId][sectorId][padrowId][padId] = histogram;
        } //End pad loop.
      } //End padrow loop.
    } //End sector loop.
  } //End TPC loop.
  
  //Variables to fill.
  Float16_t fCharge = 0;
  UShort_t fMaxADC = 0;
  UShort_t fTimeSlice = 0;
  UShort_t fNPixels = 0;
  UChar_t fNTimeSlices = 0;
  UChar_t fNPads = 0;
  UChar_t fPadrow = 0;
  UChar_t fPad = 0;

  //Loop over input files. Give progress percentage.
  double filesProcessed = 0;
  double previousPercentage = 0;
  for (auto fileIt = filenamesVector.begin(), fileEnd = filenamesVector.end();
       fileIt != fileEnd; ++fileIt) {
    ++filesProcessed;
    const double progressPercentage = round(1000*(filesProcessed - 1)/filenamesVector.size()/10.);
    if (previousPercentage != progressPercentage && fmod(progressPercentage,5) == 0)
      cout << "[INFO] Processing file " << filesProcessed
	   << " / " << filenamesVector.size() 
	   << " (" << progressPercentage << "% complete)." << endl;
    previousPercentage = progressPercentage;
    
    //Open the file for filling and get the TTree.
    const string& filename = *fileIt;
    TFile* inputFile = new TFile(filename.c_str(),"READ");
    if (inputFile == NULL || inputFile->IsZombie()) {
      cout << "[WARNING] Error opening input file! Skipping." << endl;
      continue;
    }

    if (inputFile->GetNkeys() == 0) {
      cout << "[WARNING] " << filename << " has no keys. Skipping." << endl;
      continue;
    }

    //Loop through keys in input file.
    TIter fileIter(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)fileIter())) {
      const TObject* object = key->ReadObj();
      if (object->InheritsFrom("TTree")) {
	TTree* tree = dynamic_cast<TTree*>(key->ReadObj());

	//Ignore empty trees.
	if (tree->GetEntries() == 0)
	  continue;
      
	//Identify TPC and sector.
	//Format: TTree name = [TPCName]Sector[SectorId]Clusters
	const string& treeName = tree->GetName();
	const string& tpcName = treeName.substr(0,treeName.find("Sector"));
	const size_t sectorIdStart = treeName.find("Sector") + 6;
	const size_t sectorIdStop = treeName.find("Clusters");

	const string& sectorIdString = treeName.substr(sectorIdStart,sectorIdStop - sectorIdStart);
	const det::TPCConst::EId tpcId = det::TPCConst::GetId(tpcName);
	const unsigned int sectorId = stoi(sectorIdString);
      
	const det::TPCSector& sector = tpc.GetChamber(tpcId).GetSector(sectorId);

	//Skip entries for TPCs we do not wish to calibrate.
	if (fTPCIdList.find(tpcId) == fTPCIdList.end())
	  continue;
	
      
	tree->SetBranchAddress("fCharge",&fCharge);
	tree->SetBranchAddress("fMaxADC",&fMaxADC);
	tree->SetBranchAddress("fTimeSlice",&fTimeSlice);
	tree->SetBranchAddress("fNPixels",&fNPixels);
	tree->SetBranchAddress("fNTimeSlices",&fNTimeSlices);
	tree->SetBranchAddress("fNPads",&fNPads);
	tree->SetBranchAddress("fPadrow",&fPadrow);
	tree->SetBranchAddress("fPad",&fPad);

	//Prepare QA plots.
	TString nameString = tpcName.data() +
	  TString("Sector") + Form("%i",(unsigned int)sector.GetId());
	TString titleString = tpcName.data() +
	  TString(" Sector ") + Form("%i",(unsigned int)sector.GetId());
      
	const double minADCPeakSearch = 
	  (tpcId == det::TPCConst::eVTPC1 && (sectorId == 1 || sectorId == 4)) ? 
	  fMinADCPeakSearchVTPC1Upstream : fMinADCPeakSearch;
      
	const double histogramMax = minADCPeakSearch*fHistogramPadding;
	const unsigned int nPadrows = sector.GetNPadrows();
	const unsigned int nPads = sector.GetPadrow(sector.GetNPadrows()).GetNPads();
      
	if (sectorSpectraHistograms.find(tpcId) == sectorSpectraHistograms.end() ||
	    sectorSpectraHistograms.at(tpcId).find(sectorId) == 
	    sectorSpectraHistograms.at(tpcId).end()) {

	  TH1D* spectrumNoCuts =
	    new TH1D(Form("ChargeNoCuts%s",nameString.Data()),
		     Form("%s Krypton Cluster Charges (No cuts);"
			  "Cluster Charge [ADC];Entries",
			  titleString.Data()),
		     2*fHistogramBins,0,histogramMax);
	  TH1D* spectrumAllCuts =
	    new TH1D(Form("ChargeAllCuts%s",nameString.Data()),
		     Form("%s Krypton Cluster Charges (All cuts Applied);"
			  "Cluster Charge [ADC];Entries",
			  titleString.Data()),
		     2*fHistogramBins,0,histogramMax);

	  TH2D* padEntriesNoCuts =
	    new TH2D(Form("padEntriesNoCuts%s",nameString.Data()),
		     Form("%s Entries Per Pad (No cuts);Pad Number;Padrow Number",
			  titleString.Data()),
		     nPads+2,0,nPads+2,
		     nPadrows+2,0,nPadrows+2);
	  TH2D* padEntriesAllCuts =
	    new TH2D(Form("padEntriesAllCuts%s",nameString.Data()),
		     Form("%s Entries Per Pad (All cuts Applied);Pad Number;Padrow Number",
			  titleString.Data()),
		     nPads+2,0,nPads+2,
		     nPadrows+2,0,nPadrows+2);

	  TH1D* timeSlicesNoCuts =
	    new TH1D(Form("timeSlicesNoCuts%s",nameString.Data()),
		     Form("%s Time Slices (No cuts);Time Slice;Entries",
			  titleString.Data()),
		     260,0,260);
	  TH1D* timeSlicesAllCuts =
	    new TH1D(Form("timeSlicesAllCuts%s",nameString.Data()),
		     Form("%s Time Slices (All cuts Applied);Time Slice;Entries",
			  titleString.Data()),
		     260,0,260);
      
	  TH2D* chargeVsMaxADCNoCuts =
	    new TH2D(Form("chargeVsMaxADCNoCuts%s",nameString.Data()),
		     Form("%s Charge vs. MaxADC (No cuts);Charge [ADC];MaxADC [ADC]",
			  titleString.Data()),
		     histogramMax*2,0,histogramMax*2,
		     512,0,512);
	  TH2D* chargeVsMaxADCAllCuts =
	    new TH2D(Form("chargeVsMaxADCAllCuts%s",nameString.Data()),
		     Form("%s Charge vs. MaxADC (All cuts Applied);Charge [ADC];MaxADC [ADC]",
			  titleString.Data()),
		     histogramMax*2,0,histogramMax*2,
		     512,0,512);

	  TH2D* nPadsVsNTimeSlicesNoCuts =
	    new TH2D(Form("nPadsVsNTimeSlicesNoCuts%s",nameString.Data()),
		     Form("%s nPads vs. nTimeSlices (No cuts);nPads;nTimeSlices",
			  titleString.Data()),
		     fMaxPads*5,0,fMaxPads*5,
		     fMaxTimeSlices*5,0,fMaxTimeSlices*5);
	  TH2D* nPadsVsNTimeSlicesAllCuts =
	    new TH2D(Form("nPadsVsNTimeSlicesAllCuts%s",nameString.Data()),
		     Form("%s nPads vs. nTimeSlices (All cuts Applied);nPads;nTimeSlices",
			  titleString.Data()),
		     fMaxPads*5,0,fMaxPads*5,
		     fMaxTimeSlices*5,0,fMaxTimeSlices*5);
      
	  sectorSpectraHistograms[tpcId][sectorId] = make_pair(spectrumNoCuts,spectrumAllCuts);
	  sectorPadEntries[tpcId][sectorId] = make_pair(padEntriesNoCuts,padEntriesAllCuts);
	  sectorTimeSlices[tpcId][sectorId] = make_pair(timeSlicesNoCuts,timeSlicesAllCuts);
	  sectorChargeVsMaxADC[tpcId][sectorId] = 
	    make_pair(chargeVsMaxADCNoCuts,chargeVsMaxADCAllCuts);
	  sectorNPadsVsNTimeSlices[tpcId][sectorId] = 
	    make_pair(nPadsVsNTimeSlicesNoCuts,nPadsVsNTimeSlicesAllCuts);
	}
      
	//Loop through all data to calculate total average for sectors.
	for(long double i = 0; i < tree->GetEntries(); ++i) {
	  tree->GetEntry(i);

	  const unsigned int pad = (unsigned int)fPad;
	  const unsigned int padrow = (unsigned int)fPadrow;
	  const unsigned int nPads = (unsigned int)fNPads;
	  const unsigned int nTimeSlices = (unsigned int)fNTimeSlices;
	
	  //Fill sector QA histograms (no cuts).
	  sectorSpectraHistograms[tpcId][sectorId].first->Fill(fCharge);
	  sectorPadEntries[tpcId][sectorId].first->Fill(pad,padrow);
	  sectorTimeSlices[tpcId][sectorId].first->Fill(fTimeSlice);
	  sectorChargeVsMaxADC[tpcId][sectorId].first->Fill(fCharge,fMaxADC);
	  sectorNPadsVsNTimeSlices[tpcId][sectorId].first->Fill(fNPads,fNTimeSlices);

	  //Ignore zero charge bins.
	  if (fCharge == 0)
	    continue;
	  //Cluster cuts.
	  if (nPads < fMinPads)
	    continue;
	  if (nPads > fMaxPads)
	    continue;
	  if (nTimeSlices < fMinTimeSlices)
	    continue;
	  if (nTimeSlices > fMaxTimeSlices)
	    continue;
	
	  if (fTimeSlice < fMinTimeSliceNumber)
	    continue;    
	  if (fCharge < fChargeCut && fMaxADC < fMaxADCCut)
	    continue;
	
	  if (updateGains) {
	    const det::TPCPadrow& detPadrow = sector.GetPadrow(padrow);
	    fCharge *= detPadrow.GetPadGain(pad);
	  }
	
	  //Fill pad histogram.
	  fSpectraHistograms[tpcId][sectorId][padrow][pad]->Fill(fCharge);

	  //Fill sector QA histograms (all cuts).
	  sectorSpectraHistograms[tpcId][sectorId].second->Fill(fCharge);
	  sectorPadEntries[tpcId][sectorId].second->Fill(pad,padrow);
	  sectorTimeSlices[tpcId][sectorId].second->Fill(fTimeSlice);
	  sectorChargeVsMaxADC[tpcId][sectorId].second->Fill(fCharge,fMaxADC);
	  sectorNPadsVsNTimeSlices[tpcId][sectorId].second->Fill(fNPads,fNTimeSlices);
	} //End TTree loop.

      } //End inherets from TTree.
    } //End key iteration.
    inputFile->Close();
  } //End filename loop.
  
  //Calculate peak positions.
  unordered_map<unsigned int, 
		unordered_map<unsigned int,
			      unordered_map<unsigned int,
					    unordered_map<unsigned int,double> > > > spectrumADCs;
  //Container for calculating and holding total sector averages.
  DEDXTools::SectorAveragers totalAccumulators;
  
  outputFile->cd();
  for (auto chamberIt = fSpectraHistograms.begin(), chamberEnd = fSpectraHistograms.end();
       chamberIt != chamberEnd; ++chamberIt) {
    const unsigned int tpcId = chamberIt->first;
    const SectorHistograms sectorHistograms = chamberIt->second;
    for (auto sectorIt = sectorHistograms.begin(), sectorEnd = sectorHistograms.end();
         sectorIt != sectorEnd; ++sectorIt) {
      const unsigned int sectorId = sectorIt->first;
      const PadrowHistograms padrowHistograms = sectorIt->second;
      for (auto padrowIt = padrowHistograms.begin(), padrowEnd = padrowHistograms.end();
           padrowIt != padrowEnd; ++padrowIt) {
        const unsigned int padrowId = padrowIt->first;
        const PadHistograms padHistograms = padrowIt->second;
        for (auto padIt = padHistograms.begin(), padEnd = padHistograms.end();
             padIt != padEnd; ++padIt) {
          const unsigned int padId = padIt->first;

          //Get histogram.
          TH1D* padHistogram = padIt->second;

          const int lastBin = padHistogram->GetXaxis()->GetNbins() - 1;
          const double maxCharge = padHistogram->GetXaxis()->GetBinCenter(lastBin);

          //Search for peak above minimum acceptable Krypton peak value.
          int maxBin = 0;
          double chargePeak = 0;
          double chargePeakValue = 0;
          for (int i = 0; i < lastBin; ++i) {
            const double binCenter = padHistogram->GetXaxis()->GetBinCenter(i);
	    const double minADCPeakSearch = 
	      ((det::TPCConst::EId)tpcId == 
	       det::TPCConst::eVTPC1 && (sectorId == 1 || sectorId == 4)) ? 
	      fMinADCPeakSearchVTPC1Upstream : fMinADCPeakSearch;
	    if (binCenter < minADCPeakSearch)
              continue;
            const double value = padHistogram->GetBinContent(i);
            if (value > chargePeakValue) {
              maxBin = i;
              chargePeak = binCenter;
              chargePeakValue = value;
            }
          }

          //Find where peak drops by a factor of 2 above and below.
          //Calculate peak nearest end of distribution. Count backwards from end of histogram.
          double minChargeForFit = 0;
          double maxChargeForFit = 0;
          for (int bin = maxBin; bin > 0; --bin) {
            const double binContent = padHistogram->GetBinContent(bin);
            if (binContent < 0.5*chargePeakValue) {
              minChargeForFit = padHistogram->GetXaxis()->GetBinCenter(bin);
              break;
            }
          }
          for (int bin = maxBin; bin <= lastBin; ++bin) {
            const double binContent = padHistogram->GetBinContent(bin);
            if (binContent < 0.5*chargePeakValue) {
              maxChargeForFit = padHistogram->GetXaxis()->GetBinCenter(bin);
              break;
            }
          }

          //Define fit functions.
          TF1* gausFit = new TF1("gausFit","gaus",
                                 minChargeForFit,maxChargeForFit);
          TF1* fermiFit = new TF1("fermiFit","[0]/(1+TMath::Exp([1]*(x-[2])))",
                                  chargePeak,maxCharge);
          fermiFit->FixParameter(0,chargePeakValue);
          fermiFit->SetParameter(1,0.01);
          fermiFit->SetParLimits(1,0.0001,1);
          fermiFit->SetParameter(2,chargePeak);
          
          //Don't do anything for pads with too few entries.
          if (padHistogram->GetEntries() >= fMinHistogramEntries) {
            //Perform desired fit. Store results.
            if (fFitFunction == "Gaussian") {
              padHistogram->Fit(gausFit,"R Q");
              spectrumADCs[tpcId][sectorId][padrowId][padId] = gausFit->GetParameter(1);
              totalAccumulators.AddValue(tpcId,sectorId,gausFit->GetParameter(1));
            }
            else if (fFitFunction == "Fermi") {
              padHistogram->Fit(fermiFit,"R Q");
              spectrumADCs[tpcId][sectorId][padrowId][padId] = fermiFit->GetParameter(2);
              totalAccumulators.AddValue(tpcId,sectorId,fermiFit->GetParameter(2));
            }
          }
          //Write to QA file.
	  if (padHistogram->GetEntries() > 0)
	    padHistogram->Write();
        } // Pad loop.
      } // Padrow loop.
    } // Sector loop.
  } // TPC loop.

  if (fSpectraHistograms.begin() == fSpectraHistograms.end())
    cout << "[WARNING] No histograms were filled. "
         << "Was your TPC included in the configuration file list?" << endl;
  
  //TTree for storing results.
  TTree* fResultTree = new TTree("fResultTree","Krypton Analysis Results");
  unsigned int fTPCId;
  fResultTree->Branch("fTPCId",&fTPCId);
  unsigned int fSectorId;
  fResultTree->Branch("fSectorId",&fSectorId);
  unsigned int fPadrowId;
  fResultTree->Branch("fPadrowId",&fPadrowId);
  unsigned int fPadId;
  fResultTree->Branch("fPadId",&fPadId);
  double fSpectrumADC;
  fResultTree->Branch("fSpectrumADC",&fSpectrumADC);
  double fGain;
  fResultTree->Branch("fGain",&fGain);

  //XML file writing infrastructure.
  string gainsFilename = currentWorkingDirectory + outputPrefix + "-KryptonPadGains.xml";
  ofstream gainsFileStream;
  gainsFileStream.open(gainsFilename);
  
  gainsFileStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
                  << "\n"
                  << "<PadByPadGain\n"
                  << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                  << "  xsi:noNamespaceSchemaLocation=\"[SCHEMAPATH]/TPCPadGain_DataFormat.xsd\">\n"
                  << "\n";
  
  //Calculate gains. Normalize spectrum ADC to average sector ADCs.
  for (auto chamberIt = tpc.ChambersBegin(), chamberEnd = tpc.ChambersEnd();
       chamberIt != chamberEnd; ++chamberIt) {
    const det::TPCChamber& chamber = *chamberIt;
    const unsigned int tpcId = (unsigned int)chamber.GetId();
    gainsFileStream << "  <TPC name=\"" << det::TPCConst::GetName(chamber.GetId())<< "\">\n";
    for (auto sectorIt = chamber.SectorsBegin(), sectorEnd = chamber.SectorsEnd();
         sectorIt != sectorEnd; ++sectorIt) {
      const det::TPCSector& sector = *sectorIt;
      const unsigned int sectorId = (unsigned int)sector.GetId();
      const double sectorADC = totalAccumulators.GetAverage(tpcId,sectorId);
      gainsFileStream << "    <Sector id=\"" << (unsigned int)sectorId << "\">\n";
      for (auto padrowIt = sector.PadrowsBegin(), padrowEnd = sector.PadrowsEnd();
           padrowIt != padrowEnd; ++padrowIt) {
        const det::TPCPadrow& padrow = *padrowIt;
        const unsigned int padrowId = padrow.GetId();
        gainsFileStream << "      <Padrow id=\"" << (unsigned int)padrowId << "\">\n";
        gainsFileStream << "        <PadGains> ";
        for (unsigned int padId = 1; padId <= padrow.GetNPads(); ++padId) {
          
          const double padADC = spectrumADCs[tpcId][sectorId][padrowId][padId];
          const double gain = (updateGains) ?
            padrow.GetPadGain(padId)*sectorADC/padADC : sectorADC/padADC;

          if (fTPCIdList.find(chamber.GetId()) == fTPCIdList.end())
            continue;
          
          if (updateGains) {
            cout << "TPC " << tpcId
                 << ", sector " << sectorId
                 << ", padrow " << padrowId
                 << ", pad " << padId
                 << ": Previous gain = " << padrow.GetPadGain(padId)
                 << ". Pad ADC = " << padADC
                 << ". Sector ADC = " << sectorADC
                 << ". SectorADC/PadADC = " << sectorADC/padADC
                 << ". New gain = " << gain 
                 << ". updateGains? " << updateGains << endl;
            
          }
          
          //Record information in output file.
          fTPCId = tpcId;
          fSectorId = sectorId;
          fPadrowId = padrowId;
          fPadId = padId;
          fSpectrumADC = padADC;
          fGain = (isnan(gain) || isinf(gain) ) ? 
	    0 : gain;
          fResultTree->Fill();
          
          if (gain > fMinAcceptableGain &&
              gain < fMaxAcceptableGain)
            gainsFileStream << gain << " ";
          else
            gainsFileStream << -1.0 << " ";
        } // Pad loop.
        gainsFileStream << "</PadGains>\n";
        gainsFileStream << "      </Padrow>\n";
      } // Padrow loop.
      gainsFileStream << "    </Sector>\n";
    } // Sector loop.
    gainsFileStream << "  </TPC>\n";
  } // TPC loop.
  gainsFileStream << "</PadByPadGain>\n";
  
  cout << "[INFO] Pad gains written to file " << gainsFilename << " . Thanks!" << endl;

  //Make QA plots.
  for (auto chamberIt = tpc.ChambersBegin(), chamberEnd = tpc.ChambersEnd();
       chamberIt != chamberEnd; ++chamberIt) {
    const det::TPCChamber& chamber = *chamberIt;
    const unsigned int tpcId = (unsigned int)chamber.GetId();
    if (fTPCIdList.find(chamber.GetId()) == fTPCIdList.end())
      continue;
    for (auto sectorIt = chamber.SectorsBegin(), sectorEnd = chamber.SectorsEnd();
         sectorIt != sectorEnd; ++sectorIt) {
      const det::TPCSector& sector = *sectorIt;
      unsigned int maxPadsPerPadrow = 0;
      for (auto padrowIt = sector.PadrowsBegin(), padrowEnd = sector.PadrowsEnd();
           padrowIt != padrowEnd; ++padrowIt) {
        const det::TPCPadrow& padrow = *padrowIt;
        if (maxPadsPerPadrow < padrow.GetNPads())
          maxPadsPerPadrow = padrow.GetNPads();
      }
      const unsigned int sectorId = (unsigned int)sector.GetId();
      const double sectorADC = totalAccumulators.GetAverage(tpcId,sectorId);

      //Make sector histogram.
      TString nameString =
        TString("tpc") + Form("%i",(unsigned int)chamber.GetId()) +
        TString("Sector") + Form("%i",(unsigned int)sector.GetId());
      TString titleString = "Pad Gains, " + 
        det::TPCConst::GetName(chamber.GetId()) +
        TString(" Sector ") + Form("%i",(unsigned int)sector.GetId()) +
        TString(";Pad;Padrow");
      TH2D sectorGains(nameString,titleString,maxPadsPerPadrow+1,0,maxPadsPerPadrow+1,
                       sector.GetNPadrows()+1,0,sector.GetNPadrows()+1);
      for (auto padrowIt = sector.PadrowsBegin(), padrowEnd = sector.PadrowsEnd();
           padrowIt != padrowEnd; ++padrowIt) {
        const det::TPCPadrow& padrow = *padrowIt;
        const unsigned int padrowId = padrow.GetId();
        for (unsigned int padId = 1; padId <= padrow.GetNPads(); ++padId) {
          
          const double padADC = spectrumADCs[tpcId][sectorId][padrowId][padId];
          const double gain = sectorADC/padADC;

          sectorGains.SetBinContent(padId,padrowId,gain);
          
        } // Pad loop.
      } // Padrow loop.



      //Make a canvas for holding all plots.
      TCanvas canvas;
      gStyle->SetOptStat(0);
      sectorGains.SetMinimum(0.6);
      sectorGains.SetMaximum(1.4);
      sectorGains.Draw("COLZ");
      canvas.SaveAs(gainsPDFName);
      outputFile->cd();
      sectorGains.Write();
    } // Sector loop.
  } // TPC loop.

  ///Save QA PDFs.
  TGaxis::SetMaxDigits(3);

  const double axisTitleOffset = 1.7;
  const double bottomMargin = 0.15;
  const double leftMargin = 0.14;
  const double rightMargin = 0.14;

  gStyle->SetOptStat(0);
  
  for (auto it = sectorSpectraHistograms.begin(),
  	 itEnd = sectorSpectraHistograms.end(); it != itEnd; ++it) {
    const det::TPCConst::EId tpcId = (det::TPCConst::EId)it->first;
    auto sectorMap = it->second;
    for (auto sectorIt = sectorMap.begin(), sectorEnd = sectorMap.end();
	 sectorIt != sectorEnd; ++sectorIt) {
      const int sectorId = sectorIt->first;
      auto histogramPair = sectorIt->second;


      TCanvas canvas;
      canvas.Divide(2,1);
      canvas.cd(1);
      histogramPair.first->Draw();
      histogramPair.first->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.cd(2);
      TH1D* spectrum = histogramPair.second;
      spectrum->Draw();
      //Fit around peak.
      const double minADCPeakSearch = 
	(tpcId == det::TPCConst::eVTPC1 && (sectorId == 1 || sectorId == 4)) ? 
	fMinADCPeakSearchVTPC1Upstream : fMinADCPeakSearch;
      
      int maxBin = spectrum->FindFixBin(minADCPeakSearch);
      double max = spectrum->GetBinContent(maxBin);
      for (int bin = maxBin; bin < spectrum->GetNbinsX(); ++bin) {
	if (spectrum->GetBinContent(bin) > max) {
	  maxBin = bin;
	  max = spectrum->GetBinContent(bin);
	}
      }
      int halfWidthHalfMaxBin = maxBin;
      for (int i = 0; i < 50; ++i) {
	if (spectrum->GetBinContent(maxBin+ i) < 0.7*max) {
	  halfWidthHalfMaxBin = i;
	  break;
	}
      }
      const double fitMin = spectrum->GetXaxis()->GetBinCenter(maxBin - halfWidthHalfMaxBin);
      const double fitMax = spectrum->GetXaxis()->GetBinCenter(maxBin + halfWidthHalfMaxBin);
      TFitResultPtr r = spectrum->Fit("gaus", "QSIR", "",fitMin,fitMax);
      int fitStatus = r;
      double mean = 0;
      double sigma = 0;
      if (fitStatus >= 0) {
	mean = r->Parameter(1);
	sigma = r->Parameter(2);
      }
      TLatex latexX;
      latexX.SetTextSize(0.035);
      latexX.DrawLatexNDC(.6,.82,Form("#mu = %1.4f",mean));
      latexX.DrawLatexNDC(.6,.8,Form("#sigma = %1.3f",sigma));
    
      spectrum->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.SaveAs(gainsPDFName);
    }
  }

  for (auto it = sectorSpectraHistograms.begin(),
  	 itEnd = sectorSpectraHistograms.end(); it != itEnd; ++it) {
    const det::TPCConst::EId tpcId = (det::TPCConst::EId)it->first;
    const string& tpcName = det::TPCConst::GetName(tpcId);
    auto sectorMap= it->second;
    for (auto sectorIt = sectorMap.begin(), sectorEnd = sectorMap.end();
	 sectorIt != sectorEnd; ++sectorIt) {
      const unsigned int sectorId = sectorIt->first;
      const det::TPCSector& sector = tpc.GetChamber(tpcId).GetSector(sectorId);
      const unsigned int nPadrows = sector.GetNPadrows();
      
      TH1D* gains = new TH1D(Form("%sSector%iGains",tpcName.data(),sectorId),
			     Form("%s Sector %i Gains;Gain;Entries",
				  tpcName.data(),sectorId),
			     200,0.5,1.5);
      
      fResultTree->SetBranchAddress("fTPCId",&fTPCId);
      fResultTree->SetBranchAddress("fSectorId",&fSectorId);
      fResultTree->SetBranchAddress("fPadrowId",&fPadrowId);
      fResultTree->SetBranchAddress("fPadId",&fPadId);
      fResultTree->SetBranchAddress("fGain",&fGain);

      //Create padrow color pallette.
      const map<int,int> colorsByPadrow = GetPadrowColorMap(nPadrows);

      //Container for holding gains indexed by pad and padrow.
      //Indices: gainsContainer[padrowId][padId] = gain
      map<int,map<int,double> > gainsContainer;
      
      for (double i = 0; i < fResultTree->GetEntries(); ++i) {
	fResultTree->GetEntry(i);
	if (fTPCId == (unsigned int)tpcId && fSectorId == sectorId) {
	  gains->Fill(fGain);
	  gainsContainer[fPadrowId][fPadId] = fGain;
	}
      }
      //Create TGraphs and TMultiGraph.
      TMultiGraph multigraph;
      TString gainsByPadName = Form("%sSector%iGainsByPad",tpcName.data(),sectorId);
      multigraph.SetNameTitle(gainsByPadName,
			      Form("%s Sector %i Gains Vs. Pad;Pad Id;Gain;Padrow Id",
				   tpcName.data(),sectorId));
      for (unsigned int padrowId = 1; padrowId <= nPadrows; ++padrowId) {
	vector<double> padIds;
	vector<double> gains;
	for (auto it = gainsContainer[padrowId].begin(), itEnd = gainsContainer[padrowId].end();
	     it != itEnd; ++it) {
	  const int padId = it->first;
	  const double gain = it->second;
	  padIds.push_back(padId);
	  gains.push_back(gain);
	}
	TGraph* gainGraph = new TGraph(padIds.size(),padIds.data(),gains.data());
	gainGraph->SetMarkerColor(colorsByPadrow.at(padrowId));
	gainGraph->SetMarkerStyle(8);
	gainGraph->SetMarkerSize(0.6);
	multigraph.Add(gainGraph);
      }
      
      //Create z-scale palette using dummy TH2D.
      TH2D* dummy2D = new TH2D("dummy","dummy",100,0,1,100,0,1);
      dummy2D->Fill(0.1,0.1,1);
      dummy2D->Fill(0.9,0.9,nPadrows);
      dummy2D->GetZaxis()->SetLabelSize(0.02);
      TCanvas dummy;
      dummy2D->Draw("COLZ");
      dummy.Update();
      TPaletteAxis* palette =
      	(TPaletteAxis*)dummy2D->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.9);
      palette->SetX2NDC(0.925);
      palette->SetY1NDC(0.1);
      palette->SetY2NDC(0.9);
      TLatex label;
      label.SetTextSize(0.035);
      label.SetTextAngle(90);
      
      TCanvas canvas;
      canvas.cd();
      multigraph.SetMinimum(0.5);
      multigraph.SetMaximum(1.5);
      multigraph.Draw("AP");
      palette->Draw();
      label.DrawLatexNDC(0.975,0.45,"Padrow Id");
      canvas.SaveAs(gainsPDFName);
      gains->Draw();
      canvas.SaveAs(gainsPDFName);
      // gStyle->SetPalette(55);
    }
  }  

  for (auto it = sectorPadEntries.begin(),
  	 itEnd = sectorPadEntries.end(); it != itEnd; ++it) {
    auto sectorMap = it->second;
    for (auto sectorIt = sectorMap.begin(), sectorEnd = sectorMap.end();
	 sectorIt != sectorEnd; ++sectorIt) {
      auto histogramPair = sectorIt->second;
      
      TCanvas canvas;
      canvas.Divide(2,1);
      canvas.cd(1);
      histogramPair.first->Draw("COLZ");
      histogramPair.first->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.cd(2);
      histogramPair.second->Draw("COLZ");
      histogramPair.second->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.SaveAs(gainsPDFName);
    }
  }
  
  for (auto it = sectorTimeSlices.begin(),
  	 itEnd = sectorTimeSlices.end(); it != itEnd; ++it) {
    auto sectorMap = it->second;
    for (auto sectorIt = sectorMap.begin(), sectorEnd = sectorMap.end();
	 sectorIt != sectorEnd; ++sectorIt) {
      auto histogramPair = sectorIt->second;
      TCanvas canvas;
      canvas.Divide(2,1);
      canvas.cd(1);
      histogramPair.first->Draw();
      histogramPair.first->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.cd(2);
      histogramPair.second->Draw();
      histogramPair.second->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.SaveAs(gainsPDFName);
    }
  }

  for (auto it = sectorChargeVsMaxADC.begin(),
  	 itEnd = sectorChargeVsMaxADC.end(); it != itEnd; ++it) {
    auto sectorMap = it->second;
    for (auto sectorIt = sectorMap.begin(), sectorEnd = sectorMap.end();
	 sectorIt != sectorEnd; ++sectorIt) {
      auto histogramPair = sectorIt->second;
      TCanvas canvas;
      canvas.Divide(2,1);
      canvas.cd(1);
      histogramPair.first->Draw("COLZ");
      histogramPair.first->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.cd(2);
      histogramPair.second->Draw("COLZ");
      histogramPair.second->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.SaveAs(gainsPDFName);
    }
  }

  for (auto it = sectorNPadsVsNTimeSlices.begin(),
  	 itEnd = sectorNPadsVsNTimeSlices.end(); it != itEnd; ++it) {
    auto sectorMap = it->second;
    for (auto sectorIt = sectorMap.begin(), sectorEnd = sectorMap.end();
	 sectorIt != sectorEnd; ++sectorIt) {
      auto histogramPair = sectorIt->second;
      TCanvas canvas;
      canvas.Divide(2,1);
      canvas.cd(1);
      histogramPair.first->Draw("COLZ");
      histogramPair.first->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.cd(2);
      histogramPair.second->Draw("COLZ");
      histogramPair.second->GetYaxis()->SetTitleOffset(axisTitleOffset);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      gPad->SetLogz();
      canvas.SaveAs(gainsPDFName);
    }
  } 
  
  //Close PDF.
  dummy.SaveAs(gainsCloseString);
  
  //Clean up and finish.
  gainsFileStream.close();
  fResultTree->Write();
  outputFile->Close();
  
  return exitCode; 
}

void ParseConfigFile(const std::string& configFile) {
  //Open file.  
  ifstream file(configFile);
  //Parse lines in file.
  std::string line;
  cout << "[INFO] Parsing config file:" << endl;
  while (std::getline(file, line))
    {
      std::istringstream lineString(line);
      //Ignore lines beginning with a "#".
      if (lineString.str().front() == '#')
        continue;
      //Parse everything else.
      string variableName = "";
      int foundPosition = 0;
      if (lineString.str().find("tpcList",foundPosition) != string::npos) {
        string tpcName = "";
        while (line.find("tpcListEnd") == string::npos) {
          std::getline(file,line);
          std::istringstream nameString(line);
          //Ignore lines beginning with a "#".
          if (nameString.str().front() == '#')
            continue;
          if (!(nameString >> tpcName)) {
            cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
          }
          //Translate to TPC Id.
          const det::TPCConst::EId& tpcId = det::TPCConst::GetId(tpcName);
          if (tpcId != det::TPCConst::eUnknown) {
            fTPCIdList.insert(tpcId);
            cout << "[INFO] Added TPC " << tpcName << " (ID = " << tpcId << ")" << endl;
          }
        }
      }
      if (lineString.str().find("fitFunction",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fFitFunction)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Fit function: " << fFitFunction << endl;
      }
      else if (lineString.str().find("minAcceptableGain",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinAcceptableGain)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] minAcceptableGain: " << fMinAcceptableGain << endl;
      }
      else if (lineString.str().find("maxAcceptableGain",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxAcceptableGain)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] maxAcceptableGain: " << fMaxAcceptableGain << endl;
      }
      else if (lineString.str().find("minADCPeakSearch",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinADCPeakSearch)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] minADCPeakSearch: " << fMinADCPeakSearch << endl;
      }
      else if (lineString.str().find("vtpc1UpstreamSectorsMinADCPeakSearch",
				     foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinADCPeakSearchVTPC1Upstream)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] MinADCPeakSearchVTPC1Upstream: " << fMinADCPeakSearchVTPC1Upstream << endl;
      }
      else if (lineString.str().find("minHistogramEntries",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinHistogramEntries)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] minHistogramEntries: " << fMinHistogramEntries << endl;
      }
      else if (lineString.str().find("histogramBins",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fHistogramBins)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] histogramBins: " << fHistogramBins << endl;
      }
      else if (lineString.str().find("histogramPadding",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fHistogramPadding)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] histogramPadding: " << fHistogramPadding << endl;
      }
      else if (lineString.str().find("minPads",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinPads)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] minPads: " << fMinPads << endl;
      }
      else if (lineString.str().find("maxPads",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxPads)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] maxPads: " << fMaxPads << endl;
      }
      else if (lineString.str().find("minTimeSliceNumber",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinTimeSliceNumber)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] minTimeSliceNumber: " << fMinTimeSliceNumber << endl;
      }
      else if (lineString.str().find("minTimeSlices",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinTimeSlices)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] minTimeSlices: " << fMinTimeSlices << endl;
      }
      else if (lineString.str().find("maxTimeSlices",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxTimeSlices)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] maxTimeSlices: " << fMaxTimeSlices << endl;
      }
      else if (lineString.str().find("maxADCCut",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxADCCut)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] maxADCCut: " << fMaxADCCut << endl;
      }
      else if (lineString.str().find("chargeCut",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fChargeCut)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] chargeCut: " << fChargeCut << endl;
      }
      
    } //End parsing.
  return;  
}

void ReplacePadGainPath(const std::string& bootstrap,
                        const std::string& newPadGainXML) 
{
  //Definitions.
  string defaultPadGainXML = "&configDir;/TPCPadGainFixedManager.xml";
  string tempBootstrapName = "bootstrap-temp.xml";
  //Open file.  
  ifstream file(bootstrap);
  //Create output file.
  ofstream outputFile(tempBootstrapName);
  //Parse lines in file.
  std::string line;
  cout << "[INFO] Replacing default pad gain path with new path: " << newPadGainXML << endl;
  while (std::getline(file, line)) {
    //Search for path to default TPC pad gain XML.
    size_t foundPosition = line.find(defaultPadGainXML);
    if (foundPosition != string::npos) {
      line.replace(foundPosition,
                   defaultPadGainXML.size(),
                   newPadGainXML);
    }
    //Copy string to temp file.
    outputFile << line << endl;
  }
  //Move the tmp file to original bootstrap location.
  namespace fs = boost::filesystem;
  try {
    fs::rename(tempBootstrapName,bootstrap);
  }
  catch (fs::filesystem_error& e) {
    cout << e.what() << '\n';
  }
  return;
}

map<int,int> GetPadrowColorMap(const int maxPadrows) 
{
  const int nColors = maxPadrows;
  map<int,int> colorsByPadrow;
    
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};

  // Double_t red[3]    = { 1.00, 0.00, 0.00};
  // Double_t green[3]  = { 0.00, 1.00, 0.00};
  // Double_t blue[3]   = { 1.00, 0.00, 1.00};
  // Double_t length[3] = { 0.00, 0.50, 1.00 };
 
  int color = TColor::CreateGradientColorTable(9,stops,red,green,blue,nColors);
  for (int i = 0, padrow = 1; i < nColors; ++i, ++padrow) {
    int currentColor = color + i;
    colorsByPadrow[padrow] = currentColor;
  }
  
  return colorsByPadrow;
}
