/Michel Electron Analysis Program. It  accepts multiple files. 
//* Purpose:
//* This code is designed to process PMT (Photomultiplier Tube) data to identify and analyze
//* Michel electrons resulting from muon decays. It performs three main functions:
//* 1. Calibrates PMT responses to determine single photoelectron (SPE) spectra
//* 2. Applies event selection cuts to identify good candidate events
//* 3. Analyzes muon decay time distributions and Michel electron energy spectra
//* Key Features:
//* - Automatic SPE calibration using a multi-peak fit function
//* - Advanced event selection with multiple criteria:
//*   • Pulse height requirements
//*   • Peak position consistency checks
//*   • Signal-to-noise ratio cuts
//* - Muon lifetime measurement via exponential fitting
//* - Michel electron energy spectrum analysis
//* - Automatic output organization with timestamped directories
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>
#include <ctime>
#include <cstdlib>

using namespace std;

// SPE fitting function
Double_t SPEfit(Double_t *x, Double_t *par) {
    Double_t term1 = par[0] * exp(-0.5 * pow((x[0]-par[1])/par[2], 2));
    Double_t term2 = par[3] * exp(-0.5 * pow((x[0]-par[4])/par[5], 2));
    Double_t term3 = par[6] * exp(-0.5 * pow((x[0]-sqrt(2)*par[4])/sqrt(2*pow(par[5],2)-pow(par[2],2)), 2));
    Double_t term4 = par[7] * exp(-0.5 * pow((x[0]-sqrt(3)*par[4])/sqrt(3*pow(par[5],2)-2*pow(par[2],2)), 2));
    return term1 + term2 + term3 + term4;
}

// Muon decay exponential function
Double_t DecayFit(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0]/par[1]);
}

// Function to calculate mean and RMS of a dataset
void CalculateMeanAndRMS(const vector<Double_t> &data, Double_t &mean, Double_t &rms) {
    mean = 0.0;
    for (const auto &value : data) mean += value;
    mean /= data.size();
    
    rms = 0.0;
    for (const auto &value : data) rms += pow(value - mean, 2);
    rms = sqrt(rms / data.size());
}

// Function to generate a random directory name
string generateRandomDirName(const string &baseName) {
    time_t now = time(0);
    tm *ltm = localtime(&now);
    int randomNum = rand() % 10000;
    
    string dirName = baseName + "_Output_" + to_string(1900 + ltm->tm_year) + 
                   to_string(1 + ltm->tm_mon) + to_string(ltm->tm_mday) + "_" +
                   to_string(ltm->tm_hour) + to_string(ltm->tm_min) + to_string(ltm->tm_sec) +
                   "_" + to_string(randomNum);
    return dirName;
}

// Function to analyze muon/Michel events
void analyzeMuonMichel(TTree *tree, const Double_t *mu1, const string &outputDir) {
    gErrorIgnoreLevel = kError;

    // Variables to read from the tree
    Double_t area[23], pulseH[23];
    Int_t triggerBits;
    Long64_t nsTime;

    // Set branch addresses
    tree->SetBranchAddress("area", area);
    tree->SetBranchAddress("pulseH", pulseH);
    tree->SetBranchAddress("triggerBits", &triggerBits);
    tree->SetBranchAddress("nsTime", &nsTime);

    const int nPMTs = 12;
    int pmtChannelMap[nPMTs] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

    // Create analysis subdirectory
    string analysisDir = outputDir + "/AnalysisResults";
    mkdir(analysisDir.c_str(), 0777);

    // Histograms
    TH1F *histDeltaT = new TH1F("DeltaT", "Time difference;(Time difference #mus);Events/0.1 #mus", 45, 1, 10);
    TH1F *histMichelSpectrum = new TH1F("MichelSpectrum", "Michel Electron Spectrum;Photoelectrons;Events", 90, 100, 1000);

    int muonCount = 0, michelCount = 0;
    const double muonThreshold = 50;
    const double michelThreshold = 50;

    Long64_t nEntries = tree->GetEntries();
    cout << "Analyzing " << nEntries << " events for muon/Michel decays..." << endl;

    for(Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);
        
        // Check trigger condition for muon candidates
        if(triggerBits != 34 && triggerBits != 2) continue;

        // Check PMT hits for muon
        int nPMTsHit = 0;
        for(int pmt=0; pmt<nPMTs; pmt++) {
            if(area[pmtChannelMap[pmt]] > 500) nPMTsHit++;
        }
        if(nPMTsHit < 3) continue;

        // Convert area to p.e. and sum total energy
        double totalEnergy = 0;
        for(int pmt=0; pmt<nPMTs; pmt++) {
            if(area[pmtChannelMap[pmt]] > 500 && mu1[pmt] > 0) {
                totalEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
            }
        }

        if (totalEnergy >= muonThreshold) {
            muonCount++;
            Long64_t muonTime = nsTime;
            
            // Search for Michel electrons
            for(Long64_t nextEntry=entry+1; nextEntry<nEntries; nextEntry++) {
                tree->GetEntry(nextEntry);
                double deltaT = (nsTime - muonTime)*1e-3;
                
                if(deltaT > 10) break;
                if(deltaT < 1) continue;
                
                int nPMTsHitMichel = 0;
                for(int pmt=0; pmt<nPMTs; pmt++) {
                    if(area[pmtChannelMap[pmt]] > 500) nPMTsHitMichel++;
                }
                if(nPMTsHitMichel < 3) continue;
                
                double michelEnergy = 0;
                for(int pmt=0; pmt<nPMTs; pmt++) {
                    if(area[pmtChannelMap[pmt]] > 500 && mu1[pmt] > 0) {
                        michelEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
                    }
                }
                
                if(michelEnergy >= michelThreshold) {
                    michelCount++;
                    histDeltaT->Fill(deltaT);
                    histMichelSpectrum->Fill(michelEnergy);
                    break;
                }
            }
        }
    }

    cout << "Found " << muonCount << " muon candidates and " 
         << michelCount << " Michel electrons." << endl;

    // Time difference fit and plot
    TCanvas *c1 = new TCanvas("c1", "Time Difference", 800, 600);
    TF1 *decayFit = new TF1("decayFit", DecayFit, 1, 10, 2);
    decayFit->SetParameters(histDeltaT->GetMaximum(), 2.2);
    decayFit->SetParLimits(0, 0, histDeltaT->GetMaximum() * 2);
    decayFit->SetParLimits(1, 0.1, 10);
    
    histDeltaT->Fit(decayFit, "L", "", 1, 10);
    decayFit->SetRange(2, 10);
    
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.85, 0.8, "NDC");
    pt->AddText(Form("#tau = %.2f #pm %.2f #mus", 
                    decayFit->GetParameter(1), decayFit->GetParError(1)));
    histDeltaT->GetListOfFunctions()->Add(pt);
    
    histDeltaT->Draw();
    decayFit->Draw("same");
    c1->SaveAs((analysisDir + "/time_difference.png").c_str());

    // Michel spectrum plot
    TCanvas *c2 = new TCanvas("c2", "Michel Electron Spectrum", 800, 600);
    histMichelSpectrum->Draw();
    c2->SaveAs((analysisDir + "/michel_spectrum.png").c_str());

    // Save histograms to file
    TFile *analysisResultsFile = new TFile((analysisDir + "/analysis_results.root").c_str(), "RECREATE");
    histDeltaT->Write();
    histMichelSpectrum->Write();
    analysisResultsFile->Close();

    // Cleanup
    delete histDeltaT;
    delete histMichelSpectrum;
    delete c1;
    delete c2;
}

// Main function to process events
void processEvents(const char *fileName) {
    string inputFileName(fileName);
    size_t lastDot = inputFileName.find_last_of(".");
    string baseName = inputFileName.substr(0, lastDot);
    string outputDir = generateRandomDirName(baseName);
    cout << "Creating output directory: " << outputDir << endl;

    if (mkdir(outputDir.c_str(), 0777)) {
        cerr << "Error creating directory: " << outputDir << endl;
        return;
    }

    // Open input file
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }
    cout << "Opened input file: " << fileName << endl;

    // Access TTree
    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree!" << endl;
        file->Close();
        return;
    }

    // Variables to read from the tree
    Int_t eventID;
    Int_t nSamples[23];
    Short_t adcVal[23][45];
    Double_t baselineMean[23], baselineRMS[23], pulseH[23], area[23];
    Int_t peakPosition[23];
    Long64_t nsTime;
    Int_t triggerBits;

    // Set branch addresses
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("nSamples", nSamples);
    tree->SetBranchAddress("adcVal", adcVal);
    tree->SetBranchAddress("baselineMean", baselineMean);
    tree->SetBranchAddress("baselineRMS", baselineRMS);
    tree->SetBranchAddress("pulseH", pulseH);
    tree->SetBranchAddress("area", area);
    tree->SetBranchAddress("peakPosition", peakPosition);
    tree->SetBranchAddress("nsTime", &nsTime);
    tree->SetBranchAddress("triggerBits", &triggerBits);

    // 1. CALIBRATION PHASE
    TH1F *histArea[12];
    int pmtChannelMap[12] = {0,10,7,2,6,3,8,9,11,4,5,1};
    
    for (int i=0; i<12; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area",i+1), 
                             Form("PMT %d;ADC Counts;Events",i+1), 150, -50, 400);
    }

    Long64_t nEntries = tree->GetEntries();
    cout << "Processing " << nEntries << " events for calibration..." << endl;
    for (Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);
        if (triggerBits != 16) continue;
        
        for (int pmt=0; pmt<12; pmt++) {
            histArea[pmt]->Fill(area[pmtChannelMap[pmt]]);
        }
    }

    Double_t mu1[12] = {0};
    for (int i=0; i<12; i++) {
        if (histArea[i]->GetEntries() == 0) {
            cerr << "Empty histogram for PMT " << i+1 << endl;
            continue;
        }

        TF1 *fitFunc = new TF1("fitFunc", SPEfit, -50, 400, 8);
        Double_t histMean = histArea[i]->GetMean();
        Double_t histRMS = histArea[i]->GetRMS();

        fitFunc->SetParameters(
            1000,               // par[0]
            histMean - histRMS, // par[1]
            histRMS / 2,        // par[2]
            1000,              // par[3]
            histMean,          // par[4]
            histRMS,           // par[5]
            500,               // par[6]
            500                // par[7]
        );
        histArea[i]->Fit(fitFunc, "Q0");
        mu1[i] = fitFunc->GetParameter(4);
        delete fitFunc;
    }

    cout << "\nCALIBRATION RESULTS (1PE peak positions):\n";
    cout << "PMT#  HardwareCh  mu1 [ADC]\n";
    cout << "---------------------------\n";
    for (int i=0; i<12; i++) {
        printf("PMT%02d     %2d       %6.2f\n", 
              i+1, pmtChannelMap[i], mu1[i]);
    }
    cout << endl;

    // 2. EVENT SELECTION
    vector<Long64_t> goodEvents, badEvents;
    vector<Double_t> goodRMS, badRMS;

    cout << "Processing " << nEntries << " events for selection..." << endl;
    for (Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);
        bool isGood = false;
        Double_t currentRMS;

        // Condition A: Pulse Height > 2 p.e. for at least 3 PMTs
        int countAbove2PE = 0;
        for (int pmt=0; pmt<12; pmt++) {
            if (pulseH[pmtChannelMap[pmt]] > 2 * mu1[pmt]) {
                countAbove2PE++;
            }
        }

        if (countAbove2PE >= 3) {
            vector<Double_t> peakPositions;
            for (int pmt=0; pmt<12; pmt++) {
                peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
            }
            if (peakPositions.size() > 0) {
                Double_t dummyMean;
                CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                if (currentRMS < 2.5) isGood = true;
            }
        } 
        else {
            // Condition B: Pulse Height > 3 * baseline RMS and area/height > 1.2
            int countConditionB = 0;
            for (int pmt=0; pmt<12; pmt++) {
                int ch = pmtChannelMap[pmt];
                if (pulseH[ch] > 3 * baselineRMS[ch] && (area[ch] / pulseH[ch]) > 1.2) {
                    countConditionB++;
                }
            }

            if (countConditionB >= 3) {
                vector<Double_t> peakPositions;
                for (int pmt=0; pmt<12; pmt++) {
                    peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
                }
                if (peakPositions.size() > 0) {
                    Double_t dummyMean;
                    CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                    if (currentRMS < 2.5) isGood = true;
                }
            }
        }

        if (isGood) {
            goodEvents.push_back(entry);
            goodRMS.push_back(currentRMS);
        } else {
            badEvents.push_back(entry);
            badRMS.push_back(currentRMS);
        }
    }

    // Save good events
    string goodFilePath = outputDir + "/GoodEvents_" + baseName + ".root";
    cout << "Saving good events to: " << goodFilePath << endl;
    TFile *goodFile = new TFile(goodFilePath.c_str(), "RECREATE");
    if (!goodFile || goodFile->IsZombie()) {
        cerr << "Error creating output file: " << goodFilePath << endl;
        return;
    }

    TTree *goodTree = tree->CloneTree(0);
    Double_t peakPosition_rms;
    goodTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<goodEvents.size(); i++) {
        tree->GetEntry(goodEvents[i]);
        peakPosition_rms = goodRMS[i];
        goodTree->Fill();
    }
    goodTree->Write();
    goodFile->Close();
    cout << "Saved " << goodEvents.size() << " good events." << endl;

    // Save bad events
    string badFilePath = outputDir + "/BadEvents_" + baseName + ".root";
    cout << "Saving bad events to: " << badFilePath << endl;
    TFile *badFile = new TFile(badFilePath.c_str(), "RECREATE");
    if (!badFile || badFile->IsZombie()) {
        cerr << "Error creating output file: " << badFilePath << endl;
        return;
    }

    TTree *badTree = tree->CloneTree(0);
    badTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<badEvents.size(); i++) {
        tree->GetEntry(badEvents[i]);
        peakPosition_rms = badRMS[i];
        badTree->Fill();
    }
    badTree->Write();
    badFile->Close();
    cout << "Saved " << badEvents.size() << " bad events." << endl;

    // 3. MUON/MICHEL ANALYSIS
    cout << "\nStarting muon/Michel analysis on good events..." << endl;
    TFile *goodFileForAnalysis = new TFile(goodFilePath.c_str(), "READ");
    if (!goodFileForAnalysis || goodFileForAnalysis->IsZombie()) {
        cerr << "Error reopening good events file for analysis" << endl;
        return;
    }
    
    TTree *goodTreeForAnalysis = (TTree*)goodFileForAnalysis->Get("tree");
    if (!goodTreeForAnalysis) {
        cerr << "Error accessing tree in good events file" << endl;
        goodFileForAnalysis->Close();
        return;
    }

    analyzeMuonMichel(goodTreeForAnalysis, mu1, outputDir);
    goodFileForAnalysis->Close();

    // Cleanup
    for (int i=0; i<12; i++) delete histArea[i];
    file->Close();
    cout << "Analysis complete. Results saved in " << outputDir << endl;
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file.root>" << endl;
        return 1;
    }
    processEvents(argv[1]);
    return 0;
}
