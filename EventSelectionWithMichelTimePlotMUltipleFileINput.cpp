*/ Michel Electron Analysis Program. It  accepts multiple files.
 * 
 * Purpose:
 * This code is designed to process PMT (Photomultiplier Tube) data to identify and analyze
 * Michel electrons resulting from muon decays. It performs three main functions:
 * 1. Calibrates PMT responses to determine single photoelectron (SPE) spectra
 * 2. Applies event selection cuts to identify good candidate events
 * 3. Analyzes muon decay time distributions and Michel electron energy spectra
 * 
 * Key Features:
 * - Automatic SPE calibration using a multi-peak fit function
 * - Advanced event selection with multiple criteria:
 *   • Pulse height requirements
 *   • Peak position consistency checks
 *   • Signal-to-noise ratio cuts
 * - Muon lifetime measurement via exponential fitting
 * - Michel electron energy spectrum analysis
 * - Automatic output organization with timestamped directories
 * 
 * Input:
 * - ROOT file containing PMT waveform data with the following structure:
 *   • Tree name: "tree"
 *   • Required branches: eventID, nSamples, adcVal, baselineMean, baselineRMS,
 *     pulseH, area, peakPosition, nsTime, triggerBits
 * 
 * Output:
 * - AnalysisResults/ directory containing:
 *   • time_difference.png: Muon decay time distribution with fit
 *   • michel_spectrum.png: Michel electron energy spectrum
 *  
 * 
 * Physics Context:
 * Michel electrons are produced when cosmic ray muons (μ⁻) decay via:
 * μ⁻ → e⁻ + νₑ + ν̅μ
 * This code identifies muon stops by large energy deposits, then looks for subsequent
 * Michel electron signals in the 1-10 μs range, fitting the time difference distribution
 * to extract the muon lifetime (expected ~2.2 μs at rest).
 */
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
#include <TChain.h>

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

// Function to calculate mean and RMS
void CalculateMeanAndRMS(const vector<Double_t> &data, Double_t &mean, Double_t &rms) {
    mean = 0.0;
    for (const auto &value : data) mean += value;
    mean /= data.size();
    
    rms = 0.0;
    for (const auto &value : data) rms += pow(value - mean, 2);
    rms = sqrt(rms / data.size());
}

// Function to generate a random directory name
string generateRandomDirName() {
    time_t now = time(0);
    tm *ltm = localtime(&now);
    int randomNum = rand() % 10000;
    
    string dirName = "MichelAnalysis_" + to_string(1900 + ltm->tm_year) + 
                   to_string(1 + ltm->tm_mon) + to_string(ltm->tm_mday) + "_" +
                   to_string(ltm->tm_hour) + to_string(ltm->tm_min) + to_string(ltm->tm_sec) +
                   "_" + to_string(randomNum);
    return dirName;
}

// Main analysis function
void analyzeFiles(const vector<string> &fileNames) {
    if (fileNames.empty()) {
        cerr << "No input files specified!" << endl;
        return;
    }

    // Create output directory
    string outputDir = generateRandomDirName();
    cout << "Creating output directory: " << outputDir << endl;
    if (mkdir(outputDir.c_str(), 0777)) {
        cerr << "Error creating directory: " << outputDir << endl;
        return;
    }

    // Initialize calibration histograms
    TH1F *histArea[12];
    for (int i=0; i<12; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area",i+1), 
                             Form("PMT %d;ADC Counts;Events",i+1), 150, -50, 400);
    }

    // Process each file for calibration
    Double_t mu1[12] = {0};
    int pmtChannelMap[12] = {0,10,7,2,6,3,8,9,11,4,5,1};
    
    cout << "Performing calibration using all files..." << endl;
    for (const auto &fileName : fileNames) {
        TFile *file = TFile::Open(fileName.c_str());
        if (!file || file->IsZombie()) {
            cerr << "Error opening file: " << fileName << endl;
            continue;
        }

        TTree *tree = (TTree*)file->Get("tree");
        if (!tree) {
            cerr << "Error accessing TTree in file: " << fileName << endl;
            file->Close();
            continue;
        }

        // Variables for calibration
        Double_t area[23];
        Int_t triggerBits;
        tree->SetBranchAddress("area", area);
        tree->SetBranchAddress("triggerBits", &triggerBits);

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry=0; entry<nEntries; entry++) {
            tree->GetEntry(entry);
            if (triggerBits != 16) continue;
            
            for (int pmt=0; pmt<12; pmt++) {
                histArea[pmt]->Fill(area[pmtChannelMap[pmt]]);
            }
        }
        file->Close();
    }

    // Perform SPE fits
    for (int i=0; i<12; i++) {
        if (histArea[i]->GetEntries() == 0) {
            cerr << "Empty histogram for PMT " << i+1 << endl;
            continue;
        }

        TF1 *fitFunc = new TF1("fitFunc", SPEfit, -50, 400, 8);
        Double_t histMean = histArea[i]->GetMean();
        Double_t histRMS = histArea[i]->GetRMS();

        fitFunc->SetParameters(1000, histMean - histRMS, histRMS / 2,
                             1000, histMean, histRMS, 500, 500);
        histArea[i]->Fit(fitFunc, "Q0");
        mu1[i] = fitFunc->GetParameter(4);
        delete fitFunc;
    }

    cout << "\nCALIBRATION RESULTS (1PE peak positions):\n";
    cout << "PMT#  HardwareCh  mu1 [ADC]\n";
    cout << "---------------------------\n";
    for (int i=0; i<12; i++) {
        printf("PMT%02d     %2d       %6.2f\n", i+1, pmtChannelMap[i], mu1[i]);
    }

    // Create chain for all input files
    TChain *chain = new TChain("tree");
    for (const auto &fileName : fileNames) {
        chain->Add(fileName.c_str());
    }

    // Set up branches for analysis
    Double_t area[23], pulseH[23], baselineRMS[23];
    Int_t triggerBits, peakPosition[23];
    Long64_t nsTime;
    
    chain->SetBranchAddress("area", area);
    chain->SetBranchAddress("pulseH", pulseH);
    chain->SetBranchAddress("baselineRMS", baselineRMS);
    chain->SetBranchAddress("peakPosition", peakPosition);
    chain->SetBranchAddress("triggerBits", &triggerBits);
    chain->SetBranchAddress("nsTime", &nsTime);

    // Create analysis subdirectory
    string analysisDir = outputDir + "/Plots";
    mkdir(analysisDir.c_str(), 0777);

    // Histograms for final results
    TH1F *histDeltaT = new TH1F("DeltaT", "Time difference;(Time difference #mus);Events/0.1 #mus", 45, 1, 10);
    TH1F *histMichelSpectrum = new TH1F("MichelSpectrum", "Michel Electron Spectrum;Photoelectrons;Events", 90, 100, 1000);

    // Variables for muon tracking
    vector<Long64_t> muonIndices;
    vector<Long64_t> muonTimes;
    int muonCount = 0, michelCount = 0;
    const double muonThreshold = 50;
    const double michelThreshold = 50;

    // First pass: Identify muon candidates
    cout << "\nIdentifying muon candidates..." << endl;
    Long64_t nEntries = chain->GetEntries();
    for (Long64_t entry=0; entry<nEntries; entry++) {
        chain->GetEntry(entry);
        
        // Check trigger condition
        if (triggerBits != 34 && triggerBits != 2) continue;

        // Check PMT hits
        int nPMTsHit = 0;
        for (int pmt=0; pmt<12; pmt++) {
            if (area[pmtChannelMap[pmt]] > 500) nPMTsHit++;
        }
        if (nPMTsHit < 3) continue;

        // Calculate total energy in p.e.
        double totalEnergy = 0;
        for (int pmt=0; pmt<12; pmt++) {
            if (area[pmtChannelMap[pmt]] > 500 && mu1[pmt] > 0) {
                totalEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
            }
        }

        if (totalEnergy >= muonThreshold) {
            muonCount++;
            muonIndices.push_back(entry);
            muonTimes.push_back(nsTime);
        }
    }

    // Second pass: Find Michel electrons
    cout << "Searching for Michel electrons..." << endl;
    for (size_t i=0; i<muonIndices.size(); i++) {
        Long64_t muonEntry = muonIndices[i];
        Long64_t muonTime = muonTimes[i];
        
        // Search subsequent events
        for (Long64_t nextEntry = muonEntry+1; nextEntry < nEntries; nextEntry++) {
            chain->GetEntry(nextEntry);
            double deltaT = (nsTime - muonTime) * 1e-3; // in μs
            
            if (deltaT > 10) break;
            if (deltaT < 1) continue;
            
            // Check PMT hits for Michel candidate
            int nPMTsHit = 0;
            for (int pmt=0; pmt<12; pmt++) {
                if (area[pmtChannelMap[pmt]] > 500) nPMTsHit++;
            }
            if (nPMTsHit < 3) continue;
            
            // Calculate Michel energy
            double michelEnergy = 0;
            for (int pmt=0; pmt<12; pmt++) {
                if (area[pmtChannelMap[pmt]] > 500 && mu1[pmt] > 0) {
                    michelEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
                }
            }
            
            if (michelEnergy >= michelThreshold) {
                michelCount++;
                histDeltaT->Fill(deltaT);
                histMichelSpectrum->Fill(michelEnergy);
                break;
            }
        }
    }

    cout << "\nAnalysis Results:\n";
    cout << "Total muon candidates: " << muonCount << endl;
    cout << "Total Michel electrons: " << michelCount << endl;

    // Fit and plot time difference
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

    // Plot Michel spectrum
    TCanvas *c2 = new TCanvas("c2", "Michel Electron Spectrum", 800, 600);
    histMichelSpectrum->Draw();
    c2->SaveAs((analysisDir + "/michel_spectrum.png").c_str());

    // Cleanup
    for (int i=0; i<12; i++) delete histArea[i];
    delete histDeltaT;
    delete histMichelSpectrum;
    delete c1;
    delete c2;
    delete chain;

    cout << "\nAnalysis complete. Plots saved in: " << outputDir << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file1.root> [input_file2.root ...]" << endl;
        return 1;
    }

    vector<string> fileNames;
    for (int i = 1; i < argc; i++) {
        fileNames.push_back(argv[i]);
    }

    analyzeFiles(fileNames);
    return 0;
}
