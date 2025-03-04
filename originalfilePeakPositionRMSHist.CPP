//This code calculates each event's peak Position RMS and plots the histogram of peak position RMS (of PMTs only). 
// This code is needed to see the distribution of peak position before the cut, i.e, the distribution of original data.
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <algorithm> // For max_element

using namespace std;

// PMT to ADC channel mapping from your code
const int pmtChannelMap[12] = {0,10,7,2,6,3,8,9,11,4,5,1};

void CalculateRMS(const vector<double>& positions, double& rms) {
    if(positions.size() < 2) {
        rms = -1;
        return;
    }
    
    double sum = 0.0, sumSq = 0.0;
    for(const auto& pos : positions) {
        sum += pos;
        sumSq += pos * pos;
    }
    
    const double mean = sum/positions.size();
    rms = sqrt((sumSq/positions.size()) - (mean*mean));
}

void plotPMTPeakRMS(const char* fileName) {
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree!" << endl;
        file->Close();
        return;
    }

    // Setup input branches
    Int_t peakPosition[23];
    tree->SetBranchAddress("peakPosition", peakPosition);

    vector<double> eventRMS;
    const Long64_t nEntries = tree->GetEntries();
    
    // Event processing
    for(Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);
        
        vector<double> validPositions;
        for(int pmt = 0; pmt < 12; pmt++) {
            const int ch = pmtChannelMap[pmt];
            
            // Basic validity check (20-44 bins from original code)
            if(peakPosition[ch] >= 20 && peakPosition[ch] <= 44) {
                validPositions.push_back(peakPosition[ch]);
            }
        }
        
        double rms;
        CalculateRMS(validPositions, rms);
        if(rms > 0) eventRMS.push_back(rms);
    }

    // Create histogram
    if(eventRMS.empty()) {
        cerr << "No valid RMS values found!" << endl;
        return;
    }
    
    const double maxRMS = *max_element(eventRMS.begin(), eventRMS.end());
    TH1F* hRMS = new TH1F("hRMS", "PMT Peak Position RMS Distribution;RMS (bins);Events",
                          int(maxRMS*2), 0, ceil(maxRMS*1.1));
    
    for(const auto& rms : eventRMS) {
        hRMS->Fill(rms);
    }

    // Draw and save
    TCanvas* c = new TCanvas("c", "RMS Distribution", 800, 600);
    hRMS->SetFillColor(kBlue);
    hRMS->SetFillStyle(3003);
    hRMS->Draw("HIST");
    
    TLatex tex;
    tex.SetNDC(true);
    tex.DrawLatex(0.15, 0.85, Form("Mean RMS: %.2f", hRMS->GetMean()));
    tex.DrawLatex(0.15, 0.80, Form("Entries: %d", hRMS->GetEntries()));
    
    c->SaveAs("pmt_peak_rms.png");

    // Cleanup
    delete hRMS;
    delete c;
    file->Close();
}

int main(int argc, char* argv[]) {
    if(argc != 2) {
        cerr << "Usage: " << argv[0] << " <input.root>" << endl;
        return 1;
    }
    plotPMTPeakRMS(argv[1]);
    return 0;
}
