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
#include <cmath>     // For sqrt, ceil, isnan, isinf

using namespace std;

// PMT to ADC channel mapping
const int pmtChannelMap[12] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

// Function to calculate RMS
void CalculateRMS(const vector<double>& positions, double& rms) {
    if (positions.size() < 2) {
        rms = -1;
        return;
    }

    double sum = 0.0, sumSq = 0.0;
    for (const auto& pos : positions) {
        sum += pos;
        sumSq += pos * pos;
    }

    const double mean = sum / positions.size();
    double variance = (sumSq / positions.size()) - (mean * mean);
    if (variance < 0) variance = 0; // Handle precision issues
    rms = sqrt(variance);
}

// Main function to plot PMT peak RMS
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

    // Debugging counters
    int excludedRange = 0;      // Excluded due to peakPosition outside [20, 44]
    int excludedPMTs = 0;       // Excluded due to fewer than 2 valid PMTs
    int excludedRMS = 0;        // Excluded due to invalid RMS (< 0)
    int excludedNaN = 0;        // Excluded due to NaN peakPosition
    int excludedInf = 0;        // Excluded due to Inf peakPosition

    // Event processing
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);

        vector<double> validPositions;
        for (int pmt = 0; pmt < 12; pmt++) {
            const int ch = pmtChannelMap[pmt];

            // Check for NaN or Inf
            if (isnan(peakPosition[ch])) {
                excludedNaN++;
                continue;
            }
            if (isinf(peakPosition[ch])) {
                excludedInf++;
                continue;
            }

            // Check for valid range
            if (peakPosition[ch] >= 20 && peakPosition[ch] <= 44) {
                validPositions.push_back(peakPosition[ch]);
            } else {
                excludedRange++;
            }
        }

        // Check for sufficient valid PMTs
        if (validPositions.size() < 2) {
            excludedPMTs++;
            continue;
        }

        // Calculate RMS
        double rms;
        CalculateRMS(validPositions, rms);
        if (rms < 0) { // Only exclude if RMS is negative
            excludedRMS++;
            continue;
        }

        // Store valid RMS (including RMS = 0)
        eventRMS.push_back(rms);
    }

    // Print debugging information
    cout << "Total entries: " << nEntries << endl;
    cout << "Excluded due to range [20, 44]: " << excludedRange << endl;
    cout << "Excluded due to NaN peakPosition: " << excludedNaN << endl;
    cout << "Excluded due to Inf peakPosition: " << excludedInf << endl;
    cout << "Excluded due to fewer than 2 valid PMTs: " << excludedPMTs << endl;
    cout << "Excluded due to invalid RMS (< 0): " << excludedRMS << endl;
    cout << "Valid events with RMS (including RMS = 0): " << eventRMS.size() << endl;

    // Create histogram
    if (eventRMS.empty()) {
        cerr << "No valid RMS values found!" << endl;
        return;
    }

    const double maxRMS = *max_element(eventRMS.begin(), eventRMS.end());
    const int numBins = 100; // Fixed number of bins
    TH1F* hRMS = new TH1F("hRMS", "Peak Position RMS Distribution before Cut; peakPosition RMS; Events",
                          numBins, 0, 10);

    for (const auto& rms : eventRMS) {
        hRMS->Fill(rms);
    }

    // Draw and save
    TCanvas* c = new TCanvas("c", "RMS Distribution", 800, 600);
    hRMS->SetFillColor(kBlue);
    hRMS->SetFillStyle(3003);
    hRMS->Draw("HIST");


    c->SaveAs("pmt_peak_rms.png");

    // Cleanup
    delete hRMS;
    delete c;
    file->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input.root>" << endl;
        return 1;
    }
    plotPMTPeakRMS(argv[1]);
    return 0;
}
