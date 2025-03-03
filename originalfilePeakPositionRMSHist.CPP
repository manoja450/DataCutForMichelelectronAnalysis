#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <cmath>
#include <vector>
#include "TStyle.h"
#include "TLatex.h"

using namespace std;

void plotPMTPeakRMS(const char* fileName) {
    // Open the ROOT file
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

    // Access the TTree
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree 'tree'!" << endl;
        file->Close();
        return;
    }

    // Declare tree variables
    Int_t peakPosition[23];
    
    // Set branch addresses
    tree->SetBranchAddress("peakPosition", peakPosition);

    // PMT channel mapping
    int pmtChannelMap[12] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

    // Vector to store RMS values for each event
    vector<double> eventRMS;

    // Process all events
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);

        // Collect peak positions for all PMTs in this event
        vector<double> peakPositions;
        for (int pmt = 0; pmt < 12; pmt++) {
            int adcChannel = pmtChannelMap[pmt];
            int pos = peakPosition[adcChannel];
            
            // Skip invalid positions (20-44 bins valid range)
            if (pos < 20 || pos > 44) continue;

            peakPositions.push_back(pos);
        }

        // Calculate RMS for this event
        if (peakPositions.size() > 1) {  // Require at least 2 PMTs for RMS calculation
            double sum = 0, sumSq = 0;
            for (double pos : peakPositions) {
                sum += pos;
                sumSq += pos * pos;
            }
            double mean = sum / peakPositions.size();
            double rms = sqrt((sumSq / peakPositions.size()) - (mean * mean));
            eventRMS.push_back(rms);
        }
    }

    // Create histogram for event RMS values
    TH1F* hRMS = new TH1F("hRMS", "Peak Position RMS Distribution;RMS;Events", 
                           100, 0, 10); // 100 bins, range 0 to 10

    // Fill histogram with RMS values
    for (double rms : eventRMS) {
        hRMS->Fill(rms);
    }

    // Customize ROOT style
    gStyle->SetOptStat(1111); // Enable default statistics box
    gStyle->SetTitleFont(42, "XYZ"); // Set font for axis titles
    gStyle->SetLabelFont(42, "XYZ"); // Set font for axis labels
    gStyle->SetTitleSize(0.05, "XYZ"); // Set title size
    gStyle->SetLabelSize(0.04, "XYZ"); // Set label size
    gStyle->SetPadLeftMargin(0.12); // Increase left margin
    gStyle->SetPadRightMargin(0.05); // Adjust right margin
    gStyle->SetPadBottomMargin(0.12); // Increase bottom margin
    gStyle->SetPadTopMargin(0.08); // Adjust top margin
    gStyle->SetGridColor(kGray); // Set grid color
    gStyle->SetGridStyle(2); // Set grid style (dashed lines)

    // Create and configure canvas
    TCanvas* canvas = new TCanvas("canvas", "Peak Position RMS", 1200, 800);
    canvas->SetGridx(); // Enable vertical grid lines
    canvas->SetGridy(); // Enable horizontal grid lines

    // Customize histogram
    hRMS->SetLineColor(kBlue);
    hRMS->SetLineWidth(2);
    hRMS->SetFillColor(kBlue);
    hRMS->SetFillStyle(3003); // Semi-transparent fill
    hRMS->GetXaxis()->CenterTitle(true);
    hRMS->GetYaxis()->CenterTitle(true);
    hRMS->GetXaxis()->SetTitleOffset(1.1); // Adjust title offset
    hRMS->GetYaxis()->SetTitleOffset(1.3); // Adjust title offset

    // Draw histogram
    hRMS->Draw("HIST");


    // Save the plot
    canvas->SaveAs("PeakPosition_RMS_Distribution.png");

    // Cleanup
    delete hRMS;
    delete canvas;
    file->Close();

    cout << "Peak position RMS distribution plot saved as PeakPosition_RMS_Distribution.png" << endl;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <root_file>" << endl;
        return 1;
    }

    plotPMTPeakRMS(argv[1]);
    return 0;
}
