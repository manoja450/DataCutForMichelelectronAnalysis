//This code gives the histogram of peakPosition_rms of the data after the applied cut. The data after cut applied is stored in a new root file, which contains an additional branch  peakPosition_rms.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <cmath> // For std::max

void plot_peakPosition_rms_curve(const char* fileName) {
    // Open the ROOT file
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return;
    }

    // Access the TTree
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Could not access TTree 'tree'!" << std::endl;
        file->Close();
        return;
    }

    // Set up branch access
    Double_t peakPosition_rms;
    tree->SetBranchAddress("peakPosition_rms", &peakPosition_rms);

    // Variables to track maximum value
    Double_t maxPPRMS = -1.0;

    // Create a histogram with fixed x-axis range (0 to 10)
    const int nBins = 100;
    TH1D* hist = new TH1D("hist", "Peak Position RMS Distribution After Cut;Peak Position RMS;Events/0.1 RMS", nBins, 0, 10);

    // Fill the histogram and track the maximum value
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        hist->Fill(peakPosition_rms); // Fill the histogram
        maxPPRMS = std::max(maxPPRMS, peakPosition_rms); // Track maximum value
    }

    // Check if histogram has entries
    if (hist->GetEntries() == 0) {
        std::cerr << "Error: No entries found in the histogram!" << std::endl;
        delete hist;
        file->Close();
        return;
    }

    // Print the maximum value of peakPosition_rms
    std::cout << "Maximum peakPosition_rms: " << maxPPRMS << std::endl;

    // Create and configure canvas
    TCanvas* canvas = new TCanvas("canvas", "Peak Position RMS", 800, 600);
    canvas->SetGrid();

    // Enable statistics box
    gStyle->SetOptStat(1111); // Show entries, mean, and RMS
    hist->SetStats(1);        // Ensure statistics box is enabled

    // Increase axis text size
    hist->GetXaxis()->SetLabelSize(0.04); // Increase x-axis label size
    hist->GetXaxis()->SetTitleSize(0.04); // Increase x-axis title size
    hist->GetYaxis()->SetLabelSize(0.04); // Increase y-axis label size
    hist->GetYaxis()->SetTitleSize(0.05); // Increase y-axis title size

    // Customize histogram appearance
    hist->SetFillColor(kBlue);
    hist->SetFillStyle(3003); // Semi-transparent fill

    // Explicitly set the x-axis range to ensure all data is visible
    hist->GetXaxis()->SetRangeUser(0, 10); // Ensure the range is 0 to 10

    // Draw the histogram
    hist->Draw("HIST");

    // Save the plot
    canvas->SaveAs("peakPosition_rms_histogram.png");

    // Clean up
    delete canvas;
    delete hist;
    file->Close();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    plot_peakPosition_rms_curve(argv[1]);
    return 0;
}
