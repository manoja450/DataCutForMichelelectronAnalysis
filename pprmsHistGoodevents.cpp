//This code gives the histogram of peakPosition_rms of the data after the applied cut. The data after cut applied is stored in a new root file, which contains an additional branch  peakPosition_rms.
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <climits>

void plotPPRMS(const char* fileName) {
    // Open ROOT file
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    // Access TTree
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error accessing TTree!" << std::endl;
        file->Close();
        return;
    }

    // Set up branch access
    Double_t pprms;
    tree->SetBranchAddress("peakPosition_rms", &pprms);

    // Find data range and maximum value
    Double_t maxVal = -DBL_MAX;
    Double_t minVal = DBL_MAX;
    const Long64_t nEntries = tree->GetEntries();
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (pprms > maxVal) maxVal = pprms;
        if (pprms < minVal) minVal = pprms;
    }

    // Handle edge cases
    if (maxVal == -DBL_MAX) {
        std::cerr << "No valid ppRMS values found!" << std::endl;
        file->Close();
        return;
    }

    // Create histogram with dynamic range
    const int bins = 100;
    TH1F* h = new TH1F("h", "Peak Position RMS Distribution;ppRMS [bins];Events",
                       bins, minVal, maxVal * 1.05);  // 5% headroom

    // Fill histogram
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        h->Fill(pprms);
    }

    // Create and configure canvas
    TCanvas* c = new TCanvas("c", "ppRMS Distribution", 1200, 800);
    c->SetGrid();
    
    // Customize histogram appearance
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->SetFillStyle(3003);
    h->SetFillColor(kBlue);
    h->Draw("HIST");  // Histogram step drawing

    // Add annotations
    TLatex tex;
    tex.SetNDC(true);
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.15, 0.88, Form("File: %s", fileName));
    tex.DrawLatex(0.15, 0.83, Form("Entries: %lld", nEntries));
    tex.DrawLatex(0.15, 0.78, Form("Maximum ppRMS: %.2f", maxVal));
    tex.DrawLatex(0.15, 0.73, Form("Mean: %.2f", h->GetMean()));

    // Save and clean up
    c->SaveAs("pprms_distribution.png");
    std::cout << "Maximum ppRMS value: " << maxVal << std::endl;

    delete h;
    delete c;
    file->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input.root>" << std::endl;
        return 1;
    }
    plotPPRMS(argv[1]);
    return 0;
}
