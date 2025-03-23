//This program reads a ROOT file containing event data and generates a histogram  
// showing the number of PMTs with PulseRatios greater than 1 for good events.
//The program performs the following steps:
 /// 1. Opens the specified ROOT file and accesses the TTree named "tree".
 // 2. Reads the "PulseRatios" branch, which stores pulse height ratios for different PMTs.
 // 3. Counts how many PMTs in each event have a PulseRatio greater than 1.
 //4. Creates a stacked histogram where each bin represents the number of PMTs exceeding this threshold.
 //5. Assigns different colors to different bins using multiple histograms in a THStack.
 // 6. Saves the resulting histogram as a PNG file named "PulseRatioHistogram_Pass.png".
 
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <THStack.h>

using namespace std;

void createPulseRatioHistogram(const char *passFileName) {
    TFile *passFile = TFile::Open(passFileName);
    if (!passFile || passFile->IsZombie()) {
        cerr << "Error: Cannot open file " << passFileName << endl;
        return;
    }

    TTree *passTree = (TTree*)passFile->Get("tree");
    if (!passTree) {
        cerr << "Error: Cannot access TTree!" << endl;
        passFile->Close();
        return;
    }

    Double_t PulseRatios[23];
    const int pmtChannelMap[12] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};
    passTree->SetBranchAddress("PulseRatios", PulseRatios);

    // Create multiple histograms, one for each bin
    THStack *stack = new THStack("stack", "Number of PMTs with PulseRatios > 1 (Good Events)");
    TH1F *histBins[12];
    const int colors[12] = {kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kPink, kSpring, kTeal, kAzure};

    for (int i = 0; i < 12; i++) {
        histBins[i] = new TH1F(Form("bin_%d", i+1), "", 13, 0, 13);
        histBins[i]->SetFillColor(colors[i]);
        stack->Add(histBins[i]);
    }

    Long64_t passEntries = passTree->GetEntries();
    for (Long64_t entry = 0; entry < passEntries; entry++) {
        passTree->GetEntry(entry);
        Int_t countGreaterThan1 = 0;
        for (int pmt = 0; pmt < 12; pmt++) {
            int channel = pmtChannelMap[pmt];
            if (PulseRatios[channel] > 1) {
                countGreaterThan1++;
            }
        }
        if (countGreaterThan1 >= 1 && countGreaterThan1 <= 12) {
            histBins[countGreaterThan1 - 1]->Fill(countGreaterThan1);
        }
    }

    TCanvas *canvas = new TCanvas("canvas", "PulseRatio Histogram (Good Events)", 800, 600);
    stack->Draw("hist");
    stack->GetXaxis()->SetTitle("Number of PMTs");
    stack->GetYaxis()->SetTitle("Number of Events");

    canvas->SaveAs("PulseRatioHistogram_Pass.png");

    for (int i = 0; i < 12; i++) {
        delete histBins[i];
    }
    delete stack;
    delete canvas;
    passFile->Close();
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <pass_file.root>" << endl;
        return 1;
    }
    createPulseRatioHistogram(argv[1]);
    return 0;
}
