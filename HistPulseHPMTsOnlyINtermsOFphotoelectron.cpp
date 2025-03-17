//This code gives the Histogram of pulseH of PMTs in terms of photoelectrons. We need to provide mu1 values in the code.
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TLatex.h"

using namespace std;

// Function to process the ROOT file and generate pulse height distributions for all events
void processPulseHDistributions(const char *fileName) {
    // Open the ROOT file
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

    // Access the TTree named "tree" from the ROOT file
    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree 'tree'!" << endl;
        file->Close();
        return;
    }

    // Declare variables to store data from the TTree
    Short_t adcVal[23][45]; // ADC values for 23 channels and 45 time bins
    Double_t pulseH[23];    // Pulse height for each channel

    // Set branch addresses to read data from the TTree
    tree->SetBranchAddress("adcVal", adcVal);
    tree->SetBranchAddress("pulseH", pulseH);

    // Get the total number of entries in the TTree
    Long64_t nEntries = tree->GetEntries();

    // Calibration data: 1PE peak positions for each PMT (index 0-11)
    double mu1[12] = {112.83, 112.59, 109.91, 112.83, 111.71, 109.66, 
                      101.39, 105.92, 108.13, 115.50, 101.05, 107.13};

    // Create histograms to store the pulse height distributions for each PMT in p.e. units
    TH1F *histPulseH[12];
    for (int i = 0; i < 12; i++) {
        histPulseH[i] = new TH1F(Form("PMT%d_PulseH", i + 1), 
                                 Form("PMT %d Pulse Height Distributions in terms of Photoelectrons; Photoelectrons; Events per 0.1 p.e.", i + 1), 
                                 100, 0, 10);
        histPulseH[i]->SetLineColor(kRed); // Set histogram line color to red
    }

    // Mapping of PMT channels
    int pmtChannelMap[12] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

    // Loop over all events in the TTree
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);

        // Loop through the 12 PMTs and fill their pulse height distributions in p.e.
        for (int pmt = 0; pmt < 12; pmt++) {
            int adcIndex = pmtChannelMap[pmt]; // Map PMT channels
            double pe = pulseH[adcIndex] / mu1[pmt]; // Convert ADC to p.e.
            histPulseH[pmt]->Fill(pe); // Fill the histogram with pulse height in p.e.
        }
    }

    // Create a canvas to draw the histograms
    TCanvas *canvas = new TCanvas("canvas", "PMT Pulse Height Distributions in terms of Photoelectrons", 800, 600);

    // Adjust margins and text size for better visibility
    canvas->SetLeftMargin(0.15);   // Increase left margin
    canvas->SetRightMargin(0.05);  // Adjust right margin
    canvas->SetBottomMargin(0.15); // Increase bottom margin
    canvas->SetTopMargin(0.15);    // Adjust top margin

    // Save each histogram as a PNG file
    for (int i = 0; i < 12; i++) {
        canvas->Clear(); // Clear the canvas for the next histogram

        // Adjust histogram text size
        histPulseH[i]->GetXaxis()->SetTitleSize(0.06); // Increase x-axis title size
        histPulseH[i]->GetYaxis()->SetTitleSize(0.06); // Increase y-axis title size
        histPulseH[i]->GetXaxis()->SetLabelSize(0.04); // Increase x-axis label size
        histPulseH[i]->GetYaxis()->SetLabelSize(0.04); // Increase y-axis label size

        // Move Y-axis title closer to the axis
        histPulseH[i]->GetYaxis()->SetTitleOffset(1.1); // Adjust this value as needed

        // Set title size for individual plots
        histPulseH[i]->SetTitleSize(0.13, "t"); // "t" for title, size 0.06

        histPulseH[i]->Draw(); // Draw the histogram
        canvas->SaveAs(Form("PMT%d_PulseH_Distribution.png", i + 1)); // Save as PNG
    }

    // Create a master canvas for the combined plot
    TCanvas *masterCanvas = new TCanvas("MasterCanvas", "Combined PMT Pulse Height Distributions", 3600, 3000);
    masterCanvas->Divide(3, 4, 0.01, 0.01); // Increase spacing between subplots

    // Define the layout of PMT channels on the canvas
    int layout[4][3] = {
        {9, 3, 7},  // Row 1: PMT 10, PMT 4, PMT 8
        {5, 4, 8},  // Row 2: PMT 6, PMT 5, PMT 9
        {0, 6, 1},  // Row 3: PMT 1, PMT 7, PMT 2
        {10, 11, 2} // Row 4: PMT 11, PMT 12, PMT 3
    };

    // Loop through the layout to plot histograms on the master canvas
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 3; col++) {
            int padPosition = row * 3 + col + 1; // Calculate pad position (1-12)
            masterCanvas->cd(padPosition); // Switch to the specific pad

            int pmtIndex = layout[row][col]; // Get PMT index from layout

            // Adjust histogram title and size
            histPulseH[pmtIndex]->SetTitle(""); // Clear default title
            histPulseH[pmtIndex]->GetXaxis()->SetTitleSize(0.06); // Increase x-axis title size
            histPulseH[pmtIndex]->GetYaxis()->SetTitleSize(0.06); // Increase y-axis title size
            histPulseH[pmtIndex]->GetXaxis()->SetLabelSize(0.05); // Increase x-axis label size
            histPulseH[pmtIndex]->GetYaxis()->SetLabelSize(0.05); // Increase y-axis label size

            // Move Y-axis title closer to the axis
            histPulseH[pmtIndex]->GetYaxis()->SetTitleOffset(1.1); // Adjust this value as needed

            // Ensure Y-axis label is visible
            histPulseH[pmtIndex]->GetYaxis()->SetTitle("Events per 0.1 p.e.");

            // Ensure X-axis label is visible
            histPulseH[pmtIndex]->GetXaxis()->SetTitle("Photoelectrons");

            // Adjust margins for each subplot
            gPad->SetLeftMargin(0.12);   // Increase left margin
            gPad->SetRightMargin(0.05); // Adjust right margin
            gPad->SetBottomMargin(0.12); // Increase bottom margin
            gPad->SetTopMargin(0.12);    // Increase top margin to accommodate larger title

            // Draw the histogram
            histPulseH[pmtIndex]->Draw();

            // Add custom title using TLatex
            TLatex *title = new TLatex();
            title->SetTextSize(0.14); // Set title size
            title->SetTextAlign(22);  // Center align
            title->SetNDC(true);      // Use normalized coordinates
            title->DrawLatex(0.5, 0.92, Form("PMT %d", pmtIndex + 1)); // Draw title
        }
    }

    // Save the combined canvas as a PNG file
    masterCanvas->SaveAs("Combined_PMT_PulseH_Distributions.png");

    // Clean up
    for (int i = 0; i < 12; i++) {
        delete histPulseH[i];
    }
    delete canvas;
    delete masterCanvas;
    file->Close();

    cout << "Pulse height distributions for all events saved as PNG files." << endl;
    cout << "Combined image saved as Combined_PMT_PulseH_Distributions.png" << endl;
}

// Main function to handle command-line arguments
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <root_file>" << endl;
        return 1;
    }

    const char* fileName = argv[1];
    processPulseHDistributions(fileName);

    return 0;
}
