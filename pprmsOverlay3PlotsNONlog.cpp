//This code accepts three input files.
//  And can make an overlayed histogram of three along with individual ones. All events, Good Events, Bad events.
//non-log plot
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <algorithm>
#include <cmath>

using namespace std;

const int pmtChannelMap[12] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

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
    double mean = sum / positions.size();
    double variance = (sumSq / positions.size()) - (mean * mean);
    rms = sqrt(variance > 0 ? variance : 0);
}

TH1F* processFile(const char* fileName, int color, int lineStyle, const char* histName) {
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return nullptr;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree in file: " << fileName << endl;
        file->Close();
        return nullptr;
    }

    Int_t peakPosition[23];
    if (!tree->GetBranch("peakPosition")) {
        cerr << "Branch 'peakPosition' not found in file: " << fileName << endl;
        file->Close();
        return nullptr;
    }
    tree->SetBranchAddress("peakPosition", peakPosition);

    vector<double> eventRMS;
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);
        vector<double> validPositions;
        for (int pmt = 0; pmt < 12; pmt++) {
            int ch = pmtChannelMap[pmt];
            if (ch < 0 || ch >= 23) continue;
            if (isnan(peakPosition[ch]) || isinf(peakPosition[ch])) continue;
            if (peakPosition[ch] >= 20 && peakPosition[ch] <= 44) {
                validPositions.push_back(peakPosition[ch]);
            }
        }
        if (validPositions.size() < 2) continue;
        double rms;
        CalculateRMS(validPositions, rms);
        if (rms >= 0) eventRMS.push_back(rms);
    }

    if (eventRMS.empty()) {
        file->Close();
        return nullptr;
    }

    TH1F* hRMS = new TH1F(histName, "Peak Position RMS Distribution; RMS; Events/0.1 RMS", 100, 0, 10);
    hRMS->SetDirectory(0);
    for (const auto& rms : eventRMS) hRMS->Fill(rms);
    hRMS->SetLineColor(color);
    hRMS->SetLineStyle(lineStyle);
    hRMS->SetLineWidth(2);

    file->Close();
    return hRMS;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <file3.root>" << endl;
        return 1;
    }

    // Process files with different line styles and colors
    TH1F* h1 = processFile(argv[1], kBlue, 1, "h1");   // Solid blue line (style=1)
    TH1F* h2 = processFile(argv[2], kRed, 2, "h2");    // Dashed red line (style=2)
    TH1F* h3 = processFile(argv[3], kBlack, 2, "h3");  // Dashed black line (style=2)

    if (h1 && h2 && h3) {
        // Create individual plots for each histogram
        TCanvas* c1 = new TCanvas("c1", "ALL EVENTS", 800, 600);
        //c1->SetLogy();
        h1->Draw("HIST");

        // Add legend to individual plot (ALL EVENTS)
        TLegend* legend1 = new TLegend(0.5, 0.5, 0.60, 0.60); // Just below top-right corner
        legend1->AddEntry(h1, "All Events", "l");
        legend1->SetBorderSize(0); // Remove box around the legend
        legend1->SetTextSize(0.05); // Increase font size
        legend1->Draw();

        c1->Update();
        c1->SaveAs("ALL_EVENTS.png");

        TCanvas* c2 = new TCanvas("c2", "GOOD EVENTS", 800, 600);
        //c2->SetLogy();
        h2->Draw("HIST");

        // Add legend to individual plot (GOOD EVENTS)
        TLegend* legend2 = new TLegend(0.5, 0.5, 0.60, 0.60); // Just below top-right corner
        legend2->AddEntry(h2, "Good Events", "l");
        legend2->SetBorderSize(0); // Remove box around the legend
        legend2->SetTextSize(0.04); // Increase font size
        legend2->Draw();

        c2->Update();
        c2->SaveAs("GOOD_EVENTS.png");

        TCanvas* c3 = new TCanvas("c3", "BAD EVENTS", 800, 600);
        //c3->SetLogy();
        h3->Draw("HIST");

        // Add legend to individual plot (BAD EVENTS)
        TLegend* legend3 = new TLegend(0.5, 0.5, 0.60, 0.60); // Just below top-right corner
        legend3->AddEntry(h3, "Bad Events", "l");
        legend3->SetBorderSize(0); // Remove box around the legend
        legend3->SetTextSize(0.04); // Increase font size
        legend3->Draw();

        c3->Update();
        c3->SaveAs("BAD_EVENTS.png");

        // Create combined plot
        TCanvas* c_combined = new TCanvas("c_combined", "RMS Distribution", 800, 600);
       // c_combined->SetLogy();

       //Set margins (left, right, bottom, top) - values are fractions of the canvas
c_combined->SetLeftMargin(0.15);    // Increase left margin for y-axis label
c_combined->SetRightMargin(0.08);   // Right margin
c_combined->SetBottomMargin(0.12);  // Bottom margin for x-axis label
c_combined->SetTopMargin(0.08);     // Top margin

        // Find maximum y-axis value to set the range
        double maxY = max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
        //h1->SetMaximum(maxY); // Increase range for log scale

        // Draw histograms
        h1->Draw("HIST");
        h2->Draw("HIST SAME");
        h3->Draw("HIST SAME");

        // Add legend with custom labels and position
        TLegend* legend = new TLegend(0.5, 0.5, 0.60, 0.60); // Just below top-right corner
        legend->AddEntry(h1, "All Events", "l");   // Custom label for h1
        legend->AddEntry(h2, "Good Events", "l");  // Custom label for h2
        legend->AddEntry(h3, "Bad Events", "l");   // Custom label for h3
        legend->SetBorderSize(0); // Remove box around the legend
        legend->SetTextSize(0.03); // Increase font size
        legend->Draw();

        c_combined->Update();
        c_combined->SaveAs("combined_plot.png");

        // Cleanup
        delete c1;
        delete c2;
        delete c3;
        delete c_combined;
    } else {
        cerr << "Failed to create one or more histograms." << endl;
    }

    // Cleanup
    delete h1;
    delete h2;
    delete h3;

    return 0;
}
