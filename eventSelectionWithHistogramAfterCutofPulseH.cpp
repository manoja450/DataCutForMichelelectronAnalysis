#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <vector>
#include <cmath>
#include <unistd.h>

using namespace std;

// SPE fitting function
Double_t SPEfit(Double_t *x, Double_t *par) {
    Double_t term1 = par[0] * exp(-0.5 * pow((x[0]-par[1])/par[2], 2));
    Double_t term2 = par[3] * exp(-0.5 * pow((x[0]-par[4])/par[5], 2));
    Double_t term3 = par[6] * exp(-0.5 * pow((x[0]-sqrt(2)*par[4])/sqrt(2*pow(par[5],2)-pow(par[2],2)), 2));
    Double_t term4 = par[7] * exp(-0.5 * pow((x[0]-sqrt(3)*par[4])/sqrt(3*pow(par[5],2)-2*pow(par[2],2)), 2));
    return term1 + term2 + term3 + term4;
}

void CalculateMeanAndRMS(const vector<Double_t> &data, Double_t &mean, Double_t &rms) {
    mean = 0.0;
    for (const auto &value : data) mean += value;
    mean /= data.size();
    
    rms = 0.0;
    for (const auto &value : data) rms += pow(value - mean, 2);
    rms = sqrt(rms / data.size());
}

void processEvents(const char *fileName) {
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

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
            1000,               // par[0]: Amplitude of the pedestal peak
            histMean - histRMS, // par[1]: Mean of the pedestal peak
            histRMS / 2,        // par[2]: Sigma of the pedestal peak
            1000,               // par[3]: Amplitude of the 1PE peak
            histMean,           // par[4]: Mean of the 1PE peak (mu1)
            histRMS,            // par[5]: Sigma of the 1PE peak
            500,                // par[6]: Amplitude of the 2PE peak
            500                 // par[7]: Amplitude of the 3PE peak
        );
        histArea[i]->Fit(fitFunc, "Q0");
        mu1[i] = fitFunc->GetParameter(4);
        delete fitFunc;
    }

    // Print calibration results
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
            // Calculate peak position mean and RMS
            vector<Double_t> peakPositions;
            for (int pmt=0; pmt<12; pmt++) {
                if (pulseH[pmtChannelMap[pmt]] > 5.25) { // PMT Hit > 3 condition
                    peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
                }
            }
            if (peakPositions.size() > 0) {
                Double_t dummyMean;
                CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                if (currentRMS < 2.5) isGood = true;
            }
        } 
        else {
            // Condition B: Pulse Height > 3 * baseline RMS and area/height > 1.2 for all PMTs
            bool allPassConditionB = true;
            for (int pmt=0; pmt<12; pmt++) {
                int ch = pmtChannelMap[pmt];
                if (pulseH[ch] <= 3 * baselineRMS[ch] || (area[ch] / pulseH[ch]) <= 1.2) {
                    allPassConditionB = false;
                    break;
                }
            }

            if (allPassConditionB) {
                // Calculate peak position mean and RMS
                vector<Double_t> peakPositions;
                for (int pmt=0; pmt<12; pmt++) {
                    if (pulseH[pmtChannelMap[pmt]] > 5.25) { // PMT Hit > 3 condition
                        peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
                    }
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

    // Save results
    TFile *goodFile = new TFile(Form("./GoodEvents_%d.root", getpid()), "RECREATE");
    TTree *goodTree = tree->CloneTree(0); // Clone the entire input tree structure

    // Add a new branch for peakPosition_rms
    Double_t peakPosition_rms;
    goodTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<goodEvents.size(); i++) {
        tree->GetEntry(goodEvents[i]); // Load all channel data for the event
        peakPosition_rms = goodRMS[i]; // Add the calculated RMS value
        goodTree->Fill();              // Save the event with all channel data
    }
    goodTree->Write();
    delete goodFile;

    TFile *badFile = new TFile(Form("./BadEvents_%d.root", getpid()), "RECREATE");
    TTree *badTree = tree->CloneTree(0); // Clone the entire input tree structure
    badTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<badEvents.size(); i++) {
        tree->GetEntry(badEvents[i]); // Load all channel data for the event
        peakPosition_rms = badRMS[i]; // Add the calculated RMS value
        badTree->Fill();              // Save the event with all channel data
    }
    badTree->Write();
    delete badFile;

    // 3. PLOT HISTOGRAMS
    // Create histograms for pulseH before and after the cut
    TH1F *histPulseHBefore[12];
    TH1F *histPulseHAfter[12];

    for (int i=0; i<12; i++) {
        histPulseHBefore[i] = new TH1F(Form("PMT%d_Before",i+1), Form("PMT %d (Before Cut);pulseH [ADC];Events",i+1), 100, 0, 100);
        histPulseHAfter[i] = new TH1F(Form("PMT%d_After",i+1), Form("PMT %d (After Cut);pulseH [ADC];Events",i+1), 100, 0, 100);
    }

    // Fill histograms for pulseH before the cut
    for (Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);
        for (int pmt=0; pmt<12; pmt++) {
            histPulseHBefore[pmt]->Fill(pulseH[pmtChannelMap[pmt]]); // Fill with ADC counts
        }
    }

    // Fill histograms for pulseH after the cut
    TFile *goodFileRead = TFile::Open(Form("./GoodEvents_%d.root", getpid()));
    TTree *goodTreeRead = (TTree*)goodFileRead->Get("tree");
    goodTreeRead->SetBranchAddress("pulseH", pulseH);

    for (size_t i=0; i<goodEvents.size(); i++) {
        goodTreeRead->GetEntry(i);
        for (int pmt=0; pmt<12; pmt++) {
            histPulseHAfter[pmt]->Fill(pulseH[pmtChannelMap[pmt]]); // Fill with ADC counts
        }
    }

    // Create a master canvas for the combined plot
    TCanvas *masterCanvas = new TCanvas("MasterCanvas", "Combined PMT Energy Distributions", 3600, 3000);
    masterCanvas->Divide(3, 4, 0, 0); // Adjust spacing between subplots

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

            // Adjust histogram text size
            histPulseHBefore[pmtIndex]->GetXaxis()->SetTitleSize(0.07); // Increase x-axis title size
            histPulseHBefore[pmtIndex]->GetYaxis()->SetTitleSize(0.09); // Increase y-axis title size
            histPulseHBefore[pmtIndex]->GetXaxis()->SetLabelSize(0.04); // Increase x-axis label size
            histPulseHBefore[pmtIndex]->GetYaxis()->SetLabelSize(0.04); // Increase y-axis label size

            // Ensure Y-axis label is visible
            histPulseHBefore[pmtIndex]->GetYaxis()->SetTitle("Events per 3 ADCs");
            histPulseHBefore[pmtIndex]->GetYaxis()->SetTitleOffset(0.8); // Move Y-axis title closer to the axis line

            // Ensure X-axis label is visible
            histPulseHBefore[pmtIndex]->GetXaxis()->SetTitle("pulseH [ADC]");

            // Adjust margins for each subplot
            gPad->SetLeftMargin(0.15);   // Increase left margin
            gPad->SetRightMargin(0.00); // Adjust right margin
            gPad->SetBottomMargin(0.15); // Increase bottom margin
            gPad->SetTopMargin(0.01);    // Adjust top margin

            // Draw histograms
            histPulseHBefore[pmtIndex]->SetLineColor(kBlue);
            histPulseHAfter[pmtIndex]->SetLineColor(kRed);
            histPulseHBefore[pmtIndex]->Draw();
            histPulseHAfter[pmtIndex]->Draw("SAME");

            // Draw a vertical line at the cut threshold (2 * mu1)
            TLine *cutLine = new TLine(2 * mu1[pmtIndex], 0, 2 * mu1[pmtIndex], histPulseHBefore[pmtIndex]->GetMaximum());
            cutLine->SetLineColor(kBlack);
            cutLine->SetLineStyle(2);
            cutLine->Draw();

            // Add a legend
            TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
            legend->AddEntry(histPulseHBefore[pmtIndex], "Before Cut", "l");
            legend->AddEntry(histPulseHAfter[pmtIndex], "After Cut", "l");
            legend->Draw();
        }
    }

    // Save the master canvas
    masterCanvas->SaveAs("CombinedPMTEnergyDistributions.png");

    // Cleanup
    for (int i=0; i<12; i++) {
        delete histPulseHBefore[i];
        delete histPulseHAfter[i];
    }
    delete masterCanvas;
    file->Close();
    goodFileRead->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file.root>" << endl;
        return 1;
    }
    processEvents(argv[1]);
    return 0;
}