#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
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

    // Branch variables
    Int_t eventID;
    Int_t nSamples[23];
    Short_t adcVal[23][45];
    Double_t baselineMean[23], baselineRMS[23], pulseH[23], area[23];
    Int_t peakPosition[23];
    Long64_t nsTime;
    Int_t triggerBits;

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
            1000,               // par[0]
            histMean - histRMS,  // par[1]
            histRMS / 2,         // par[2]
            1000,               // par[3]
            histMean,           // par[4] (mu1)
            histRMS,            // par[5]
            500,                // par[6]
            500                 // par[7]
        );
        histArea[i]->Fit(fitFunc, "Q0");
        mu1[i] = fitFunc->GetParameter(4);
        delete fitFunc;
    }

    // Print calibration results with headers
    cout << "\nCALIBRATION RESULTS (1PE peak positions):\n";
    cout << "==========================================\n";
    cout << "PMT#  HardwareCh  mu1 [ADC]\n";
    cout << "---------------------------\n";
    for (int i=0; i<12; i++) {
        printf("PMT%02d     %2d       %6.2f\n", 
              i+1, pmtChannelMap[i], mu1[i]);
    }
    cout << "==========================================\n\n";

    // 2. EVENT SELECTION (same as before)

    // 3. PLOT PULSEH DISTRIBUTIONS
    TFile *originalFile = TFile::Open(fileName);
    TTree *originalTree = (TTree*)originalFile->Get("tree");
    TFile *goodFile = TFile::Open(Form("./GoodEvents_%d.root", getpid()));
    TTree *goodTree = (TTree*)goodFile->Get("tree");

    TH1F *histOriginal[12];
    TH1F *histGood[12];

    // Initialize histograms for pulseH
    for(int pmt=0; pmt<12; pmt++) {
        histOriginal[pmt] = new TH1F(Form("PMT%d_Original",pmt+1),
                                    Form("PMT %d;pulseH [ADC];Events",pmt+1), 
                                    40, -10, 30);
        
        histGood[pmt] = new TH1F(Form("PMT%d_Good",pmt+1),
                                Form("PMT %d;pulseH [ADC];Events",pmt+1), 
                                40, -10, 30);
        
        // Fill original histograms
        originalTree->Draw(Form("pulseH[%d] >> PMT%d_Original", 
                              pmtChannelMap[pmt], pmt+1), "");
        // Fill good histograms
        goodTree->Draw(Form("pulseH[%d] >> PMT%d_Good", 
                          pmtChannelMap[pmt], pmt+1), "");
        
        // Style settings
        histOriginal[pmt]->SetLineColor(kBlue);
        histGood[pmt]->SetLineColor(kRed);
        histOriginal[pmt]->SetLineWidth(2);
        histGood[pmt]->SetLineWidth(2);
    }

    // Create combined canvas
    TCanvas *masterCanvas = new TCanvas("MasterCanvas", 
                                       "PMT Pulse Height Distributions", 
                                       3600, 3000);
    masterCanvas->Divide(3, 4, 0, 0);

    int layout[4][3] = {
        {9, 3, 7},  // Row 1: PMT 10, 4, 8
        {5, 4, 8},  // Row 2: PMT 6, 5, 9
        {0, 6, 1},  // Row 3: PMT 1, 7, 2
        {10, 11, 2} // Row 4: PMT 11, 12, 3
    };

    for (int row=0; row<4; row++) {
        for (int col=0; col<3; col++) {
            int padIndex = row*3 + col + 1;
            masterCanvas->cd(padIndex);
            
            int pmtIndex = layout[row][col];
            TH1F *hOrig = histOriginal[pmtIndex];
            TH1F *hGood = histGood[pmtIndex];
            
            // Set axis ranges
            double ymax = std::max(hOrig->GetMaximum(), hGood->GetMaximum()) * 1.2;
            hOrig->SetMaximum(ymax);
            
            // Draw histograms
            hOrig->Draw("HIST");
            hGood->Draw("HIST SAME");
            
            // Add legend
            TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
            leg->AddEntry(hOrig, "Original", "l");
            leg->AddEntry(hGood, "Good Events", "l");
            leg->SetBorderSize(0);
            leg->Draw();
            
            // Style adjustments
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.15);
            gPad->SetTopMargin(0.1);
            
            hOrig->GetXaxis()->SetTitleSize(0.07);
            hOrig->GetYaxis()->SetTitleSize(0.07);
            hOrig->GetXaxis()->SetLabelSize(0.05);
            hOrig->GetYaxis()->SetLabelSize(0.05);
            hOrig->GetYaxis()->SetTitleOffset(1.2);
        }
    }

    masterCanvas->SaveAs("PMT_PulseH_Comparison.png");

    // Cleanup
    for(int pmt=0; pmt<12; pmt++) {
        delete histOriginal[pmt];
        delete histGood[pmt];
    }
    for (int i=0; i<12; i++) delete histArea[i];
    
    delete masterCanvas;
    originalFile->Close();
    goodFile->Close();
    file->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file.root>" << endl;
        return 1;
    }
    processEvents(argv[1]);
    return 0;
}
