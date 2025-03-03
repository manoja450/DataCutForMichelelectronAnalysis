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
    // 1. FILE INPUT AND TREE SETUP
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

    // 2. BRANCH VARIABLES SETUP
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

    // 3. CALIBRATION PHASE
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

    // 4. EVENT SELECTION
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
            vector<Double_t> peakPositions;
            for (int pmt=0; pmt<12; pmt++) {
                if (pulseH[pmtChannelMap[pmt]] > 5.25) {
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
            bool allPassConditionB = true;
            for (int pmt=0; pmt<12; pmt++) {
                int ch = pmtChannelMap[pmt];
                if (pulseH[ch] <= 3 * baselineRMS[ch] || (area[ch] / pulseH[ch]) <= 1.2) {
                    allPassConditionB = false;
                    break;
                }
            }

            if (allPassConditionB) {
                vector<Double_t> peakPositions;
                for (int pmt=0; pmt<12; pmt++) {
                    if (pulseH[pmtChannelMap[pmt]] > 5.25) {
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

    // 5. SAVE RESULTS
    TFile *goodFile = new TFile(Form("./GoodEvents_%d.root", getpid()), "RECREATE");
    TTree *goodTree = tree->CloneTree(0);
    Double_t peakPosition_rms;
    goodTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<goodEvents.size(); i++) {
        tree->GetEntry(goodEvents[i]);
        peakPosition_rms = goodRMS[i];
        goodTree->Fill();
    }
    goodTree->Write();
    delete goodFile;

    TFile *badFile = new TFile(Form("./BadEvents_%d.root", getpid()), "RECREATE");
    TTree *badTree = tree->CloneTree(0);
    badTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<badEvents.size(); i++) {
        tree->GetEntry(badEvents[i]);
        peakPosition_rms = badRMS[i];
        badTree->Fill();
    }
    badTree->Write();
    delete badFile;

    // 6. CREATE COMPARISON PLOTS
    TFile *originalFile = TFile::Open(fileName);
    TTree *originalTree = (TTree*)originalFile->Get("tree");
    goodFile = TFile::Open(Form("./GoodEvents_%d.root", getpid()));
    goodTree = (TTree*)goodFile->Get("tree");

    TH1F *histOriginal[12];
    TH1F *histGood[12];

    for(int pmt=0; pmt<12; pmt++) {
        histOriginal[pmt] = new TH1F(Form("PMT%d_Original",pmt+1),
                                    Form("PMT %d;ADC Counts;Events",pmt+1), 
                                    200, -100, 900);
        histGood[pmt] = new TH1F(Form("PMT%d_Good",pmt+1),
                                Form("PMT %d;ADC Counts;Events",pmt+1), 
                                200, -100, 900);
        
        originalTree->Draw(Form("pulseH> PMT%d_Original", 
                              pmtChannelMap[pmt], pmt+1), "");
        goodTree->Draw(Form("pulseH>> PMT%d_Good", 
                          pmtChannelMap[pmt], pmt+1), "");
        
        histOriginal[pmt]->SetLineColor(kBlue);
        histGood[pmt]->SetLineColor(kRed);
        histOriginal[pmt]->SetLineWidth(1);
        histGood[pmt]->SetLineWidth(1);
    }

    TCanvas *masterCanvas = new TCanvas("MasterCanvas", 
                                       " pulseH ", 
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
            
            double ymax = std::max(hOrig->GetMaximum(), hGood->GetMaximum()) * 1.2;
            hOrig->SetMaximum(ymax);
            
            hOrig->Draw("HIST");
            hGood->Draw("HIST SAME");
            
            TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg->AddEntry(hOrig, "Original", "l");
            leg->AddEntry(hGood, "Good Events", "l");
            leg->SetBorderSize(0);
            leg->Draw();
            
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

    masterCanvas->SaveAs("PMT_ADC_Comparison.png");

    // 7. CLEANUP
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
