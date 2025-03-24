//working .donot edit this
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// SPE fitting function
Double_t SPEfit(Double_t *x, Double_t *par) {
    Double_t term1 = par[0] * exp(-0.5 * pow((x[0]-par[1])/par[2], 2));
    Double_t term2 = par[3] * exp(-0.5 * pow((x[0]-par[4])/par[5], 2));
    Double_t term3 = par[6] * exp(-0.5 * pow((x[0]-sqrt(2)*par[4])/sqrt(2*pow(par[5],2)-pow(par[2],2)), 2));
    Double_t term4 = par[7] * exp(-0.5 * pow((x[0]-sqrt(3)*par[4])/sqrt(3*pow(par[5],2)-2*pow(par[2],2)), 2));
    return term1 + term2 + term3 + term4;
}

// Muon decay exponential function
Double_t DecayFit(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0]/par[1]);
}

// Function to perform SPE calibration
void performCalibration(TTree *tree, Double_t *mu1) {
    gErrorIgnoreLevel = kError;

    Double_t area[23], pulseH[23];
    Int_t triggerBits;

    tree->SetBranchAddress("area", area);
    tree->SetBranchAddress("pulseH", pulseH);
    tree->SetBranchAddress("triggerBits", &triggerBits);

    const int nPMTs = 12;
    int pmtChannelMap[nPMTs] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

    // SPE Calibration Phase
    TH1F *histArea[nPMTs];
    for (int i=0; i<nPMTs; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area",i+1), 
                              Form("PMT %d;ADC Counts;Events",i+1), 150, -50, 400);
    }

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);
        if (triggerBits != 16) continue;
        
        for (int pmt=0; pmt<nPMTs; pmt++) {
            histArea[pmt]->Fill(area[pmtChannelMap[pmt]]);
        }
    }

    for (int i=0; i<nPMTs; i++) {
        if (histArea[i]->GetEntries() == 0) {
            cerr << "Empty histogram for PMT " << i+1 << endl;
            continue;
        }

        TF1 *fitFunc = new TF1("fitFunc", SPEfit, -50, 400, 8);
        Double_t histMean = histArea[i]->GetMean();
        Double_t histRMS = histArea[i]->GetRMS();

        fitFunc->SetParameters(1000, histMean - histRMS, histRMS / 2,
                               1000, histMean, histRMS, 500, 500);
        histArea[i]->Fit(fitFunc, "Q0");
        mu1[i] = fitFunc->GetParameter(4);
        delete fitFunc;
    }

    // Cleanup
    for (int i=0; i<nPMTs; i++) {
        delete histArea[i];
    }
}

// Function to analyze muon/Michel events
void analyzeMuonMichel(TTree *tree, const Double_t *mu1) {
    gErrorIgnoreLevel = kError;

    Double_t area[23], pulseH[23];
    Int_t triggerBits;
    Long64_t nsTime;

    tree->SetBranchAddress("area", area);
    tree->SetBranchAddress("pulseH", pulseH);
    tree->SetBranchAddress("triggerBits", &triggerBits);
    tree->SetBranchAddress("nsTime", &nsTime);

    const int nPMTs = 12;
    int pmtChannelMap[nPMTs] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

    // Muon/Michel Analysis
    TH1F *histDeltaT = new TH1F("DeltaT", "Time difference;Time Difference(#mus);Events/0.1 #mus", 100, 1, 10);
    TH1F *histMichelSpectrum = new TH1F("MichelSpectrum", 
                                       "Michel Electron Spectrum;Photoelectrons;Events/11 Photoelectrons", 
                                       100, 100, 1000); // Adjusted x-axis

    int muonCount = 0, michelCount = 0;
    const double muonThreshold = 100;
    const double michelThreshold = 50;

    Long64_t nEntries = tree->GetEntries();
    for(Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);
        
        if(triggerBits != 34 && triggerBits != 2 && triggerBits !=32) continue;

        // Check PMT hits for muon
        int nPMTsHit = 0;
        for(int pmt=0; pmt<nPMTs; pmt++) {
            if(area[pmtChannelMap[pmt]] > 0) nPMTsHit++;
        }
        if(nPMTsHit < 3) continue;

        // Calculate energy with individual SPE factors
        double currentEnergy = 0;
        for(int pmt=0; pmt<nPMTs; pmt++) {
            if(mu1[pmt] > 0) currentEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
        }

        if (currentEnergy >= muonThreshold) {
            muonCount++;
            Long64_t muonTime = nsTime;

            for(Long64_t nextEntry=entry+1; nextEntry<nEntries; nextEntry++) {
                tree->GetEntry(nextEntry);
                double deltaT = (nsTime - muonTime)*1e-3;
                
                if(deltaT > 10) break;
                if(deltaT < 1) continue;
                
                // Check PMT hits for Michel
                int nPMTsHitMichel = 0;
                for(int pmt=0; pmt<nPMTs; pmt++) {
                    if(area[pmtChannelMap[pmt]] > 0) nPMTsHitMichel++;
                }
                if(nPMTsHitMichel < 3) continue;

                // Calculate Michel energy
                double michelEnergy = 0;
                for(int pmt=0; pmt<nPMTs; pmt++) {
                    if(mu1[pmt] > 0) michelEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
                }

                if(michelEnergy >= michelThreshold) {
                    michelCount++;
                    histDeltaT->Fill(deltaT);
                    histMichelSpectrum->Fill(michelEnergy);
                    break;
                }
            }
        }
    }

    // Time difference fit
    TCanvas *c1 = new TCanvas("c1", "Time Difference", 800, 600);
    TF1 *decayFit = new TF1("decayFit", DecayFit, 0, 10, 2);
    decayFit->SetParameters(histDeltaT->GetMaximum(), 2.2);
    histDeltaT->Fit(decayFit, "R");
    
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.85, 0.8, "NDC");
    pt->AddText(Form("#tau = %.2f #pm %.2f #mus", 
                   decayFit->GetParameter(1), decayFit->GetParError(1)));
    histDeltaT->GetListOfFunctions()->Add(pt);
    histDeltaT->Draw();
    c1->SaveAs("time_difference.png");

    // Michel spectrum
    TCanvas *c2 = new TCanvas("c2", "Michel Electron Spectrum", 800, 600);
    histMichelSpectrum->Draw();
    c2->SaveAs("michel_spectrum.png");

    // Cleanup
    delete histDeltaT;
    delete histMichelSpectrum;
    delete c1;
    delete c2;
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file.root>" << endl;
        return 1;
    }

    const int nPMTs = 12;
    Double_t mu1[nPMTs] = {0};

    TFile *inputFile = TFile::Open(argv[1]);
    if (!inputFile) return 1;

    TTree *tree = (TTree*)inputFile->Get("tree");
    if (!tree) {
        inputFile->Close();
        return 1;
    }

    performCalibration(tree, mu1);
    analyzeMuonMichel(tree, mu1);

    inputFile->Close();
    return 0;
}
