//This code is designed to select data for Michel electron analysis. It applies a series of cuts to remove afterpulsing events. 
//After applying the cut it stores the good data in a new root file, with an additional branch for peakPosition_rms.
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
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

    Short_t adcVal[23][45];
    Double_t area[23], pulseH[23], baselineRMS[23];
    Int_t peakPosition[23], triggerBits;

    tree->SetBranchAddress("adcVal", adcVal);
    tree->SetBranchAddress("area", area);
    tree->SetBranchAddress("pulseH", pulseH);
    tree->SetBranchAddress("peakPosition", peakPosition);
    tree->SetBranchAddress("baselineRMS", baselineRMS);
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

        // Condition A Check
        bool allAbove2PE = true;
        for (int pmt=0; pmt<12; pmt++) {
            if (pulseH[pmtChannelMap[pmt]] <= 2*mu1[pmt]) {
                allAbove2PE = false;
                break;
            }
        }

        vector<Double_t> peakPositions;
        if (allAbove2PE) {
            for (int pmt=0; pmt<12; pmt++) {
                peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
            }
            Double_t dummyMean;
            CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
            if (currentRMS < 2.5) isGood = true;
        } 
        else {  // Condition B Check
            bool allPassConditionB = true;
            for (int pmt=0; pmt<12; pmt++) {
                int ch = pmtChannelMap[pmt];
                if (pulseH[ch] <= 3*baselineRMS[ch] || 
                   (area[ch]/pulseH[ch]) <= 1.0) {
                    allPassConditionB = false;
                    break;
                }
            }

            if (allPassConditionB) {
                peakPositions.clear();
                for (int pmt=0; pmt<12; pmt++) {
                    peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
                }
                Double_t dummyMean;
                CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                if (currentRMS < 2.5) isGood = true;
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
    TTree *goodTree = tree->CloneTree(0);
    Double_t rmsValue;
    goodTree->Branch("peakPosition_rms", &rmsValue, "peakPosition_rms/D");
    for (size_t i=0; i<goodEvents.size(); i++) {
        tree->GetEntry(goodEvents[i]);
        rmsValue = goodRMS[i];
        goodTree->Fill();
    }
    goodTree->Write();
    delete goodFile;

    TFile *badFile = new TFile(Form("./BadEvents_%d.root", getpid()), "RECREATE");
    TTree *badTree = tree->CloneTree(0);
    badTree->Branch("peakPosition_rms", &rmsValue, "peakPosition_rms/D");
    for (size_t i=0; i<badEvents.size(); i++) {
        tree->GetEntry(badEvents[i]);
        rmsValue = badRMS[i];
        badTree->Fill();
    }
    badTree->Write();
    delete badFile;

    // Cleanup
    for (int i=0; i<12; i++) delete histArea[i];
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
       

  
   
