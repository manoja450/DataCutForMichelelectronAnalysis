//This code apply first cut criteria and saves events in Good and Bad root files. 
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h> // For mkdir

using namespace std;

// SPE fitting function
Double_t SPEfit(Double_t *x, Double_t *par) {
    Double_t term1 = par[0] * exp(-0.5 * pow((x[0]-par[1])/par[2], 2));
    Double_t term2 = par[3] * exp(-0.5 * pow((x[0]-par[4])/par[5], 2));
    Double_t term3 = par[6] * exp(-0.5 * pow((x[0]-sqrt(2)*par[4])/sqrt(2*pow(par[5],2)-pow(par[2],2)), 2));
    Double_t term4 = par[7] * exp(-0.5 * pow((x[0]-sqrt(3)*par[4])/sqrt(3*pow(par[5],2)-2*pow(par[2],2)), 2));
    return term1 + term2 + term3 + term4;
}

// Main function to process events
void processEvents(const char *fileName) {
    // Extract the base name of the input file (without the .root extension)
    string inputFileName(fileName);
    size_t lastDot = inputFileName.find_last_of(".");
    string baseName = inputFileName.substr(0, lastDot);

    // Create a directory for output files
    string outputDir = baseName + "_Output";
    if (mkdir(outputDir.c_str(), 0777) != 0 && errno != EEXIST) {
        cerr << "Error creating directory: " << outputDir << endl;
        return;
    }

    // Open the input ROOT file
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

    // Access the TTree
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
    
    // Initialize histograms for each PMT
    for (int i=0; i<12; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area",i+1), 
                              Form("PMT %d;ADC Counts;Events",i+1), 150, -50, 400);
    }

    // Loop over all events for calibration
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);
        if (triggerBits != 16) continue; // Apply trigger condition
        
        // Fill histograms with area values
        for (int pmt=0; pmt<12; pmt++) {
            histArea[pmt]->Fill(area[pmtChannelMap[pmt]]);
        }
    }

    // Fit histograms to extract 1PE peak positions
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
        histArea[i]->Fit(fitFunc, "Q0"); // Perform the fit
        mu1[i] = fitFunc->GetParameter(4); // Extract 1PE peak position
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

    // Loop over all events for selection
    for (Long64_t entry=0; entry<nEntries; entry++) {
        tree->GetEntry(entry);

        // Check if at least 3 PMTs have pulse height > 2PE
        int countAbove2PE = 0;
        for (int pmt=0; pmt<12; pmt++) {
            if (pulseH[pmtChannelMap[pmt]] > 2 * mu1[pmt]) {
                countAbove2PE++;
            }
        }

        // Categorize events
        if (countAbove2PE >= 3) {
            goodEvents.push_back(entry); // Good event
        } else {
            badEvents.push_back(entry); // Bad event
        }
    }

    // Save results
    string goodFilePath = outputDir + "/GoodEvents_" + baseName + ".root";
    TFile *goodFile = new TFile(goodFilePath.c_str(), "RECREATE");
    TTree *goodTree = tree->CloneTree(0); // Clone the input tree structure

    for (size_t i=0; i<goodEvents.size(); i++) {
        tree->GetEntry(goodEvents[i]); // Load event data
        goodTree->Fill();              // Save event
    }
    goodTree->Write();
    delete goodFile;

    string badFilePath = outputDir + "/BadEvents_" + baseName + ".root";
    TFile *badFile = new TFile(badFilePath.c_str(), "RECREATE");
    TTree *badTree = tree->CloneTree(0); // Clone the input tree structure

    for (size_t i=0; i<badEvents.size(); i++) {
        tree->GetEntry(badEvents[i]); // Load event data
        badTree->Fill();              // Save event
    }
    badTree->Write();
    delete badFile;

    // Cleanup
    for (int i=0; i<12; i++) delete histArea[i];
    file->Close();
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file.root>" << endl;
        return 1;
    }
    processEvents(argv[1]);
    return 0;
}
