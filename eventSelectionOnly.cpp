//This code is designed to select data for Michel electron analysis. It applies a series of cuts to remove afterpulsing events. 
//After applying the cut it stores the good data in a new root file inside the unique directory, with an additional branch for peakPosition_rms.
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h> // For mkdir
#include <ctime>      // For timestamp
#include <cstdlib>    // For rand()

using namespace std;

// SPE fitting function
Double_t SPEfit(Double_t *x, Double_t *par) {
    Double_t term1 = par[0] * exp(-0.5 * pow((x[0]-par[1])/par[2], 2));
    Double_t term2 = par[3] * exp(-0.5 * pow((x[0]-par[4])/par[5], 2));
    Double_t term3 = par[6] * exp(-0.5 * pow((x[0]-sqrt(2)*par[4])/sqrt(2*pow(par[5],2)-pow(par[2],2)), 2));
    Double_t term4 = par[7] * exp(-0.5 * pow((x[0]-sqrt(3)*par[4])/sqrt(3*pow(par[5],2)-2*pow(par[2],2)), 2));
    return term1 + term2 + term3 + term4;
}

// Function to calculate mean and RMS of a dataset
void CalculateMeanAndRMS(const vector<Double_t> &data, Double_t &mean, Double_t &rms) {
    mean = 0.0;
    for (const auto &value : data) mean += value;
    mean /= data.size();
    
    rms = 0.0;
    for (const auto &value : data) rms += pow(value - mean, 2);
    rms = sqrt(rms / data.size());
}

// Function to generate a random directory name
string generateRandomDirName(const string &baseName) {
    // Get current timestamp
    time_t now = time(0);
    tm *ltm = localtime(&now);

    // Generate a random number
    int randomNum = rand() % 10000; // Random number between 0 and 9999

    // Create directory name
    string dirName = baseName + "_Output_" + to_string(1900 + ltm->tm_year) + 
                     to_string(1 + ltm->tm_mon) + to_string(ltm->tm_mday) + "_" +
                     to_string(ltm->tm_hour) + to_string(ltm->tm_min) + to_string(ltm->tm_sec) +
                     "_" + to_string(randomNum);
    return dirName;
}

// Main function to process events
void processEvents(const char *fileName) {
    // Extract the base name of the input file (without the .root extension)
    string inputFileName(fileName);
    size_t lastDot = inputFileName.find_last_of(".");
    string baseName = inputFileName.substr(0, lastDot);

    // Generate a unique directory name
    string outputDir = generateRandomDirName(baseName);
    cout << "Creating output directory: " << outputDir << endl;

    // Create the directory
    if (mkdir(outputDir.c_str(), 0777)) {
        cerr << "Error creating directory: " << outputDir << endl;
        return;
    }

    // Open the input ROOT file
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }
    cout << "Opened input file: " << fileName << endl;

    // Access the TTree
    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree!" << endl;
        file->Close();
        return;
    }
    cout << "Accessed TTree: tree" << endl;

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
    cout << "Processing " << nEntries << " events for calibration..." << endl;
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
    vector<Double_t> goodRMS, badRMS;

    // Loop over all events for selection
    cout << "Processing " << nEntries << " events for selection..." << endl;
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
                peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]); // Include all PMTs
            }
            if (peakPositions.size() > 0) {
                Double_t dummyMean;
                CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                if (currentRMS < 2.5) isGood = true;
            }
        } 
        else {
            // New Condition B: Pulse Height > 3 * baseline RMS and area/height > 1.2 for at least 3 PMTs
            int countConditionB = 0;
            for (int pmt=0; pmt<12; pmt++) {
                int ch = pmtChannelMap[pmt];
                if (pulseH[ch] > 3 * baselineRMS[ch] && (area[ch] / pulseH[ch]) > 1.2) {
                    countConditionB++;
                }
            }

            if (countConditionB >= 3) {
                // Calculate peak position mean and RMS
                vector<Double_t> peakPositions;
                for (int pmt=0; pmt<12; pmt++) {
                    peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]); // Include all PMTs
                }
                if (peakPositions.size() > 0) {
                    Double_t dummyMean;
                    CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                    if (currentRMS < 2.5) isGood = true;
                }
            }
        }

        // Categorize events
        if (isGood) {
            goodEvents.push_back(entry);
            goodRMS.push_back(currentRMS);
        } else {
            badEvents.push_back(entry);
            badRMS.push_back(currentRMS);
        }
    }

    // Save results
    string goodFilePath = outputDir + "/GoodEvents_" + baseName + ".root";
    cout << "Saving good events to: " << goodFilePath << endl;
    TFile *goodFile = new TFile(goodFilePath.c_str(), "RECREATE");
    if (!goodFile || goodFile->IsZombie()) {
        cerr << "Error creating output file: " << goodFilePath << endl;
        return;
    }

    TTree *goodTree = tree->CloneTree(0); // Clone the input tree structure
    Double_t peakPosition_rms;
    goodTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<goodEvents.size(); i++) {
        tree->GetEntry(goodEvents[i]); // Load event data
        peakPosition_rms = goodRMS[i]; // Add RMS value
        goodTree->Fill();              // Save event
    }
    goodTree->Write();
    goodFile->Close();
    cout << "Saved " << goodEvents.size() << " good events." << endl;

    string badFilePath = outputDir + "/BadEvents_" + baseName + ".root";
    cout << "Saving bad events to: " << badFilePath << endl;
    TFile *badFile = new TFile(badFilePath.c_str(), "RECREATE");
    if (!badFile || badFile->IsZombie()) {
        cerr << "Error creating output file: " << badFilePath << endl;
        return;
    }

    TTree *badTree = tree->CloneTree(0); // Clone the input tree structure
    badTree->Branch("peakPosition_rms", &peakPosition_rms, "peakPosition_rms/D");

    for (size_t i=0; i<badEvents.size(); i++) {
        tree->GetEntry(badEvents[i]); // Load event data
        peakPosition_rms = badRMS[i]; // Add RMS value
        badTree->Fill();              // Save event
    }
    badTree->Write();
    badFile->Close();
    cout << "Saved " << badEvents.size() << " bad events." << endl;

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
