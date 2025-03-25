#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TChain.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>
#include <ctime>
#include <cstdlib>

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

void CalculateMeanAndRMS(const vector<Double_t> &data, Double_t &mean, Double_t &rms) {
    mean = 0.0;
    for (const auto &value : data) mean += value;
    mean /= data.size();
    
    rms = 0.0;
    for (const auto &value : data) rms += pow(value - mean, 2);
    rms = sqrt(rms / data.size());
}

string generateRandomDirName() {
    time_t now = time(0);
    tm *ltm = localtime(&now);
    int randomNum = rand() % 10000;
    
    string dirName = "MichelAnalysis_" + to_string(1900 + ltm->tm_year) + 
                   to_string(1 + ltm->tm_mon) + to_string(ltm->tm_mday) + "_" +
                   to_string(ltm->tm_hour) + to_string(ltm->tm_min) + to_string(ltm->tm_sec) +
                   "_" + to_string(randomNum);
    return dirName;
}

void performCalibration(const string &calibFileName, Double_t *mu1) {
    TFile *calibFile = TFile::Open(calibFileName.c_str());
    if (!calibFile || calibFile->IsZombie()) {
        cerr << "Error opening calibration file: " << calibFileName << endl;
        exit(1);
    }

    TTree *calibTree = (TTree*)calibFile->Get("tree");
    if (!calibTree) {
        cerr << "Error accessing tree in calibration file" << endl;
        calibFile->Close();
        exit(1);
    }

    TH1F *histArea[12];
    int pmtChannelMap[12] = {0,10,7,2,6,3,8,9,11,4,5,1};
    
    for (int i=0; i<12; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area",i+1), 
                             Form("PMT %d;ADC Counts;Events",i+1), 150, -50, 400);
    }

    // Variables for calibration
    Int_t triggerBits;
    Double_t area[23];
    calibTree->SetBranchAddress("triggerBits", &triggerBits);
    calibTree->SetBranchAddress("area", area);

    Long64_t nEntries = calibTree->GetEntries();
    cout << "Processing " << nEntries << " events from " << calibFileName << " for calibration..." << endl;
    
    for (Long64_t entry=0; entry<nEntries; entry++) {
        calibTree->GetEntry(entry);
        if (triggerBits != 16) continue;
        
        for (int pmt=0; pmt<12; pmt++) {
            histArea[pmt]->Fill(area[pmtChannelMap[pmt]]);
        }
    }

    for (int i=0; i<12; i++) {
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
        delete histArea[i];
    }

    cout << "\nCALIBRATION RESULTS (1PE peak positions):\n";
    cout << "PMT#  HardwareCh  mu1 [ADC]\n";
    cout << "---------------------------\n";
    for (int i=0; i<12; i++) {
        printf("PMT%02d     %2d       %6.2f\n", i+1, pmtChannelMap[i], mu1[i]);
    }
    cout << endl;

    calibFile->Close();
}

void analyzeMuonMichel(TChain *analysisChain, const Double_t *mu1, const string &outputDir) {
    gErrorIgnoreLevel = kError;

    // Variables to read from the tree
    Double_t area[23], pulseH[23], baselineRMS[23];
    Int_t triggerBits, peakPosition[23];
    Long64_t nsTime;

    // Set branch addresses
    analysisChain->SetBranchAddress("triggerBits", &triggerBits);
    analysisChain->SetBranchAddress("area", area);
    analysisChain->SetBranchAddress("pulseH", pulseH);
    analysisChain->SetBranchAddress("baselineRMS", baselineRMS);
    analysisChain->SetBranchAddress("peakPosition", peakPosition);
    analysisChain->SetBranchAddress("nsTime", &nsTime);

    const int nPMTs = 12;
    int pmtChannelMap[nPMTs] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

    // Create output directory
    mkdir(outputDir.c_str(), 0777);

    TH1F *histDeltaT = new TH1F("DeltaT", "Time difference;(Time difference #mus);Events/0.1 #mus", 45, 1, 10);
    TH1F *histMichelSpectrum = new TH1F("MichelSpectrum", "Michel Electron Spectrum;Photoelectrons;Events", 90, 100, 1000);

    int muonCount = 0, michelCount = 0;
    const double muonThreshold = 50;
    const double michelThreshold = 50;

    Long64_t nEntries = analysisChain->GetEntries();
    cout << "Analyzing " << nEntries << " events for muon/Michel decays..." << endl;

    for(Long64_t entry=0; entry<nEntries; entry++) {
        analysisChain->GetEntry(entry);
        
        // Event selection cuts
        bool isGood = false;
        Double_t currentRMS;

        // Condition A: Pulse Height > 2 p.e. for at least 3 PMTs
        int countAbove2PE = 0;
        for (int pmt=0; pmt<nPMTs; pmt++) {
            if (pulseH[pmtChannelMap[pmt]] > 2 * mu1[pmt]) {
                countAbove2PE++;
            }
        }

        if (countAbove2PE >= 3) {
            vector<Double_t> peakPositions;
            for (int pmt=0; pmt<nPMTs; pmt++) {
                peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
            }
            if (peakPositions.size() > 0) {
                Double_t dummyMean;
                CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                if (currentRMS < 2.5) isGood = true;
            }
        } 
        else {
            // Condition B: Pulse Height > 3 * baseline RMS and area/height > 1.2
            int countConditionB = 0;
            for (int pmt=0; pmt<nPMTs; pmt++) {
                int ch = pmtChannelMap[pmt];
                if (pulseH[ch] > 3 * baselineRMS[ch] && (area[ch] / pulseH[ch]) > 1.2) {
                    countConditionB++;
                }
            }

            if (countConditionB >= 3) {
                vector<Double_t> peakPositions;
                for (int pmt=0; pmt<nPMTs; pmt++) {
                    peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
                }
                if (peakPositions.size() > 0) {
                    Double_t dummyMean;
                    CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
                    if (currentRMS < 2.5) isGood = true;
                }
            }
        }

        if (!isGood) continue;

        // Muon/Michel analysis
        if(triggerBits != 34 && triggerBits != 2) continue;

        int nPMTsHit = 0;
        for(int pmt=0; pmt<nPMTs; pmt++) {
            if(area[pmtChannelMap[pmt]] > 500) nPMTsHit++;
        }
        if(nPMTsHit < 3) continue;

        double totalEnergy = 0;
        for(int pmt=0; pmt<nPMTs; pmt++) {
            if(area[pmtChannelMap[pmt]] > 500 && mu1[pmt] > 0) {
                totalEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
            }
        }

        if (totalEnergy >= muonThreshold) {
            muonCount++;
            Long64_t muonTime = nsTime;
            
            for(Long64_t nextEntry=entry+1; nextEntry<nEntries; nextEntry++) {
                analysisChain->GetEntry(nextEntry);
                double deltaT = (nsTime - muonTime)*1e-3;
                
                if(deltaT > 10) break;
                if(deltaT < 1) continue;
                
                // Michel selection cuts
                bool isMichelGood = false;
                Double_t michelRMS;

                int michelCountAbove2PE = 0;
                for (int pmt=0; pmt<nPMTs; pmt++) {
                    if (pulseH[pmtChannelMap[pmt]] > 2 * mu1[pmt]) {
                        michelCountAbove2PE++;
                    }
                }

                if (michelCountAbove2PE >= 3) {
                    vector<Double_t> michelPeakPositions;
                    for (int pmt=0; pmt<nPMTs; pmt++) {
                        michelPeakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
                    }
                    if (michelPeakPositions.size() > 0) {
                        Double_t dummyMean;
                        CalculateMeanAndRMS(michelPeakPositions, dummyMean, michelRMS);
                        if (michelRMS < 2.5) isMichelGood = true;
                    }
                } 
                else {
                    int michelCountConditionB = 0;
                    for (int pmt=0; pmt<nPMTs; pmt++) {
                        int ch = pmtChannelMap[pmt];
                        if (pulseH[ch] > 3 * baselineRMS[ch] && (area[ch] / pulseH[ch]) > 1.2) {
                            michelCountConditionB++;
                        }
                    }

                    if (michelCountConditionB >= 3) {
                        vector<Double_t> michelPeakPositions;
                        for (int pmt=0; pmt<nPMTs; pmt++) {
                            michelPeakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
                        }
                        if (michelPeakPositions.size() > 0) {
                            Double_t dummyMean;
                            CalculateMeanAndRMS(michelPeakPositions, dummyMean, michelRMS);
                            if (michelRMS < 2.5) isMichelGood = true;
                        }
                    }
                }

                if (!isMichelGood) continue;
                
                int nPMTsHitMichel = 0;
                for(int pmt=0; pmt<nPMTs; pmt++) {
                    if(area[pmtChannelMap[pmt]] > 500) nPMTsHitMichel++;
                }
                if(nPMTsHitMichel < 3) continue;
                
                double michelEnergy = 0;
                for(int pmt=0; pmt<nPMTs; pmt++) {
                    if(area[pmtChannelMap[pmt]] > 500 && mu1[pmt] > 0) {
                        michelEnergy += area[pmtChannelMap[pmt]] / mu1[pmt];
                    }
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

    cout << "Found " << muonCount << " muon candidates and " 
         << michelCount << " Michel electrons." << endl;

    // Time difference fit and plot
    TCanvas *c1 = new TCanvas("c1", "Time Difference", 800, 600);
    TF1 *decayFit = new TF1("decayFit", DecayFit, 1, 10, 2);
    decayFit->SetParameters(histDeltaT->GetMaximum(), 2.2);
    decayFit->SetParLimits(0, 0, histDeltaT->GetMaximum() * 2);
    decayFit->SetParLimits(1, 0.1, 10);
    
    histDeltaT->Fit(decayFit, "L", "", 1, 10);
    decayFit->SetRange(2, 10);
    
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.85, 0.8, "NDC");
    pt->AddText(Form("#tau = %.2f #pm %.2f #mus", 
                    decayFit->GetParameter(1), decayFit->GetParError(1)));
    histDeltaT->GetListOfFunctions()->Add(pt);
    
    histDeltaT->Draw();
    decayFit->Draw("same");
    c1->SaveAs((outputDir + "/time_difference.png").c_str());

    // Michel spectrum plot
    TCanvas *c2 = new TCanvas("c2", "Michel Electron Spectrum", 800, 600);
    histMichelSpectrum->Draw();
    c2->SaveAs((outputDir + "/michel_spectrum.png").c_str());

    // Cleanup
    delete histDeltaT;
    delete histMichelSpectrum;
    delete c1;
    delete c2;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file1.root> [input_file2.root ...]" << endl;
        cerr << "Note: First file will be used for calibration, all files for analysis" << endl;
        return 1;
    }

    string outputDir = generateRandomDirName();
    cout << "Creating output directory: " << outputDir << endl;
    mkdir(outputDir.c_str(), 0777);

    // Perform calibration using only the first file
    Double_t mu1[12] = {0};
    performCalibration(argv[1], mu1);

    // Create chain for analysis using all files
    TChain *analysisChain = new TChain("tree");
    for (int i = 1; i < argc; i++) {
        analysisChain->Add(argv[i]);
        cout << "Added file for analysis: " << argv[i] << endl;
    }

    // Perform analysis using the calibration from first file
    analyzeMuonMichel(analysisChain, mu1, outputDir);
    delete analysisChain;

    cout << "Analysis complete. Results saved in " << outputDir << endl;
    return 0;
}
