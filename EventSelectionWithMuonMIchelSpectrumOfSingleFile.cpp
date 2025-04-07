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

// 4-Gaussian SPE fitting function
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

void performCalibration(const string &calibFileName, Double_t *mu1, Double_t *mu1_err) {
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
    Long64_t nLEDFlashes[12] = {0};
    
    for (int i=0; i<12; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area",i+1), 
                             Form("PMT %d;ADC Counts;Events",i+1), 150, -50, 400);
    }

    Int_t triggerBits;
    Double_t area[23];
    calibTree->SetBranchAddress("triggerBits", &triggerBits);
    calibTree->SetBranchAddress("area", area);

    Long64_t nEntries = calibTree->GetEntries();
    cout << "Processing " << nEntries << " calibration events (triggerBits=16)..." << endl;
    
    for (Long64_t entry=0; entry<nEntries; entry++) {
        calibTree->GetEntry(entry);
        if (triggerBits != 16) continue;
        
        for (int pmt=0; pmt<12; pmt++) {
            histArea[pmt]->Fill(area[pmtChannelMap[pmt]]);
            nLEDFlashes[pmt]++;
        }
    }

    for (int i=0; i<12; i++) {
        if (histArea[i]->GetEntries() < 1000) {
            cerr << "Warning: Insufficient calibration data for PMT " << i+1 << endl;
            mu1[i] = 0;
            mu1_err[i] = 0;
            continue;
        }

        TF1 *fitFunc = new TF1("fitFunc", SPEfit, -50, 400, 8);
        Double_t histMean = histArea[i]->GetMean();
        Double_t histRMS = histArea[i]->GetRMS();

        fitFunc->SetParameters(1000, histMean-histRMS, histRMS/2,
                             1000, histMean, histRMS,
                             500, 200);
        
        histArea[i]->Fit(fitFunc, "Q0", "", -50, 400);
        mu1[i] = fitFunc->GetParameter(4); // 1PE mean position
        Double_t sigma_mu1 = fitFunc->GetParError(4); // Fit error on 1PE mean
        Double_t sigma1 = fitFunc->GetParameter(5); // Width of SPE Gaussian
        
        // Calculate calibration error: σ_cal = √(σ²_μ₁ + (σ₁/√n)²)
        mu1_err[i] = sqrt(pow(sigma_mu1, 2) + pow(sigma1/sqrt(nLEDFlashes[i]), 2));
        
        delete fitFunc;
        delete histArea[i];
    }

    cout << "\nCALIBRATION RESULTS (1PE peak positions with errors):\n";
    cout << "PMT#  HardwareCh  mu1 [ADC]  Error [ADC]  N_events\n";
    cout << "--------------------------------------------------\n";
    for (int i=0; i<12; i++) {
        printf("PMT%02d     %2d       %6.2f      %5.2f      %6lld\n", 
              i+1, pmtChannelMap[i], mu1[i], mu1_err[i], nLEDFlashes[i]);
    }
    cout << endl;

    calibFile->Close();
}

bool passesMainEventSelection(const Double_t *pulseH, const Double_t *baselineRMS, 
                            const Double_t *area, const Int_t *peakPosition,
                            const Double_t *mu1, const int *pmtChannelMap) {
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
            peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
        }
        Double_t dummyMean, currentRMS;
        CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
        if (currentRMS < 2.5) return true;
    } 
    else {
        // Condition B: Pulse Height > 3 * baseline RMS and area/height > 1.2
        int countConditionB = 0;
        for (int pmt=0; pmt<12; pmt++) {
            int ch = pmtChannelMap[pmt];
            if (pulseH[ch] > 3 * baselineRMS[ch] && (area[ch] / pulseH[ch]) > 1.2) {
                countConditionB++;
            }
        }

        if (countConditionB >= 3) {
            vector<Double_t> peakPositions;
            for (int pmt=0; pmt<12; pmt++) {
                peakPositions.push_back(peakPosition[pmtChannelMap[pmt]]);
            }
            Double_t dummyMean, currentRMS;
            CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
            if (currentRMS < 2.5) return true;
        }
    }
    return false;
}

void analyzeMuonMichel(TChain *analysisChain, const Double_t *mu1, const string &outputDir) {
    gErrorIgnoreLevel = kError;

    // Variables to read from the tree
    Double_t area[23], pulseH[23], baselineRMS[23];
    Int_t triggerBits, peakPosition[23];
    Long64_t nsTime;

    analysisChain->SetBranchAddress("triggerBits", &triggerBits);
    analysisChain->SetBranchAddress("area", area);
    analysisChain->SetBranchAddress("pulseH", pulseH);
    analysisChain->SetBranchAddress("baselineRMS", baselineRMS);
    analysisChain->SetBranchAddress("peakPosition", &peakPosition);
    analysisChain->SetBranchAddress("nsTime", &nsTime);

    // PMT and SiPM configuration
    const int nPMTs = 12;
    int pmtChannelMap[nPMTs] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};
    
    // SiPM configuration as per requirements
    const int nSiPMs = 10;
    int sipmChannelMap[nSiPMs] = {12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
    double sipmThresholds[nSiPMs] = {800, 800, 1100, 1200, 550, 600, 650, 450, 600, 650};
    double pmtThresholds[nPMTs] = {4800, 6000, 5000, 6000, 6000, 4700, 4500, 3000, 2000, 5000, 4500, 4800};

    mkdir(outputDir.c_str(), 0777);

    // Create histograms
    TH1F *histDeltaT = new TH1F("DeltaT", "Time difference between Muon and Michel;Time difference (#mus);Events/0.1 #mus", 90, 0, 10);
    TH1F *histMichelSpectrum = new TH1F("MichelSpectrum", "Michel Electron Energy Spectrum;Photoelectrons;Events", 100, 0, 1000);
    TH1F *histMuonPMTHits = new TH1F("MuonPMTHits", "PMT Multiplicity for Muon Events;Number of PMTs hit;Events", 13, 0, 13);
    TH1F *histSiPMMultiplicity = new TH1F("SiPMMultiplicity", "SiPM Multiplicity for Muon Events;Number of SiPMs hit;Events", 11, 0, 11);
    TH1F *histTriggerBits = new TH1F("TriggerBits", "Trigger Bits Distribution;Trigger Bits;Events", 64, 0, 64);

    int muonCount = 0, michelCount = 0;
    int totalTrigger34Events = 0;

    Long64_t nEntries = analysisChain->GetEntries();
    cout << "Analyzing " << nEntries << " events for muon/Michel analysis..." << endl;

    for(Long64_t entry=0; entry<nEntries; entry++) {
        analysisChain->GetEntry(entry);
        
        // Fill trigger bits histogram for monitoring
        histTriggerBits->Fill(triggerBits);
        
        // =============================================
        // ONLY PROCESS EVENTS WITH triggerBits == 34 or 2 
        // =============================================
        if (triggerBits != 34 && triggerBits != 2) continue;
        totalTrigger34Events++;
        
        // Apply main event selection (only on triggerBits=34 or 2 events)
        if (!passesMainEventSelection(pulseH, baselineRMS, area, peakPosition, mu1, pmtChannelMap)) {
            continue;
        }

        // =============================================
        // (i) Muon-like signal selection with SiPMs
        // =============================================
        
        // Count number of SiPMs above threshold (shower rejection cut)
        int sipmHitCount = 0;
        for (int i = 0; i < nSiPMs; i++) {
            int ch = sipmChannelMap[i];
            if (area[ch] >= sipmThresholds[i]) {
                sipmHitCount++;
            }
        }
        
        // Reject events with >2 hit SiPMs (multiplicity ≤ 2)
        if (sipmHitCount > 2) continue;
        histSiPMMultiplicity->Fill(sipmHitCount);

        // =============================================
        // (ii) Cherenkov signal requirement (PMTs)
        // =============================================
        
        // Count number of PMTs above threshold (must be >90% = at least 11/12 PMTs)
        int pmtHitCount = 0;
        for (int pmt = 0; pmt < nPMTs; pmt++) {
            int ch = pmtChannelMap[pmt];
            if (area[ch] >= pmtThresholds[pmt]) {
                pmtHitCount++;
            }
        }
        
        if (pmtHitCount < 11) continue;
        histMuonPMTHits->Fill(pmtHitCount);

        // =============================================
        // Both conditions satisfied - Cosmic Muon Candidate
        // =============================================
        muonCount++;
        Long64_t muonTime = nsTime;
        
        // =============================================
        // (iii) Michel electron search in subsequent events
        // =============================================
        for(Long64_t nextEntry = entry + 1; nextEntry < nEntries; nextEntry++) {
            analysisChain->GetEntry(nextEntry);
            
            // Only consider subsequent events with triggerBits==34 or 2
            if (triggerBits != 34 && triggerBits != 2) continue;
            
            double deltaT = (nsTime - muonTime) * 1e-3; // Convert to µs

            // Only look within 10 µs window
            if (deltaT > 10) break;
            // Skip events too close in time (2 µs minimum)
            if (deltaT < 2) continue;

            // Apply main event selection to Michel candidate
            if (!passesMainEventSelection(pulseH, baselineRMS, area, peakPosition, mu1, pmtChannelMap)) {
                continue;
            }

            // Michel condition: at least 11 PMTs with >=2 PE
            int michelPMTCount = 0;
            double michelEnergy = 0.0;
            for (int pmt = 0; pmt < nPMTs; pmt++) {
                double pe = area[pmtChannelMap[pmt]] / mu1[pmt];
                if (pe >= 2) {
                    michelPMTCount++;
                    michelEnergy += pe;
                }
            }

            if (michelPMTCount >= 11) {
                michelCount++;
                histDeltaT->Fill(deltaT);
                histMichelSpectrum->Fill(michelEnergy);
                break; // Take first valid Michel in the time window
            }
        }
    }

    cout << "\nANALYSIS SUMMARY:" << endl;
    cout << "Total events processed: " << nEntries << endl;
    cout << "Events with triggerBits==34: " << totalTrigger34Events << endl;
    cout << "Muon candidates found: " << muonCount << endl;
    cout << "Michel electrons found: " << michelCount << endl;

    // =============================================
    // Save results and create plots
    // =============================================
    
    // 1. Time difference between muon and Michel
    TCanvas *c1 = new TCanvas("c1", "Time Difference", 800, 600);
    TF1 *decayFit = new TF1("decayFit", DecayFit, 2, 10, 2);
    decayFit->SetParameters(histDeltaT->GetMaximum(), 2.2);
    decayFit->SetParLimits(0, 0, histDeltaT->GetMaximum() * 2);
    decayFit->SetParLimits(2, 0.1, 10);
    
    histDeltaT->Fit(decayFit, "L", "", 1, 10);
    
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.85, 0.8, "NDC");
    pt->AddText(Form("#tau = %.2f #pm %.2f #mus", 
                    decayFit->GetParameter(1), decayFit->GetParError(1)));
    pt->AddText(Form("N_{Michel} = %d", michelCount));
    histDeltaT->GetListOfFunctions()->Add(pt);
    
    histDeltaT->Draw();
    c1->SaveAs((outputDir + "/time_difference.png").c_str());

    // 2. Michel electron energy spectrum
    TCanvas *c2 = new TCanvas("c2", "Michel Spectrum", 800, 600);
    histMichelSpectrum->Draw();
    c2->SaveAs((outputDir + "/michel_spectrum.png").c_str());

    // 3. PMT multiplicity for muon events
    TCanvas *c3 = new TCanvas("c3", "PMT Multiplicity", 800, 600);
    histMuonPMTHits->Draw();
    c3->SaveAs((outputDir + "/muon_pmt_hits.png").c_str());

    // 4. SiPM multiplicity for muon events
    TCanvas *c4 = new TCanvas("c4", "SiPM Multiplicity", 800, 600);
    histSiPMMultiplicity->Draw();
    c4->SaveAs((outputDir + "/sipm_multiplicity.png").c_str());

    // 5. Trigger bits distribution
    TCanvas *c5 = new TCanvas("c5", "Trigger Bits", 800, 600);
    histTriggerBits->Draw();
    c5->SaveAs((outputDir + "/trigger_bits.png").c_str());

    // Save histograms to root file
    TFile *outFile = new TFile((outputDir + "/results.root").c_str(), "RECREATE");
    histDeltaT->Write();
    histMichelSpectrum->Write();
    histMuonPMTHits->Write();
    histSiPMMultiplicity->Write();
    histTriggerBits->Write();
    outFile->Close();

    // Cleanup
    delete histDeltaT;
    delete histMichelSpectrum;
    delete histMuonPMTHits;
    delete histSiPMMultiplicity;
    delete histTriggerBits;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete outFile;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file1.root> [input_file2.root ...]" << endl;
        cerr << "Note: First file used for calibration (triggerBits=16 only), all files for analysis" << endl;
        return 1;
    }

    string outputDir = generateRandomDirName();
    cout << "Creating output directory: " << outputDir << endl;
    mkdir(outputDir.c_str(), 0777);

    // Perform calibration using first file (triggerBits=16 only)
    Double_t mu1[12] = {0};
    Double_t mu1_err[12] = {0};
    performCalibration(argv[1], mu1, mu1_err);

    // Create chain for analysis using all files
    TChain *analysisChain = new TChain("tree");
    for (int i = 1; i < argc; i++) {
        analysisChain->Add(argv[i]);
        cout << "Added file for analysis: " << argv[i] << endl;
    }

    // Analyze only events with triggerBits==34
    analyzeMuonMichel(analysisChain, mu1, outputDir);
    delete analysisChain;

    cout << "Analysis complete. Results saved in " << outputDir << endl;
    return 0;
}
