#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include <fstream>
#include <cassert>
#include <string>
#include <iomanip>
#include <set>
#include "helpers.h"

using namespace fastjet;

/**
 * Describe parameters taken and print errors if parameters are incorrectly
 * given.
 */
void usage(std::ostream &out, const char *msg) {
    out << msg << std::endl;
    out << std::endl;
    out << "    Usage:" << std::endl;
    out << "            makeJets.out data histos txt" << std::endl;
    out << "    data        - root file from JetNtuple analyzer to use" << std::endl;
    out << "    histos      - root file to write histograms to" << std::endl;
//    out << "    txt         - text file to write jet data to" << std::endl;
    exit(1);
}

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        usage(std::cerr, "Incorrect number of parameters given.");
    }

    std::string  dataFile(argv[1]);
    std::string  histosFile(argv[2]);
    //std::string  txtFile(argv[3]);
    std::string dataPath = "./data/root/" + dataFile;
    std::string histosPath = "./data/root/" + histosFile;
    //std::string txtPath = "./data/txt/" + txtFile;
    // open data file and get trees

    TFile f(dataPath.c_str());
    TDirectoryFile* df = (TDirectoryFile*) f.Get("demo");
    TTree* tree = (TTree*) df->Get("eventTree");

    std::ofstream write_out("./data/txt/partHadE.txt");
    assert(write_out.is_open());


    TFile h(histosPath.c_str(), "RECREATE");
    // set up histograms
    TH1* numMatchesPartonRecoHist = new TH1I("numMatchesPartonRecoHist", "Number of parton Jets with reco jet match (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* numMatchesPartonGenHist = new TH1I("numMatchesPartonGenHist", "Number of parton jets with gen jet match", 5, 0, 5);
    TH1* numMatchesGenPartonHist = new TH1I("numMatchesGenPartonHist", "Number of gen jets with parton jet match", 5, 0, 5);
    TH1* numMatchesGenRecoHist = new TH1I("numMatchesGenRecoHist", "Number of gen jets with reco jet match", 5, 0, 5);
    TH1* numMatchesRecoPartonHist = new TH1I("numMatchesRecoPartonHist", "Number of reco Jets with parton jet match (dR < 0.35) for each reco jet with Pt > 30GeV", 5, 0, 5);
    TH1* numMatchesRecoGenHist = new TH1I("numMatchesRecoGenHist", "Number of reco jes with gen jet match", 5, 0, 5);
    

    TH1* partonPtNoRecoMatchHist = new TH1F("partonPtNoRecoMatchHist", "Pt distribution of parton jets with no reco match", 200, 0, 200);
    TH1* recoPtNoPartonMatchHist = new TH1F("recoPtNoPartonMatchHist", "Pt distribution of reco jets with no parton match", 200, 0, 200);

    TH1* minDRRecoPartonNoMatchHist = new TH1F("minDRRecoPartonNoMatchHist", "minDR for each reco jet with no parton match", 50, 0, 5);

    TH1* emFractionRecoPartonNoMatchHist = new TH1F("emFractionRecoPartonNoMatchHist", "EM fraction of reco jets with no parton Match", 50, 0, 1);
    TH1* emFractionRecoPartonMatchHist = new TH1F("emFractionRecoPartonMatchHist", "EM fraction of reco jets with parton match", 50, 0, 1);

    TH1* partonPtHist = new TH1F("partonPtHist", "parton pt", 200, 0, 200);



    TH1* recoEtaDistributionPartonMatchHist = new TH1F("recoEtaDistributionPartonMatchHist", "Eta distribution of reco jets with a parton match", 50, -5, 5);
    TH1* recoEtaDistributionPtDiscrepencyHist = new TH1F("recoEtaDistributionPtDiscrepencyHist", "Eta distribution of reco jet with recoPt/partonPt < 0.75", 50, -5, 5);

    TH1* emFractionRecoPartonMatchPtDiscrepencyHist = new TH1F("emFractionRecoPartonMatchPtDiscrepencyHist", 
            "EM Fraction of reco jets with recoPt/partonPt < 0.75", 50, 0, 1);



    // file to write data as txt
    //std::ofstream write_out("./data/txt/" + txtFile);
    //assert(write_out.is_open());
    std::vector<Float_t>* pfJetPt = 0;
    std::vector<Float_t>* pfJetEta = 0;
    std::vector<Float_t>* pfJetPhi = 0;
    std::vector<Float_t>* pfJetE = 0;
    std::vector<Float_t>* pfJetPhotonEnergy = 0;
    std::vector<Float_t>* pfJetElectronEnergy = 0;
    std::vector<Float_t>* pfJetMuonEnergy = 0;

    std::vector<Float_t>* genJetPt = 0;
    std::vector<Float_t>* genJetEta = 0;
    std::vector<Float_t>* genJetPhi = 0;

    
    std::vector<Float_t>* partonPx = 0;
    std::vector<Float_t>* partonPy = 0;
    std::vector<Float_t>* partonPz = 0;
    std::vector<Float_t>* partonE = 0;
    std::vector<Float_t>* partonPt = 0;
    std::vector<Int_t>* partonPdgId = 0;
    std::vector<Int_t>* partonStatus = 0;

    std::vector<Float_t>* hadronPx = 0;
    std::vector<Float_t>* hadronPy = 0;
    std::vector<Float_t>* hadronPz = 0;
    std::vector<Float_t>* hadronE = 0;

    tree->SetBranchAddress("pfJetPt", &pfJetPt);
    tree->SetBranchAddress("pfJetEta", &pfJetEta);
    tree->SetBranchAddress("pfJetPhi", &pfJetPhi);
    tree->SetBranchAddress("pfJetE", &pfJetE);
    tree->SetBranchAddress("pfJetPhotonEnergy", &pfJetPhotonEnergy);
    tree->SetBranchAddress("pfJetElectronEnergy", &pfJetElectronEnergy);
    tree->SetBranchAddress("pfJetMuonEnergy", &pfJetMuonEnergy);

    tree->SetBranchAddress("genJetPt", &genJetPt);
    tree->SetBranchAddress("genJetEta", &genJetEta);
    tree->SetBranchAddress("genJetPhi", &genJetPhi);


    tree->SetBranchAddress("partonPx", &partonPx);
    tree->SetBranchAddress("partonPy", &partonPy);
    tree->SetBranchAddress("partonPz", &partonPz);
    tree->SetBranchAddress("partonE", &partonE);
    tree->SetBranchAddress("partonPt", &partonPt);
    tree->SetBranchAddress("partonPdgId", &partonPdgId); 
    tree->SetBranchAddress("partonStatus", &partonStatus);

    tree->SetBranchAddress("hadronPx", &hadronPx);
    tree->SetBranchAddress("hadronPy", &hadronPy);
    tree->SetBranchAddress("hadronPz", &hadronPz);
    tree->SetBranchAddress("hadronE", &hadronE);


    //pfCandTree variables
    int numEvents = tree->GetEntries();

    int numPartonJets = 0;
    int numPartonJetswMatch = 0;
    int numRecoJets = 0;
    int numRecoJetswMatch = 0;

    // loop through events
    for (int i=0; i < numEvents; i++) {
        //std::cout << "Event: " << i << std::endl;
        //std::cout << "=======================================" << std::endl;
        tree->GetEntry(i);

        // create vector with all parton 4-momenta
        int numPartons = partonPx->size();
        float partonPxTot = 0;
        float partonPyTot = 0;
        float partonPzTot = 0;
        float partonETot = 0;
        std::vector<PseudoJet> particles;
        for (int j=0; j < numPartons; j++) {
            if ((*partonStatus)[j] != 72) {
                particles.push_back(PseudoJet((*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]));
                partonPxTot += (*partonPx)[j];
                partonPyTot += (*partonPy)[j];
                partonPzTot += (*partonPz)[j];
                partonETot += (*partonE)[j];

            }

            partonPtHist->Fill((*partonPt)[j]);
        }

        std::cout << "Parton Px Tot: " << partonPxTot << std::endl;
        std::cout << "Parton Py Tot: " << partonPyTot << std::endl;
        std::cout << "Parton Pz Tot: " << partonPzTot << std::endl;
        std::cout << "Parton E Tot: " << partonETot << std::endl;

        int numHadrons = hadronPx->size();
        float hadronPxTot = 0;
        float hadronPyTot = 0;
        float hadronPzTot = 0;
        float hadronETot = 0;
        for (int j=0; j < numHadrons; j++) {
            hadronPxTot += (*hadronPx)[j];
            hadronPyTot += (*hadronPy)[j];
            hadronPzTot += (*hadronPz)[j];
            hadronETot += (*hadronE)[j];
        }

        std::cout << "Hadron Px Tot: " << hadronPxTot << std::endl;
        std::cout << "Hadron Py Tot: " << hadronPyTot << std::endl;
        std::cout << "Hadron Pz Tot: " << hadronPzTot << std::endl;
        std::cout << "Hadron E Tot: " << hadronETot << std::endl;

        write_out << partonETot << " " << hadronETot << std::endl;


        // cluster partons into jet with anti-kt algorithm
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);

        std::vector<PseudoJet> partonJets = cs.inclusive_jets(); // get new parton jets
        
        // go through each parton jet and find how many reco and gen jet
        // matches there are (and find min dR)
        if (partonJets.size() > 0) {
            for (size_t j=0; j < partonJets.size(); j++) {
                if (partonJets[j].pt() >= 20.0) {
                    if (i % 1000 == 0) {
                        std::cout << "Event: " << i << std::endl;
                        std::cout << "Parton Jet Pt: " << partonJets[j].pt() << std::endl;

                    }
                    numPartonJets++;

                    int genJetMatches = 0;
                    float minDR = 10.0;
                    for (size_t k=0; k < genJetPt->size(); k++) {
                        if ((*genJetPt)[k] >= 30.0) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], partonJets[j].rap(), partonJets[j].phi_std());
                            if (dR < 0.35) genJetMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesPartonGenHist->Fill(genJetMatches);
                    int pfJetMatches = 0;
                    minDR = 10.0;
                    for (size_t k=0; k < pfJetPt->size(); k++) {
                        if ((*pfJetPt)[k] >= 30.0) {
                            float dR = deltaR((*pfJetEta)[k], (*pfJetPhi)[k], partonJets[j].rap(), partonJets[j].phi_std());
                            if (dR < 0.35) pfJetMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesPartonRecoHist->Fill(pfJetMatches);
                    if (pfJetMatches == 1) numPartonJetswMatch++;
                    if (pfJetMatches == 0) partonPtNoRecoMatchHist->Fill(partonJets[j].pt());
                }
            }
        }

        // go through each gen  jet and find how many parton and reco jet
        // matches there are (and find min dR)
        if (genJetPt->size() > 0) {
            for (size_t j=0; j < genJetPt->size(); j++) {
                if ((*genJetPt)[j] >= 30.0) {
                    int partonJetMatches = 0;
                    float minDR = 10.0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if (partonJets[k].pt() >= 20.0) {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), (*genJetEta)[j], (*genJetPhi)[j]);
                            if (dR < 0.35) partonJetMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesGenPartonHist->Fill(partonJetMatches);
                    int pfJetMatches = 0;
                    minDR = 10.0;
                    for (size_t k=0; k < pfJetPt->size(); k++) {
                        if ((*pfJetPt)[k] >= 30.0) {
                            float dR = deltaR((*pfJetEta)[k], (*pfJetPhi)[k], (*genJetEta)[j], (*genJetPhi)[j]);
                            if (dR < 0.35) pfJetMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesGenRecoHist->Fill(pfJetMatches);
                }
            }
        }


        if (pfJetPt->size() > 0) {
            for (size_t j=0; j < pfJetPt->size(); j++) {
                if ((*pfJetPt)[j] > 30.0) {
                    numRecoJets++;
                    int partonJetMatches = 0;
                    float minDR = 10.0;
                    int matchIndex = -1;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if (partonJets[k].pt() >= 20.0) {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), (*pfJetEta)[j], (*pfJetPhi)[j]);
                            if (dR < 0.35) partonJetMatches++;
                            if (dR < minDR) {
                                minDR = dR;
                                matchIndex = k;
                            }
                        }
                    }
                    numMatchesRecoPartonHist->Fill(partonJetMatches);
                    if (partonJetMatches > 0) {
                        numRecoJetswMatch++;
                        float emFraction = ( (*pfJetPhotonEnergy)[j] + (*pfJetElectronEnergy)[j] + (*pfJetMuonEnergy)[j] ) / (*pfJetE)[j];
                        emFractionRecoPartonMatchHist->Fill(emFraction);
                        recoEtaDistributionPartonMatchHist->Fill((*pfJetEta)[j]);
                        float ptRatio = (*pfJetEta)[j] / partonJets[matchIndex].rap();
                        if (ptRatio < 0.75) {
                            recoEtaDistributionPtDiscrepencyHist->Fill((*pfJetEta)[j]);
                            emFractionRecoPartonMatchPtDiscrepencyHist->Fill(emFraction);
                        }
                    }
                    if (partonJetMatches == 0) {
                        recoPtNoPartonMatchHist->Fill((*pfJetPt)[j]);
                        minDRRecoPartonNoMatchHist->Fill(minDR);
                        /*
                        std::cout << "Event: " << i << " has a reco jet with no parton match." << std::endl;
                        std::cout << "  Reco Jet Pt: " << (*pfJetPt)[j] << std::endl;
                        std::cout << "  Reco Jet Eta: " << (*pfJetEta)[j] << std::endl;
                        std::cout << "  Reco Jet Phi: " << (*pfJetPhi)[j] << std::endl;
                        std::cout << "  Reco Jet Photon Energy: " << (*pfJetPhotonEnergy)[j] << std::endl;
                        std::cout << "  Reco Jet Electron Energy: " << (*pfJetElectronEnergy)[j] << std::endl;
                        std::cout << "  Reco Jet Muon Energy: " << (*pfJetMuonEnergy)[j] << std::endl;
                        std::cout << "  Reco Jet Energy: " << (*pfJetE)[j] << std::endl;
                        */
                        float emFraction = ( (*pfJetPhotonEnergy)[j] + (*pfJetElectronEnergy)[j] + (*pfJetMuonEnergy)[j] ) / (*pfJetE)[j];
                        emFractionRecoPartonNoMatchHist->Fill(emFraction);
                    }

                    int genJetMatches = 0;
                    minDR = 10.0;
                    for (size_t k=0; k < genJetPt->size(); k++) {
                        if ((*genJetPt)[k] >= 20.0) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], (*pfJetEta)[j], (*pfJetPhi)[j]);
                            if (dR < 0.35) genJetMatches++;
                        }
                    }
                    numMatchesRecoGenHist->Fill(genJetMatches);
                }
            }
        }
    }

    numMatchesPartonRecoHist->Write();
    numMatchesPartonGenHist->Write();
    numMatchesGenPartonHist->Write();
    numMatchesGenRecoHist->Write();
    numMatchesRecoPartonHist->Write();
    numMatchesRecoGenHist->Write();
    

    partonPtNoRecoMatchHist->Write();
    recoPtNoPartonMatchHist->Write();

    minDRRecoPartonNoMatchHist->Write();

    emFractionRecoPartonNoMatchHist->Write();
    emFractionRecoPartonMatchHist->Write();

    partonPtHist->Write();

    recoEtaDistributionPartonMatchHist->Write();
    recoEtaDistributionPtDiscrepencyHist->Write();

    emFractionRecoPartonMatchPtDiscrepencyHist->Write();

    delete numMatchesPartonRecoHist;
    delete numMatchesPartonGenHist;
    delete numMatchesGenPartonHist;
    delete numMatchesGenRecoHist;
    delete numMatchesRecoPartonHist;
    delete numMatchesRecoGenHist;

    delete partonPtNoRecoMatchHist;
    delete recoPtNoPartonMatchHist;

    delete minDRRecoPartonNoMatchHist;

    delete emFractionRecoPartonNoMatchHist;
    delete emFractionRecoPartonMatchHist;

    delete partonPtHist;

    delete recoEtaDistributionPartonMatchHist;
    delete recoEtaDistributionPtDiscrepencyHist;

    delete emFractionRecoPartonMatchPtDiscrepencyHist;

    write_out.close();


    return 0;
}
