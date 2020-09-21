#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include <fstream>
#include <cassert>
#include <string>
#include <iomanip>
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
    out << "    txt         - text file to write jet data to" << std::endl;
    exit(1);
}

int main(int argc, char const *argv[]) {
    if (argc != 4) {
        usage(std::cerr, "Incorrect number of parameters given.");
    }

    std::string  dataFile(argv[1]);
    std::string  histosFile(argv[2]);
    std::string  txtFile(argv[3]);
    std::string dataPath = "./data/root/" + dataFile;
    std::string histosPath = "./data/root/" + histosFile;
    std::string txtPath = "./data/txt/" + txtFile;
    // open data file and get trees

    TFile f(dataPath.c_str());
    TDirectoryFile* df = (TDirectoryFile*) f.Get("AK4jets");
    TTree* tree2 = (TTree*) df->Get("jetTree;3");
    TTree* tree1 = (TTree*) df->Get("genPartTree");
    TTree* tree3 = (TTree*) df->Get("genJetTree");
    TTree* tree4 = (TTree*) df->Get("pfCandTree");


    TFile h(histosPath.c_str(), "new");
    // set up histograms
    TH1* numMatchesPartonRecoHist = new TH1F("numMatchesPartonRecoHist", "Number of reco jet matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* numMatchesPartonGenHist = new TH1F("numMatchesPartonGenHist", "Number of gen jet matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* minDRPartonsRecoHist = new TH1F("minDRPartonsRecoHist", "Distribution of min dR for each parton jet (matching parton to reco)", 100, 0, 1.5);
    TH1* minDRPartonsGenHist = new TH1F("minDRPartonsGenHist", "Distribution of min dR for each parton jet (matching parton to gen)", 100, 0, 1.5);
    TH1* numPartonJetsHist = new TH1I("numPartonJetsHist", "Number of parton jets per event (Pt > 20)", 10, 0, 10);
    TH1* numRecoJetsHist = new TH1I("numRecoJetsHist", "Number of reco jets per event", 10, 0, 10);
    TH1* partonJetPtNoGenMatchHist = new TH1F("partonJetPtNoGenMatchHist", "Distribution of Pt for parton jets with no genJet match", 100, 0, 200);
    TH1* partonJetPtNoRecoMatchHist = new TH1F("partonJetPtNoRecoMatchHist", "Distribution of Pt for parton jets with no recoJet match", 200, 0, 400);
    TH1* partonJetPhiHist = new TH1F("partonJetPhiHist", "Distribution of Phi for parton jets", 100, -3.5, 3.5);
    TH1* partonJetEtaHist = new TH1F("partonJetEtaHist", "Distribution of Eta for parton jets", 100, -5, 5);
    TH1* partonJetPhiNoRecoMatchHist = new TH1F("partonJetPhiNoRecoMatchHist", "Distribution of phi for parton jets with no reco match", 100, -3.5, 3.5);
    TH1* partonJetEtaNoRecoMatchHist = new TH1F("partonJetEtaNoRecoMatchHist", "Distribution of eta for parton jets with no reco match", 100, -5, 5);
    TH1* highestPtPartonJetNoRecoMatchHist = new TH1I("highestPtPartonJetNoRecoMatchHist)", "PGDIDs of highest pt parton in each parton jet with no reco match", 33, -6, 27);
    TH1* partonNoRecoMatchPdgIdHist = new TH1I("partonNoRecoMatchPDgIdHist", "PGDIDs of partons in parton jets with no reco match", 33, -6, 27);
    TH1* numMatchesPartonPfCandHist = new TH1F("numMatchesPartonPfCandHist", "Number of pf Cand Jet matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);

    // file to write data as txt
    std::ofstream write_out("./data/txt/" + txtFile);
    assert(write_out.is_open());

    //pfCandTree variables
    std::vector<Float_t>* pfCandPx = 0;
    std::vector<Float_t>* pfCandPy = 0;
    std::vector<Float_t>* pfCandPz = 0;
    std::vector<Float_t>* pfCandE = 0;

    // genPartTree variables
    std::vector<Int_t>* pdgId = 0;
    std::vector<Float_t>* partonPx = 0;
    std::vector<Float_t>* partonPy = 0;
    std::vector<Float_t>* partonPz = 0;
    std::vector<Float_t>* partonE = 0;
    std::vector<Float_t>* Pt = 0;
    std::vector<Float_t>* Vx = 0;
    std::vector<Float_t>* Vy = 0;
    std::vector<Float_t>* Vz = 0;
    ULong64_t partonEvent = 0;

    // jetTree variables
    ULong64_t recoJetEvent = 0;
    Float_t recoJetEta = 0;
    Float_t recoJetPhi = 0;
    Float_t recoJetPt = 0;
    Int_t recoJetGenMatch = 0;

    // genJetTree variables
    std::vector<Float_t>* genJetPt = 0;
    std::vector<Float_t>* genJetEta = 0;
    std::vector<Float_t>* genJetPhi = 0;
    ULong64_t genJetEvent = 0;

    tree1->SetBranchAddress("genPartPdgId", &pdgId);
    tree1->SetBranchAddress("genPartPx", &partonPx);
    tree1->SetBranchAddress("genPartPy", &partonPy);
    tree1->SetBranchAddress("genPartPz", &partonPz);
    tree1->SetBranchAddress("genPartE", &partonE);
    tree1->SetBranchAddress("event", &partonEvent);
    tree1->SetBranchAddress("genPartPt", &Pt);
    tree1->SetBranchAddress("genPartVx", &Vx);
    tree1->SetBranchAddress("genPartVy", &Vy);
    tree1->SetBranchAddress("genPartVz", &Vz);

    tree2->SetBranchAddress("event", &recoJetEvent);
    tree2->SetBranchAddress("jetEta", &recoJetEta);
    tree2->SetBranchAddress("jetPhi", &recoJetPhi);
    tree2->SetBranchAddress("jetPt", &recoJetPt);
    tree2->SetBranchAddress("jetGenMatch", &recoJetGenMatch);


    tree3->SetBranchAddress("genJetPt", &genJetPt);
    tree3->SetBranchAddress("genJetEta", &genJetEta);
    tree3->SetBranchAddress("genJetPhi", &genJetPhi);
    tree3->SetBranchAddress("genJetEvent", &genJetEvent);

    tree4->SetBranchAddress("pfCandPx", &pfCandPx);
    tree4->SetBranchAddress("pfCandPy", &pfCandPy);
    tree4->SetBranchAddress("pfCandPz", &pfCandPz);
    tree4->SetBranchAddress("pfCandE", &pfCandE);



    int numEvents = tree1->GetEntries();

    // used to break out of while loop for getting reco jet data
    int numRecoJetsTot = tree2->GetEntries();
    // used for getting reco jet data
    int jetIndex = 0;

    // need to get first reco jet event number
    tree2->GetEntry(jetIndex);

    // loop through events
    for (size_t i=0; i < numEvents; i++) {
        std::vector<float> recoJetEtaVec;
        std::vector<float> recoJetPhiVec;
        std::vector<float> recoJetPtVec;
        std::vector<int> recoJetGenMatchVec;
        std::vector<float> partonVx;
        std::vector<float> partonVy;
        std::vector<float> partonVz;

        tree1->GetEntry(i);
        tree3->GetEntry(i);
        tree4->GetEntry(i);
        assert(genJetEvent == partonEvent);
        // save info on genJets
        for (int j=0; j < genJetPt->size(); j++) {
            write_out << 1 << " " << genJetEvent << " " << 0 << " " << (*genJetPt)[j] << " " << (*genJetEta)[j] << " " << (*genJetPhi)[j] << " " << 0 << " " << 0 << " " << 0 << "\n";
        }

        // get eta and phi from reco jets
        int numRecoJets = 0;
        while (true) {
            jetIndex++;
            tree2->GetEntry(jetIndex);
            if (recoJetEvent != partonEvent || jetIndex > numRecoJetsTot) break; //TODO: This seems a little hacky, figure out better logic
            numRecoJets++;
            recoJetEtaVec.push_back(recoJetEta);
            recoJetPhiVec.push_back(recoJetPhi);
            recoJetPtVec.push_back(recoJetPt);
            recoJetGenMatchVec.push_back(recoJetGenMatch);
            write_out << 2 << " " << recoJetEvent << " " << recoJetGenMatch << " " << recoJetPt << " " << recoJetEta << " " << recoJetPhi << " " << 0 << " " << 0 << " " << 0 << "\n";
        }

        // create vector with all parton 4-momenta
        int numPartons = pdgId->size();
        std::vector<PseudoJet> particles;
        for (size_t j=0; j < numPartons; j++) {
            particles.push_back( PseudoJet( (*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]) );
            particles[j].set_user_index(j); // set index to be able to identify constituent particles later
        }

        // cluster partons into jet with anti-kt algorithm
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R); // define jet algorithm as anti-kt with R=0.4
        ClusterSequence cs(particles, jet_def); // run the clustering

        std::vector<PseudoJet> partonJets = cs.inclusive_jets(); // get new jets

        // create vector with pf cand 4 momenta
        int numPfCands = pfCandPx->size();
        std::vector<PseudoJet> pfCands;
        for (int j=0; j < numPfCands; j++) {
            pfCands.push_back( PseudoJet((*pfCandPx)[j], (*pfCandPy)[j], (*pfCandPz)[j], (*pfCandE)[j]));
        }

        ClusterSequence pfseq(pfCands, jet_def);
        std::vector<PseudoJet> pfJets = pfseq.inclusive_jets();

        for (int j=0; j < pfJets.size(); j++) {
            write_out << std::setprecision(10) << 3 << " " << partonEvent << " " << -1 << " " << pfJets[j].pt() <<
                " " << pfJets[j].rap() << " " << pfJets[j].phi_std() << " " << 0 << " " << 0 << " " << 0 << "\n";
        }

        // constituent parton)
        for (size_t j=0; j < partonJets.size(); j++) {
            std::vector<PseudoJet> constituents = sorted_by_pt(partonJets[j].constituents());
            int index = constituents[0].user_index();
            write_out << std::setprecision(10) << 0 << " " << partonEvent << " " << -1 << " " << partonJets[j].pt() <<
                " " << partonJets[j].rap() << " " << partonJets[j].phi_std() <<  " " << (*Vx)[index]
                << " " <<(*Vy)[index] << " " << (*Vz)[index] << "\n";
        }


        // go through each parton jet and find how many reco and gen Jet
        // matches there are (and find min dR)
        int numPartonJets = 0;
        if (partonJets.size() > 0) {
            for (int j=0; j < partonJets.size(); j++) {
                if (partonJets[j].pt() > 20) {
                    numPartonJets++;
                    float partonEta = partonJets[j].rap();
                    float partonPhi = partonJets[j].phi_std();
                    partonJetPhiHist->Fill(partonPhi);
                    partonJetEtaHist->Fill(partonEta);

                    int recoMatches = 0;
                    float minDR = 10.0;
                    for (int k=0; k < recoJetEtaVec.size(); k++) {
                            float dR = deltaR(recoJetEtaVec[k], recoJetPhiVec[k], partonEta, partonPhi);
                            if (dR < 0.35) recoMatches++;
                            if (dR < minDR) minDR = dR;
                    }
                    numMatchesPartonRecoHist->Fill(recoMatches);
                    minDRPartonsRecoHist->Fill(minDR);
                    if (recoMatches == 0) {
                        partonJetPtNoRecoMatchHist->Fill(partonJets[j].pt());
                        partonJetPhiNoRecoMatchHist->Fill(partonPhi);
                        partonJetEtaNoRecoMatchHist->Fill(partonEta);
                        std::vector<PseudoJet> constituents = sorted_by_pt(partonJets[j].constituents()); // get consituent partons in the jet
                        for (int k=0; k<constituents.size(); k++) {
                            int index = constituents[k].user_index(); // get index of parton in jet with highest pt
                            int partonPdgId = (*pdgId)[index];
                            if (partonPdgId > 21 || partonPdgId < -6) {
                                std::cout << "Parton with pdgid: " << partonPdgId << " did not have reco match" << std::endl;
                            }
                            if (k == 0) highestPtPartonJetNoRecoMatchHist->Fill(partonPdgId);
                            partonNoRecoMatchPdgIdHist->Fill(partonPdgId);
                        }
                    }
                    int genMatches = 0;
                    minDR = 10.0;
                    for (int k=0; k < genJetPt->size(); k++) {
                        if ((*genJetPt)[k] > 30) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], partonEta, partonPhi);
                            if (dR < 0.35) genMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesPartonGenHist->Fill(genMatches);
                    minDRPartonsGenHist->Fill(minDR);
                    if (genMatches == 0) partonJetPtNoGenMatchHist->Fill(partonJets[j].pt());

                    int pfCandMatches = 0;
                    minDR = 10.9;
                    for (int k=0; k < pfJets.size(); k++) {
                        if (pfJets[k].pt() > 30) {
                            float dR = deltaR(pfJets[k].rap(), pfJets[k].phi_std(), partonEta, partonPhi);
                            if (dR < 0.35) pfCandMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesPartonPfCandHist->Fill(pfCandMatches++);
                }
            }
        }
        numPartonJetsHist->Fill(numPartonJets);
        numRecoJetsHist->Fill(numRecoJets);
    }

    numMatchesPartonRecoHist->Write();
    numMatchesPartonGenHist->Write();
    minDRPartonsRecoHist->Write();
    minDRPartonsGenHist->Write();
    numPartonJetsHist->Write();
    numRecoJetsHist->Write();
    partonJetPtNoGenMatchHist->Write();
    partonJetPtNoRecoMatchHist->Write();
    partonJetPhiHist->Write();
    partonJetEtaHist->Write();
    partonJetPhiNoRecoMatchHist->Write();
    partonJetEtaNoRecoMatchHist->Write();
    highestPtPartonJetNoRecoMatchHist->Write();
    partonNoRecoMatchPdgIdHist->Write();
    numMatchesPartonPfCandHist->Write();

    delete numMatchesPartonRecoHist;
    delete numMatchesPartonGenHist;
    delete minDRPartonsRecoHist;
    delete minDRPartonsGenHist;
    delete numPartonJetsHist;
    delete numRecoJetsHist;
    delete partonJetPtNoGenMatchHist;
    delete partonJetPtNoRecoMatchHist;
    delete partonJetPhiHist;
    delete partonJetEtaHist;
    delete partonJetPhiNoRecoMatchHist;
    delete partonJetEtaNoRecoMatchHist;
    delete highestPtPartonJetNoRecoMatchHist;
    delete partonNoRecoMatchPdgIdHist;
    delete numMatchesPartonPfCandHist;

    write_out.close();

    return 0;
}
