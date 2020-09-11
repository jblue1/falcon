#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include <cmath>
#include <fstream>
#include <cassert>

using namespace fastjet;

float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return sqrt( pow(eta2 - eta1, 2) + pow(phi2 - phi1, 2));
}

int main(int argc, char const *argv[]) {
    // open data file and get trees
    TFile f("./data/JetNtuple_RunIISummer16_13TeV_MC.root");
    TDirectoryFile* df = (TDirectoryFile*) f.Get("AK4jets");
    TTree* tree1 = (TTree*) df->Get("genPartTree");
    TTree* tree2 = (TTree*) df->Get("jetTree;3");
    // set up histagrams
    TH1* numPartonJetsHist = new TH1I("numPartonJetsHist", "Distribution of parton jets per event", 31, 0, 30);
    TH1* numPartonJetsCutHist = new TH1I("numPartonJetsCutHist", "Distribution of parton jets with Pt > 20GeV per event", 31, 0, 30);
    TH1* numRecoJetsHist = new TH1I("numRecoJetsHist", "Distribution of reco jets per event", 31, 0, 30);
    TH1* minDeltaRPartonRecoHist = new TH1F("minDeltaRPartonRecoHist", "Distribution of min delta R for each reco jet", 100, 0.0, 10.0);
    TH1* minDeltaRPartonRecoHistZoom = new TH1F("minDeltaRPartonRecoHistZoom", "Distribution of min delta R for each reco jet", 100, 0.0, 2.0);
    TH1* partonJetPtCutHist = new TH1F("partonJetPtCutHist", "Parton jet Pt distribution (Pt > 20GeV)", 201, 0.0, 1000.0);
    TH1* recoJetPtHist = new TH1F("recoJetPtHist", "Reco jet Pt distribution", 200, 0.0, 1000.0);
    TH1* partonJetPtHist = new TH1F("partonJetPtHist", "Parton jet Pt disribution", 201, 0.0, 1000.0);
    TH1* minDeltaRPartonGenHist = new TH1F("minDeltaRPartonGenHist", "Distribution of min delta R for each gen jet", 100, 0.0, 10.0);
    TH1* minDeltaRPartonRecoGenFoundHist = new TH1F("minDeltaRPartonRecoGenFoundHist", "Distribution of min delta R for each reco jet with gen jet match", 100, 0.0, 10.0);
    TFile h("./data/histos.root", "new");

    // file to write data as txt
    std::ofstream write_out("./data/jetInfo.txt");
    assert(write_out.is_open());

    // genPartTree variables
    std::vector<Int_t>* pdgId = 0;
    std::vector<Float_t>* Px = 0;
    std::vector<Float_t>* Py = 0;
    std::vector<Float_t>* Pz = 0;
    std::vector<Float_t>* E = 0;
    std::vector<Float_t>* Pt = 0;

    // jetTree variables
    ULong64_t event = 0;
    ULong64_t recoJetEvent = 0;
    Float_t recoJetEta = 0;
    Float_t recoJetPhi = 0;
    Float_t recoJetPt = 0;
    Int_t recoJetGenMatch = 0;
    Float_t genJetPt = 0;
    Float_t genJetEta = 0;
    Float_t genJetPhi = 0;

    // hold the number of jets in each event for comparison
    //std::vector<int> numGenJetsVec;
    //std::vector<int> numRecoJetsVec;

    tree1->SetBranchAddress("genPartPdgId", &pdgId);
    tree1->SetBranchAddress("genPartPx", &Px);
    tree1->SetBranchAddress("genPartPy", &Py);
    tree1->SetBranchAddress("genPartPz", &Pz);
    tree1->SetBranchAddress("genPartE", &E);
    tree1->SetBranchAddress("event", &event);
    tree1->SetBranchAddress("genPartPt", &Pt);

    tree2->SetBranchAddress("event", &recoJetEvent);
    tree2->SetBranchAddress("jetEta", &recoJetEta);
    tree2->SetBranchAddress("jetPhi", &recoJetPhi);
    tree2->SetBranchAddress("jetPt", &recoJetPt);
    tree2->SetBranchAddress("jetGenMatch", &recoJetGenMatch);
    tree2->SetBranchAddress("genJetPt", &genJetPt);
    tree2->SetBranchAddress("genJetEta", &genJetEta);
    tree2->SetBranchAddress("genJetPhi", &genJetPhi);

    int numEvents = tree1->GetEntries();
    int numRecoJetsTot = tree2->GetEntries();
    int jetIndex = 0;

    tree2->GetEntry(jetIndex);

    // loop through events
    for (size_t i=0; i < numEvents; i++) {
        std::vector<float> recoJetEtaVec;
        std::vector<float> recoJetPhiVec;
        std::vector<float> partonJetEtaVec;
        std::vector<float> partonJetPhiVec;
        std::vector<int> recoJetGenMatchVec;
        std::vector<int> genJetEtaVec;
        std::vector<int> genJetPhiVec;

        tree1->GetEntry(i);
         
        // get eta and phi from reco jets
        int numRecoJets = 0;
        while (true) {
            jetIndex++;
            tree2->GetEntry(jetIndex);
            if (recoJetEvent != event || jetIndex > numRecoJetsTot) break; //TODO: This seems a little hacky, figure out better logic
            numRecoJets++;
            recoJetEtaVec.push_back(recoJetEta);
            recoJetPhiVec.push_back(recoJetPhi);
            recoJetPtHist->Fill(recoJetPt);
            recoJetGenMatchVec.push_back(recoJetGenMatch);
            if (recoJetGenMatch == 1) {
                genJetEtaVec.push_back(genJetEta);
                genJetPhiVec.push_back(genJetPhi);
            }
            else {
                genJetEtaVec.push_back(0.0/0.0);
                genJetPhiVec.push_back(0.0/0.0);

            }
            write_out << 2 << " " << recoJetEvent << " " << recoJetPt << " " << recoJetEta << " " << recoJetPhi << "\n";

        }

        numRecoJetsHist->Fill(numRecoJets);

        int numPartons = pdgId->size();
        std::vector<PseudoJet> particles;
        for (size_t j=0; j < numPartons; j++) {
            particles.push_back( PseudoJet( (*Px)[j], (*Py)[j], (*Pz)[j], (*E)[j]) );
        }

        // cluster partons into jet with anti-kt algorithm
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> partonJets = sorted_by_pt(cs.inclusive_jets());

        //numGenJetsVec.push_back(jets.size());
        numPartonJetsHist->Fill(partonJets.size());

        // get eta and phi from parton gets
        int numPartonJetsCut = 0;
        for (size_t j=0; j < partonJets.size(); j++) {
            partonJetEtaVec.push_back(partonJets[j].rap());
            partonJetPhiVec.push_back(partonJets[j].phi());
            partonJetPtHist->Fill(partonJets[j].pt());
            write_out << 0 << " " << event << " " << partonJets[j].pt() << " " << partonJets[j].rap() << " " << partonJets[j].phi() << "\n";
            if (partonJets[j].pt() > 20.0){
                numPartonJetsCut++;
                partonJetPtCutHist->Fill(partonJets[j].pt());
            }
        }
        numPartonJetsCutHist->Fill(numPartonJetsCut);

        // go through each reco Jet and see how many matches there are
        if (recoJetEtaVec.size() > 0) {
            for (int j=0; j < recoJetEtaVec.size(); j++) {
                float recoEta = recoJetEtaVec[j];
                float recoPhi = recoJetPhiVec[j];
                float genEta = genJetEtaVec[j];
                float genPhi = genJetPhiVec[j];
                int matches = 0;
                float minDRPartonsReco = 10.0;
                float minDRPartonsGen = 10.0;
                for (size_t k=0; k < partonJetEtaVec.size(); k++) {
                    float dRPartonsReco = deltaR(recoEta, recoPhi, partonJetEtaVec[k], partonJetPhiVec[k]);
                    float dRPartonsGen = deltaR(genEta, genPhi, partonJetEtaVec[k], partonJetPhiVec[k]);
                    if (dRPartonsReco < minDRPartonsReco) minDRPartonsReco = dRPartonsReco;
                    if (dRPartonsGen < minDRPartonsGen) minDRPartonsGen = dRPartonsGen;
                    if (dRPartonsReco < 0.4) matches++;
                }
                minDeltaRPartonRecoHist->Fill(minDRPartonsReco);
                minDeltaRPartonRecoHistZoom->Fill(minDRPartonsReco);
                if (!isnan(genEta)) minDeltaRPartonGenHist->Fill(minDRPartonsGen);
                if (recoJetGenMatchVec[j] == 1) minDeltaRPartonRecoGenFoundHist->Fill(minDRPartonsReco);

            }
        }
    }

    numPartonJetsHist->Write();
    numPartonJetsCutHist->Write();
    numRecoJetsHist->Write();
    minDeltaRPartonRecoHist->Write();
    minDeltaRPartonRecoHistZoom->Write();
    recoJetPtHist->Write();
    partonJetPtHist->Write();
    partonJetPtCutHist->Write();
    minDeltaRPartonGenHist->Write();
    minDeltaRPartonRecoGenFoundHist->Write();

    delete numPartonJetsHist;
    delete numPartonJetsCutHist;
    delete numRecoJetsHist;
    delete minDeltaRPartonRecoHist;
    delete minDeltaRPartonRecoHistZoom;
    delete recoJetPtHist;
    delete partonJetPtHist;
    delete partonJetPtCutHist;
    delete minDeltaRPartonGenHist;
    delete minDeltaRPartonRecoGenFoundHist;

    write_out.close();

    return 0;
}
