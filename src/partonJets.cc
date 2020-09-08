#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include <cmath>
#include <assert.h>


using namespace fastjet;

float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return sqrt( pow(eta2 - eta1, 2) + pow(phi2 - phi1, 2));
}

int main(int argc, char const *argv[]) {
    // open data file and get trees
    TFile f("../data/JetNtuple_RunIISummer16_13TeV_MC.root");
    TDirectoryFile* df = (TDirectoryFile*) f.Get("AK4jets");
    TTree* tree1 = (TTree*) df->Get("genPartTree");
    TTree* tree2 = (TTree*) df->Get("jetTree");
    // set up histagrams
    TH1* numGenJetsHist = new TH1I("numGenJets", "Distribution of parton jets per event", 31, 0, 30);
    TH1* numRecoJetsHist = new TH1I("numRecoJets", "Distribution of reco jets per event", 31, 0, 30);
    TH1* minDeltaRHist = new TH1F("minDeltaRHist", "Distribution of min delta R for each reco jet", 100, 0.0, 15.0);
    TFile h("../data/histos.root", "new");

    // genPartTree variables
    std::vector<Int_t>* pdgId = 0;
    std::vector<Float_t>* Px = 0;
    std::vector<Float_t>* Py = 0;
    std::vector<Float_t>* Pz = 0;
    std::vector<Float_t>* E = 0;

    // jetTree variables
    ULong64_t event = 0;
    ULong64_t jetEvent = 0;
    Float_t jetEta = 0;
    Float_t jetPhi = 0;

    // hold the number of jets in each event for comparison
    std::vector<int> numGenJetsVec;
    std::vector<int> numRecoJetsVec;

    tree1->SetBranchAddress("genPartPdgId", &pdgId);
    tree1->SetBranchAddress("genPartPx", &Px);
    tree1->SetBranchAddress("genPartPy", &Py);
    tree1->SetBranchAddress("genPartPz", &Pz);
    tree1->SetBranchAddress("genPartE", &E);
    tree1->SetBranchAddress("event", &event);

    tree2->SetBranchAddress("event", &jetEvent);
    tree2->SetBranchAddress("jetEta", &jetEta);
    tree2->SetBranchAddress("jetPhi", &jetPhi);

    int numEvents = tree1->GetEntries();
    int jetIndex = 0;

    tree2->GetEntry(jetIndex);

    // loop through events
    for (int i=0; i < numEvents; i++) {
        std::vector<float> recoJetEta;
        std::vector<float> recoJetPhi;
        std::vector<float> genJetEta;
        std::vector<float> genJetPhi;

        // get eta and phi fro reco jets
        int numRecoJets = 0;
        while (jetEvent == event) {
            jetIndex++;
            numRecoJets++;
            tree2->GetEntry(jetIndex);
            recoJetEta.push_back(jetEta);
            recoJetPhi.push_back(jetPhi);
        }

        numRecoJetsVec.push_back(numRecoJets);
        numRecoJetsHist->Fill(numRecoJets);

        tree1->GetEntry(i);
        int numPartons = pdgId->size();
        std::vector<PseudoJet> particles;
        for (int j=0; j < numPartons; j++) {
            particles.push_back( PseudoJet( (*Px)[j], (*Py)[j], (*Pz)[j], (*E)[j]) );
        }

        // cluster partons into jet with anti-kt algorithm
        double R = 0.7;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

        numGenJetsVec.push_back(jets.size());
        numGenJetsHist->Fill(jets.size());

        // get eta and phi from parton gets
        for (int j=0; j < jets.size(); j++) {
            genJetEta.push_back(jets[j].rap());
            genJetPhi.push_back(jets[j].phi());
        }

        // go through each reco Jet and see how many matches there are
        if (recoJetEta.size() > 0) {
            for (int j=0; j < recoJetEta.size(); j++) {
                float recoEta = recoJetEta[j];
                float recoPhi = recoJetPhi[j];
                int matches = 0;
                float mindR = 10.0;
                for (int k=0; k < genJetEta.size(); k++) {
                    float dR = deltaR(recoEta, recoPhi, genJetEta[k], genJetPhi[k]);
                    if (dR < mindR) mindR = dR;
                    if (dR < 0.4) matches++;
                }
                minDeltaRHist->Fill(mindR);

            }
        }
    }
    numRecoJetsHist->Write();
    numGenJetsHist->Write();
    minDeltaRHist->Write();
    delete numGenJetsHist;
    delete numRecoJetsHist;
    delete minDeltaRHist;

    return 0;
}
