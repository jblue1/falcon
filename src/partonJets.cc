#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"


using namespace fastjet;

int main(int argc, char const *argv[]) {
    TFile f("../data/JetNtuple_RunIISummer16_13TeV_MC.root");
    TDirectoryFile* df = (TDirectoryFile*) f.Get("AK4jets");
    TTree* tree1 = (TTree*) df->Get("genPartTree");
    TTree* tree2 = (TTree*) df->Get("jetTree");
    TH1* numGenJetsHist = new TH1I("numGenJets", "Distribution of parton jets per event", 31, 0, 30);
    TH1* numRecoJetsHist = new TH1I("numRecoJets", "Distribution of reco jets per event", 31, 0, 30);
    TFile h("histos.root", "new");

    // genPartTree variables
    std::vector<Int_t>* pdgId = 0;
    std::vector<Float_t>* Px = 0;
    std::vector<Float_t>* Py = 0;
    std::vector<Float_t>* Pz = 0;
    std::vector<Float_t>* E = 0;

    // jetTree variables
    ULong64_t* event = 0;
    ULong64_t* jetEvent = 0;
    Float_t* jetEta = 0;
    Float_t* jetPhi = 0;

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

    int numEvents = tree1->GetEntries();
    int jetIndex = 0;
    tree2->GetEntry(jetIndex);
    for (int i=0; i < numEvents; i++) {
        tree1->GetEntry(i);
        int numRecoJets = 0;
        while (jetEvent == event) {
            jetIndex++;
            numRecoJets++;
            tree2->GetEntry(jetIndex);
        }
        numRecoJetsVec.push_back(numRecoJets);
        numRecoJetsHist->Fill(numRecoJets);

        int numPartons = pdgId->size();
        std::vector<PseudoJet> particles;
        for (int j=0; j < numPartons; j++) {
            particles.push_back( PseudoJet( (*Px)[j], (*Py)[j], (*Pz)[j], (*E)[j]) );
        }

        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);

        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
        numGenJetsVec.push_back(jets.size());
        numGenJetsHist->Fill(jets.size());
    }
    numRecoJetsHist->Write();
    numGenJetsHist->Write();
    delete numGenJetsHist;
    delete numRecoJetsHist;

    return 0;
}
