#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include "helpers.h"
#include "string.h"

using namespace fastjet;


int main(int arc, char const *argv[]) {
    std::string dataFile(argv[1]);
    std::string txtFile(argv[2]);
    std::string dataPath = "./data/root/" + dataFile;
    std::string txtPath = "./data/txt/" + txtFile;

    TFile f(dataPath.c_str());
    TDirectoryFile* df = (TDirectoryFile*) f.Get("AK4jets");
    TTree* tree1 = (TTree*) df->Get("genPartTree");
    TTree* tree3 = (TTree*) df->Get("genJetTree");


    std::ofstream write_out("./data/txt/" + txtFile);
    assert(write_out.is_open());

    write_out << "The file contains 4 vectors of matched parton and gen jets. The columns are" << std::endl;
    write_out << "Px_p Py_p Pz_p E_p Px_g Py_g Pz_g E_g" << std::endl;

    std::vector<Float_t>* partonPx = 0;
    std::vector<Float_t>* partonPy = 0;
    std::vector<Float_t>* partonPz = 0;
    std::vector<Float_t>* partonE = 0;
    ULong64_t partonEvent = 0;

    std::vector<Float_t> *genJetPt = 0;
    std::vector<Float_t> *genJetEta = 0;
    std::vector<Float_t> *genJetPhi = 0;
    std::vector<Float_t> *genJetPx = 0;
    std::vector<Float_t> *genJetPy = 0;
    std::vector<Float_t> *genJetPz = 0;
    std::vector<Float_t> *genJetE = 0;
    ULong64_t genJetEvent = 0;

    tree1->SetBranchAddress("genPartPx", &partonPx);
    tree1->SetBranchAddress("genPartPy", &partonPy);
    tree1->SetBranchAddress("genPartPz", &partonPz);
    tree1->SetBranchAddress("genPartE", &partonE);
    tree1->SetBranchAddress("event", &partonEvent);

    tree3->SetBranchAddress("genJetPt", &genJetPt);
    tree3->SetBranchAddress("genJetEta", &genJetEta);
    tree3->SetBranchAddress("genJetPhi", &genJetPhi);
    tree3->SetBranchAddress("genJetPx", &genJetPx);
    tree3->SetBranchAddress("genJetPy", &genJetPy);
    tree3->SetBranchAddress("genJetPz", &genJetPz);
    tree3->SetBranchAddress("genJetE", &genJetE);
    tree3->SetBranchAddress("genJetEvent", &genJetEvent);

    int numEvents = tree1->GetEntries();

    for (int i=0; i < numEvents; i++) {
        tree1->GetEntry(i);
        tree3->GetEntry(i);

        assert(genJetEvent == partonEvent);

        int numPartons = partonPx->size();
        std::vector<PseudoJet> partons;
        for (int j=0; j < numPartons; j++) {
            partons.push_back(PseudoJet( (*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]));
        }

        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(partons, jet_def);

        std::vector<PseudoJet> partonJets = cs.inclusive_jets();

        if (partonJets.size() > 0) {
            for (int j=0; j< partonJets.size(); j++) {
                if (partonJets[j].pt() > 20) {
                    float partonEta = partonJets[j].rap();
                    float partonPhi = partonJets[j].phi_std();
                    float minDR = 10.0;
                    int jetIndex = 0;
                    for (int k=0; k < genJetPt->size(); k++) {
                        if ((*genJetPt)[k] > 30) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], partonEta, partonPhi);
                            if (dR < minDR) {
                                minDR = dR;
                                jetIndex = k;
                            }
                        }
                    }
                    if (minDR < 0.35) {
                        write_out << partonJets[j].px() << " " << partonJets[j].py() << " " << partonJets[j].pz()
                            << " " <<  partonJets[j].e() << " " << (*genJetPx)[jetIndex] << " "
                            << (*genJetPy)[jetIndex] << " " << (*genJetPz)[jetIndex] << " "
                            << (*genJetE)[jetIndex] << std::endl;
                        std::cout << "Parton Pt: " << partonJets[j].pt() << std::endl;
                        std::cout << "Gen Pt: " << (*genJetPt)[jetIndex] << std::endl;
                        std::cout << "Parton Px: " << partonJets[j].px() << std::endl;
                        std::cout << "Gen Px: " << (*genJetPx)[jetIndex] << std::endl;
                        std::cout << "Parton Py: " << partonJets[j].py() << std::endl;
                        std::cout << "Gen Py: " << (*genJetPy)[jetIndex] << std::endl;
                        }
                    }
                }
            }
        }


    write_out.close();

    return 0;
}
