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
    write_out << "Px_p Py_p Pz_p E_p tbar bbar cba sbar ubar dbar d u s c b t g Px_g Py_g Pz_g E_g" << std::endl;

    std::vector<Float_t>* partonPx = 0;
    std::vector<Float_t>* partonPy = 0;
    std::vector<Float_t>* partonPz = 0;
    std::vector<Float_t>* partonE = 0;
    std::vector<Int_t>* partonPdgId = 0;
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
    tree1->SetBranchAddress("genPartPdgId", &partonPdgId);

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
            partons[j].set_user_index(j);
        }

        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(partons, jet_def);

        std::vector<PseudoJet> partonJets = cs.inclusive_jets();

        if (partonJets.size() > 0) {
            for (size_t j=0; j< partonJets.size(); j++) {
                if (partonJets[j].pt() > 20) {
                    int tbar = 0;
                    int bbar = 0;
                    int cbar = 0;
                    int sbar = 0;
                    int ubar = 0;
                    int dbar = 0;
                    int d = 0;
                    int u = 0;
                    int s = 0;
                    int c = 0;
                    int b = 0;
                    int t = 0;
                    int g = 0;

                    std::vector<PseudoJet> constituents = sorted_by_pt(partonJets[j].constituents());
                    int index = constituents[0].user_index();
                    std::cout << "Highest Pt PdgId: " << (*partonPdgId)[index] << std::endl;
                    int pdgId = (*partonPdgId)[index];
                    
                    if (pdgId == -6) tbar = 1;
                    else if (pdgId == -5) bbar = 1;
                    else if (pdgId == -4) cbar = 1;
                    else if (pdgId == -3) sbar = 1;
                    else if (pdgId == -2) ubar = 1;
                    else if (pdgId == -1) dbar = 1;
                    else if (pdgId == 1) d = 1;
                    else if (pdgId == 2) u = 1;
                    else if (pdgId == 3) s = 1;
                    else if (pdgId == 4) c = 1;
                    else if (pdgId == 5) b = 1;
                    else if (pdgId == 6) t = 1;
                    else if (pdgId == 21) g = 1;
                    else {
                        std::cout << "Unknown Pdgid: " << pdgId << std::endl;
                        exit(1);
                    }

                    float partonEta = partonJets[j].rap();
                    float partonPhi = partonJets[j].phi_std();
                    float minDR = 10.0;
                    int jetIndex = 0;
                    for (size_t k=0; k < genJetPt->size(); k++) {
                        if ((*genJetPt)[k] > 30) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], partonEta, partonPhi);
                            if (dR < minDR) {
                                minDR = dR;
                                jetIndex = k;
                            }
                        }
                    }
                    if (minDR < 0.35) {
                        int sum = tbar + bbar + cbar + sbar + ubar + dbar + d + u + s + c + b + t + g;
                        assert(sum == 1);
                        write_out 
                            << partonJets[j].pt() << " " 
                            << partonJets[j].rap() << " " 
                            << partonJets[j].phi_std() << " " 
                            << partonJets[j].e() << " " 
                            << tbar << " "
                            << bbar << " "
                            << cbar << " "
                            << sbar << " "
                            << ubar << " "
                            << dbar << " "
                            << d << " "
                            << u << " "
                            << s << " "
                            << c << " "
                            << b << " "
                            << t << " "
                            << g << " "
                            << (*genJetPt)[jetIndex] << " "
                            << (*genJetEta)[jetIndex] << " " 
                            << (*genJetPhi)[jetIndex] << " "
                            << (*genJetE)[jetIndex] << std::endl;
                        //std::cout << "Parton Pt: " << partonJets[j].pt() << std::endl;
                        //std::cout << "Gen Pt: " << (*genJetPt)[jetIndex] << std::endl;
                        //std::cout << "Parton Eta: " << partonJets[j].rap() << std::endl;
                        //std::cout << "Gen Eta: " << (*genJetEta)[jetIndex] << std::endl;
                        //std::cout << "Parton Phi: " << partonJets[j].phi() << std::endl;
                        //std::cout << "Gen Phi: " << (*genJetPhi)[jetIndex] << std::endl;
                        }
                    }
                }
            }
        }


    write_out.close();

    return 0;
}
