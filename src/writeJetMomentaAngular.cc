/*
 * This program makes jets out of partons present in the dataset, runs matching
 * between the parton and reco level jets, and then writes out the 4-momenta of
 * the matched jets in a text file
 */

#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include "helpers.h"
#include "string.h"

using namespace fastjet;

int main(int arc, char const *argv[])
{
    std::string dataFile(argv[1]);
    std::string txtFile(argv[2]);
    std::string dataPath = "./data/raw/" + dataFile;
    std::string txtPath = "./data/processed/" + txtFile;

    TFile f(dataPath.c_str());
    TDirectoryFile *df = (TDirectoryFile *)f.Get("demo");
    TTree *tree = (TTree *)df->Get("eventTree");

    std::ofstream write_out(txtPath);
    assert(write_out.is_open());

    write_out << "The file contains 4 vectors of matched parton and reco jets. The columns are" << std::endl;
    write_out << "Pt_p Eta_p Phi_p E_p Pt_r Eta_r Phi_r E_r" << std::endl;

    std::vector<Float_t> *partonPx = 0;
    std::vector<Float_t> *partonPy = 0;
    std::vector<Float_t> *partonPz = 0;
    std::vector<Float_t> *partonE = 0;
    std::vector<Int_t> *partonStatus = 0;

    std::vector<Float_t> *genJetPt = 0;
    std::vector<Float_t> *genJetEta = 0;
    std::vector<Float_t> *genJetPhi = 0;
    std::vector<Float_t> *genJetE = 0;

    std::vector<Float_t> *pfJetPt = 0;
    std::vector<Float_t> *pfJetEta = 0;
    std::vector<Float_t> *pfJetPhi = 0;
    std::vector<Float_t> *pfJetPx = 0;
    std::vector<Float_t> *pfJetPy = 0;
    std::vector<Float_t> *pfJetPz = 0;
    std::vector<Float_t> *pfJetE = 0;

    tree->SetBranchAddress("pfJetPt", &pfJetPt);
    tree->SetBranchAddress("pfJetEta", &pfJetEta);
    tree->SetBranchAddress("pfJetPhi", &pfJetPhi);
    tree->SetBranchAddress("pfJetPx", &pfJetPx);
    tree->SetBranchAddress("pfJetPy", &pfJetPy);
    tree->SetBranchAddress("pfJetPz", &pfJetPz);
    tree->SetBranchAddress("pfJetE", &pfJetE);

    tree->SetBranchAddress("genJetPt", &genJetPt);
    tree->SetBranchAddress("genJetEta", &genJetEta);
    tree->SetBranchAddress("genJetPhi", &genJetPhi);
    tree->SetBranchAddress("genJetE", &genJetE);

    tree->SetBranchAddress("partonPx", &partonPx);
    tree->SetBranchAddress("partonPy", &partonPy);
    tree->SetBranchAddress("partonPz", &partonPz);
    tree->SetBranchAddress("partonE", &partonE);
    tree->SetBranchAddress("partonStatus", &partonStatus);

    int numEvents = tree->GetEntries();
    int numJets = 0;

    for (int i = 0; i < numEvents; i++)
    {
        tree->GetEntry(i);

        int numPartons = partonPx->size();
        std::vector<PseudoJet> partons;
        for (int j = 0; j < numPartons; j++)
        {
            partons.push_back(PseudoJet((*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]));
            partons[j].set_user_index(j);
        }

        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(partons, jet_def);

        std::vector<PseudoJet> partonJets = cs.inclusive_jets();

        if (partonJets.size() > 0)
        {
            for (size_t j = 0; j < partonJets.size(); j++)
            {
                if (partonJets[j].pt() > 20)
                {
                    float minDRPfJet = 10.0;
                    int pfJetIndex = 0;
                    for (size_t k = 0; k < pfJetPt->size(); k++)
                    {
                        float dR = deltaR((*pfJetEta)[k], (*pfJetPhi)[k], partonJets[j].rap(),
                                          partonJets[j].phi_std());
                        if (dR < minDRPfJet)
                        {
                            minDRPfJet = dR;
                            pfJetIndex = k;
                        }
                    }

                    if (minDRPfJet < 0.35)
                    {
                        float partonJetPSquared = pow(partonJets[j].px(), 2) + pow(partonJets[j].py(), 2) + pow(partonJets[j].pz(), 2);
                        float pfJetPSquared = pow((*pfJetPx)[pfJetIndex], 2) + pow((*pfJetPy)[pfJetIndex], 2) + pow((*pfJetPz)[pfJetIndex], 2);
                        float partonJetMSquared = pow(partonJets[j].e(), 2) - partonJetPSquared;
                        float pfJetMSquared = pow((*pfJetE)[pfJetIndex], 2) - pfJetPSquared;
                        if (partonJetMSquared > 0 && pfJetMSquared > 0)
                        {
                            numJets++;
                            write_out
                                << partonJets[j].pt() << " "
                                << partonJets[j].eta() << " "
                                << partonJets[j].phi_std() << " "
                                << sqrt(partonJetMSquared) << " "
                                << (*pfJetPt)[pfJetIndex] << " "
                                << (*pfJetEta)[pfJetIndex] << " "
                                << (*pfJetPhi)[pfJetIndex] << " "
                                << sqrt(pfJetMSquared) << std::endl;
                        }
                    }
                }
            }
        }
    }
    std::cout << "Wrote " << numJets << " to " << txtPath << std::endl;
    write_out.close();

    return 0;
}
