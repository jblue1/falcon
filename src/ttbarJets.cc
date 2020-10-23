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


    //TFile h(histosPath.c_str(), "RECREATE");
    // set up histograms

    // file to write data as txt
    //std::ofstream write_out("./data/txt/" + txtFile);
    //assert(write_out.is_open());
    std::vector<Float_t>* pfJetPt = 0;
    std::vector<Float_t>* pfJetEta = 0;
    std::vector<Float_t>* pfJetPhi = 0;

    std::vector<Float_t>* partonPx = 0;
    std::vector<Float_t>* partonPy = 0;
    std::vector<Float_t>* partonPz = 0;
    std::vector<Float_t>* partonE = 0;

    tree->SetBranchAddress("pfJetPt", &pfJetPt);
    tree->SetBranchAddress("pfJetEta", &pfJetEta);
    tree->SetBranchAddress("pfJetPhi", &pfJetPhi);

    tree->SetBranchAddress("partonPx", &partonPx);
    tree->SetBranchAddress("partonPy", &partonPy);
    tree->SetBranchAddress("partonPz", &partonPz);
    tree->SetBranchAddress("partonE", &partonE);


    //pfCandTree variables
    int numEvents = tree->GetEntries();

    int numPartonJets = 0;
    int numPartonJetswMatch = 0;
    int numRecoJets = 0;
    int numRecoJetswMatch = 0;

    // loop through events
    for (int i=0; i < numEvents; i++) {
        std::cout << "Event: " << i << std::endl;
        std::cout << "=======================================" << std::endl;
        tree->GetEntry(i);

        // create vector with all parton 4-momenta
        int numPartons = partonPx->size();
        std::vector<PseudoJet> particles;
        for (int j=0; j < numPartons; j++) {
            particles.push_back(PseudoJet((*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]));
        }

        // cluster partons into jet with anti-kt algorithm
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);

        std::vector<PseudoJet> partonJets = cs.inclusive_jets(); // get new parton jets
        
        // go through each parton jet and find how many reco and gen jet
        // matches there are (and find min dR)
        if (partonJets.size() > 0) {
            for (size_t j=0; j < partonJets.size(); j++) {
                
                if (partonJets[j].pt() > 20) {
                    numPartonJets++;
                    float partonEta = partonJets[j].rap();
                    float partonPhi = partonJets[j].phi_std();
                    float partonPt = partonJets[j].pt();
                    //std::cout << "PartonJet: " << j << std::endl;
                    //std::cout << "Eta: " << partonEta << std::endl;
                    //std::cout << "Phi: " << partonPhi << std::endl;
                    //std::cout << "Pt: " << partonPt << std::endl;

                    int pfJetMatches = 0;
                    float minDR = 10.0;
                    for (size_t k=0; k < pfJetPt->size(); k++) {
                        if ((*pfJetPt)[k] >= 30.0) {
                            //std::cout << "  RecoJet: " << k << std::endl;
                            //std::cout << "  Eta: " << (*pfJetEta)[k] << std::endl;
                            //std::cout << "  Phi: " << (*pfJetPhi)[k] << std::endl;
                            float dR = deltaR((*pfJetEta)[k], (*pfJetPhi)[k], partonEta, partonPhi);
                            //std::cout << "  DR: " << dR << std::endl;
                            if (dR < 0.35) pfJetMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    if (pfJetMatches == 1) numPartonJetswMatch++;
                    //std::cout << "Num Matches: " << pfJetMatches << std::endl;
                    //std::cout << std::endl << std::endl;
                }
            }
        }

        if (pfJetPt->size() > 0) {
            for (size_t j=0; j < pfJetPt->size(); j++) {
                if ((*pfJetPt)[j] > 30.0) {
                    numRecoJets++;
                    //std::cout << "RecoJet: " << j << std::endl;
                    //std::cout << "Eta: " << (*pfJetEta)[j] << std::endl;
                    //std::cout << "Phi: " << (*pfJetPhi)[j] << std::endl;
                    //std::cout << "Pt: " << (*pfJetPt)[j] << std::endl; 

                    int partonJetMatches = 0;
                    float minDR = 10.0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if (partonJets[k].pt() > 20.0) {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), (*pfJetEta)[j], (*pfJetPhi)[j]);
                            if (j == 3 && i == 3) {
                                std::cout << "  PartonJet: " << j << std::endl;
                                std::cout << "  Eta: " << partonJets[k].rap() << std::endl;
                                std::cout << "  Phi: " << partonJets[k].phi_std() << std::endl;
                                std::cout << "  Pt: " << partonJets[k].pt() << std::endl;
                                std::cout << "  DR: " << dR << std::endl;
                            }

                            if (dR < 0.35) partonJetMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    
                    if (partonJetMatches == 1) numRecoJetswMatch++;
                    if (partonJetMatches == 0) {
                        std::cout << "Reco Jet with no match number: " << j << std::endl;
                        std::cout << "Reco Jet with no match Pt: " << (*pfJetPt)[j] << std::endl;
                        std::cout << "Reco Jet with no match Eta: " << (*pfJetEta)[j] << std::endl;
                        std::cout << "Reco Jet with no match Phi: " << (*pfJetPhi)[j] << std::endl;
                        std::cout << "Reco Jet with no match minDR: " << minDR << std::endl;
                    }
                    //std::cout << "Num matches for this reco jet: " << partonJetMatches << std::endl << std::endl;
                }
            }

        }
    }
    std::cout << "Num Reco Jets: " << numRecoJets << std::endl;
    std::cout << "Num Reco Jets w/ Match: " << numRecoJetswMatch << std::endl;
    std::cout << "Num Parton Jets: " << numPartonJets << std::endl;
    std::cout << "Num Parton Jets w/ Match: " << numPartonJetswMatch << std::endl;
    return 0;
}
