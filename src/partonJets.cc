#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include <cmath>
#include <fstream>
#include <cassert>
#include <string>
#include <numeric>
#include <iomanip>

using namespace fastjet;

/**
 * Calculate the average value of a vector
 */
float vecAvg(std::vector<float> &vec) {
    float numerator = std::accumulate(vec.begin(), vec.end(), 0.0);
    int denominator = (float) vec.size();
    return (float) numerator/denominator;
}

/**
 * Calculate the difference between two angles in the range (-pi, pi).
 */
float angleDif(float angle1,float angle2) {
    return atan2(sin(angle1-angle2), cos(angle1-angle2));
}

/**
 * Calculate delta R for two points in eta phi space.
 */
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return sqrt( pow(eta2 - eta1, 2) + pow(angleDif(phi2, phi1), 2));
}

/**
 * Describe parameters taken and print erros if parameters are incorrectly
 * given.
 */
void usage(std::ostream &out, const char *msg) {
    out << msg << std::endl;
    out << std::endl;
    out << "    Usage:" << std::endl;
    out << "            makeJets.out data histos txt" << std::endl;
    out << "    data        - root file from JetNtuple analyzer to use" << std::endl;
    out << "    histos      - root file to write histograms to" << std::endl;
    out << "    txt         - text file to write jet eta phi data to" << std::endl;
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
    TTree* tree1 = (TTree*) df->Get("genPartTree");
    TTree* tree2 = (TTree*) df->Get("jetTree");
    TTree* tree3 = (TTree*) df->Get("genJetTree");
    // set up histograms
    TH1* numMatchesPartonRecoHist = new TH1F("numMatchesPartonJetHist", "Number of reco jet matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* numMatchesPartonGenHist = new TH1F("numMatchesPartonJetHist", "Number of gen jet matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TFile h(histosPath.c_str(), "new");

    // file to write data as txt
    std::ofstream write_out("./data/txt/" + txtFile);
    assert(write_out.is_open());
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

    int numEvents = tree1->GetEntries();
    int numRecoJetsTot = tree2->GetEntries();
    int jetIndex = 0;

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

        int numPartons = pdgId->size();
        std::vector<PseudoJet> particles;
        for (size_t j=0; j < numPartons; j++) {
            particles.push_back( PseudoJet( (*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]) );
            particles[j].set_user_index(j);
        }

        // cluster partons into jet with anti-kt algorithm
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R); // define jet algorithm as anti-kt with R=0.4
        ClusterSequence cs(particles, jet_def); // run the clustering
        //std::vector<PseudoJet> partonJets = sorted_by_pt(cs.inclusive_jets()); // sort jets by Pt

        std::vector<PseudoJet> partonJets = cs.inclusive_jets(); // sort jets by Pt

        //numGenJetsVec.push_back(jets.size());

        // get eta and phi from parton gets
        int numPartonJetsCut = 0;
        for (size_t j=0; j < partonJets.size(); j++) {
            std::vector<PseudoJet> constituents = partonJets[j].constituents();
            std::vector<float> constitVxVec;
            std::vector<float> constitVyVec;
            std::vector<float> constitVzVec;
            for (int ii=0; ii < constituents.size(); ii++) {
                int partonIndex = constituents[ii].user_index();
                constitVxVec.push_back((*Vx)[partonIndex]);
                constitVyVec.push_back((*Vy)[partonIndex]);
                constitVzVec.push_back((*Vz)[partonIndex]);
            }
            write_out << std::setprecision(10) << 0 << " " << partonEvent << " " << -1 << " " << partonJets[j].pt() <<
                " " << partonJets[j].rap() << " " << partonJets[j].phi_std() <<  " " << vecAvg(constitVxVec)
                << " " << vecAvg(constitVyVec) << " " << vecAvg(constitVzVec)<< "\n";
        }


        // go through each parton jet and find how many reco and gen Jet
        // matches there are
        if (partonJets.size() > 0) {
            for (int j=0; j < partonJets.size(); j++) {
                if (partonJets[j].pt() > 20) {
                    float partonEta = partonJets[j].rap();
                    float partonPhi = partonJets[j].phi_std();

                    int recoMatches = 0;
                    for (int k=0; k < recoJetEtaVec.size(); k++) {
                        if (recoJetPtVec[k] > 30) {
                            float dR = deltaR(recoJetEtaVec[k], recoJetPhiVec[k], partonEta, partonPhi);
                            if (dR < 0.35) recoMatches++;
                        }
                    }
                    numMatchesPartonRecoHist->Fill(recoMatches);
                    int genMatches = 0;
                    for (int k=0; k < genJetPt->size(); k++) { 
                        if ((*genJetPt)[k] > 30) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], partonEta, partonPhi);
                            if (dR < 0.35) genMatches++;
                        }
                    }
                    numMatchesPartonGenHist->Fill(genMatches);
                }
            }
        }
    }

    numMatchesPartonRecoHist->Write();
    numMatchesPartonGenHist->Write();

    delete numMatchesPartonRecoHist;
    delete numMatchesPartonGenHist;

    write_out.close();

    return 0;
}
