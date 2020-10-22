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
    TTree* tree2 = (TTree*) df->Get("jetTree");
    TTree* tree1 = (TTree*) df->Get("genPartTree");
    TTree* tree3 = (TTree*) df->Get("genJetTree");
    TTree* tree4 = (TTree*) df->Get("pfCandTree");


    TFile h(histosPath.c_str(), "RECREATE");
    // set up histograms
    TH1* numMatchesPartonRecoCHS_20_30_Hist = new TH1F("numMatchesPartonRecoCHS_20_30_Hist", "Number of reco jet (CHS) matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* numMatchesPartonGen_20_30_Hist = new TH1F("numMatchesPartonGen_20_30_Hist", "Number of gen jet matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* numMatchesPartonGen_20_15_Hist = new TH1F("numMatchesPartonGen_20_15_Hist", "Number of gen jet matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* numMatchesPartonRecoNoCHSAllHist = new TH1F("numMatchesPartonRecoNoCHSAllHist", "Number of reco jet (no CHS, all particles) matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1* numMatchesPartonRecoNoCHSHadHist = new TH1F("numMatchesPartonRecoNoCHSHadHist", "Number of reco jet (no CHS, just hadrons) matches (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);

    TH1* minDRPartonsRecoCHSHist = new TH1F("minDRPartonsRecoCHSHist", "Distribution of min dR for each parton jet (matching parton to reco)", 100, 0, 1.5);
    TH1* minDRPartonsGenHist = new TH1F("minDRPartonsGenHist", "Distribution of min dR for each parton jet (matching parton to gen)", 100, 0, 1.5);
    TH1* minDRPartonsRecoNoCHSAllHist = new TH1F("minDRPartonsRecoNoCHSAllHist", "Distribution of min dR for each parton jet (matching parton to reco (no CHS, all particles)", 100, 0, 1.5);
    TH1* minDRPartonsRecoNoCHSHadHist = new TH1F("minDRPartonsRecoNoCHSHadHist", "Distribution of min dR for each parton jet (matching parton to reco (no CHS, just hadrons)", 100, 0, 1.5);

    TH1* partonPtNoGenMatchHist = new TH1F("partonPtNoGenMatchHist", "Distribution of Pt for parton jets with no genJet match", 100, 0, 200);
    TH1* partonPtNoRecoCHSMatchHist = new TH1F("partonPtNoRecoCHSMatchHist", "Distribution of Pt for parton jets with no reco jet (CHS) match", 200, 0, 400);
    TH1* partonPtNoRecoNoCHSAllMatchHist = new TH1F("partonPtNoRecoCHSAllMatchHist", "Distribution of Pt for parton jets with no reco jet (no CHS, all particles) match", 200, 0, 400);
    TH1* partonPtNoRecoNoCHSHadMatchHist = new TH1F("partonPtNoRecoCHSHadMatchHist", "Distribution of Pt for parton jets with no reco jet (no CHS, just hadrons) match", 200, 0, 400);

    TH1* numLowPtPartonNoGenMatchHist = new TH1I("numLowPtPartonNoGenMatchHist", "Number of partons with Pt<30Gev with no Gen match", 1, 0, 1);
    TH1* numLowPtPartonNoRecoCHSMatchHist = new TH1I("numLowPtPartonNoRecoCHSMatchHist", "Number of partons with Pt<30GeV with no reco match (CHS)", 1, 0, 1);
    TH1* numLowPtPartonNoRecoNoCHSAllMatchHist = new TH1I("numLowPtPartonNoRecoNoCHSAllMatchHist", "Number of partons with Pt<30GeV with no reco match (no CHS, all particles)", 1, 0, 1);
    TH1* numLowPtPartonNoRecoNoCHSHadMatchHist = new TH1I("numLowPtPartonNoRecoNoCHSHadMatchHist", "Number of partons with Pt<30GeV with no reco match (no CHS, just hadrons)", 1, 0, 1);

    TH1* numMatchesRecoCHSPartonHist = new TH1I("numMatchesRecoCHSPartonHist", "Number of parton jet matches (dR < 0.35) for each reco jet (CHS)", 10, 0, 10);
    TH1* numMatchesGenPartonHist = new TH1I("numMatchesGenPartonHist", "Number of parton jet matches (dR < 0.35) for each gen jet", 10, 0, 10);
    TH1* numMatchesRecoNoCHSAllPartonHist = new TH1I("numMatchesRecoNoCHSAllPartonHist", "Number of parton jet matches (dR < 0.35) for each reco jet (no CHS, all particles)", 10, 0, 10);
    TH1* numMatchesRecoNoCHSHadPartonHist = new TH1I("numMatchesRecoNoCHSHadPartonHist", "Number of parton jet matches (dR < 0.35) for each reco jet (no CHS, just hadrons)", 10, 0, 10);

    TH1* numMatchesGenRecoNoCHSAllHist = new TH1I("numMatchesGenRecoNoCHSAllHist", "number of gen jets with Reco NoCHSALL matches", 10, 0, 10);
    TH1* numMatchesRecoNoCHSAllGen_30_20_Hist = new TH1I("numMatchesRecoNoCHSAllGen_30_20_Hist", "Number of reco jets (No CHS, all particles) with no gen matches", 10, 0, 10);
    TH1* numMatchesRecoNoCHSHadGen_30_20_Hist = new TH1I("numMatchesRecoNoCHSHadGen_30_20_Hist", "Number of reco jets (No CHS, just hadrons) with no gen matches", 10, 0, 10);

    TH1* partonJetPhiHist = new TH1F("partonJetPhiHist", "Distribution of Phi for parton jets", 100, -3.5, 3.5);
    TH1* partonJetEtaHist = new TH1F("partonJetEtaHist", "Distribution of Eta for parton jets", 100, -5, 5);
    TH1* partonJetPhiNoRecoMatchHist = new TH1F("partonJetPhiNoRecoMatchHist", "Distribution of phi for parton jets with no reco match", 100, -3.5, 3.5);
    TH1* partonJetEtaNoRecoMatchHist = new TH1F("partonJetEtaNoRecoMatchHist", "Distribution of eta for parton jets with no reco match", 100, -5, 5);

    TH1* highestPtPartonJetNoRecoMatchHist = new TH1I("highestPtPartonJetNoRecoMatchHist)", "PGDIDs of highest pt parton in each parton jet with no reco match", 33, -6, 27);
    TH1* partonNoRecoMatchPdgIdHist = new TH1I("partonNoRecoMatchPDgIdHist", "PGDIDs of partons in parton jets with no reco match", 33, -6, 27);
    TH1* recoNoCHSAllNoPartonMatchPdgIdHist = new TH1I("recoNoCHSAllNoPartonmatchPdgIdHist", "PGDIDs of pf candidates in reco jets with no parton match", 424, -212, 212);

    TH1* recoNoCHSAllPtHist = new TH1F("recoNoCHSAllPtHist", "Pt Distribution of reco jets (no CHS, all particles)", 1000, 0, 1000);
    TH1* recoNoCHSAllNoPartonMatchPtHist = new TH1F("recoNoCHSAllNoPartonMatchPtHist", "Pt Distribution of reco jets (no CHS, all particles) with no parton jet match", 1000, 0, 1000);

    TH1* numMatchesRecoCHSGenHist = new TH1I("numMatchesRecoCHSGenHist", "Number of reco jets (CHS) with a gen jet match", 10, 0, 10);

    TH1* partonPtHist = new TH1F("partonPtHist", "Parton Pt", 200, 0, 200);
    TH1* pfCandPtHist = new TH1F("pfCandPtHist", "Pf Cand Pt", 200, 0, 200);

    TH1* numPartonVerticesHist = new TH1I("numPartonVerticesHist", "Number of parton vertices in each event", 50, 0, 50);
    TH1* numPfCandVerticesHist = new TH1I("numPfCandVerticesHist", "Number of PF cand vertices in each event", 1500, 0, 1500);
    TH1* numRecoNoCHSAllVerticesHist = new TH1I("numRecoNoCHSAllVerticesHist", "Number of vertices of reco jets (no CHS, all particles) in each event", 50, 0, 50);

    TH1* numRecoVerticesHist = new TH1I("numRecoVerticesHist", "Number of good primary vertices from the reco (CHS) clustering", 50, 0, 50);


    // file to write data as txt
    std::ofstream write_out("./data/txt/" + txtFile);
    assert(write_out.is_open());

    //pfCandTree variables
    std::vector<Float_t>* pfCandPx = 0;
    std::vector<Float_t>* pfCandPy = 0;
    std::vector<Float_t>* pfCandPz = 0;
    std::vector<Float_t>* pfCandE = 0;
    std::vector<Int_t>* pfCandPdgId = 0;
    std::vector<Float_t>* pfCandVx = 0;
    std::vector<Float_t>* pfCandVy = 0;
    std::vector<Float_t>* pfCandVz = 0;
    std::vector<Float_t>* pfCandPt = 0;

    // genPartTree variables
    std::vector<Int_t>* partonPdgId = 0;
    std::vector<Float_t>* partonPx = 0;
    std::vector<Float_t>* partonPy = 0;
    std::vector<Float_t>* partonPz = 0;
    std::vector<Float_t>* partonE = 0;
    std::vector<Float_t>* partonPt = 0;
    std::vector<Float_t>* partonVx = 0;
    std::vector<Float_t>* partonVy = 0;
    std::vector<Float_t>* partonVz = 0;
    ULong64_t partonEvent = 0;

    // jetTree variables
    ULong64_t recoJetEvent = 0;
    Float_t recoJetEta = 0;
    Float_t recoJetPhi = 0;
    Float_t recoJetPt = 0;
    Int_t recoJetGenMatch = 0;
    Int_t numRecoVertices = 0;

    // genJetTree variables
    std::vector<Float_t>* genJetPt = 0;
    std::vector<Float_t>* genJetEta = 0;
    std::vector<Float_t>* genJetPhi = 0;
    ULong64_t genJetEvent = 0;

    tree1->SetBranchAddress("genPartPdgId", &partonPdgId);
    tree1->SetBranchAddress("genPartPx", &partonPx);
    tree1->SetBranchAddress("genPartPy", &partonPy);
    tree1->SetBranchAddress("genPartPz", &partonPz);
    tree1->SetBranchAddress("genPartE", &partonE);
    tree1->SetBranchAddress("event", &partonEvent);
    tree1->SetBranchAddress("genPartPt", &partonPt);
    tree1->SetBranchAddress("genPartVx", &partonVx);
    tree1->SetBranchAddress("genPartVy", &partonVy);
    tree1->SetBranchAddress("genPartVz", &partonVz);

    tree2->SetBranchAddress("event", &recoJetEvent);
    tree2->SetBranchAddress("jetEta", &recoJetEta);
    tree2->SetBranchAddress("jetPhi", &recoJetPhi);
    tree2->SetBranchAddress("jetPt", &recoJetPt);
    tree2->SetBranchAddress("jetGenMatch", &recoJetGenMatch);
    tree2->SetBranchAddress("PV_npvsGood", &numRecoVertices);

    tree3->SetBranchAddress("genJetPt", &genJetPt);
    tree3->SetBranchAddress("genJetEta", &genJetEta);
    tree3->SetBranchAddress("genJetPhi", &genJetPhi);
    tree3->SetBranchAddress("genJetEvent", &genJetEvent);

    tree4->SetBranchAddress("pfCandPx", &pfCandPx);
    tree4->SetBranchAddress("pfCandPy", &pfCandPy);
    tree4->SetBranchAddress("pfCandPz", &pfCandPz);
    tree4->SetBranchAddress("pfCandE", &pfCandE);
    tree4->SetBranchAddress("pfCandPdgId", &pfCandPdgId);
    tree4->SetBranchAddress("pfCandVx", &pfCandVx);
    tree4->SetBranchAddress("pfCandVy", &pfCandVy);
    tree4->SetBranchAddress("pfCandVz", &pfCandVz);
    tree4->SetBranchAddress("pfCandPt", &pfCandPt);

    int numEvents = tree1->GetEntries();

    // used to break out of while loop for getting reco jet data
    int numRecoJetsTot = tree2->GetEntries();
    // used for getting reco jet data
    int jetIndex = 0;

    // need to get first reco jet event number
    tree2->GetEntry(jetIndex);

    // loop through events
    for (int i=0; i < numEvents; i++) {
        //if (i % 500 == 0) std::cout << "Event: " << i << std::endl;
        std::cout << "Event: " << i << std::endl;
        //std::cout << "=======================================" << std::endl;
        std::vector<float> recoJetEtaVec;
        std::vector<float> recoJetPhiVec;
        std::vector<float> recoJetPtVec;
        std::vector<int> recoJetGenMatchVec;

        std::set<std::vector<float>> partonVertices;
        std::set<std::vector<float>> pfCandVertices;
        std::set<std::vector<float>> recoNoCHSAllVertices;

        tree1->GetEntry(i);
        tree3->GetEntry(i);
        tree4->GetEntry(i);
        assert(genJetEvent == partonEvent);

        // save info on genJets
        for (size_t j=0; j < genJetPt->size(); j++) {
            if ((*genJetPt)[j] > 30) {
            write_out
                << 1 << " "
                << genJetEvent << " "
                << 0 << " "
                << (*genJetPt)[j] << " "
                << (*genJetEta)[j] << " "
                << (*genJetPhi)[j] << " "
                << 0 << " "
                << 0 << " "
                << 0 << " "
                << 0 << " "
                << 0 << " "
                << 0 << " "
                << 0 << " "
                << 0 << " "
                << 0 << "\n";
            }
        }

        // get eta and phi from reco jets
        int numRecoJets = 0;
        int count = 0;
        while (true) {
            jetIndex++;
            tree2->GetEntry(jetIndex);
            if (recoJetEvent != partonEvent || jetIndex > numRecoJetsTot) break; //TODO: This seems a little hacky, figure out better logic
            if (count == 0) numRecoVerticesHist->Fill(numRecoVertices);
            count++;

            if (recoJetPt >= 30.0) {
                numRecoJets++;
                recoJetEtaVec.push_back(recoJetEta);
                recoJetPhiVec.push_back(recoJetPhi);
                recoJetPtVec.push_back(recoJetPt);
                recoJetGenMatchVec.push_back(recoJetGenMatch);
                numMatchesRecoCHSGenHist->Fill(recoJetGenMatch);
                write_out
                    << 2 << " "
                    << recoJetEvent << " "
                    << recoJetGenMatch << " "
                    << recoJetPt << " "
                    << recoJetEta << " "
                    << recoJetPhi << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << "\n";
            }
        }

        // create vector with all parton 4-momenta
        int numPartons = partonPdgId->size();
        std::vector<PseudoJet> particles;
        for (int j=0; j < numPartons; j++) {
            partonPtHist->Fill((*partonPt)[j]);
            particles.push_back( PseudoJet( (*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]) );
            particles[j].set_user_index(j); // set index to be able to identify constituent particles later
            std::vector<float> vertex;
            vertex.push_back((*partonVx)[j]);
            vertex.push_back((*partonVy)[j]);
            vertex.push_back((*partonVz)[j]);
            partonVertices.insert(vertex);
        }

        // cluster partons into jet with anti-kt algorithm
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R); // define jet algorithm as anti-kt with R=0.4
        ClusterSequence cs(particles, jet_def); // run the clustering

        std::vector<PseudoJet> partonJets = cs.inclusive_jets(); // get new jets

        // create vector with pf cand 4 momenta
        int numPfCands = pfCandPx->size();
        std::vector<PseudoJet> pfCandsAll;
        std::vector<PseudoJet> pfCandsHad;
        for (int j=0; j < numPfCands; j++) {
            std::vector<float> vertex;
            vertex.push_back((*pfCandVx)[j]);
            vertex.push_back((*pfCandVy)[j]);
            vertex.push_back((*pfCandVz)[j]);
            pfCandVertices.insert(vertex);

            pfCandPtHist->Fill((*pfCandPt)[j]);
            pfCandsAll.push_back( PseudoJet((*pfCandPx)[j], (*pfCandPy)[j], (*pfCandPz)[j], (*pfCandE)[j]));
            pfCandsAll[j].set_user_index(j);
            int pdgId = (*pfCandPdgId)[j];
            if (pdgId == 211 || pdgId == -211 || pdgId == 130 || pdgId == 1 || pdgId == -1) {
                pfCandsHad.push_back( PseudoJet((*pfCandPx)[j], (*pfCandPy)[j], (*pfCandPz)[j], (*pfCandE)[j]));
                // TODO: figure out why this line causes a seg fault
                // pfCandsHad[j].set_user_index(j);
            }
        }
        numPfCandVerticesHist->Fill(pfCandVertices.size());
        assert(pfCandsHad.size() <= pfCandsAll.size());

        ClusterSequence pfSeqAll(pfCandsAll, jet_def);
        std::vector<PseudoJet> pfJetsAll = pfSeqAll.inclusive_jets();

        ClusterSequence pfSeqHad(pfCandsHad, jet_def);
        std::vector<PseudoJet> pfJetsHad = pfSeqHad.inclusive_jets();
        //std::cout << "PF Candidates:" << std::endl;

        for (size_t j=0; j < pfJetsAll.size(); j++) {
            std::vector<PseudoJet> constituents = sorted_by_pt(pfJetsAll[j].constituents());
            int index0 = 0;
            int index1 = 0;
            int index2 = 0;
            int index3 = 0;
            int index4 = 0;
            if (constituents.size() > 0) index0 = constituents[0].user_index();
            if (constituents.size() > 1) index1 = constituents[1].user_index();
            if (constituents.size() > 2) index2 = constituents[2].user_index();
            if (constituents.size() > 3) index3 = constituents[3].user_index();
            if (constituents.size() > 4) index4 = constituents[4].user_index();
            if (pfJetsAll[j].pt() >= 30.0) {
                std::vector<float> vertex;
                vertex.push_back((*pfCandVx)[index0]);
                vertex.push_back((*pfCandVy)[index0]);
                vertex.push_back((*pfCandVz)[index0]);
                recoNoCHSAllVertices.insert(vertex);


                                //std::cout << "  Jet Pt: " << pfJetsAll[j].pt() << std::endl;
                //std::cout << std::setprecision(12) << "  PF cand Vx: " << (*pfCandVx)[index0] << std::endl;
                //std::cout << std::setprecision(12) << "  PF cand Vy: " << (*pfCandVy)[index0] << std::endl;
                //std::cout << std::setprecision(12) << "  PF cand Vz: " << (*pfCandVz)[index0] << std::endl;
                //std::cout << std::endl;

                write_out
                    << std::setprecision(10)
                    << 3 << " "
                    << partonEvent << " "
                    << -1 << " "
                    << pfJetsAll[j].pt() << " "
                    << pfJetsAll[j].rap() << " "
                    << pfJetsAll[j].phi_std() << " "
                    << 0 << " "
                    << 0 << " "
                    << 0 << " "
                    << constituents.size() << " "
                    << (*pfCandPdgId)[index0] << " "
                    << (*pfCandPdgId)[index1] << " "
                    << (*pfCandPdgId)[index2] << " "
                    << (*pfCandPdgId)[index3] << " "
                    << (*pfCandPdgId)[index4] << "\n";
            }
        }
        numRecoNoCHSAllVerticesHist->Fill(recoNoCHSAllVertices.size());


        int numPartonJets = 0;
        // constituent parton)

        //std::cout << "Partons" << std::endl;
        for (size_t j=0; j < partonJets.size(); j++) {
            std::vector<PseudoJet> constituents = sorted_by_pt(partonJets[j].constituents());
            int index0 = 0;
            int index1 = 0;
            int index2 = 0;
            int index3 = 0;
            int index4 = 0;
            if (constituents.size() > 0) index0 = constituents[0].user_index();
            if (constituents.size() > 1) index1 = constituents[1].user_index();
            if (constituents.size() > 2) index2 = constituents[2].user_index();
            if (constituents.size() > 3) index3 = constituents[3].user_index();
            if (constituents.size() > 4) index4 = constituents[4].user_index();
            if (partonJets[j].pt() >= 20.0) {
                numPartonJets++;
                
                //std::cout << "  Jet Pt: " << partonJets[j].pt() << std::endl;
                //std::cout << std::setprecision(12) <<  "  Parton Vx: " << (*partonVx)[index0] << std::endl;
                //std::cout << std::setprecision(12) <<  "  Parton Vy: " << (*partonVy)[index0] << std::endl;
                //std::cout << std::setprecision(12) <<  "  Parton Vz: " << (*partonVz)[index0] << std::endl;
               // std::cout << std::endl;
            write_out
                << std::setprecision(10) 
                << 0 << " " 
                << partonEvent << " " 
                << -1 << " " 
                << partonJets[j].pt() << " " 
                << partonJets[j].rap() << " " 
                << partonJets[j].phi_std() <<  " " 
                << (*partonVx)[index0] << " "
                << (*partonVy)[index0] << " " 
                << (*partonVz)[index0] << " "
                << constituents.size() << " " 
                << (*partonPdgId)[index0]<< " " 
                << (*partonPdgId)[index1]<< " " 
                << (*partonPdgId)[index2]<< " " 
                << (*partonPdgId)[index3]<< " " 
                << (*partonPdgId)[index4]<< "\n";
            }
        }
        numPartonVerticesHist->Fill(partonVertices.size());


        // go through each parton jet and find how many reco and gen Jet
        // matches there are (and find min dR)
        //int numPartonJets = 0;
        if (partonJets.size() > 0) {
            for (size_t j=0; j < partonJets.size(); j++) {
                float partonEta = partonJets[j].rap();
                float partonPhi = partonJets[j].phi_std();

                // before instituting a Pt cutoff on the parton jets, match
                // parton to gen jets with no Pt cut on either
                int matches = 0;
                float minDR = 10.0;
                if (partonJets[j].pt() > 20.0) {
                    for (size_t kk=0; kk < genJetPt->size(); kk++) {
                        if ((*genJetPt)[kk] > 15.0) {
                            float dR = deltaR((*genJetEta)[kk], (*genJetPhi)[kk], partonEta, partonPhi);
                            
                            if (dR < 0.35) {
                                matches++;
                            }
                            if (dR < minDR) minDR = dR;
                        }
                    }
                numMatchesPartonGen_20_15_Hist->Fill(matches);
                }

                // institue Pt cutoff on parton jets
                if (partonJets[j].pt() > 20) {
                    partonJetPhiHist->Fill(partonPhi);
                    partonJetEtaHist->Fill(partonEta);
                    
                    //find parton to recoCHS matches
                    int recoCHSMatches = 0;
                    float minDR = 10.0;
                    for (size_t k=0; k < recoJetEtaVec.size(); k++) {
                        if (recoJetPtVec[k] >= 30.0) {
                            float dR = deltaR(recoJetEtaVec[k], recoJetPhiVec[k], partonEta, partonPhi);
                            if (dR < 0.35) recoCHSMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesPartonRecoCHS_20_30_Hist->Fill(recoCHSMatches);
                    minDRPartonsRecoCHSHist->Fill(minDR);
                    if (recoCHSMatches == 0) {
                        partonPtNoRecoCHSMatchHist->Fill(partonJets[j].pt());
                        if (partonJets[j].pt() < 30) numLowPtPartonNoRecoCHSMatchHist->Fill(0);
                        partonJetPhiNoRecoMatchHist->Fill(partonPhi);
                        partonJetEtaNoRecoMatchHist->Fill(partonEta);
                        // check the pdgIds of partons in jet with no matches
                        std::vector<PseudoJet> constituents = sorted_by_pt(partonJets[j].constituents()); // get consituent partons in the jet
                        for (size_t k=0; k<constituents.size(); k++) {
                            int index = constituents[k].user_index(); // get index of parton in jet with highest pt
                            int pdgId = (*partonPdgId)[index];
                            if (pdgId > 21 || pdgId < -6) {
                            }
                            if (k == 0) highestPtPartonJetNoRecoMatchHist->Fill(pdgId);
                            partonNoRecoMatchPdgIdHist->Fill(pdgId);
                        }
                    }

                    // find parton to gen matches normal Pt cutoffs
                    int genMatches = 0;
                    minDR = 10.0;
                    for (size_t k=0; k < genJetPt->size(); k++) {
                        if ((*genJetPt)[k] >= 30.0) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], partonEta, partonPhi);
                            if (dR < 0.35) genMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesPartonGen_20_30_Hist->Fill(genMatches);
                    minDRPartonsGenHist->Fill(minDR);
                    if (genMatches == 0) {
                        partonPtNoGenMatchHist->Fill(partonJets[j].pt());
                        if (partonJets[j].pt() < 30) numLowPtPartonNoGenMatchHist->Fill(0);
                    }

                    int recoNoCHSAllMatches = 0;
                    minDR = 10.0;
                    for (size_t k=0; k < pfJetsAll.size(); k++) {
                        if (pfJetsAll[k].pt() >= 30.0) {
                            float dR = deltaR(pfJetsAll[k].rap(), pfJetsAll[k].phi_std(), partonEta, partonPhi);
                            if (dR < 0.35) recoNoCHSAllMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    if ( (float) recoNoCHSAllMatches < 0.5 ) {
                        partonPtNoRecoNoCHSAllMatchHist->Fill(partonJets[j].pt());
                        if (partonJets[j].pt() < 30) numLowPtPartonNoRecoNoCHSAllMatchHist->Fill(0);
                    }
                    numMatchesPartonRecoNoCHSAllHist->Fill(recoNoCHSAllMatches++);
                    minDRPartonsRecoNoCHSAllHist->Fill(minDR);
             

                    int recoNoCHSHadMatches = 0;
                    minDR = 10.0;
                    for (size_t k=0; k < pfJetsHad.size(); k++) {
                        if (pfJetsHad[k].pt() >= 30.0) {
                            float dR = deltaR(pfJetsHad[k].rap(), pfJetsHad[k].phi_std(), partonEta, partonPhi);
                            if (dR < 0.35) recoNoCHSHadMatches++;
                            if (dR < minDR) minDR = dR;
                        }
                    }
                    numMatchesPartonRecoNoCHSHadHist->Fill(recoNoCHSHadMatches);
                    minDRPartonsRecoNoCHSHadHist->Fill(minDR);
                    if (recoNoCHSHadMatches == 0) {
                        partonPtNoRecoNoCHSHadMatchHist->Fill(partonJets[j].pt());
                        if (partonJets[j].pt() < 30) numLowPtPartonNoRecoNoCHSHadMatchHist->Fill(0);
                    }
                }
            }
        }

        // go through each gen Jet and see how many parton matches there are
        if (genJetPt->size() > 0) {
            for (size_t j=0; j < genJetPt->size(); j++) {
                if ((*genJetPt)[j] >= 30.0) {
                    int genPartonMatches = 0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if (partonJets[k].pt() > 20) {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), (*genJetEta)[j], (*genJetPhi)[j]);
                            if (dR < 0.35) genPartonMatches++;
                        }
                    }
                    numMatchesGenPartonHist->Fill(genPartonMatches);
                    int genRecoMatches = 0;
                    for (size_t k=0; k < pfJetsAll.size(); k++) {
                       if (pfJetsAll[k].pt() >= 30) {
                           float dR = deltaR((*genJetEta)[j], (*genJetPhi)[j], pfJetsAll[k].rap(), pfJetsAll[k].phi_std());
                           if (dR < 0.35) genRecoMatches++;
                       }
                    }
                    numMatchesGenRecoNoCHSAllHist->Fill(genRecoMatches);
                }
            }
        }

        // go through each reco jet (CHS) and see how many parton matches there
        // are
        if (recoJetEtaVec.size() > 0) {
            for (size_t j=0; j < recoJetEtaVec.size(); j++) {
                if (recoJetPtVec[j] >= 30.0) {
                    int recoCHSPartonMatches = 0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if (partonJets[k].pt() > 20) {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), recoJetEtaVec[j], recoJetPhiVec[j]);
                            if (dR < 0.35) recoCHSPartonMatches++;
                        }
                    }
                    numMatchesRecoCHSPartonHist->Fill(recoCHSPartonMatches);
                }
            }
        }

        // go through each reco jet (NoCHSAll) and see how many parton and gen matches
        // there are
        if (pfJetsAll.size() > 0) {
            for (size_t j=0; j < pfJetsAll.size(); j++) {
                if (pfJetsAll[j].pt() >= 30.0) {
                    recoNoCHSAllPtHist->Fill(pfJetsAll[j].pt());
                    int recoNoCHSAllPartonMatches = 0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if (partonJets[k].pt() > 20) {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), pfJetsAll[j].rap(), pfJetsAll[j].phi_std());
                            if (dR < 0.35) recoNoCHSAllPartonMatches++;
                        }
                    }
                    numMatchesRecoNoCHSAllPartonHist->Fill(recoNoCHSAllPartonMatches);
                    if (recoNoCHSAllPartonMatches == 0) {
                        std::vector<PseudoJet> constituents = sorted_by_pt(pfJetsAll[j].constituents()); // get pf cands in the jet
                        recoNoCHSAllNoPartonMatchPtHist->Fill(pfJetsAll[j].pt());
                        for (size_t k=0; k < constituents.size(); k++) {
                            int index = constituents[k].user_index();
                            int pdgId = (*pfCandPdgId)[index];
                            recoNoCHSAllNoPartonMatchPdgIdHist->Fill(pdgId);
                        }
                    }
                    int recoNoCHSAllGenMatches = 0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if ((*genJetPt)[k] > 20.0) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], pfJetsAll[j].rap(), pfJetsAll[j].phi_std());
                            if (dR < 0.35) recoNoCHSAllGenMatches++;
                        }
                    }
                    numMatchesRecoNoCHSAllGen_30_20_Hist->Fill(recoNoCHSAllGenMatches);
                }
            }
        }

        // go through each reco jet (NoCHSHad) and see how many parton and gen matches
        // there are
        if (pfJetsHad.size() > 0) {
            for (size_t j=0; j < pfJetsHad.size(); j++) {
                if (pfJetsHad[j].pt() >= 30.0) {
                    int recoNoCHSHadPartonMatches = 0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if (partonJets[k].pt() > 20) {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), pfJetsHad[j].rap(), pfJetsHad[j].phi_std());
                            if (dR < 0.35) recoNoCHSHadPartonMatches++;
                        }
                    }
                    numMatchesRecoNoCHSHadPartonHist->Fill(recoNoCHSHadPartonMatches);
                    int recoNoCHSHadGenMatches = 0;
                    for (size_t k=0; k < partonJets.size(); k++) {
                        if ((*genJetPt)[k] > 20.0) {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], pfJetsHad[j].rap(), pfJetsHad[j].phi_std());
                            if (dR < 0.35) recoNoCHSHadGenMatches++;
                        }
                    }
                    numMatchesRecoNoCHSHadGen_30_20_Hist->Fill(recoNoCHSHadGenMatches);
                }
            }
        }
    }

    numMatchesPartonRecoCHS_20_30_Hist->Write();
    numMatchesPartonGen_20_30_Hist->Write();
    numMatchesPartonGen_20_15_Hist->Write();
    numMatchesPartonRecoNoCHSAllHist->Write();
    numMatchesPartonRecoNoCHSHadHist->Write();

    minDRPartonsRecoCHSHist->Write();
    minDRPartonsGenHist->Write();
    minDRPartonsRecoNoCHSAllHist->Write();
    minDRPartonsRecoNoCHSHadHist->Write();

    partonPtNoGenMatchHist->Write();
    partonPtNoRecoCHSMatchHist->Write();
    partonPtNoRecoNoCHSAllMatchHist->Write();
    partonPtNoRecoNoCHSHadMatchHist->Write();

    numLowPtPartonNoGenMatchHist->Write();
    numLowPtPartonNoRecoCHSMatchHist->Write();
    numLowPtPartonNoRecoNoCHSAllMatchHist->Write();
    numLowPtPartonNoRecoNoCHSHadMatchHist->Write();

    numMatchesRecoCHSPartonHist->Write();
    numMatchesGenPartonHist->Write();
    numMatchesRecoNoCHSAllPartonHist->Write();
    numMatchesRecoNoCHSHadPartonHist->Write();

    partonJetPhiHist->Write();
    partonJetEtaHist->Write();
    partonJetPhiNoRecoMatchHist->Write();
    partonJetEtaNoRecoMatchHist->Write();

    highestPtPartonJetNoRecoMatchHist->Write();
    partonNoRecoMatchPdgIdHist->Write();
    recoNoCHSAllNoPartonMatchPdgIdHist->Write();

    numMatchesGenRecoNoCHSAllHist->Write();
    numMatchesRecoNoCHSAllGen_30_20_Hist->Write();
    numMatchesRecoNoCHSHadGen_30_20_Hist->Write();

    recoNoCHSAllPtHist->Write();
    recoNoCHSAllNoPartonMatchPtHist->Write();

    numMatchesRecoCHSGenHist->Write();

    partonPtHist->Write();
    pfCandPtHist->Write();

    numPartonVerticesHist->Write();
    numPfCandVerticesHist->Write();
    numRecoNoCHSAllVerticesHist->Write();
    numRecoVerticesHist->Write();

    delete numMatchesPartonRecoCHS_20_30_Hist;
    delete numMatchesPartonGen_20_30_Hist;
    delete numMatchesPartonGen_20_15_Hist;
    delete numMatchesPartonRecoNoCHSAllHist;
    delete numMatchesPartonRecoNoCHSHadHist;

    delete minDRPartonsRecoCHSHist;
    delete minDRPartonsGenHist;
    delete minDRPartonsRecoNoCHSAllHist;
    delete minDRPartonsRecoNoCHSHadHist;

    delete partonPtNoGenMatchHist;
    delete partonPtNoRecoCHSMatchHist;
    delete partonPtNoRecoNoCHSAllMatchHist;
    delete partonPtNoRecoNoCHSHadMatchHist;

    delete numLowPtPartonNoGenMatchHist;
    delete numLowPtPartonNoRecoCHSMatchHist;
    delete numLowPtPartonNoRecoNoCHSAllMatchHist;
    delete numLowPtPartonNoRecoNoCHSHadMatchHist;

    delete numMatchesRecoCHSPartonHist;
    delete numMatchesGenPartonHist;
    delete numMatchesRecoNoCHSAllPartonHist;
    delete numMatchesRecoNoCHSHadPartonHist;

    delete partonJetPhiHist;
    delete partonJetEtaHist;
    delete partonJetPhiNoRecoMatchHist;
    delete partonJetEtaNoRecoMatchHist;

    delete highestPtPartonJetNoRecoMatchHist;
    delete partonNoRecoMatchPdgIdHist;
    delete recoNoCHSAllNoPartonMatchPdgIdHist;

    delete numMatchesGenRecoNoCHSAllHist;
    delete numMatchesRecoNoCHSAllGen_30_20_Hist;
    delete numMatchesRecoNoCHSHadGen_30_20_Hist;

    delete recoNoCHSAllPtHist;
    delete recoNoCHSAllNoPartonMatchPtHist;
    
    delete numMatchesRecoCHSGenHist;

    delete partonPtHist;
    delete pfCandPtHist;

    delete numPartonVerticesHist;
    delete numPfCandVerticesHist;
    delete numRecoNoCHSAllVerticesHist;
    delete numRecoVerticesHist;

    write_out.close();

    return 0;
}
