/**
 * This program clusters partons in dataset into jets, then runs jet matching
 * between the parton, generated, and reco level jets. It also writes out the
 * 4-momenta of each particle and jet to a specified location. Finally, it also produces the
 * following histograms concering the the jets and the matching
 * 
 * numMatchesPartonRecoHist - the number of reco jet matches for each parton jet
 * numMatchesPartonGenHist - the number of gen jet matches for each parton jet
 * numMatchesGenPartonHist - the number of parton jet matches for each gen jet
 * numMatchesGenRecoHist - the number of reco jet matches for each gen jet
 * numMatchesRecoPartonHist - the number of parton jet matches for each reco jet
 * numMatchesRecoGenHist - the number of gen jet matches for each reco jet
 * partonPtNoRecoMatchHist - the distribution of transverse momentum for parton
 * jets that did not have a reco jet match
 *
 * recoPtNoPartonMatchHist - the distribution of transverse momentum for reco
 * jets that did not have a parton jet match
 *
 * minDRRecoPartonNoMatchHist - the distribution of the smallest deltaR value
 * found matcheing reco jets to parton jets when no match (deltaR < 0.35) was
 * found
 *
 * emFractionRecoPartonNoMatchHist - distribution of the fraction of EM energy
 * in each reco jet with no parton jet match
 *
 * emFractionRecoPartonMatchHist - distribution of the fraction of EM energy in
 * each reco jet with a parton jet match
 *
 * partonPtHist - distribution of the transverse momentum of each parton
 *
 * recoEtaDistributionPartonMatchHist - distribution of the pseudo-rapidity of
 * reco jets with a parton jet match
 *
 * recoEtaDistributionPtDiscrepencyHist - distriubtion of the pseudo-rapidity
 * of reco jets with parton jet match, but when the ratio of recoJetPt /
 * partonJetPt is less than 0.75
 *
 * emFractionRecoPartonMatchPtDiscrepency - distribution of the fraction of 
 * EM energy in each recoJet with a parton jet match when the ratio of
 * recoJetPt / parotnJetPt is less than 0.75
 */

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
void usage(std::ostream &out, const char *msg)
{
    out << msg << std::endl;
    out << std::endl;
    out << "    Usage:" << std::endl;
    out << "            makeJets.out data histos txt" << std::endl;
    out << "    data        - root file from JetNtuple analyzer to use" << std::endl;
    out << "    histos      - root file to write histograms to" << std::endl;
    exit(1);
}

int main(int argc, char const *argv[])
{
    if (argc != 4)
    {
        usage(std::cerr, "Incorrect number of parameters given.");
    }

    std::string dataFile(argv[1]);
    std::string histosFile(argv[2]);
    std::string txtFile(argv[3]);
    std::string dataPath = "./data/raw/" + dataFile;
    std::string histosPath = "./data/plots/" + histosFile;
    std::string txtPath1 = "./data/processed/" + txtFile + "particles.txt";
    std::string txtPath2 = "./data/processed/" + txtFile + "jets.txt";

    // open data file and get trees
    TFile f(dataPath.c_str());
    TDirectoryFile *df = (TDirectoryFile *)f.Get("demo");
    TTree *tree = (TTree *)df->Get("eventTree");

    std::ofstream write_out_particles(txtPath1);
    std::ofstream write_out_jets(txtPath2);
    assert(write_out_particles.is_open());
    assert(write_out_jets.is_open());

    write_out_particles << "event status pdgId Pt Eta Phi E" << std::endl;
    write_out_jets << "event part/gen/reco(0/1/2) Pt Eta Phi E" << std::endl;

    TFile h(histosPath.c_str(), "RECREATE");
    // set up histograms
    TH1 *numMatchesPartonRecoHist = new TH1I("numMatchesPartonRecoHist", "Number of parton Jets with reco jet match (dR < 0.35) for each parton jet with Pt > 20GeV", 5, 0, 5);
    TH1 *numMatchesPartonGenHist = new TH1I("numMatchesPartonGenHist", "Number of parton jets with gen jet match", 5, 0, 5);
    TH1 *numMatchesGenPartonHist = new TH1I("numMatchesGenPartonHist", "Number of gen jets with parton jet match", 5, 0, 5);
    TH1 *numMatchesGenRecoHist = new TH1I("numMatchesGenRecoHist", "Number of gen jets with reco jet match", 5, 0, 5);
    TH1 *numMatchesRecoPartonHist = new TH1I("numMatchesRecoPartonHist", "Number of reco Jets with parton jet match (dR < 0.35) for each reco jet with Pt > 30GeV", 5, 0, 5);
    TH1 *numMatchesRecoGenHist = new TH1I("numMatchesRecoGenHist", "Number of reco jes with gen jet match", 5, 0, 5);

    TH1 *partonPtNoRecoMatchHist = new TH1F("partonPtNoRecoMatchHist", "Pt distribution of parton jets with no reco match", 200, 0, 200);
    TH1 *recoPtNoPartonMatchHist = new TH1F("recoPtNoPartonMatchHist", "Pt distribution of reco jets with no parton match", 200, 0, 200);

    TH1 *minDRRecoPartonNoMatchHist = new TH1F("minDRRecoPartonNoMatchHist", "minDR for each reco jet with no parton match", 50, 0, 5);

    TH1 *emFractionRecoPartonNoMatchHist = new TH1F("emFractionRecoPartonNoMatchHist", "EM fraction of reco jets with no parton Match", 50, 0, 1);
    TH1 *emFractionRecoPartonMatchHist = new TH1F("emFractionRecoPartonMatchHist", "EM fraction of reco jets with parton match", 50, 0, 1);

    TH1 *partonPtHist = new TH1F("partonPtHist", "parton pt", 200, 0, 200);

    TH1 *recoEtaDistributionPartonMatchHist = new TH1F("recoEtaDistributionPartonMatchHist", "Eta distribution of reco jets with a parton match", 50, -5, 5);
    TH1 *recoEtaDistributionPtDiscrepencyHist = new TH1F("recoEtaDistributionPtDiscrepencyHist", "Eta distribution of reco jet with recoPt/partonPt < 0.75", 50, -5, 5);

    TH1 *emFractionRecoPartonMatchPtDiscrepencyHist = new TH1F("emFractionRecoPartonMatchPtDiscrepencyHist",
                                                               "EM Fraction of reco jets with recoPt/partonPt < 0.75", 50, 0, 1);

    std::vector<Float_t> *pfJetPt = 0;
    std::vector<Float_t> *pfJetEta = 0;
    std::vector<Float_t> *pfJetPhi = 0;
    std::vector<Float_t> *pfJetE = 0;
    std::vector<Float_t> *pfJetPhotonEnergy = 0;
    std::vector<Float_t> *pfJetElectronEnergy = 0;
    std::vector<Float_t> *pfJetMuonEnergy = 0;

    std::vector<Float_t> *genJetPt = 0;
    std::vector<Float_t> *genJetEta = 0;
    std::vector<Float_t> *genJetPhi = 0;
    std::vector<Float_t> *genJetE = 0;

    std::vector<Float_t> *partonPx = 0;
    std::vector<Float_t> *partonPy = 0;
    std::vector<Float_t> *partonPz = 0;
    std::vector<Float_t> *partonE = 0;
    std::vector<Float_t> *partonPt = 0;
    std::vector<Float_t> *partonEta = 0;
    std::vector<Float_t> *partonPhi = 0;
    std::vector<Int_t> *partonPdgId = 0;
    std::vector<Int_t> *partonStatus = 0;

    std::vector<Float_t> *hadronPx = 0;
    std::vector<Float_t> *hadronPy = 0;
    std::vector<Float_t> *hadronPz = 0;
    std::vector<Float_t> *hadronPt = 0;
    std::vector<Float_t> *hadronEta = 0;
    std::vector<Float_t> *hadronPhi = 0;
    std::vector<Float_t> *hadronE = 0;
    std::vector<Int_t> *hadronStatus = 0;
    std::vector<Int_t> *hadronPdgId = 0;

    tree->SetBranchAddress("pfJetPt", &pfJetPt);
    tree->SetBranchAddress("pfJetEta", &pfJetEta);
    tree->SetBranchAddress("pfJetPhi", &pfJetPhi);
    tree->SetBranchAddress("pfJetE", &pfJetE);
    tree->SetBranchAddress("pfJetPhotonEnergy", &pfJetPhotonEnergy);
    tree->SetBranchAddress("pfJetElectronEnergy", &pfJetElectronEnergy);
    tree->SetBranchAddress("pfJetMuonEnergy", &pfJetMuonEnergy);

    tree->SetBranchAddress("genJetPt", &genJetPt);
    tree->SetBranchAddress("genJetEta", &genJetEta);
    tree->SetBranchAddress("genJetPhi", &genJetPhi);
    tree->SetBranchAddress("genJetE", &genJetE);

    tree->SetBranchAddress("partonPx", &partonPx);
    tree->SetBranchAddress("partonPy", &partonPy);
    tree->SetBranchAddress("partonPz", &partonPz);
    tree->SetBranchAddress("partonE", &partonE);
    tree->SetBranchAddress("partonPt", &partonPt);
    tree->SetBranchAddress("partonEta", &partonEta);
    tree->SetBranchAddress("partonPhi", &partonPhi);
    tree->SetBranchAddress("partonPdgId", &partonPdgId);
    tree->SetBranchAddress("partonStatus", &partonStatus);

    tree->SetBranchAddress("hadronPx", &hadronPx);
    tree->SetBranchAddress("hadronPy", &hadronPy);
    tree->SetBranchAddress("hadronPz", &hadronPz);
    tree->SetBranchAddress("hadronE", &hadronE);
    tree->SetBranchAddress("hadronEta", &hadronEta);
    tree->SetBranchAddress("hadronPhi", &hadronPhi);
    tree->SetBranchAddress("hadronStatus", &hadronStatus);
    tree->SetBranchAddress("hadronPdgId", &hadronPdgId);
    tree->SetBranchAddress("hadronPt", &hadronPt);

    //pfCandTree variables
    int numEvents = tree->GetEntries();

    std::cout << "Num Events: " << numEvents << std::endl;

    int numPartonJets = 0;
    int numPartonJetswMatch = 0;
    int numRecoJets = 0;
    int numRecoJetswMatch = 0;

    // loop through events
    for (int i = 0; i < numEvents; i++)
    {
        tree->GetEntry(i);

        // create vector with all parton 4-momenta
        int numPartons = partonPx->size();
        float partonPxTot = 0;
        float partonPyTot = 0;
        float partonPzTot = 0;
        float partonETot = 0;
        std::vector<PseudoJet> particles;
        for (int j = 0; j < numPartons; j++)
        {
                particles.push_back(PseudoJet((*partonPx)[j], (*partonPy)[j], (*partonPz)[j], (*partonE)[j]));
                partonPxTot += (*partonPx)[j];
                partonPyTot += (*partonPy)[j];
                partonPzTot += (*partonPz)[j];
                partonETot += (*partonE)[j];
                partonPtHist->Fill((*partonPt)[j]);
                write_out_particles << i << " "
                                    << (*partonStatus)[j] << " "
                                    << (*partonPdgId)[j] << " "
                                    << (*partonPt)[j] << " "
                                    << (*partonEta)[j] << " "
                                    << (*partonPhi)[j] << " "
                                    << (*partonE)[j] << std::endl;
        }
        std::cout << partonETot << std::endl;

        int numHadrons = hadronPx->size();
        float hadronPxTot = 0;
        float hadronPyTot = 0;
        float hadronPzTot = 0;
        float hadronETot = 0;
        for (int j = 0; j < numHadrons; j++)
        {
            hadronPxTot += (*hadronPx)[j];
            hadronPyTot += (*hadronPy)[j];
            hadronPzTot += (*hadronPz)[j];
            hadronETot += (*hadronE)[j];
            write_out_particles << i << " "
                                << (*hadronStatus)[j] << " "
                                << (*hadronPdgId)[j] << " "
                                << (*hadronPt)[j] << " "
                                << (*hadronEta)[j] << " "
                                << (*hadronPhi)[j] << " "
                                << (*hadronE)[j] << std::endl;
        }

        // cluster partons into jet with anti-kt algorithm
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);

        std::vector<PseudoJet> partonJets = cs.inclusive_jets(); // get new parton jets

        // go through each parton jet and find how many reco and gen jet
        // matches there are (and find min dR)
        if (partonJets.size() > 0)
        {
            for (size_t j = 0; j < partonJets.size(); j++)
            {
                write_out_jets << i << " "
                               << 0 << " "
                               << partonJets[j].pt() << " "
                               << partonJets[j].eta() << " "
                               << partonJets[j].phi_std() << " "
                               << partonJets[j].E() << std::endl;
                if (partonJets[j].pt() >= 20.0)
                {
                    numPartonJets++;

                    int genJetMatches = 0;
                    float minDR = 10.0;
                    for (size_t k = 0; k < genJetPt->size(); k++)
                    {
                        if ((*genJetPt)[k] >= 30.0)
                        {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], partonJets[j].rap(), partonJets[j].phi_std());
                            if (dR < 0.35)
                                genJetMatches++;
                            if (dR < minDR)
                                minDR = dR;
                        }
                    }
                    numMatchesPartonGenHist->Fill(genJetMatches);
                    int pfJetMatches = 0;
                    minDR = 10.0;
                    for (size_t k = 0; k < pfJetPt->size(); k++)
                    {
                        if ((*pfJetPt)[k] >= 30.0)
                        {
                            float dR = deltaR((*pfJetEta)[k], (*pfJetPhi)[k], partonJets[j].rap(), partonJets[j].phi_std());
                            if (dR < 0.35)
                                pfJetMatches++;
                            if (dR < minDR)
                                minDR = dR;
                        }
                    }
                    numMatchesPartonRecoHist->Fill(pfJetMatches);
                    if (pfJetMatches == 1)
                        numPartonJetswMatch++;
                    if (pfJetMatches == 0)
                        partonPtNoRecoMatchHist->Fill(partonJets[j].pt());
                }
            }
        }

        // go through each gen  jet and find how many parton and reco jet
        // matches there are (and find min dR)
        if (genJetPt->size() > 0)
        {
            for (size_t j = 0; j < genJetPt->size(); j++)
            {
                write_out_jets << i << " "
                               << 1 << " "
                               << (*genJetPt)[j] << " "
                               << (*genJetEta)[j] << " "
                               << (*genJetPhi)[j] << " "
                               << (*genJetE)[j] << std::endl;
                if ((*genJetPt)[j] >= 30.0)
                {
                    int partonJetMatches = 0;
                    float minDR = 10.0;
                    for (size_t k = 0; k < partonJets.size(); k++)
                    {
                        if (partonJets[k].pt() >= 20.0)
                        {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), (*genJetEta)[j], (*genJetPhi)[j]);
                            if (dR < 0.35)
                                partonJetMatches++;
                            if (dR < minDR)
                                minDR = dR;
                        }
                    }
                    numMatchesGenPartonHist->Fill(partonJetMatches);
                    int pfJetMatches = 0;
                    minDR = 10.0;
                    for (size_t k = 0; k < pfJetPt->size(); k++)
                    {
                        if ((*pfJetPt)[k] >= 30.0)
                        {
                            float dR = deltaR((*pfJetEta)[k], (*pfJetPhi)[k], (*genJetEta)[j], (*genJetPhi)[j]);
                            if (dR < 0.35)
                                pfJetMatches++;
                            if (dR < minDR)
                                minDR = dR;
                        }
                    }
                    numMatchesGenRecoHist->Fill(pfJetMatches);
                }
            }
        }

        if (pfJetPt->size() > 0)
        {
            for (size_t j = 0; j < pfJetPt->size(); j++)
            {
                write_out_jets << i << " "
                               << 2 << " "
                               << (*pfJetPt)[j] << " "
                               << (*pfJetEta)[j] << " "
                               << (*pfJetPhi)[j] << " "
                               << (*pfJetE)[j] << std::endl;
                if ((*pfJetPt)[j] > 30.0)
                {
                    numRecoJets++;
                    int partonJetMatches = 0;
                    float minDR = 10.0;
                    int matchIndex = -1;
                    for (size_t k = 0; k < partonJets.size(); k++)
                    {
                        if (partonJets[k].pt() >= 20.0)
                        {
                            float dR = deltaR(partonJets[k].rap(), partonJets[k].phi_std(), (*pfJetEta)[j], (*pfJetPhi)[j]);
                            if (dR < 0.35)
                                partonJetMatches++;
                            if (dR < minDR)
                            {
                                minDR = dR;
                                matchIndex = k;
                            }
                        }
                    }
                    numMatchesRecoPartonHist->Fill(partonJetMatches);
                    if (partonJetMatches > 0)
                    {
                        numRecoJetswMatch++;
                        float emFraction = ((*pfJetPhotonEnergy)[j] + (*pfJetElectronEnergy)[j] + (*pfJetMuonEnergy)[j]) / (*pfJetE)[j];
                        emFractionRecoPartonMatchHist->Fill(emFraction);
                        recoEtaDistributionPartonMatchHist->Fill((*pfJetEta)[j]);
                        float ptRatio = (*pfJetEta)[j] / partonJets[matchIndex].rap();
                        if (ptRatio < 0.75)
                        {
                            recoEtaDistributionPtDiscrepencyHist->Fill((*pfJetEta)[j]);
                            emFractionRecoPartonMatchPtDiscrepencyHist->Fill(emFraction);
                        }
                    }
                    if (partonJetMatches == 0)
                    {
                        recoPtNoPartonMatchHist->Fill((*pfJetPt)[j]);
                        minDRRecoPartonNoMatchHist->Fill(minDR);
                        float emFraction = ((*pfJetPhotonEnergy)[j] + (*pfJetElectronEnergy)[j] + (*pfJetMuonEnergy)[j]) / (*pfJetE)[j];
                        emFractionRecoPartonNoMatchHist->Fill(emFraction);
                    }

                    int genJetMatches = 0;
                    minDR = 10.0;
                    for (size_t k = 0; k < genJetPt->size(); k++)
                    {
                        if ((*genJetPt)[k] >= 20.0)
                        {
                            float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k], (*pfJetEta)[j], (*pfJetPhi)[j]);
                            if (dR < 0.35)
                                genJetMatches++;
                        }
                    }
                    numMatchesRecoGenHist->Fill(genJetMatches);
                }
            }
        }
    }

    numMatchesPartonRecoHist->Write();
    numMatchesPartonGenHist->Write();
    numMatchesGenPartonHist->Write();
    numMatchesGenRecoHist->Write();
    numMatchesRecoPartonHist->Write();
    numMatchesRecoGenHist->Write();

    partonPtNoRecoMatchHist->Write();
    recoPtNoPartonMatchHist->Write();

    minDRRecoPartonNoMatchHist->Write();

    emFractionRecoPartonNoMatchHist->Write();
    emFractionRecoPartonMatchHist->Write();

    partonPtHist->Write();

    recoEtaDistributionPartonMatchHist->Write();
    recoEtaDistributionPtDiscrepencyHist->Write();

    emFractionRecoPartonMatchPtDiscrepencyHist->Write();

    delete numMatchesPartonRecoHist;
    delete numMatchesPartonGenHist;
    delete numMatchesGenPartonHist;
    delete numMatchesGenRecoHist;
    delete numMatchesRecoPartonHist;
    delete numMatchesRecoGenHist;

    delete partonPtNoRecoMatchHist;
    delete recoPtNoPartonMatchHist;

    delete minDRRecoPartonNoMatchHist;

    delete emFractionRecoPartonNoMatchHist;
    delete emFractionRecoPartonMatchHist;

    delete partonPtHist;

    delete recoEtaDistributionPartonMatchHist;
    delete recoEtaDistributionPtDiscrepencyHist;

    delete emFractionRecoPartonMatchPtDiscrepencyHist;

    write_out_particles.close();
    write_out_jets.close();

    return 0;
}
