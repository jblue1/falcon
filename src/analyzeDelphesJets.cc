#include "ExRootAnalysis/ExRootTreeReader.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "fastjet/ClusterSequence.hh"
#include "helpers.h"
#include <cassert>
#include <fstream>
#include <iostream>

using namespace fastjet;

/** Cluster a given set of particles into jets using the anti-kt algorithm with
 * R = 0.4.
 */
std::vector<PseudoJet> cluster_jets(std::vector<Float_t> *px,
                                    std::vector<Float_t> *py,
                                    std::vector<Float_t> *pz,
                                    std::vector<Float_t> *energy) {
  int numberOfParticles = px->size();
  std::vector<PseudoJet> particles;
  for (int j = 0; j < numberOfParticles; j++) {
    particles.push_back(PseudoJet((*px)[j], (*py)[j], (*pz)[j], (*energy)[j]));
    particles[j].set_user_index(j);
  }

  // make jets
  double R = 0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  ClusterSequence cs(particles, jet_def);

  std::vector<PseudoJet> jets = cs.inclusive_jets();

  return jets;
}

int main(int argc, char const *argv[]) {
  std::ofstream write_out("test.txt");
  assert(write_out.is_open());

  write_out
      << "The file contains 4 vectors of matched parton, reco, and gen jets. "
         "The columns are"
      << std::endl;
  write_out
      << "Parton Pt | Parton Eta | Parton Phi | Delphes Reco Pt | Delphes Reco "
         "Eta | Delphes Reco Phi | Gen Pt | Gen Eta | Gen Phi | Reco Pt | Reco "
         "Eta |"
         " Reco Phi | PfCand Pt | PfCand Eta | PfCand Phi"
      << std::endl;

  /* Set up variables for extracting the Delphes-reco jets */
  TChain chain("Delphes");
  chain.Add(argv[1]);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfDelphesEvents = treeReader->GetEntries();
  TClonesArray *branchRecoJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

  /* Set up variables for extracting the parton jets */
  std::string dataFile(argv[2]);
  // open data file and get trees
  TFile f(dataFile.c_str());
  TDirectoryFile *df = (TDirectoryFile *)f.Get("demo");
  TTree *partonTree = (TTree *)df->Get("eventTree");

  std::vector<Float_t> *partonPx = 0;
  std::vector<Float_t> *partonPy = 0;
  std::vector<Float_t> *partonPz = 0;
  std::vector<Float_t> *partonE = 0;

  std::vector<Float_t> *pfCandPx = 0;
  std::vector<Float_t> *pfCandPy = 0;
  std::vector<Float_t> *pfCandPz = 0;
  std::vector<Float_t> *pfCandE = 0;

  std::vector<Float_t> *genJetPt = 0;
  std::vector<Float_t> *genJetEta = 0;
  std::vector<Float_t> *genJetPhi = 0;
  std::vector<Float_t> *genJetE = 0;

  std::vector<Float_t> *recoJetPt = 0;
  std::vector<Float_t> *recoJetEta = 0;
  std::vector<Float_t> *recoJetPhi = 0;
  std::vector<Float_t> *recoJetE = 0;

  partonTree->SetBranchAddress("partonPx", &partonPx);
  partonTree->SetBranchAddress("partonPy", &partonPy);
  partonTree->SetBranchAddress("partonPz", &partonPz);
  partonTree->SetBranchAddress("partonE", &partonE);

  partonTree->SetBranchAddress("pfCandPx", &pfCandPx);
  partonTree->SetBranchAddress("pfCandPy", &pfCandPy);
  partonTree->SetBranchAddress("pfCandPz", &pfCandPz);
  partonTree->SetBranchAddress("pfCandE", &pfCandE);

  partonTree->SetBranchAddress("genJetPt", &genJetPt);
  partonTree->SetBranchAddress("genJetEta", &genJetEta);
  partonTree->SetBranchAddress("genJetPhi", &genJetPhi);
  partonTree->SetBranchAddress("genJetE", &genJetE);

  partonTree->SetBranchAddress("pfJetPt", &recoJetPt);
  partonTree->SetBranchAddress("pfJetEta", &recoJetEta);
  partonTree->SetBranchAddress("pfJetPhi", &recoJetPhi);
  partonTree->SetBranchAddress("pfJetE", &recoJetE);

  int numberOfCMSSWEvents = partonTree->GetEntries();

  assert(numberOfCMSSWEvents == numberOfDelphesEvents);

  for (int i = 0; i < numberOfDelphesEvents; i++) {
    partonTree->GetEntry(i);
    treeReader->ReadEntry(i);
    int delphesRecoJetEntries = branchRecoJet->GetEntries();
    int delphesGenJetEntries = branchGenJet->GetEntries();

    std::vector<PseudoJet> partonJets =
        cluster_jets(partonPx, partonPy, partonPz, partonE);

    std::vector<PseudoJet> pfCandJets =
        cluster_jets(pfCandPx, pfCandPy, pfCandPz, pfCandE);

    if (partonJets.size() > 0) {
      // loop through parton jets looking for matches

      for (int j = 0; j < partonJets.size(); j++) {

        if (partonJets[j].pt() > 20) {
          float minDRDelphesRecoJet = 10.0;
          int delphesRecoJetIndex = 0;

          float minDRGenJet = 10.0;
          int genJetIndex = 0;

          float minDRDelphesGenJet = 10.0;
          int delphesGenJetIndex = 0;

          float minDRRecoJet = 10.0;
          int recoJetIndex = 0;

          float minDRPfCandJet = 10.0;
          int pfCandJetIndex = 0;

          for (int k = 0; k < delphesRecoJetEntries; k++) {
            Jet *jet = (Jet *)branchRecoJet->At(k);
            float dR = deltaR(jet->Eta, jet->Phi, partonJets[j].rap(),
                              partonJets[j].phi_std());

            if (dR < minDRDelphesRecoJet) {
              minDRDelphesRecoJet = dR;
              delphesRecoJetIndex = k;
            }
          }

          for (int k = 0; k < recoJetPt->size(); k++) {
            if ((*recoJetPt)[k] > 10.0) {
              float dR = deltaR((*recoJetEta)[k], (*recoJetPhi)[k],
                                partonJets[j].rap(), partonJets[j].phi_std());

              if (dR < minDRRecoJet) {
                minDRRecoJet = dR;
                recoJetIndex = k;
              }
            }
          }

          if (minDRDelphesRecoJet < 0.35 && minDRRecoJet < 0.35) {
            Jet *matchedDelphesRecoJet =
                (Jet *)branchRecoJet->At(delphesRecoJetIndex);
            Jet *matchedDelphesGenJet =
                (Jet *)branchGenJet->At(delphesGenJetIndex);
            write_out << partonJets[j].pt() << " " << partonJets[j].rap() << " "
                      << partonJets[j].phi_std() << " "
                      << partonJets[j].e() << " "
                      << matchedDelphesRecoJet->PT << " "
                      << matchedDelphesRecoJet->Eta << " "
                      << matchedDelphesRecoJet->Phi << " "
                      << (*recoJetPt)[recoJetIndex] << " "
                      << (*recoJetEta)[recoJetIndex] << " "
                      << (*recoJetPhi)[recoJetIndex] << " " << std::endl;
          }
        }
      }
    }
  }

  return 0;
}
