#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "fastjet/ClusterSequence.hh"
#include "helpers.h"

using namespace fastjet;

int main(int argc, char const *argv[]) {
  std::ofstream write_out("test_CMS_Phase_II.txt");
  assert(write_out.is_open());

  write_out << "The file contains 4 vectors of matched parton and reco jets. "
               "The columns are"
            << std::endl;
  write_out << "Pt_p Eta_p Phi_p E_p Pt_r Eta_r Phi_r E_r" << std::endl;

  /* Set up variables for extracting the Delphes-reco jets */
  TChain chain("Delphes");
  chain.Add(argv[1]);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfDelphesEvents = treeReader->GetEntries();
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

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

  std::vector<Float_t> *genJetPt = 0;
  std::vector<Float_t> *genJetEta = 0;
  std::vector<Float_t> *genJetPhi = 0;
  std::vector<Float_t> *genJetE = 0;

  partonTree->SetBranchAddress("partonPx", &partonPx);
  partonTree->SetBranchAddress("partonPy", &partonPy);
  partonTree->SetBranchAddress("partonPz", &partonPz);
  partonTree->SetBranchAddress("partonE", &partonE);

  partonTree->SetBranchAddress("genJetPt", &genJetPt);
  partonTree->SetBranchAddress("genJetEta", &genJetEta);
  partonTree->SetBranchAddress("genJetPhi", &genJetPhi);
  partonTree->SetBranchAddress("genJetE", &genJetE);

  int numberOfCMSSWEvents = partonTree->GetEntries();

  assert(numberOfCMSSWEvents == numberOfDelphesEvents);

  for (int i = 0; i < numberOfCMSSWEvents; i++) {
    partonTree->GetEntry(i);
    int numPartons = partonPx->size();
    std::vector<PseudoJet> partons;
    for (int j = 0; j < numPartons; j++) {
      partons.push_back(PseudoJet((*partonPx)[j], (*partonPy)[j],
                                  (*partonPz)[j], (*partonE)[j]));
      partons[j].set_user_index(j);
    }

    // make parton jets
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(partons, jet_def);

    std::vector<PseudoJet> partonJets = cs.inclusive_jets();

    if (partonJets.size() > 0) {
      // loop through parton jets looking for matches

      for (int j = 0; j < partonJets.size(); j++) {
        treeReader->ReadEntry(i);
        int jetEntries = branchJet->GetEntries();

        if (partonJets[j].pt() > 20) {
          float minDRDelphesJet = 10.0;
          int delphesJetIndex = 0;

          float minDRGenJet = 10.0;
          int genJetIndex = 0;
          for (int k = 0; k < jetEntries; k++) {
            Jet *jet = (Jet *)branchJet->At(k);
            float dR = deltaR(jet->Eta, jet->Phi, partonJets[j].rap(),
                              partonJets[j].phi_std());

            if (dR < minDRDelphesJet) {
              minDRDelphesJet = dR;
              delphesJetIndex = k;
            }
          }

          for (int k = 0; k < genJetPt->size(); k++) {
            if ((*genJetPt)[k] > 10.0) {
              float dR = deltaR((*genJetEta)[k], (*genJetPhi)[k],
                                partonJets[j].rap(), partonJets[j].phi_std());

              if (dR < minDRGenJet) {
                minDRGenJet = dR;
                genJetIndex = k;
              }
            }
          }
          if (minDRDelphesJet < 0.35 && minDRGenJet < 0.35) {
            Jet *matchedJet = (Jet *)branchJet->At(delphesJetIndex);
            write_out << partonJets[j].pt() << " " << partonJets[j].rap() << " "
                      << partonJets[j].phi_std() << " "
                      << " " << matchedJet->PT << " " << matchedJet->Eta << " "
                      << matchedJet->Phi << " " << (*genJetPt)[genJetIndex]
                      << " " << (*genJetEta)[genJetIndex] << " "
                      << (*genJetPhi)[genJetIndex] << std::endl;
          }
        }
      }
    }
  }

  return 0;
}