/** C++ code to check the momenta and masses of particles and jets in CMS ttbar
 * events.
 */

#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "fastjet/ClusterSequence.hh"
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
  out << "            makeJets.out data" << std::endl;
  out << "    data        - root file from JetNtuple analyzer to use"
      << std::endl;
  exit(1);
}

int main(int argc, char const *argv[])
{
  if (argc != 2)
  {
    usage(std::cerr, "Incorrect number of parameters given.");
  }

  std::string dataFile(argv[1]);
  std::string dataPath = "./data/raw/" + dataFile;

  // open data file and get trees
  TFile f(dataPath.c_str());
  TDirectoryFile *df = (TDirectoryFile *)f.Get("demo");
  TTree *tree = (TTree *)df->Get("eventTree");

  std::vector<Float_t> *pfJetPt = 0;
  std::vector<Float_t> *pfJetEta = 0;
  std::vector<Float_t> *pfJetPhi = 0;
  std::vector<Float_t> *pfJetPx = 0;
  std::vector<Float_t> *pfJetPy = 0;
  std::vector<Float_t> *pfJetPz = 0;
  std::vector<Float_t> *pfJetE = 0;

  std::vector<Float_t> *partonPx = 0;
  std::vector<Float_t> *partonPy = 0;
  std::vector<Float_t> *partonPz = 0;
  std::vector<Float_t> *partonE = 0;

  tree->SetBranchAddress("pfJetPt", &pfJetPt);
  tree->SetBranchAddress("pfJetEta", &pfJetEta);
  tree->SetBranchAddress("pfJetPhi", &pfJetPhi);
  tree->SetBranchAddress("pfJetPx", &pfJetPx);
  tree->SetBranchAddress("pfJetPy", &pfJetPy);
  tree->SetBranchAddress("pfJetPz", &pfJetPz);
  tree->SetBranchAddress("pfJetE", &pfJetE);

  tree->SetBranchAddress("partonPx", &partonPx);
  tree->SetBranchAddress("partonPy", &partonPy);
  tree->SetBranchAddress("partonPz", &partonPz);
  tree->SetBranchAddress("partonE", &partonE);

  int numEvents = tree->GetEntries();

  std::ofstream write_out("./data/processed/JetMasses.txt");
  write_out.precision(std::numeric_limits<double>::digits10);

  write_out << "p1 = px^2 + py^2 + pz^2" << std::endl;
  write_out << "partonPt partonEta partonPhi partonE partonM^2 partonP1 recoPt recoEta recoPhi recoE recoP1" << std::endl;

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
            write_out
                << partonJets[j].pt() << " "
                << partonJets[j].rap() << " "
                << partonJets[j].phi_std() << " "
                << partonJets[j].e() << " "
                << partonJets[j].m2() << " "
                << sqrt(pow(partonJets[j].px(), 2) + pow(partonJets[j].py(), 2) + pow(partonJets[j].pz(), 2)) << " "
                << (*pfJetPt)[pfJetIndex] << " "
                << (*pfJetEta)[pfJetIndex] << " "
                << (*pfJetPhi)[pfJetIndex] << " "
                << (*pfJetE)[pfJetIndex] << " "
                << sqrt(pow((*pfJetPx)[pfJetIndex], 2) + pow((*pfJetPy)[pfJetIndex], 2) + pow((*pfJetPz)[pfJetIndex], 2)) << std::endl;
          }
        }
      }
    }
  }
  write_out.close();

  exit(0);
}
