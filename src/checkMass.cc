/** C++ code to check the momenta and masses of particles and jets in CMS ttbar
 * events.
 */

#include "TFile.h"
#include "TTree.h"
#include "fastjet/ClusterSequence.hh"
#include <cmath>
#include <fstream>

using namespace fastjet;

/**
 * Describe parameters taken and print errors if parameters are incorrectly
 * given.
 */
void usage(std::ostream &out, const char *msg) {
  out << msg << std::endl;
  out << std::endl;
  out << "    Usage:" << std::endl;
  out << "            makeJets.out data" << std::endl;
  out << "    data        - root file from JetNtuple analyzer to use"
      << std::endl;
  exit(1);
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
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

  tree->SetBranchAddress("pfJetPt", &pfJetPt);
  tree->SetBranchAddress("pfJetEta", &pfJetEta);
  tree->SetBranchAddress("pfJetPhi", &pfJetPhi);
  tree->SetBranchAddress("pfJetPx", &pfJetPx);
  tree->SetBranchAddress("pfJetPy", &pfJetPy);
  tree->SetBranchAddress("pfJetPz", &pfJetPz);
  tree->SetBranchAddress("pfJetE", &pfJetE);

  for (int i = 0; i < 1; i++) {
    tree->GetEntry(i);
    int numJets = pfJetPx->size();

    for (int j = 0; j < numJets; j++) {
      std::cout << "Jet Mass: "
                << pow((*pfJetE)[j], 2) - pow((*pfJetPx)[j], 2) -
                       pow((*pfJetPy)[j], 2) - pow((*pfJetPz)[j], 2)
                << std::endl;
    }
  }

  exit(0);
}